from __future__ import annotations

from dataclasses import dataclass
from datetime import date as dt_date
from datetime import datetime

import astropy.units as u
import numpy as np
from astroplan import moon
from astropy.coordinates import AltAz, SkyCoord, get_body, get_sun
from astropy.table import Table
from astropy.time import Time

from .astro_utils import apparent_coords, build_site, local_night_window_utc, moon_clearance_ok, stable_visibility_mask
from .config import SchedulerConfig, SiteConfig


@dataclass(frozen=True)
class CandidateCluster:
    center_index: int
    field_indices: tuple[int, ...]
    score: float


class StripScheduler:
    def __init__(self, footprints: Table, history: Table, site_config: SiteConfig, scheduler_config: SchedulerConfig):
        self.footprints = footprints.copy()
        self.history = history.copy()
        self.site_config = site_config
        self.config = scheduler_config
        self.site = build_site(site_config.lat_deg, site_config.lon_deg, site_config.height_m)
        self.coords = SkyCoord(
            ra=np.asarray(self.footprints["center_ra"], dtype=float) * u.deg,
            dec=np.asarray(self.footprints["center_dec"], dtype=float) * u.deg,
        )
        self._active_night_tag: str | None = None

    def _wrapped_ra_diff_deg(self, idx_a: int, idx_b: int) -> float:
        ra_a = float(self.footprints[idx_a]["center_ra"])
        ra_b = float(self.footprints[idx_b]["center_ra"])
        return abs(((ra_b - ra_a + 180.0) % 360.0) - 180.0)

    def _slew_overhead_s(self, from_idx: int | None, to_idx: int) -> float:
        return float(self.config.slew_overhead_s)

    def _moon_metadata(self, obstime: Time, idx: int) -> tuple[float, float]:
        moon_coord = apparent_coords(get_body("moon", obstime, location=self.site), obstime, self.site)
        field_coord = apparent_coords(self.coords[idx], obstime, self.site)
        moon_fraction = float(moon.moon_illumination(obstime))
        moon_sep_deg = float(field_coord.separation(moon_coord).deg)
        return moon_fraction, moon_sep_deg

    def _exposure_end_time(self, start_time: Time) -> Time:
        return start_time + (self.config.readout_overhead_s + self.config.exposure_time_s) * u.s

    def _simulate_sequence_end(
        self,
        sequence: list[int],
        start_time: Time,
        previous_idx: int | None,
        repeats: int,
        night_end: Time,
    ) -> tuple[Time, int | None, bool]:
        t_now = start_time
        last_idx = previous_idx
        for _ in range(repeats):
            for idx in sequence:
                slew_s = self._slew_overhead_s(last_idx, idx)
                exp_start = t_now + slew_s * u.s
                exp_end = self._exposure_end_time(exp_start)
                if exp_end > night_end:
                    return t_now, last_idx, False
                t_now = exp_end
                last_idx = idx
        return t_now, last_idx, True

    def _stable_mask(self, obstime: Time, window_minutes: int | None = None) -> np.ndarray:
        return stable_visibility_mask(
            self.coords,
            obstime,
            self.site,
            self.config.min_altitude_deg,
            self.config.visibility_window_minutes if window_minutes is None else window_minutes,
            self.config.visibility_min_fraction,
            (
                self.config.moon_sep_new_deg,
                self.config.moon_sep_half_deg,
                self.config.moon_sep_full_deg,
                self.config.moon_sep_margin_deg,
            ),
        )

    def _weights(self, obstime: Time, stable_mask: np.ndarray) -> np.ndarray:
        altaz = self.coords.transform_to(AltAz(obstime=obstime, location=self.site))
        alt_score = np.clip((altaz.alt.deg - self.config.min_altitude_deg) / 30.0, 0.0, 1.0)
        visible_days = np.clip(np.asarray(self.footprints["visible_days"], dtype=float), 1.0, 365.0)
        visible_score = 1.0 - 0.9 * (visible_days - 1.0) / 364.0
        history_count = np.asarray(self.history["exposure_count"], dtype=float)
        if history_count.size:
            hist_score = 1.0 - history_count / max(float(np.max(history_count)), 1.0)
        else:
            hist_score = np.ones(len(self.footprints), dtype=float)
        recent_score = self._recent_night_score()
        moon_ok = moon_clearance_ok(
            self.coords,
            obstime,
            self.site,
            self.config.moon_sep_new_deg,
            self.config.moon_sep_half_deg,
            self.config.moon_sep_full_deg,
            self.config.moon_sep_margin_deg,
        ).astype(float)
        weight = (
            self.config.w_visible_days * visible_score
            + self.config.w_altitude * alt_score
            + self.config.w_history * hist_score
            + self.config.w_recent_night * recent_score
            + self.config.w_moon_clearance * moon_ok
        )
        weight[~stable_mask] = -np.inf
        return weight

    def _recent_night_score(self) -> np.ndarray:
        if self._active_night_tag is None:
            return np.ones(len(self.footprints), dtype=float)
        active_night = datetime.strptime(self._active_night_tag, "%Y%m%d").date()
        cooldown = max(int(self.config.recent_cooldown_nights), 1)
        raw = np.asarray(self.history["last_observed_night"]).astype("U8")
        score = np.ones(len(raw), dtype=float)
        for idx, value in enumerate(raw):
            value = value.strip()
            if not value or value == "--":
                continue
            try:
                days_since = (active_night - datetime.strptime(value, "%Y%m%d").date()).days
            except ValueError:
                continue
            score[idx] = float(np.clip(days_since / cooldown, 0.0, 1.0))
        return score

    def _cluster_indices_near_center(self, center_index: int, stable_mask: np.ndarray) -> np.ndarray:
        stable_idx = np.where(stable_mask)[0]
        if stable_idx.size == 0:
            return stable_idx
        stable_coords = self.coords[stable_idx]
        stable_sep = self.coords[center_index].separation(stable_coords).deg
        stable_order = np.argsort(stable_sep)
        target_size = min(self.config.cluster_field_count, stable_idx.size)
        return stable_idx[stable_order[:target_size]]

    def _order_cluster_fields(self, field_indices: tuple[int, ...]) -> list[int]:
        sequence = [int(idx) for idx in field_indices]
        if len(sequence) <= 1:
            return sequence
        ra = np.asarray([float(self.footprints[idx]["center_ra"]) for idx in sequence], dtype=float)
        order = np.argsort(ra)
        sorted_idx = [sequence[i] for i in order]
        sorted_ra = ra[order]
        cyclic_gaps = np.diff(np.concatenate([sorted_ra, [sorted_ra[0] + 360.0]]))
        start = int((np.argmax(cyclic_gaps) + 1) % len(sorted_idx))
        return sorted_idx[start:] + sorted_idx[:start]

    def select_best_cluster(self, obstime: Time) -> CandidateCluster:
        stable_mask = self._stable_mask(obstime)
        weights = self._weights(obstime, stable_mask)
        valid_idx = np.where(np.isfinite(weights))[0]
        if valid_idx.size == 0:
            raise RuntimeError("No valid fields found for this time.")
        center_index = int(valid_idx[np.argmax(weights[valid_idx])])
        chosen_idx = self._cluster_indices_near_center(center_index, stable_mask)
        score = float(np.sum(weights[chosen_idx]))
        return CandidateCluster(center_index=center_index, field_indices=tuple(int(idx) for idx in chosen_idx), score=score)

    def _append_near_sun_cycle(
        self,
        rows: list[dict],
        cycle_id: int,
        cycle_start: Time,
        night_end: Time,
        previous_idx: int | None,
    ) -> tuple[Time, int | None]:
        remaining_minutes = max(0.0, (night_end - cycle_start).to_value(u.min))
        stable_mask = self._stable_mask(cycle_start, int(np.floor(remaining_minutes)))
        if not np.any(stable_mask):
            return cycle_start, previous_idx
        sun_coord = get_sun(cycle_start)
        sep_to_sun = self.coords.separation(sun_coord).deg
        candidate_idx = np.where(stable_mask)[0]
        candidate_idx = candidate_idx[np.argsort(sep_to_sun[candidate_idx])]
        ordered_candidates = np.array(sorted(candidate_idx, key=lambda idx: float(self.footprints[idx]["center_ra"])), dtype=int)
        chosen: np.ndarray | None = None
        repeat_count = self.config.near_sun_revisit_count
        for trial_repeats in range(self.config.near_sun_revisit_count, 0, -1):
            best_size = 0
            for size in range(1, len(ordered_candidates) + 1):
                trial = ordered_candidates[:size].tolist()
                _, _, ok = self._simulate_sequence_end(trial, cycle_start, previous_idx, trial_repeats, night_end)
                if not ok:
                    break
                best_size = size
            if best_size > 0:
                chosen = ordered_candidates[:best_size]
                repeat_count = trial_repeats
                break
        if chosen is None or chosen.size == 0:
            return cycle_start, previous_idx
        t_now = cycle_start
        last_idx = previous_idx
        for rep in range(1, repeat_count + 1):
            for idx in chosen:
                slew_s = self._slew_overhead_s(last_idx, int(idx))
                start_time = t_now + slew_s * u.s
                exp_end = self._exposure_end_time(start_time)
                if exp_end > night_end:
                    return t_now, last_idx
                row = self.footprints[idx]
                moon_phase, moon_sep_deg = self._moon_metadata(start_time, int(idx))
                rows.append(
                    {
                        "field_id": str(row["field_id"]),
                        "cycle_id": cycle_id,
                        "repeat": rep,
                        "mode": "near_sun",
                        "start_utc": start_time.isot,
                        "ra": float(row["center_ra"]),
                        "dec": float(row["center_dec"]),
                        "dec_band": int(row["dec_band"]),
                        "ra_band": int(row["ra_band"]),
                        "strip_score": float(-sep_to_sun[idx]),
                        "strip_ra_min_deg": float(self.footprints[chosen[0]]["center_ra"]),
                        "strip_ra_max_deg": float(self.footprints[chosen[-1]]["center_ra"]),
                        "strip_dec_count": repeat_count,
                        "strip_ra_count": int(chosen.size),
                        "exptime_s": float(self.config.exposure_time_s),
                        "moon_phase_fraction": moon_phase,
                        "moon_sep_deg": moon_sep_deg,
                    }
                )
                t_now = exp_end
                last_idx = int(idx)
        return t_now, last_idx

    def _apply_cycle_to_history(self, field_counts: dict[int, int], night_tag: str) -> None:
        history_field_ids = np.asarray(self.history["field_id"]).astype("U16")
        idx_by_field = {field_id: i for i, field_id in enumerate(history_field_ids)}
        for idx, exposure_count in field_counts.items():
            field_id = str(self.footprints[idx]["field_id"])
            hidx = idx_by_field.get(field_id)
            if hidx is None:
                continue
            self.history["exposure_count"][hidx] += exposure_count
            self.history["last_observed_night"][hidx] = night_tag

    def _append_cluster_cycle(
        self,
        rows: list[dict],
        cluster: CandidateCluster,
        cycle_id: int,
        cycle_start: Time,
        night_end: Time,
        night_tag: str,
        previous_idx: int | None,
    ) -> tuple[Time, int | None]:
        sequence = self._order_cluster_fields(cluster.field_indices)
        t_now = cycle_start
        applied_counts: dict[int, int] = {}
        last_idx = previous_idx
        for repeat in range(1, self.config.revisit_count + 1):
            target_start = cycle_start + (repeat - 1) * self.config.revisit_spacing_s * u.s
            if t_now < target_start:
                t_now = target_start
            rep_block_start: Time | None = None
            for idx in sequence:
                slew_s = self._slew_overhead_s(last_idx, int(idx))
                exp_start = t_now + slew_s * u.s
                exp_end = self._exposure_end_time(exp_start)
                if rep_block_start is None:
                    rep_block_start = exp_start
                if (exp_end - rep_block_start).sec > self.config.max_block_duration_s:
                    break
                if exp_end > night_end:
                    self._apply_cycle_to_history(applied_counts, night_tag)
                    return t_now, last_idx
                row = self.footprints[idx]
                moon_phase, moon_sep_deg = self._moon_metadata(exp_start, int(idx))
                rows.append(
                    {
                        "field_id": str(row["field_id"]),
                        "cycle_id": cycle_id,
                        "repeat": repeat,
                        "mode": "normal",
                        "start_utc": exp_start.isot,
                        "ra": float(row["center_ra"]),
                        "dec": float(row["center_dec"]),
                        "dec_band": int(row["dec_band"]),
                        "ra_band": int(row["ra_band"]),
                        "strip_score": cluster.score,
                        "strip_ra_min_deg": float(self.footprints[sequence[0]]["center_ra"]),
                        "strip_ra_max_deg": float(self.footprints[sequence[-1]]["center_ra"]),
                        "strip_dec_count": 1,
                        "strip_ra_count": len(sequence),
                        "exptime_s": float(self.config.exposure_time_s),
                        "moon_phase_fraction": moon_phase,
                        "moon_sep_deg": moon_sep_deg,
                    }
                )
                applied_counts[idx] = applied_counts.get(idx, 0) + 1
                t_now = exp_end
                last_idx = int(idx)
        self._apply_cycle_to_history(applied_counts, night_tag)
        return t_now, last_idx

    def build_plan(self, night: dt_date) -> list[dict]:
        t0, t1 = local_night_window_utc(night, self.site)
        rows: list[dict] = []
        cycle_id = 1
        t_now = t0
        last_idx: int | None = None
        night_tag = night.strftime("%Y%m%d")
        self._active_night_tag = night_tag
        while t_now < t1:
            remaining_minutes = max(0.0, (t1 - t_now).to_value(u.min))
            if remaining_minutes < self.config.near_sun_trigger_minutes:
                break
            try:
                cluster = self.select_best_cluster(t_now)
            except RuntimeError:
                break
            cycle_end, last_idx = self._append_cluster_cycle(rows, cluster, cycle_id, t_now, t1, night_tag, last_idx)
            if cycle_end <= t_now:
                break
            cycle_id += 1
            t_now = cycle_end
        if t_now < t1 and max(0.0, (t1 - t_now).to_value(u.min)) < self.config.near_sun_trigger_minutes:
            near_sun_end, last_idx = self._append_near_sun_cycle(rows, cycle_id, t_now, t1, last_idx)
            if near_sun_end > t_now:
                t_now = near_sun_end
        self._active_night_tag = None
        return rows
