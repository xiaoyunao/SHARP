from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path


@dataclass(frozen=True)
class SiteConfig:
    lat_deg: float = 40.393
    lon_deg: float = 117.575
    height_m: float = 960.0
    timezone_name: str = "Asia/Shanghai"
    mpc_code: str = "327"


@dataclass(frozen=True)
class SchedulerConfig:
    exposure_time_s: float = 30.0
    readout_overhead_s: float = 15.0
    slew_overhead_s: float = 35.0
    max_block_duration_s: float = 1800.0
    visibility_window_minutes: int = 90
    visibility_min_fraction: float = 0.5
    min_altitude_deg: float = 30.0
    cluster_field_count: int = 45
    revisit_count: int = 3
    revisit_spacing_s: float = 1800.0
    moon_sep_new_deg: float = 2.0
    moon_sep_half_deg: float = 25.0
    moon_sep_full_deg: float = 40.0
    moon_sep_margin_deg: float = 0.50
    near_sun_revisit_count: int = 3
    near_sun_trigger_minutes: int = 60
    w_visible_days: float = 1.0
    w_altitude: float = 1.2
    w_history: float = 2.5
    w_recent_night: float = 1.8
    w_moon_clearance: float = 0.4
    recent_cooldown_nights: int = 3


@dataclass(frozen=True)
class ServerPaths:
    root_dir: Path
    workspace: Path
    footprints_path: Path
    processed_root: Path
    publish_dir: Path | None = None

    @property
    def plans_dir(self) -> Path:
        return self.workspace / "plans"

    @property
    def history_dir(self) -> Path:
        return self.workspace / "history"

    @property
    def logs_dir(self) -> Path:
        return self.root_dir / "logs"

    @property
    def plots_dir(self) -> Path:
        return self.root_dir / "plots"

    def nightly_plots_dir(self, night_tag: str) -> Path:
        return self.plots_dir / night_tag

    @property
    def history_path(self) -> Path:
        return self.history_dir / "exposure_history.fits"

    @property
    def ingest_state_path(self) -> Path:
        return self.history_dir / "l2_ingest_state.json"

    def ensure(self) -> None:
        for directory in [self.root_dir, self.workspace, self.plans_dir, self.history_dir, self.logs_dir, self.plots_dir]:
            directory.mkdir(parents=True, exist_ok=True)
        if self.publish_dir is not None:
            self.publish_dir.mkdir(parents=True, exist_ok=True)


def load_ingest_state(path: Path) -> dict[str, dict[str, int]]:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def save_ingest_state(path: Path, state: dict[str, dict[str, int]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(state, indent=2, sort_keys=True) + "\n", encoding="utf-8")
