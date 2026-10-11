#!/usr/bin/env python3
"""Candidate-level morphology + empirical-PSF gate, with immutable raw catalogs.

The two majorities are independent: their abnormal exposures need not coincide.
Only catalog-suspicious candidates need image fitting. Missing evidence that can
change a majority is quarantined for review, never silently treated as a pass.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import tempfile

import numpy as np
from astropy.io import fits
from astropy.table import Table
from scipy.ndimage import map_coordinates, shift
from scipy.optimize import least_squares
from scipy.spatial import cKDTree

POLICY = {
    'version': 'joint-psf-tracklet-20261011', 'min_observations': 3,
    'fraction': '2/3', 'fwhm_ratio_gt': 1.5, 'class_ann_lt': 0.1,
    'psf_width_gt': 1.3, 'delta_chi2_ge': 9.0,
    'local_reference_min': 20, 'local_radii_pix': [1000, 1800],
    'psf_train': 20, 'psf_validate_min': 3, 'psf_validate_max': 5,
    'psf_radius_pix': 2000, 'fit_radius_pix': 10,
}
Y, X = np.mgrid[-25:26, -25:26]
RADIUS = np.hypot(X, Y)
FIT_MASK = RADIUS <= 10
SKY_MASK = (RADIUS >= 17) & (RADIUS <= 24)


def digest(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def row_hash(row: dict) -> str:
    return digest(json.dumps(row, sort_keys=True, separators=(',', ':')).encode())


def atomic_json(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode='w', dir=path.parent, delete=False) as f:
        tmp = Path(f.name)
        try:
            json.dump(value, f, indent=2, ensure_ascii=False, allow_nan=False)
            f.write('\n')
            f.flush()
            os.fsync(f.fileno())
        except BaseException:
            tmp.unlink(missing_ok=True)
            raise
    tmp.replace(path)


def signature(path: Path) -> dict:
    try:
        s = path.stat()
        return {'path': str(path), 'size': s.st_size, 'mtime_ns': s.st_mtime_ns}
    except FileNotFoundError:
        return {'path': str(path), 'missing': True}


def paths_for(night: str, processed_root: Path) -> dict[str, Path]:
    if len(night) != 8 or not night.isdigit():
        raise ValueError('Night must be YYYYMMDD')
    l4 = processed_root / night / 'L4'
    return {
        'raw': l4 / f'{night}_unknown_links.json',
        'report': l4 / f'{night}_unknown_morphology.json',
        'passed': l4 / f'{night}_unknown_links_screened.json',
        'quarantine': l4 / f'{night}_unknown_links_quarantine.json',
    }


def split_values(value) -> list[str]:
    return [s.strip() for s in str(value or '').split(';') if s.strip()]


def detection_keys(row: dict) -> list[tuple[str, int]]:
    names = split_values(row.get('image_names'))
    ids = [int(v) for v in split_values(row.get('objids'))]
    times = [float(v) for v in split_values(row.get('mjds'))]
    n = len(names)
    if n < 3 or n != len(ids) or n != len(times) or int(row.get('n_obs', -1)) != n:
        raise ValueError('Inconsistent or fewer than three observations')
    if len(set(names)) != n or len(set(times)) != n or not np.isfinite(times).all():
        raise ValueError('Duplicate exposure/epoch or nonfinite epoch')
    if any(Path(name).name != name or not name.endswith('.fits.gz') for name in names):
        raise ValueError('Invalid image name')
    return list(zip(names, ids))


def majority(values: list[bool | None]) -> str:
    if len(values) < POLICY['min_observations']:
        return 'review'
    required = max(2, (2 * len(values) + 2) // 3)
    bad = sum(v is True for v in values)
    missing = sum(v is None for v in values)
    if bad >= required:
        return 'abnormal'
    return 'review' if bad + missing >= required else 'normal'


def decide(catalog_values: list[bool | None], psf_values: list[bool | None]) -> str:
    state = majority(catalog_values)
    if state != 'abnormal':
        return 'pass' if state == 'normal' else 'review'
    if len(psf_values) != len(catalog_values):
        return 'review'
    state = majority(psf_values)
    return {'abnormal': 'reject', 'normal': 'pass', 'review': 'review'}[state]


def numeric(table: Table, name: str) -> np.ndarray:
    return np.asarray(np.ma.filled(np.ma.asarray(table[name], dtype=float), np.nan))


def catalog_frame(task) -> dict:
    root, night, image_name, objids = task
    path = Path(root) / night / 'L2' / image_name.replace('.fits.gz', '_cat.fits.gz')
    output = {str(oid): {'image_name': image_name, 'objID': oid,
                        'catalog_status': 'unavailable', 'catalog_bad': None} for oid in objids}
    before = signature(path)
    try:
        table = Table.read(path)
        ids = numeric(table, 'objID')
        xy = np.column_stack([numeric(table, 'X_Win'), numeric(table, 'Y_Win')])
        fwhm, ann = numeric(table, 'FWHM'), numeric(table, 'Class_ANN')
        good = ((numeric(table, 'Flag') == 0) & (ann > .9)
                & (numeric(table, 'Mag_PSF') > 14) & (numeric(table, 'Mag_PSF') < 19)
                & (numeric(table, 'MagErr_PSF') < .05) & (fwhm > 0)
                & (numeric(table, 'E') < .25) & np.isfinite(xy).all(axis=1))
        ref_xy, ref_fw, ref_ids = xy[good], fwhm[good], ids[good]
        tree = cKDTree(ref_xy)
        for oid in objids:
            result = output[str(oid)]
            where = np.flatnonzero(ids == oid)
            if len(where) != 1:
                result['reason'] = 'missing_or_duplicate_objID'
                continue
            i = int(where[0])
            if not np.isfinite([*xy[i], fwhm[i], ann[i]]).all() or fwhm[i] <= 0:
                result['reason'] = 'invalid_shape_measurement'
                continue
            neighbors = []
            for radius in POLICY['local_radii_pix']:
                neighbors = [j for j in tree.query_ball_point(xy[i], radius) if ref_ids[j] != oid]
                if len(neighbors) >= POLICY['local_reference_min']:
                    break
            result['local_reference_count'] = len(neighbors)
            if len(neighbors) < POLICY['local_reference_min']:
                result['reason'] = 'insufficient_local_reference_stars'
                continue
            reference = float(np.median(ref_fw[neighbors]))
            ratio = float(fwhm[i] / reference)
            result.update(catalog_status='ok', fwhm=float(fwhm[i]),
                          local_fwhm=reference, fwhm_ratio=ratio, class_ann=float(ann[i]),
                          catalog_bad=bool(ratio > POLICY['fwhm_ratio_gt'] and ann[i] < POLICY['class_ann_lt']))
    except Exception as exc:
        for result in output.values():
            result['reason'] = f'{type(exc).__name__}: {exc}'
    if signature(path) != before:
        raise RuntimeError(f'Input changed while reading {path}')
    return {'image_name': image_name, 'observations': output, 'inputs': [before]}


def subtract_sky(stamp):
    if stamp.shape != (51, 51) or not np.isfinite(stamp).all():
        raise ValueError('Invalid 51x51 stamp')
    design = np.column_stack([np.ones(stamp.size), X.ravel(), Y.ravel()])
    good = SKY_MASK.ravel().copy()
    for _ in range(4):
        if good.sum() < 30:
            raise ValueError('Insufficient background pixels')
        p = np.linalg.lstsq(design[good], stamp.ravel()[good], rcond=None)[0]
        residual = stamp.ravel() - design @ p
        noise = 1.4826 * np.median(np.abs(residual[good] - np.median(residual[good])))
        if not np.isfinite(noise) or noise <= 0:
            raise ValueError('Invalid background noise')
        good = SKY_MASK.ravel() & (np.abs(residual) < 3 * noise)
    return residual.reshape(stamp.shape), float(noise)


def fit_stamp(stamp, x_offset, y_offset, gain, template):
    data, noise = subtract_sky(stamp)
    if not np.isfinite(gain) or gain <= 0:
        raise ValueError('Invalid image gain')
    error = np.sqrt(noise**2 + np.maximum(data, 0) / gain)
    initial_flux = max(np.sum(data[RADIUS < 10]), 1)

    def model(p, broad):
        flux, dx, dy, bkg = p[:4]
        width = p[4] if broad else 1.
        shape = map_coordinates(template, [(Y-dy)/width+25, (X-dx)/width+25],
                                order=3, mode='constant', cval=0) / width**2
        return flux * shape + bkg

    def fit(broad):
        start = [initial_flux, x_offset, y_offset, 0]
        lower = [0, x_offset-2, y_offset-2, -5*noise]
        upper = [initial_flux*10, x_offset+2, y_offset+2, 5*noise]
        if broad:
            start.append(1.2); lower.append(.7); upper.append(4.)
        result = least_squares(lambda p: ((data-model(p, broad))/error)[FIT_MASK], start,
                               bounds=(lower, upper), x_scale='jac', max_nfev=300)
        chi = float(np.sum(((data-model(result.x, broad))/error)[FIT_MASK]**2))
        return result, chi

    fixed, chi_fixed = fit(False)
    extended, chi_extended = fit(True)
    bound = bool(np.any(np.abs(fixed.x[1:3] - [x_offset, y_offset]) > 1.99)
                 or np.any(np.abs(extended.x[1:3] - [x_offset, y_offset]) > 1.99)
                 or extended.x[4] > 3.99 or extended.x[4] < .701)
    width, delta = float(extended.x[4]), chi_fixed-chi_extended
    return {'fit_success': bool(fixed.success and extended.success), 'at_bound': bound,
            'width_scale': width, 'delta_chi2': delta,
            'reduced_chi2': chi_fixed/(int(FIT_MASK.sum())-4)}


def psf_frame(task) -> dict:
    root, night, image_name, objids = task
    l1 = Path(root) / night / 'L1' / image_name
    l2 = Path(root) / night / 'L2' / image_name.replace('.fits.gz', '_cat.fits.gz')
    inputs = [signature(p) for p in (l1, l2)]
    output = {str(oid): {'psf_status': 'unavailable', 'psf_bad': None} for oid in objids}
    try:
        table = Table.read(l2)
        ids = numeric(table, 'objID')
        xy = np.column_stack([numeric(table, 'X_Win'), numeric(table, 'Y_Win')])
        finite_xy = np.isfinite(xy).all(axis=1)
        nearest = np.zeros(len(table))
        nearest[finite_xy] = cKDTree(xy[finite_xy]).query(xy[finite_xy], k=2)[0][:, 1]
        with fits.open(l1, memmap=False) as hdus:
            data, header = hdus[1].data, hdus[1].header
            height, width = data.shape
            good = ((numeric(table, 'Flag') == 0) & (numeric(table, 'Class_ANN') > .95)
                    & (numeric(table, 'Mag_PSF') > 14) & (numeric(table, 'Mag_PSF') < 17.5)
                    & (numeric(table, 'MagErr_PSF') < .03) & (numeric(table, 'E') < .2)
                    & finite_xy & (nearest > 22) & (xy[:, 0] > 30) & (xy[:, 0] < width-30)
                    & (xy[:, 1] > 30) & (xy[:, 1] < height-30))
            gain = float(header['GAIN'])

            def stamp_at(i):
                xx, yy = xy[i]-1
                cx, cy = int(round(xx)), int(round(yy))
                if cx < 25 or cy < 25 or cx+25 >= width or cy+25 >= height:
                    raise ValueError('Source stamp crosses image edge')
                # Match the native float32 extraction used in the accepted trial.
                stamp = data[cy-25:cy+26, cx-25:cx+26].astype('float32').astype(float)
                return stamp, float(xx-cx), float(yy-cy)

            for oid in objids:
                result = output[str(oid)]
                try:
                    where = np.flatnonzero(ids == oid)
                    if len(where) != 1:
                        raise ValueError('Missing or duplicate objID')
                    target = int(where[0])
                    distances = np.linalg.norm(xy-xy[target], axis=1)
                    ref = np.flatnonzero(good & (ids != oid) & (distances <= POLICY['psf_radius_pix']))
                    ref = ref[np.argsort(distances[ref])][:25]
                    result['reference_count'] = len(ref)
                    if len(ref) < 23:
                        raise ValueError('Insufficient PSF references including held-out validation')
                    training, validation = ref[:20], ref[20:]
                    shapes = []
                    for i in training:
                        stamp, dx, dy = stamp_at(i)
                        clean, _ = subtract_sky(stamp)
                        centered = shift(clean, (-dy, -dx), order=3, mode='constant', cval=0)
                        flux = float(np.sum(centered[RADIUS <= 10]))
                        if not np.isfinite(flux) or flux <= 0:
                            raise ValueError('Invalid reference-star flux')
                        shapes.append(centered/flux)
                    template = np.median(shapes, axis=0)
                    template[RADIUS > 15] = 0
                    template /= np.sum(template[RADIUS <= 10])
                    heldout = [fit_stamp(*stamp_at(i), gain, template) for i in validation]
                    widths = [r['width_scale'] for r in heldout if r['fit_success'] and not r['at_bound']]
                    median = float(np.median(widths)) if widths else None
                    result.update(reference_objids=[int(ids[i]) for i in training],
                                  heldout_count=len(widths), heldout_width_median=median)
                    if len(widths) < 3 or not .85 <= median <= 1.15:
                        raise ValueError('PSF template failed held-out-star validation')
                    fit = fit_stamp(*stamp_at(target), gain, template)
                    result.update(fit)
                    if not fit['fit_success'] or fit['at_bound'] or not np.isfinite(list(fit.values())).all():
                        raise ValueError('PSF fit failed or hit a bound')
                    result.update(psf_status='ok', psf_bad=bool(
                        fit['width_scale'] > POLICY['psf_width_gt'] and fit['delta_chi2'] >= POLICY['delta_chi2_ge']))
                except Exception as exc:
                    result['psf_reason'] = f'{type(exc).__name__}: {exc}'
    except Exception as exc:
        for result in output.values():
            result['psf_reason'] = f'{type(exc).__name__}: {exc}'
    if inputs != [signature(p) for p in (l1, l2)]:
        raise RuntimeError(f'PSF inputs changed while reading {image_name}')
    return {'image_name': image_name, 'observations': output, 'inputs': inputs}


def run_frames(function, tasks, workers):
    if not tasks:
        return []
    output = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(function, task) for task in tasks]
        for i, future in enumerate(as_completed(futures), 1):
            output.append(future.result())
            if i % 25 == 0 or i == len(tasks):
                print(f'[{function.__name__}] {i}/{len(tasks)} frames', flush=True)
    return output


def read_current_report(night: str, processed_root: Path) -> dict:
    paths = paths_for(night, processed_root)
    report = json.loads(paths['report'].read_text())
    if report['policy'] != POLICY or report['engine_sha256'] != digest(Path(__file__).read_bytes()):
        raise ValueError('Morphology policy/engine changed; rerun screening')
    if report['night'] != night or report['raw_sha256'] != digest(paths['raw'].read_bytes()):
        raise ValueError('Morphology report does not match current raw catalog')
    for kind in ('passed', 'quarantine'):
        if report['outputs'][kind]['sha256'] != digest(paths[kind].read_bytes()):
            raise ValueError(f'Morphology {kind} catalog changed')
    if any(s != signature(Path(s['path'])) for s in report['inputs']):
        raise ValueError('Morphology L1/L2 inputs changed; rerun screening')
    return report


def approved_rows(night: str, processed_root: Path, rows: list[dict], strict=False) -> tuple[list[dict], dict]:
    """Validate provenance and intersect with the approved candidate set.

    Strict mode is used again immediately before transport; it cannot silently
    remove rows after the submission payload has already been prepared.
    """
    report = read_current_report(night, processed_root)
    raw = json.loads(paths_for(night, processed_root)['raw'].read_text())
    allowed = {r['row_sha256'] for r in report['candidates'] if r['status'] == 'pass'}
    originals = {row_hash(r) for r in raw}
    hashes = [row_hash(r) for r in rows]
    if len(hashes) != len(set(hashes)) or not set(hashes) <= originals:
        raise ValueError('Candidate is duplicated or does not match screened raw catalog')
    passed = [r for r, h in zip(rows, hashes) if h in allowed]
    if strict and len(passed) != len(rows):
        raise ValueError('Submission contains morphology-rejected or unresolved candidates')
    return passed, report


def screen_night(night: str, processed_root=Path('/processed1'), workers=4, force=False) -> dict:
    processed_root = Path(processed_root)
    paths = paths_for(night, processed_root)
    if not paths['raw'].exists():
        raise FileNotFoundError(paths['raw'])
    with Path(str(paths['report']) + '.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if not force:
            try:
                return read_current_report(night, processed_root)
            except (OSError, ValueError, KeyError, TypeError):
                pass
        raw_bytes = paths['raw'].read_bytes()
        rows = json.loads(raw_bytes)
        if not isinstance(rows, list):
            raise ValueError('Unknown catalog must be a list')
        candidates, frame_ids = [], defaultdict(set)
        identities = set()
        for row in rows:
            identity = str(row['linkage_id'])
            if identity in identities:
                raise ValueError(f'Duplicate linkage_id {identity}')
            identities.add(identity)
            record = {'linkage_id': row['linkage_id'], 'trk_sub': row.get('trk_sub') or row.get('trkSub', ''),
                      'row_sha256': row_hash(row), 'status': 'review', 'observations': []}
            try:
                keys = detection_keys(row)
                record['observations'] = [{'image_name': n, 'objID': oid} for n, oid in keys]
                for name, oid in keys:
                    frame_ids[name].add(oid)
            except (ValueError, TypeError) as exc:
                record['reason'] = str(exc)
            candidates.append(record)
        tasks = [(str(processed_root), night, name, sorted(ids)) for name, ids in sorted(frame_ids.items())]
        measurements, inputs = {}, {}
        for frame in run_frames(catalog_frame, tasks, workers):
            for oid, result in frame['observations'].items():
                measurements[(frame['image_name'], int(oid))] = result
            for s in frame['inputs']:
                inputs[s['path']] = s
        psf_ids = defaultdict(set)
        for record in candidates:
            for observation in record['observations']:
                observation.update(measurements[(observation['image_name'], observation['objID'])])
            state = majority([o['catalog_bad'] for o in record['observations']])
            record['catalog_vote'] = state
            if state == 'abnormal':
                for o in record['observations']:
                    psf_ids[o['image_name']].add(o['objID'])
        tasks = [(str(processed_root), night, name, sorted(ids)) for name, ids in sorted(psf_ids.items())]
        psf = {}
        for frame in run_frames(psf_frame, tasks, workers):
            for oid, result in frame['observations'].items():
                psf[(frame['image_name'], int(oid))] = result
            for s in frame['inputs']:
                if s['path'] in inputs and inputs[s['path']] != s:
                    raise RuntimeError('L2 changed between catalog and PSF stages')
                inputs[s['path']] = s
        for record in candidates:
            for o in record['observations']:
                o.update(psf.get((o['image_name'], o['objID']), {}))
            record['status'] = decide([o['catalog_bad'] for o in record['observations']],
                                      [o.get('psf_bad') for o in record['observations']])
            record['catalog_bad_count'] = sum(o['catalog_bad'] is True for o in record['observations'])
            record['psf_bad_count'] = sum(o.get('psf_bad') is True for o in record['observations'])
            record.setdefault('reason', {'pass': 'joint_rule_not_confirmed',
                                         'reject': 'catalog_and_psf_majorities_abnormal',
                                         'review': 'insufficient_evidence_for_majority'}[record['status']])
        if paths['raw'].read_bytes() != raw_bytes or any(s != signature(Path(s['path'])) for s in inputs.values()):
            raise RuntimeError('Inputs changed during screening; outputs not published')
        passed = [row for row, record in zip(rows, candidates) if record['status'] == 'pass']
        quarantine = [{**row, 'morphology_status': record['status'], 'morphology_reason': record['reason']}
                      for row, record in zip(rows, candidates) if record['status'] != 'pass']
        atomic_json(paths['passed'], passed)
        atomic_json(paths['quarantine'], quarantine)
        counts = Counter(record['status'] for record in candidates)
        report = {'night': night, 'created_utc': datetime.now(timezone.utc).isoformat(),
                  'policy': POLICY, 'engine_sha256': digest(Path(__file__).read_bytes()),
                  'raw_catalog': str(paths['raw']), 'raw_sha256': digest(raw_bytes),
                  'inputs': sorted(inputs.values(), key=lambda s: s['path']),
                  'counts': {k: counts[k] for k in ('pass', 'reject', 'review')},
                  'n_candidates': len(rows), 'n_catalog_frames': len(frame_ids), 'n_psf_frames': len(psf_ids),
                  'outputs': {k: {'path': str(paths[k]), 'sha256': digest(paths[k].read_bytes())}
                              for k in ('passed', 'quarantine')}, 'candidates': candidates}
        # The report is the last-written commit marker for both output catalogs.
        atomic_json(paths['report'], report)
        return report


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('night')
    ap.add_argument('--processed-root', type=Path, default=Path('/processed1'))
    ap.add_argument('--workers', type=int, default=4)
    ap.add_argument('--force', action='store_true')
    args = ap.parse_args()
    if not 1 <= args.workers <= 16:
        ap.error('--workers must be between 1 and 16')
    report = screen_night(args.night, args.processed_root, args.workers, args.force)
    print(json.dumps({k: report[k] for k in ('night', 'n_candidates', 'counts', 'n_psf_frames')}), flush=True)


if __name__ == '__main__':
    main()
