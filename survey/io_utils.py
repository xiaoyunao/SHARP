from __future__ import annotations

import json
from pathlib import Path

import astropy.units as u
from astropy.coordinates import SkyCoord


def write_plan_json(plan: list[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(plan, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def format_target_name(field_id: str) -> str:
    if field_id.startswith("MP") and field_id[2:].isdigit():
        return f"MP_{field_id[2:]}"
    if field_id.isdigit():
        return f"MP_{int(field_id):04d}"
    return field_id


def write_plan_txt(plan: list[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = []
    for row in plan:
        coord = SkyCoord(ra=float(row["ra"]) * u.deg, dec=float(row["dec"]) * u.deg)
        lines.append(
            "\t".join(
                [
                    format_target_name(str(row["field_id"])),
                    coord.ra.to_string(unit=u.hour, sep=":", precision=3, pad=True),
                    coord.dec.to_string(unit=u.deg, sep=":", precision=2, alwayssign=True, pad=True),
                    "0",
                    "W",
                    "1",
                    str(int(round(float(row["exptime_s"])))),
                    "0",
                ]
            )
        )
    path.write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")
