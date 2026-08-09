#!/usr/bin/env python3
"""Audit the twelve planned paper figures and build a QA contact sheet.

For every expected figure this tool records PNG dimensions/DPI/SHA256 and PDF
page count/page size/renderability/SHA256.  It produces a 3 x 4 PNG contact
sheet plus machine-readable JSON.  Existing outputs are never overwritten.

The default is strict: missing PNG/PDF inputs stop before writing.  Use
``--allow-incomplete`` deliberately to render missing tiles and write JSON with
``status: incomplete`` while the figure set is still being assembled.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import shutil
import subprocess
import tempfile
import textwrap
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import PIL
import pypdf
from PIL import Image, ImageDraw, ImageFont
from pypdf import PdfReader, PdfWriter


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_FIGURES_DIR = PROJECT_ROOT / "figures"
DEFAULT_CONTACT_SHEET = PROJECT_ROOT / "qa" / "figure_contact_sheet.png"
DEFAULT_QA_JSON = PROJECT_ROOT / "qa" / "figure_qa.json"


EXPECTED_FIGURES = [
    (1, "fig01_system_architecture", "System architecture"),
    (2, "fig02_footprint_coverage", "Footprint coverage"),
    (3, "fig03_data_accounting", "Data accounting"),
    (4, "fig04_scheduler_example", "Scheduler example"),
    (5, "fig05_known_method_and_residuals", "Known method and residuals"),
    (6, "fig06_unknown_pipeline_examples", "Unknown pipeline examples"),
    (7, "fig07_nightly_exposure_history", "Nightly exposure history"),
    (8, "fig08_known_results", "Known-object results"),
    (9, "fig09_unknown_funnel", "Unknown-object funnel"),
    (10, "fig10_high_confidence_distributions", "High-confidence distributions"),
    (11, "fig11_twilight_context", "Twilight context"),
    (12, "fig12_operations_timeline", "Operations timeline"),
]


CONTACT_COLUMNS = 3
CONTACT_ROWS = 4
CELL_WIDTH = 1000
CELL_HEIGHT = 680
SHEET_MARGIN = 36
SHEET_HEADER_HEIGHT = 124
TILE_HEADER_HEIGHT = 70
TILE_FOOTER_HEIGHT = 62
CONTACT_DPI = 150


COLORS = {
    "background": "#f4f5f7",
    "tile": "#ffffff",
    "matte": "#e6e9ed",
    "ink": "#252a31",
    "muted": "#666e78",
    "border": "#939aa4",
    "pass": "#4c78a8",
    "incomplete": "#8b8f94",
    "failed": "#e15759",
    "missing_fill": "#e3e4e6",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def font_path(*, bold: bool = False) -> str:
    """Return a bundled deterministic font without relying on UI font state."""

    import matplotlib

    name = "DejaVuSans-Bold.ttf" if bold else "DejaVuSans.ttf"
    path = Path(matplotlib.get_data_path()) / "fonts" / "ttf" / name
    if not path.is_file():
        raise FileNotFoundError(f"Bundled QA font not found: {path}")
    return str(path)


def load_font(size: int, *, bold: bool = False) -> ImageFont.FreeTypeFont:
    return ImageFont.truetype(font_path(bold=bold), size=size)


def refuse_existing(output_image: Path, output_json: Path) -> None:
    existing = [str(path) for path in [output_image, output_json] if path.exists()]
    if existing:
        raise FileExistsError("Refusing to overwrite existing QA output(s): " + ", ".join(existing))


def png_audit(
    path: Path,
    *,
    minimum_width: int,
    minimum_height: int,
    target_dpi: float,
    dpi_tolerance: float,
) -> dict[str, Any]:
    result: dict[str, Any] = {
        "path": str(path.resolve(strict=False)),
        "exists": path.is_file(),
        "size_bytes": None,
        "sha256": None,
        "width_px": None,
        "height_px": None,
        "mode": None,
        "dpi_x": None,
        "dpi_y": None,
        "openable": False,
        "minimum_dimensions_ok": False,
        "target_dpi_ok": False,
        "error": None,
    }
    if not path.is_file():
        return result
    result["size_bytes"] = path.stat().st_size
    result["sha256"] = sha256(path)
    try:
        with Image.open(path) as image:
            image.verify()
        with Image.open(path) as image:
            width, height = image.size
            dpi = image.info.get("dpi")
            result.update(
                {
                    "width_px": int(width),
                    "height_px": int(height),
                    "mode": image.mode,
                    "openable": True,
                    "minimum_dimensions_ok": width >= minimum_width and height >= minimum_height,
                }
            )
            if dpi and len(dpi) >= 2:
                result["dpi_x"] = float(dpi[0])
                result["dpi_y"] = float(dpi[1])
                result["target_dpi_ok"] = (
                    abs(float(dpi[0]) - target_dpi) <= dpi_tolerance
                    and abs(float(dpi[1]) - target_dpi) <= dpi_tolerance
                )
    except Exception as exc:
        result["error"] = f"{type(exc).__name__}: {exc}"
    return result


def pdf_audit(path: Path, render_dir: Path, pdftoppm: str | None) -> dict[str, Any]:
    result: dict[str, Any] = {
        "path": str(path.resolve(strict=False)),
        "exists": path.is_file(),
        "size_bytes": None,
        "sha256": None,
        "openable": False,
        "encrypted": None,
        "page_count": None,
        "single_page_ok": False,
        "first_page_width_pt": None,
        "first_page_height_pt": None,
        "first_page_renderable": False,
        "render_width_px": None,
        "render_height_px": None,
        "error": None,
        "render_error": None,
    }
    if not path.is_file():
        return result
    result["size_bytes"] = path.stat().st_size
    result["sha256"] = sha256(path)
    try:
        reader = PdfReader(path)
        result["encrypted"] = bool(reader.is_encrypted)
        if reader.is_encrypted:
            raise ValueError("encrypted PDF is not accepted for paper-figure QA")
        result["page_count"] = len(reader.pages)
        result["single_page_ok"] = len(reader.pages) == 1
        if reader.pages:
            box = reader.pages[0].mediabox
            result["first_page_width_pt"] = float(box.width)
            result["first_page_height_pt"] = float(box.height)
        result["openable"] = True
    except Exception as exc:
        result["error"] = f"{type(exc).__name__}: {exc}"
        return result

    if pdftoppm is None:
        result["render_error"] = "pdftoppm not found"
        return result
    prefix = render_dir / path.stem
    command = [
        pdftoppm,
        "-f",
        "1",
        "-l",
        "1",
        "-singlefile",
        "-png",
        "-r",
        "48",
        str(path),
        str(prefix),
    ]
    try:
        completed = subprocess.run(
            command,
            check=True,
            capture_output=True,
            text=True,
            timeout=60,
        )
        rendered = prefix.with_suffix(".png")
        if not rendered.is_file():
            raise FileNotFoundError(f"pdftoppm returned no image: {completed.stderr.strip()}")
        with Image.open(rendered) as image:
            image.verify()
        with Image.open(rendered) as image:
            result["render_width_px"], result["render_height_px"] = map(int, image.size)
        result["first_page_renderable"] = True
    except Exception as exc:
        result["render_error"] = f"{type(exc).__name__}: {exc}"
    return result


def figure_status(png: dict[str, Any], pdf: dict[str, Any]) -> tuple[str, list[str]]:
    issues: list[str] = []
    if not png["exists"]:
        issues.append("missing PNG")
    if not pdf["exists"]:
        issues.append("missing PDF")
    if issues:
        return "incomplete", issues
    checks = [
        (png["openable"], "PNG cannot be decoded"),
        (png["minimum_dimensions_ok"], "PNG below configured minimum dimensions"),
        (png["target_dpi_ok"], "PNG DPI metadata is not within target tolerance"),
        (pdf["openable"], "PDF cannot be parsed"),
        (pdf["single_page_ok"], "PDF is not exactly one page"),
        (pdf["first_page_renderable"], "PDF first page cannot be rendered"),
    ]
    issues.extend(message for ok, message in checks if not ok)
    return ("pass" if not issues else "failed"), issues


def audit_figure_set(
    figures_dir: Path,
    *,
    minimum_width: int,
    minimum_height: int,
    target_dpi: float,
    dpi_tolerance: float,
) -> dict[str, Any]:
    if minimum_width <= 0 or minimum_height <= 0:
        raise ValueError("minimum PNG dimensions must be positive")
    if target_dpi <= 0 or dpi_tolerance < 0:
        raise ValueError("target DPI must be positive and tolerance non-negative")
    pdftoppm = shutil.which("pdftoppm")
    rows: list[dict[str, Any]] = []
    with tempfile.TemporaryDirectory(prefix="paper_figure_pdf_render.") as temp:
        render_dir = Path(temp)
        for number, stem, title in EXPECTED_FIGURES:
            png = png_audit(
                figures_dir / f"{stem}.png",
                minimum_width=minimum_width,
                minimum_height=minimum_height,
                target_dpi=target_dpi,
                dpi_tolerance=dpi_tolerance,
            )
            pdf = pdf_audit(figures_dir / f"{stem}.pdf", render_dir, pdftoppm)
            status, issues = figure_status(png, pdf)
            rows.append(
                {
                    "figure_number": number,
                    "stem": stem,
                    "title": title,
                    "status": status,
                    "issues": issues,
                    "png": png,
                    "pdf": pdf,
                }
            )
    missing = [row for row in rows if row["status"] == "incomplete"]
    failed = [row for row in rows if row["status"] == "failed"]
    overall = "incomplete" if missing else ("failed" if failed else "complete")
    return {
        "schema_version": "1.0",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "status": overall,
        "expected_figure_count": len(EXPECTED_FIGURES),
        "present_png_count": sum(bool(row["png"]["exists"]) for row in rows),
        "present_pdf_count": sum(bool(row["pdf"]["exists"]) for row in rows),
        "pass_count": sum(row["status"] == "pass" for row in rows),
        "incomplete_count": len(missing),
        "failed_count": len(failed),
        "missing_png_stems": [row["stem"] for row in rows if not row["png"]["exists"]],
        "missing_pdf_stems": [row["stem"] for row in rows if not row["pdf"]["exists"]],
        "validation_contract": {
            "minimum_png_width_px": minimum_width,
            "minimum_png_height_px": minimum_height,
            "target_png_dpi": target_dpi,
            "png_dpi_tolerance": dpi_tolerance,
            "required_pdf_page_count": 1,
            "pdf_first_page_render_required": True,
        },
        "tool_versions": {
            "Pillow": PIL.__version__,
            "pypdf": pypdf.__version__,
            "pdftoppm": pdftoppm,
        },
        "figures": rows,
    }


def wrap_title(text: str, width: int = 38) -> list[str]:
    return textwrap.wrap(text, width=width, break_long_words=False, break_on_hyphens=False)[:2]


def draw_missing_tile(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], label: str) -> None:
    x0, y0, x1, y1 = box
    draw.rectangle(box, fill=COLORS["missing_fill"], outline=COLORS["border"], width=2)
    draw.line((x0 + 28, y0 + 28, x1 - 28, y1 - 28), fill="#b1b4b8", width=5)
    draw.line((x1 - 28, y0 + 28, x0 + 28, y1 - 28), fill="#b1b4b8", width=5)
    font = load_font(36, bold=True)
    bounds = draw.textbbox((0, 0), label, font=font)
    width = bounds[2] - bounds[0]
    height = bounds[3] - bounds[1]
    draw.rectangle(
        (
            (x0 + x1 - width) // 2 - 18,
            (y0 + y1 - height) // 2 - 12,
            (x0 + x1 + width) // 2 + 18,
            (y0 + y1 + height) // 2 + 12,
        ),
        fill="#ffffff",
        outline=COLORS["border"],
        width=2,
    )
    draw.text(((x0 + x1 - width) // 2, (y0 + y1 - height) // 2), label, font=font, fill=COLORS["ink"])


def render_contact_sheet(audit: dict[str, Any]) -> Image.Image:
    width = SHEET_MARGIN * 2 + CONTACT_COLUMNS * CELL_WIDTH
    height = SHEET_MARGIN * 2 + SHEET_HEADER_HEIGHT + CONTACT_ROWS * CELL_HEIGHT
    sheet = Image.new("RGB", (width, height), COLORS["background"])
    draw = ImageDraw.Draw(sheet)
    title_font = load_font(34, bold=True)
    subtitle_font = load_font(19)
    tile_title_font = load_font(21, bold=True)
    tile_meta_font = load_font(16)
    status_font = load_font(16, bold=True)

    draw.text(
        (SHEET_MARGIN, 24),
        "SHARP PASP figure QA contact sheet",
        font=title_font,
        fill=COLORS["ink"],
    )
    subtitle = (
        f"Status: {audit['status'].upper()} | "
        f"PNG {audit['present_png_count']}/12 | PDF {audit['present_pdf_count']}/12 | "
        f"pass {audit['pass_count']}/12"
    )
    draw.text((SHEET_MARGIN, 72), subtitle, font=subtitle_font, fill=COLORS["muted"])

    for index, row in enumerate(audit["figures"]):
        grid_row, grid_column = divmod(index, CONTACT_COLUMNS)
        x0 = SHEET_MARGIN + grid_column * CELL_WIDTH
        y0 = SHEET_MARGIN + SHEET_HEADER_HEIGHT + grid_row * CELL_HEIGHT
        x1 = x0 + CELL_WIDTH - 18
        y1 = y0 + CELL_HEIGHT - 18
        status_color = COLORS[
            "pass" if row["status"] == "pass" else ("failed" if row["status"] == "failed" else "incomplete")
        ]
        draw.rounded_rectangle(
            (x0, y0, x1, y1),
            radius=12,
            fill=COLORS["tile"],
            outline=status_color,
            width=4,
        )
        title_lines = wrap_title(f"Fig. {row['figure_number']}  {row['title']}")
        for line_index, line in enumerate(title_lines):
            draw.text(
                (x0 + 20, y0 + 13 + line_index * 25),
                line,
                font=tile_title_font,
                fill=COLORS["ink"],
            )
        status_text = row["status"].upper()
        status_bounds = draw.textbbox((0, 0), status_text, font=status_font)
        draw.text(
            (x1 - 20 - (status_bounds[2] - status_bounds[0]), y0 + 16),
            status_text,
            font=status_font,
            fill=status_color,
        )

        image_box = (
            x0 + 18,
            y0 + TILE_HEADER_HEIGHT,
            x1 - 18,
            y1 - TILE_FOOTER_HEIGHT,
        )
        png_path = Path(row["png"]["path"])
        if row["png"]["exists"] and row["png"]["openable"]:
            draw.rectangle(image_box, fill=COLORS["matte"], outline=COLORS["border"], width=1)
            with Image.open(png_path) as source:
                thumbnail = source.convert("RGB")
                target_width = image_box[2] - image_box[0] - 12
                target_height = image_box[3] - image_box[1] - 12
                thumbnail.thumbnail((target_width, target_height), Image.Resampling.LANCZOS)
                paste_x = image_box[0] + (image_box[2] - image_box[0] - thumbnail.width) // 2
                paste_y = image_box[1] + (image_box[3] - image_box[1] - thumbnail.height) // 2
                sheet.paste(thumbnail, (paste_x, paste_y))
        else:
            draw_missing_tile(draw, image_box, "MISSING PNG")

        png = row["png"]
        pdf = row["pdf"]
        if png["exists"] and png["openable"]:
            dpi_text = (
                "NA"
                if png["dpi_x"] is None
                else f"{png['dpi_x']:.1f}x{png['dpi_y']:.1f} dpi"
            )
            png_meta = f"PNG {png['width_px']}x{png['height_px']} px; {dpi_text}"
        else:
            png_meta = "PNG missing or unreadable"
        pdf_meta = (
            f"PDF {pdf['page_count']} page; render {'OK' if pdf['first_page_renderable'] else 'FAIL'}"
            if pdf["exists"] and pdf["openable"]
            else "PDF missing or unreadable"
        )
        hash_meta = (
            f"SHA PNG {png['sha256'][:10]}... | PDF {pdf['sha256'][:10]}..."
            if png["sha256"] and pdf["sha256"]
            else "SHA incomplete"
        )
        draw.text((x0 + 20, y1 - 52), f"{png_meta} | {pdf_meta}", font=tile_meta_font, fill=COLORS["ink"])
        draw.text((x0 + 20, y1 - 28), hash_meta, font=tile_meta_font, fill=COLORS["muted"])
    return sheet


def encode_contact_sheet(sheet: Image.Image) -> bytes:
    buffer = io.BytesIO()
    sheet.save(buffer, format="PNG", dpi=(CONTACT_DPI, CONTACT_DPI), optimize=True)
    return buffer.getvalue()


def write_outputs(
    audit: dict[str, Any],
    sheet: Image.Image,
    output_image: Path,
    output_json: Path,
) -> tuple[Path, Path]:
    refuse_existing(output_image, output_json)
    output_image.parent.mkdir(parents=True, exist_ok=True)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    image_bytes = encode_contact_sheet(sheet)
    audit["contact_sheet"] = {
        "path": str(output_image.resolve(strict=False)),
        "sha256": sha256_bytes(image_bytes),
        "size_bytes": len(image_bytes),
        "width_px": sheet.width,
        "height_px": sheet.height,
        "dpi": CONTACT_DPI,
    }
    with output_image.open("xb") as handle:
        handle.write(image_bytes)
    with output_json.open("x", encoding="utf-8") as handle:
        json.dump(audit, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return output_image, output_json


def create_fixture(figures_dir: Path, count: int) -> None:
    figures_dir.mkdir(parents=True, exist_ok=True)
    for number, stem, title in EXPECTED_FIGURES[:count]:
        image = Image.new("RGB", (420, 260), "white")
        draw = ImageDraw.Draw(image)
        draw.rectangle((12, 12, 408, 248), outline=COLORS["pass"], width=5)
        draw.text((30, 42), f"Fig. {number}", font=load_font(32, bold=True), fill=COLORS["ink"])
        draw.text((30, 95), title, font=load_font(18), fill=COLORS["muted"])
        image.save(figures_dir / f"{stem}.png", dpi=(300, 300))
        writer = PdfWriter()
        writer.add_blank_page(width=420, height=260)
        with (figures_dir / f"{stem}.pdf").open("wb") as handle:
            writer.write(handle)


def run_self_test() -> None:
    with tempfile.TemporaryDirectory(prefix="figure_contact_sheet_smoke.") as temp:
        root = Path(temp)
        complete_dir = root / "complete"
        create_fixture(complete_dir, len(EXPECTED_FIGURES))
        audit = audit_figure_set(
            complete_dir,
            minimum_width=400,
            minimum_height=240,
            target_dpi=300.0,
            dpi_tolerance=1.0,
        )
        if audit["status"] != "complete" or audit["pass_count"] != 12:
            raise AssertionError(f"complete fixture failed QA: {audit['status']}")
        image_path = root / "out" / "contact.png"
        json_path = root / "out" / "qa.json"
        sheet = render_contact_sheet(audit)
        write_outputs(audit, sheet, image_path, json_path)
        with Image.open(image_path) as image:
            if image.size != (sheet.width, sheet.height):
                raise AssertionError("contact-sheet size changed on write")
        saved = json.loads(json_path.read_text(encoding="utf-8"))
        if saved["contact_sheet"]["sha256"] != sha256(image_path):
            raise AssertionError("contact-sheet hash mismatch")
        try:
            write_outputs(audit, sheet, image_path, json_path)
        except FileExistsError:
            pass
        else:
            raise AssertionError("overwrite refusal did not trigger")

        incomplete_dir = root / "incomplete"
        create_fixture(incomplete_dir, 11)
        incomplete = audit_figure_set(
            incomplete_dir,
            minimum_width=400,
            minimum_height=240,
            target_dpi=300.0,
            dpi_tolerance=1.0,
        )
        if incomplete["status"] != "incomplete" or incomplete["incomplete_count"] != 1:
            raise AssertionError("incomplete fixture was not identified")
        incomplete_sheet = render_contact_sheet(incomplete)
        write_outputs(
            incomplete,
            incomplete_sheet,
            root / "incomplete_out" / "contact.png",
            root / "incomplete_out" / "qa.json",
        )
    print("synthetic self-test passed (complete, incomplete, PDF render, hash, overwrite refusal)")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--figures-dir", type=Path, default=DEFAULT_FIGURES_DIR)
    parser.add_argument("--output-image", type=Path, default=DEFAULT_CONTACT_SHEET)
    parser.add_argument("--output-json", type=Path, default=DEFAULT_QA_JSON)
    parser.add_argument("--minimum-png-width", type=int, default=1800)
    parser.add_argument("--minimum-png-height", type=int, default=1000)
    parser.add_argument("--target-png-dpi", type=float, default=300.0)
    parser.add_argument("--dpi-tolerance", type=float, default=1.0)
    parser.add_argument(
        "--allow-incomplete",
        action="store_true",
        help="Write an explicitly INCOMPLETE contact sheet/JSON when files are missing",
    )
    parser.add_argument("--self-test", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.self_test:
        run_self_test()
        return
    figures_dir = args.figures_dir.expanduser().resolve(strict=False)
    output_image = args.output_image.expanduser().resolve(strict=False)
    output_json = args.output_json.expanduser().resolve(strict=False)
    refuse_existing(output_image, output_json)
    audit = audit_figure_set(
        figures_dir,
        minimum_width=args.minimum_png_width,
        minimum_height=args.minimum_png_height,
        target_dpi=args.target_png_dpi,
        dpi_tolerance=args.dpi_tolerance,
    )
    if audit["status"] == "incomplete" and not args.allow_incomplete:
        missing = sorted(set(audit["missing_png_stems"] + audit["missing_pdf_stems"]))
        raise FileNotFoundError(
            "Figure set is incomplete; refusing to write final QA outputs. "
            f"Missing stems: {', '.join(missing)}. Use --allow-incomplete only "
            "for an explicitly provisional contact sheet."
        )
    sheet = render_contact_sheet(audit)
    image_path, json_path = write_outputs(audit, sheet, output_image, output_json)
    print(f"wrote {image_path}")
    print(f"wrote {json_path}")
    print(
        f"status={audit['status']} pass={audit['pass_count']}/12 "
        f"png={audit['present_png_count']}/12 pdf={audit['present_pdf_count']}/12"
    )


if __name__ == "__main__":
    main()
