#!/usr/bin/env python3
"""Build and execute the lightweight, frozen-data paper analysis notebook."""

from __future__ import annotations

import argparse
import os
import tempfile
from pathlib import Path

import nbformat
from nbclient import NotebookClient


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = PROJECT_ROOT / "notebooks" / "paper_analysis_summary.ipynb"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--timeout", type=int, default=120)
    return parser.parse_args()


def build_notebook() -> nbformat.NotebookNode:
    cells = [
        nbformat.v4.new_markdown_cell(
            "# SHARP PASP frozen analysis summary\n\n"
            "This executed notebook is a lightweight audit entry point for the "
            "2025-11-15 through 2026-07-15 paper interval. It reads only frozen "
            "CSV/JSON products in this analysis directory; it does not scan the "
            "production filesystem or call external services."
        ),
        nbformat.v4.new_code_cell(
            "from pathlib import Path\n"
            "import json\n"
            "import pandas as pd\n"
            "from IPython.display import display, Markdown\n\n"
            "ROOT = Path.cwd().resolve()\n"
            "if not (ROOT / 'tables').is_dir():\n"
            "    ROOT = ROOT.parent\n"
            "assert (ROOT / 'tables').is_dir(), ROOT\n"
            "print(f'Frozen analysis root: {ROOT}')"
        ),
        nbformat.v4.new_markdown_cell("## Data accounting"),
        nbformat.v4.new_code_cell(
            "accounting = pd.read_csv(ROOT / 'tables/table2_data_accounting.csv')\n"
            "wanted = [\n"
            "    'strict_raw_frames', 'l2_catalogs', 'known_predictions',\n"
            "    'known_matches_1arcsec', 'unknown_catalog_linkages',\n"
            "    'human_review_marked_real_linkages',\n"
            "    'submission_selected_linkages', 'posthoc_retained_linkages',\n"
            "]\n"
            "display(accounting.loc[accounting.metric.isin(wanted), "
            "['metric','value','unit','grain','definition']].reset_index(drop=True))"
        ),
        nbformat.v4.new_markdown_cell("## Known-object recovery and astrometry"),
        nbformat.v4.new_code_cell(
            "known = pd.read_csv(ROOT / 'tables/table3_known_recovery_astrometry.csv')\n"
            "wanted = [\n"
            "    'source_prediction_rows_including_duplicates', 'predicted_nominal_n',\n"
            "    'matched_1arcsec_n', 'match_fraction_nominal',\n"
            "    'median_radial_residual_arcsec', 'p90_radial_residual_arcsec',\n"
            "    'random_shift_match_fraction',\n"
            "]\n"
            "display(known.loc[known.metric.isin(wanted), "
            "['section','metric','value','unit','denominator','caveat']].reset_index(drop=True))"
        ),
        nbformat.v4.new_markdown_cell("## Unknown-object funnel and identity guardrails"),
        nbformat.v4.new_code_cell(
            "funnel = pd.read_csv(ROOT / 'tables/table4_unknown_funnel_retention.csv')\n"
            "wanted = [\n"
            "    'l2_source_detection', 'tracklet', 'shared_linkage', 'orbit_fit_ok_linkage',\n"
            "    'orbit_is_good_linkage', 'unknown_catalog_linkage',\n"
            "    'human_review_marked_real_linkage', 'submission_selected_linkage',\n"
            "    'posthoc_retained_linkage',\n"
            "]\n"
            "display(funnel.loc[funnel.stage.isin(wanted), "
            "['stage','value','unit','scope','denominator_stage','retention_fraction']].reset_index(drop=True))\n\n"
            "links = pd.read_csv(ROOT / 'tables/table5_retained_links.csv')\n"
            "identity = pd.DataFrame({\n"
            "    'metric': ['retained single-night linkages', 'linear-motion candidate components',\n"
            "               'JPL second-pass C/2025 Y1 candidates', 'authoritative MPC states pending'],\n"
            "    'value': [len(links), links.linear_motion_candidate_group_id.nunique(),\n"
            "              links.jpl_second_pass_numerically_confirmed_candidate.astype(str).str.lower().isin(['true','1']).sum(),\n"
            "              links.mpc_ingest_state.eq('pending').sum()],\n"
            "})\n"
            "display(identity)\n"
            "display(Markdown('**Guardrail:** candidate groups and JPL numerical associations are not independent-object, discovery, designation, or MPC-ingest counts.'))"
        ),
        nbformat.v4.new_markdown_cell("## Scheduler realization and site sensitivity"),
        nbformat.v4.new_code_cell(
            "scheduler = json.loads((ROOT / 'snapshot/scheduler/scheduler_mode_summary.json').read_text())\n"
            "cohort = scheduler['cohort_accounting']\n"
            "display(pd.DataFrame({'metric': list(cohort), 'value': list(cohort.values())}))\n\n"
            "site = json.loads((ROOT / 'snapshot/orbit_site_comparison/orbit_site_sensitivity_summary.json').read_text())\n"
            "site_rows = []\n"
            "for scope in ['all_orbit_links','formal_unknown_catalog','high_confidence_58']:\n"
            "    item = site[scope]\n"
            "    site_rows.append({\n"
            "        'scope': scope,\n"
            "        'rows': item.get('rows', item.get('n_rows')),\n"
            "        'fit_ok_flips': item.get('fit_ok_flips'),\n"
            "        'is_good_flips': item.get('is_good_flips'),\n"
            "    })\n"
            "display(pd.DataFrame(site_rows))\n"
            "display(Markdown('The 960 m sensitivity rerun caused no acceptance flips, but degenerate short-arc orbital elements remain unsuitable for population inference.'))"
        ),
        nbformat.v4.new_markdown_cell("## Reproducibility checks and remaining external inputs"),
        nbformat.v4.new_code_cell(
            "assert len(links) == 58\n"
            "assert links.fit_ok.astype(str).str.lower().eq('true').all()\n"
            "assert links.is_good.astype(str).str.lower().eq('true').all()\n"
            "assert links.mpc_ingest_state.eq('pending').all()\n"
            "assert accounting.loc[accounting.metric.eq('human_review_marked_real_linkages'), 'value'].iat[0] == 68\n"
            "assert accounting.loc[accounting.metric.eq('submission_selected_linkages'), 'value'].iat[0] == 67\n"
            "assert accounting.loc[accounting.metric.eq('posthoc_retained_linkages'), 'value'].iat[0] == 58\n"
            "print('Notebook closure assertions: PASS')\n\n"
            "display(Markdown('Remaining author/upstream inputs: canonical surveyed MPC 327 coordinates and facility metadata; signed science-night quality policy; authoritative 58+9 MPC row records and treatment of six C/2025 Y1 candidates; upstream image-level injection/detection products; weather/equipment and Slurm resource records.'))"
        ),
    ]
    notebook = nbformat.v4.new_notebook(cells=cells)
    notebook.metadata.update(
        {
            "kernelspec": {
                "display_name": "Python 3",
                "language": "python",
                "name": "python3",
            },
            "language_info": {"name": "python", "version": "3"},
        }
    )
    return notebook


def main() -> None:
    args = parse_args()
    output = args.output.resolve()
    staged = output.with_suffix(output.suffix + ".inprogress")
    if output.exists() or staged.exists():
        raise FileExistsError(f"refusing to overwrite notebook output: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    notebook = build_notebook()
    try:
        # The workstation's otherwise usable Python 3.12 kernel currently has
        # an old pkg_resources/debugpy combination that still imports the
        # removed pkgutil.ImpImporter/ImpLoader names.  Keep the compatibility
        # shim process-local instead of changing the user's Python install.
        with tempfile.TemporaryDirectory(prefix="sharp-notebook-kernel-") as temp_dir:
            compat = Path(temp_dir) / "sitecustomize.py"
            compat.write_text(
                "import pkgutil\n"
                "if not hasattr(pkgutil, 'ImpImporter'):\n"
                "    pkgutil.ImpImporter = type('ImpImporter', (), {})\n"
                "if not hasattr(pkgutil, 'ImpLoader'):\n"
                "    pkgutil.ImpLoader = type('ImpLoader', (), {})\n",
                encoding="utf-8",
            )
            kernel_env = dict(os.environ)
            old_pythonpath = kernel_env.get("PYTHONPATH", "")
            kernel_env["PYTHONPATH"] = (
                f"{temp_dir}{os.pathsep}{old_pythonpath}" if old_pythonpath else temp_dir
            )
            executed = NotebookClient(
                notebook,
                timeout=args.timeout,
                kernel_name="python3",
                resources={"metadata": {"path": str(PROJECT_ROOT)}},
            ).execute(env=kernel_env)
        nbformat.write(executed, staged)
        staged.replace(output)
    except Exception:
        staged.unlink(missing_ok=True)
        raise
    print(f"[done] executed notebook: {output}")


if __name__ == "__main__":
    main()
