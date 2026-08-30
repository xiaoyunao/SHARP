"""Shared style and I/O helpers for the SHARP PASP analysis bundle."""
from pathlib import Path
import hashlib
import json
import matplotlib as mpl
import matplotlib.pyplot as plt

BLUE = "#2A6FBB"
ORANGE = "#E17C24"
CYAN = "#4DBBD5"
GREEN = "#2CA25F"
PURPLE = "#7B61A8"
GRAY = "#555555"

def apply_style():
    mpl.rcParams.update({
        "font.family": "serif", "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
        "font.size": 9, "axes.labelsize": 10, "axes.titlesize": 10,
        "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
        "axes.linewidth": 0.8, "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": True, "ytick.right": True, "savefig.bbox": "tight",
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })

def finish(fig, stem: Path):
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem.with_suffix(".pdf"))
    fig.savefig(stem.with_suffix(".png"), dpi=300)
    plt.close(fig)

def write_json(path: Path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, sort_keys=True, default=str) + "\n")

def sha256(path: Path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()
