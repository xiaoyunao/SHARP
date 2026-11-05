"""SHARP adaptation of Pal 5 New step5_stream_model.plot_density_track_summary."""
from pathlib import Path

BLUE, ORANGE, GREEN = '#1f77b4', '#ff8c00', '#2ca02c'
MODES = dict(EV=BLUE, AM=ORANGE, TJ=GREEN)
STRATEGIES = dict(A='0.45', P0='0.20', B=BLUE, C=ORANGE, D=GREEN)
GRID_ALPHA = .15
ERRORBAR = dict(elinewidth=1.6, capsize=4, capthick=1.6, linestyle='none',
                markerfacecolor='white', markeredgewidth=1.6, markersize=8.5)


def apply_style():
    import matplotlib as mpl
    mpl.rcParams.update({
        'font.family': 'DejaVu Sans', 'mathtext.fontset': 'dejavusans',
        'font.size': 13, 'axes.labelsize': 17, 'axes.titlesize': 18,
        'xtick.labelsize': 13, 'ytick.labelsize': 13, 'legend.fontsize': 12,
        'figure.dpi': 140, 'savefig.dpi': 300, 'figure.facecolor': 'white',
        'axes.facecolor': 'white', 'axes.edgecolor': 'black', 'grid.alpha': GRID_ALPHA,
        'pdf.fonttype': 42,
    })


def save_pair(fig, path, description=''):
    path = Path(path)
    fig.savefig(path.with_suffix('.png'), dpi=300, bbox_inches='tight',
                metadata={'Description': description})
    fig.savefig(path.with_suffix('.pdf'), bbox_inches='tight',
                metadata={'Subject': description})
