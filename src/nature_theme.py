from __future__ import annotations
import matplotlib as mpl

# Restrained Nature-family inspired palette.
NAVY = '#17365D'
TEAL = '#0F6D7A'
TEAL_LIGHT = '#5FA9B7'
ORANGE = '#D45A1A'
ORANGE_LIGHT = '#EBA06B'
SLATE = '#4B6070'
GRAY = '#AEB8BF'
LIGHT_GRAY = '#E8EDF0'
PALE_BLUE = '#EAF4F6'
PALE_ORANGE = '#FBF0E8'
WHITE = '#FFFFFF'
BLACK = '#17202A'

MODEL_COLORS = {
    'MTI': ORANGE,
    'LVEF': NAVY,
    'vEprime': TEAL_LIGHT,
    'MTI+LVEF': TEAL,
    'Clinical': GRAY,
    'Clinical+MTI': SLATE,
}


def apply_nature_style():
    mpl.rcParams.update({
        'font.family': 'Liberation Sans',
        'font.size': 9.0,
        'axes.titlesize': 10.5,
        'axes.titleweight': 'bold',
        'axes.labelsize': 9.2,
        'axes.labelcolor': BLACK,
        'axes.edgecolor': BLACK,
        'axes.linewidth': 0.8,
        'xtick.labelsize': 8.2,
        'ytick.labelsize': 8.2,
        'xtick.color': BLACK,
        'ytick.color': BLACK,
        'xtick.major.width': 0.7,
        'ytick.major.width': 0.7,
        'legend.fontsize': 8.0,
        'legend.frameon': False,
        'figure.facecolor': WHITE,
        'axes.facecolor': WHITE,
        'savefig.facecolor': WHITE,
        'savefig.bbox': 'tight',
        'svg.fonttype': 'none',
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
        'lines.linewidth': 1.8,
    })
