"""figstyle.py - shared style and helpers for the final Fig. 7-9 of the JTM manuscript (WP6a). Imported by make_fig7_final.py, make_fig8_final.py and make_fig9_final.py.
Layout is set in millimeters from the top-left corner of a 170 mm wide figure; text is Arial at 6 pt or larger at print size; PDF fonts are embedded (TrueType, type 42); PNG at 300 dpi.
Colors: one validated categorical blue and orange from the reference palette of the dataviz skill (slots 1 and 2), one blue ordinal ramp, neutral grays; no color key is ever drawn
inside a figure (user rule of 4 Nov, applied again on 7 Nov): no key strip, no marker-plus-text key and no text that only names a curve, line, symbol or shade; they are explained in the legend text.
Only numbers, axis texts, panel titles and cell labels (NE, absent) stand in the image. Statistical parameters (P, q, r, E, z, T, rho) and gene symbols are set in italics (it() and bit() for mathtext;
fontstyle='italic' for gene symbols, whose digits mathtext would set upright). Usage: imported, not run."""
import os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, to_rgb

sys.stdout.reconfigure(encoding='utf8')
HERE = os.path.dirname(os.path.abspath(__file__))
V2 = os.path.dirname(HERE)                                   # strengthen_v2
W4 = os.path.join(V2, 'wp4_islets')                          # sources of Fig. 7 and 8 (read only)
W2 = os.path.join(V2, 'wp2_benchmark2')                      # sources of Fig. 9 (read only)
FIGW = 170.0                                                 # mm, JTM full page width

INK, INK2, MUTED, GRID, NULLBAR = '#0b0b0b', '#52514e', '#898781', '#e1e0d9', '#c9c9c9'
BLUE, ORANGE, AQUA, GREYLINE = '#2a78d6', '#eb6834', '#1baf7a', '#6f6e69'
BLUE_TXT = '#1c5cab'
RAMP = ['#cde2fb', '#9ec5f4', '#6da7ec', '#3987e5', '#1c5cab', '#0d366b']      # blue, light to dark (reference sequential hue, steps 100, 200, 300, 400, 550, 700)
ORD4 = ['#86b6ef', '#3987e5', '#1c5cab', '#0d366b']                         # ordinal ramp for r = 0.3, 0.5, 0.7, 0.9 (light end no lighter than step 250)
HEAT = LinearSegmentedColormap.from_list('blue_seq', RAMP)
MINUS = '−'


def setup():
    plt.rcParams.update({'font.family': 'Arial', 'font.size': 7, 'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5,
                         'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5, 'axes.labelsize': 7, 'axes.titlesize': 7, 'axes.unicode_minus': True,
                         'mathtext.fontset': 'custom', 'mathtext.rm': 'Arial', 'mathtext.it': 'Arial:italic', 'mathtext.bf': 'Arial:bold', 'mathtext.bfit': 'Arial:italic:bold', 'mathtext.default': 'regular',
                         'pdf.fonttype': 42, 'ps.fonttype': 42, 'savefig.dpi': 300, 'text.color': INK, 'axes.labelcolor': INK, 'xtick.color': INK, 'ytick.color': INK,
                         'axes.edgecolor': INK, 'lines.solid_capstyle': 'butt'})
    try:                                    # subscripts at 86% of the text size (matplotlib default 70%) so that a subscript in 7 pt text is still 6 pt (journal minimum at print size)
        import matplotlib._mathtext as mt
        mt.SHRINK_FACTOR = 0.86
    except Exception as e:                  # pragma: no cover
        print('could not change the subscript size:', e)


def new_fig(height_mm):
    fig = plt.figure(figsize=(FIGW / 25.4, height_mm / 25.4))
    fig._h_mm = float(height_mm)
    return fig


def ax_mm(fig, left, top, width, height, **kw):
    """axes placed in millimeters: left and top are measured from the top-left corner of the figure"""
    H = fig._h_mm
    return fig.add_axes([left / FIGW, (H - top - height) / H, width / FIGW, height / H], **kw)


def tag(fig, s, x=0.8, y=0.6):
    """bold panel tag at (x, y) mm from the top-left corner"""
    fig.text(x / FIGW, 1 - y / fig._h_mm, s, fontsize=9, fontweight='bold', va='top', ha='left')


def text_mm(fig, x, y, s, **kw):
    fig.text(x / FIGW, 1 - y / fig._h_mm, s, **kw)


def clean(ax):
    for sp in ('top', 'right'): ax.spines[sp].set_visible(False)


def it(s):
    """statistical parameter or symbol in italics inside ordinary text (mathtext): fs.it('P') -> italic P. Not for text with digits that must be italic (mathtext sets digits upright): use fontstyle='italic' there (gene symbols)"""
    return r'$\mathit{%s}$' % s


def bit(s):
    """the same in bold italics, for bold titles"""
    return r'$\mathbfit{%s}$' % s


def sg(v, nd=3):
    """signed number with a true minus sign"""
    s = ('%+.' + str(nd) + 'f') % v
    return s.replace('-', MINUS)


def num(v, nd=3):
    s = ('%.' + str(nd) + 'f') % v
    return s.replace('-', MINUS)


def luminance(c):
    r, g, b = [(x / 12.92) if x <= 0.03928 else ((x + 0.055) / 1.055) ** 2.4 for x in to_rgb(c)]
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def save(fig, stem):
    out = os.environ.get('WP6A_OUT', HERE)                    # check_assets.py redirects the output to a scratch folder to test that the scripts reproduce the saved PNG
    png, pdf = os.path.join(out, stem + '.png'), os.path.join(out, stem + '.pdf')
    fig.savefig(png, dpi=300)
    fig.savefig(pdf, metadata={'Creator': 'matplotlib; WP6a make script', 'Title': stem})
    plt.close(fig)
    print('saved', png, '|', pdf)
