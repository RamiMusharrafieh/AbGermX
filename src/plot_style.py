"""Shared figure style for every plot this toolkit produces.

Design rules applied here (see the notes on each helper for why):

* One validated categorical palette, assigned in a fixed order, never a
  generated or rainbow hue. The eight hues clear the OKLab colour-vision
  separation gates on adjacent pairs (worst protan dE 9.1, worst normal-vision
  dE 19.6 against a #fcfcfb surface).
* Species scatters carry far more than eight categories, so identity is encoded
  by hue AND marker shape together: eight hues x three shapes gives 24 unique
  combinations where colour alone could not. Hue varies fastest, so any two
  species adjacent in the legend always differ in colour, not just in shape.
* Magnitude (the identity heatmap) uses a single-hue light-to-dark ramp rather
  than a multi-hue one, so "darker" reads unambiguously as "higher".
* Recessive chrome: hairline grid one shade off the surface, no top/right
  spines, thin marks with a thin surface-coloured ring so overlapping points
  stay separable.
"""
from matplotlib.colors import LinearSegmentedColormap

# --- ink and surfaces ---------------------------------------------------------

SURFACE = '#fcfcfb'
TEXT_PRIMARY = '#0b0b0b'
TEXT_SECONDARY = '#52514e'
TEXT_MUTED = '#898781'
GRIDLINE = '#e1e0d9'
AXIS_LINE = '#c3c2b7'

# --- categorical palette (fixed order) ----------------------------------------

CATEGORICAL = [
    '#2a78d6',  # blue
    '#eb6834',  # orange
    '#1baf7a',  # aqua
    '#eda100',  # yellow
    '#e87ba4',  # magenta
    '#008300',  # green
    '#4a3aa7',  # violet
    '#e34948',  # red
]

# The first three slots are the ones that also clear the all-pairs gates, so
# they are what emphasis plots use when only a handful of series are coloured.
EMPHASIS = CATEGORICAL[:3]

MARKERS = ['o', '^', 's']

# Context marks: everything that is not the subject of an emphasis plot.
CONTEXT_COLOR = '#c9c8c2'
QUERY_COLOR = TEXT_PRIMARY

# --- sequential ramp (magnitude) ----------------------------------------------

BLUE_RAMP = [
    '#cde2fb', '#b7d3f6', '#9ec5f4', '#86b6ef', '#6da7ec', '#5598e7',
    '#3987e5', '#2a78d6', '#256abf', '#1c5cab', '#184f95', '#104281', '#0d366b',
]
SEQUENTIAL_CMAP = LinearSegmentedColormap.from_list('abgermx_blue', BLUE_RAMP)

BASE_RCPARAMS = {
    'figure.facecolor': SURFACE,
    'figure.dpi': 110,
    'savefig.facecolor': SURFACE,
    'savefig.dpi': 220,
    'savefig.bbox': 'tight',
    'axes.facecolor': SURFACE,
    'axes.edgecolor': AXIS_LINE,
    'axes.linewidth': 0.8,
    'axes.labelcolor': TEXT_SECONDARY,
    'axes.labelsize': 10.5,
    'axes.labelpad': 8,
    'axes.titlesize': 12,
    'axes.titlecolor': TEXT_PRIMARY,
    'axes.titleweight': 'normal',
    'axes.titlepad': 12,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'axes.grid': True,
    'axes.axisbelow': True,
    'grid.color': GRIDLINE,
    'grid.linewidth': 0.7,
    'grid.linestyle': '-',
    'xtick.color': TEXT_MUTED,
    'ytick.color': TEXT_MUTED,
    'xtick.labelcolor': TEXT_SECONDARY,
    'ytick.labelcolor': TEXT_SECONDARY,
    'xtick.labelsize': 9.5,
    'ytick.labelsize': 9.5,
    'xtick.major.width': 0.8,
    'ytick.major.width': 0.8,
    'xtick.major.size': 3.5,
    'ytick.major.size': 3.5,
    'legend.frameon': False,
    'legend.fontsize': 9,
    'legend.title_fontsize': 9.5,
    'legend.labelcolor': TEXT_SECONDARY,
    # Fallback chain: renders in the platform UI sans where available and
    # degrades to DejaVu Sans (bundled with matplotlib) everywhere else.
    'font.family': ['Helvetica Neue', 'Helvetica', 'Arial', 'DejaVu Sans'],
    'font.size': 10,
    'text.color': TEXT_PRIMARY,
    'figure.constrained_layout.use': False,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
}


def use_style():
    """Apply the toolkit style to the global matplotlib rcParams."""
    import matplotlib.pyplot as plt
    plt.rcParams.update(BASE_RCPARAMS)


def species_styles(species_list):
    """Map each species to a distinct (colour, marker) pair, in a fixed order.

    Assignment is positional over the caller's (sorted) species list, so a
    species keeps its appearance as long as the list does. Hue cycles fastest,
    which keeps legend neighbours colour-separated; the marker shape is what
    makes the pairing unique once the eight hues run out.
    """
    styles = {}
    for i, sp in enumerate(species_list):
        styles[sp] = dict(
            color=CATEGORICAL[i % len(CATEGORICAL)],
            marker=MARKERS[(i // len(CATEGORICAL)) % len(MARKERS)],
        )
    return styles


def scatter_species(ax, x, y, style, label=None, size=22, alpha=0.85, zorder=3):
    """Draw one species' points with a thin surface-coloured ring, so that
    overlapping marks stay countable instead of merging into a blob."""
    return ax.scatter(
        x, y, s=size, alpha=alpha, label=label, zorder=zorder,
        color=style['color'], marker=style['marker'],
        linewidths=0.5, edgecolors=SURFACE,
    )


def species_legend(fig, handles, labels, ncol=6, y=0.0, title=None,
                   marker_size=46, italic=None):
    """One shared legend below the plot area, laid out in columns.

    A legend below the axes (rather than squeezed into a tall right-hand
    column) lets the plotting area keep the full figure width, which matters
    when two panels sit side by side.

    Every swatch is drawn at one size regardless of how big the mark is in the
    plot, so an emphasised marker does not turn into a blob in the legend.
    ``italic`` optionally restricts italics to the entries that are species
    names; pass a collection of labels, or leave it None to italicise all.
    """
    legend = fig.legend(
        handles, labels, loc='upper center', bbox_to_anchor=(0.5, y),
        ncol=ncol, title=title, frameon=False, handletextpad=0.5,
        columnspacing=1.4, borderaxespad=0, labelspacing=0.7,
    )
    for handle in legend.legend_handles:
        if hasattr(handle, 'set_sizes'):
            handle.set_sizes([marker_size])
    for text, label in zip(legend.get_texts(), labels):
        text.set_color(TEXT_SECONDARY)
        if italic is None or label in italic:
            text.set_fontstyle('italic')
    if legend.get_title() is not None:
        legend.get_title().set_color(TEXT_PRIMARY)
    return legend


def figure_title(fig, title, subtitle=None, x=0.0, y=1.0):
    """Left-aligned title with an optional quieter subtitle underneath.

    Left alignment (rather than centred) keeps the title anchored to the same
    edge as the y-axis label, and leaves room for a subtitle carrying the
    sample size and method details that would otherwise bloat the title.
    """
    fig.text(x, y, title, ha='left', va='bottom', fontsize=15,
             fontweight='bold', color=TEXT_PRIMARY)
    if subtitle:
        fig.text(x, y - 0.018, subtitle, ha='left', va='top', fontsize=10.5,
                 color=TEXT_SECONDARY)


def style_axes(ax, grid_axis='both'):
    """Apply the recessive grid/spine treatment to one axes."""
    ax.grid(True, axis=grid_axis, color=GRIDLINE, linewidth=0.7, zorder=0)
    ax.set_axisbelow(True)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    for side in ('left', 'bottom'):
        ax.spines[side].set_color(AXIS_LINE)
        ax.spines[side].set_linewidth(0.8)
    return ax


def query_projection(ax, coords, species, query_coord, highlight, query_name,
                     point_size=16, context_size=11):
    """Draw a query point against a natural/simulated background distribution.

    The question this chart answers is "where does one sequence sit relative to
    the natural spread", so colour goes to the few species being called out and
    everything else recedes to a single neutral cloud. Colouring all two dozen
    species instead produced an unreadable confetti field in which the query
    marker was the only thing anyone could actually locate.

    ``highlight`` is an ordered sequence of species names (at most three, the
    number this palette separates under all-pairs colour-vision checks).
    """
    highlight = list(highlight)[:len(EMPHASIS)]
    rest = ~species.isin(highlight) if hasattr(species, 'isin') else None
    if rest is None:
        import numpy as _np
        rest = ~_np.isin(species, highlight)

    ax.scatter(coords[rest, 0], coords[rest, 1], s=context_size, alpha=0.5,
               color=CONTEXT_COLOR, linewidths=0, zorder=2,
               label=f'other species (n={int(rest.sum())})')
    for i, sp in enumerate(highlight):
        mask = (species == sp).to_numpy() if hasattr(species, 'to_numpy') else (species == sp)
        ax.scatter(coords[mask, 0], coords[mask, 1], s=point_size, alpha=0.85,
                   color=EMPHASIS[i], linewidths=0.4, edgecolors=SURFACE, zorder=3,
                   label=sp)
    ax.scatter(query_coord[:, 0], query_coord[:, 1], s=200, marker='D',
               color=QUERY_COLOR, edgecolor=SURFACE, linewidth=1.6, zorder=6,
               label=f'{query_name} (query)')
    style_axes(ax)
    return ax
