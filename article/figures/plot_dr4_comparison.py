
"""
Produce plots comparing the RAVE DR4 parmeters to the Cannon/RAVE labels.
"""


import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.ticker import MaxNLocator

try:
    combined_table

except NameError:
    from rave_io import (rave_kordopatis_dr4, rave_cannon_dr1)

    from astropy.table import join

    combined_table = join(rave_cannon_dr1, rave_kordopatis_dr4, keys=("Name", ))

else:
    print("Warning: Using pre-loaded data!")

QC = combined_table["OK"] * (combined_table["QK_1"] == 0)
combined_table = combined_table[QC]

N_bins = 50

all_columns = [
    # RAVE label, DR4 label
    ("TEFF", "TeffK_1"),
    ("LOGG", "loggK_1"),
    ("FE_H", "__M_H_K_1"),
]
label_limits = {
    "TEFF": (3000, 8000),
    "LOGG": (0, 5),
    "FE_H": (-2.5, 0.75)
}
latex_labels = (r"$T_{\rm eff}$", r"$\log{g}$", r"$[{\rm Fe/H}]$")


K, factor = (len(all_columns), 3.5)
lbdim = 0.2 * factor
trdim = 0.1 * factor
whspace = 0.05
yspace = factor
xspace = factor * K + factor * (K - 1) * whspace + lbdim * (K - 1)
xdim = lbdim + xspace + trdim
ydim = lbdim + yspace + trdim

fig, axes = plt.subplots(1, 3, figsize=(xdim, ydim))
fig.subplots_adjust(
    left=lbdim/xdim, bottom=lbdim/ydim, right=(xspace + lbdim)/xdim,
    top=(yspace + lbdim)/ydim, wspace=whspace, hspace=whspace)

for i, (ax, columns, latex_label) in enumerate(zip(axes, all_columns, latex_labels)):

    cannon_label, dr4_label = columns

    x = combined_table[cannon_label]._data
    y = combined_table[dr4_label]._data
    
    finite = np.isfinite(x*y)
    limits = label_limits[cannon_label]
    bins = np.linspace(limits[0], limits[1], N_bins + 1)
    ax.set_axis_bgcolor("#CCCCCC")
    ax.hist2d(x[finite], y[finite], bins=(bins, bins), norm=LogNorm(),
        cmap="viridis")

    # Common limits, labels.
    ax.set_xlim(limits)
    ax.set_ylim(limits)
    ax.xaxis.set_major_locator(MaxNLocator(6))
    ax.yaxis.set_major_locator(MaxNLocator(6))
    ax.set_xlabel(" ".join([latex_label, r"$({\rm UNRAVE})$"]))
    ax.set_ylabel(" ".join([latex_label, r"$({\rm RAVE}$ ${\rm DR4})$"]))

    [_.set_rotation(30) for _ in ax.get_xticklabels()]
    [_.set_rotation(30) for _ in ax.get_yticklabels()]

    diff = y - x
    print(latex_labels, np.nanmean(diff), np.nanstd(diff))

fig.tight_layout()

fig.savefig("dr4-comparison.pdf", dpi=300)
fig.savefig("dr4-comparison.png")
