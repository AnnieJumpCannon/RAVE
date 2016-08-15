

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.ticker import MaxNLocator

try:
    data_table

except NameError:
    from rave_io import (rave_cannon_dr1, ges_idr4)

    from astropy.table import join

    data_table = join(rave_cannon_dr1, ges_idr4, keys=("Name", ))

else:
    print("Warning: Using pre-loaded data!")


# Plot labels against literature.

label_names = ("TEFF_1", "LOGG_2", "FE_H")
ges_label_names = {
    "TEFF_1": "TEFF_2",
    "LOGG_1": "LOGG_2",
    "FE_H": "FEH",
}

label_limits = {
    "TEFF_1": (3000, 8000),
    "LOGG_1": (0, 5),
    "FE_H": (-2, 0.75)
}

latex_labels = {
    "TEFF_1": r"$T_{\rm eff}$",
    "LOGG_1": r"$\log{g}$",
    "FE_H": r"$[{\rm Fe/H}]$"
}


# Exclude ones we think are bad.
ok = data_table["OK"]

data_table = data_table[ok]

K = len(label_names)
factor = 3.5
lbdim = 0.2 * factor
trdim = 0.1 * factor
whspace = 0.05
yspace = factor
xspace = factor * K + factor * (K - 1) * whspace + lbdim * (K - 1)
xdim = lbdim + xspace + trdim
ydim = lbdim + yspace + trdim

fig, axes = plt.subplots(1, K, figsize=(xdim, ydim))
fig.subplots_adjust(
    left=lbdim/xdim, bottom=lbdim/ydim, right=(xspace + lbdim)/xdim,
    top=(yspace + lbdim)/ydim, wspace=whspace, hspace=whspace)


for i, (ax, label_name) in enumerate(zip(axes, label_names)):

    x = data_table[label_name]
    #xerr
    y = data_table[ges_label_names[label_name]]
    # yerr

    ax.scatter(x, y, c=data_table["SNRK"], alpha=0.75, s=50, 
        linewidths=0.5, edgecolors="#000000")

    ax.set_xlim(limits[label_name])
    ax.set_ylim(limits[label_name])
    ax.xaxis.set_major_locator(MaxNLocator(6))
    ax.yaxis.set_major_locator(MaxNLocator(6))
    ax.set_xlabel(" ".join([latex_labels[label_name], r"$({\rm UNRAVE})$"]))
    ax.set_ylabel(" ".join([latex_labels[label_name], r"$({\rm GES})$"]))

    [_.set_rotation(30) for _ in ax.get_xticklabels()]
    [_.set_rotation(30) for _ in ax.get_yticklabels()]

fig.tight_layout()

fig.savefig("ges-comparison.pdf", dpi=300)
fig.savefig("ges-comparison.png", dpi=300)
