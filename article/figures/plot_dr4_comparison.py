
"""
Produce plots comparing the RAVE DR4 parmeters to the Cannon/RAVE labels.
"""

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

try:
    rave_cannon_dr1

except NameError:
    from rave_io import (rave_kordopatis_dr4, rave_cannon_dr1)

else:
    print("Warning: Using pre-loaded data!")


fig, axes = plt.subplots(2, 3)

labels = [
    # RAVE label, DR4 label
    ("TEFF", "TeffK"),
    ("LOGG", "loggK"),
    ("FEH", "__M_H_K"),
]
limits = {
    "TEFF": (3500, 7000),
    "LOGG": (0, 5),
}
latex_labels = (r"$T_{\rm eff}$", r"$\log{g}$", r"$[{\rm Fe/H}]$")

for i, ((cannon_label, dr4_label), latex_label) \
in enumerate(zip(labels, latex_labels)):


    x = rave_cannon_dr1[cannon_label]._data
    #xerr = rave_cannon_dr1["E_{}".format(cannon_label)]
    y = rave_cannon_dr1[dr4_label]._data
    #yerr = rave_kordopatis_dr4["E_{}".format(dr4_label)]

    #axes[i, 0].scatter(x, y, facecolor="#000000", alpha=0.1)
    finite = np.isfinite(x*y)
    axes[i, 0].hist2d(x[finite], y[finite], bins=50, norm=LogNorm())

    # Common limits, labels.
    axes[i, 0].set_xlim(limits[cannon_label])
    axes[i, 0].set_ylim(limits[cannon_label])
    axes[i, 0].set_xlabel("{} (This study)".format(latex_label))
    axes[i, 0].set_ylabel("{} (RAVE DR4)".format(latex_label))

    axes[i, 1].scatter(x, y - x)
    axes[i, 1].set_xlim(limits[cannon_label])


    raise a