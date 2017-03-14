
"""
Plot giant abundances w.r.t. GES.
"""


import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable


try:
    rave_cannon_dr1, kordopatis_comparisons

except NameError: # Do you know who I am? That's Jeff Vader!

    from rave_io import get_cannon_dr1, get_kordopatis_comparisons

    rave_cannon_dr1 = get_cannon_dr1()

    kordopatis_comparisons = get_kordopatis_comparisons()
    from astropy.table import join

    data_table = join(rave_cannon_dr1, kordopatis_comparisons, keys=("Name", ))

else:
    print("Warning: Using pre-loaded data.")


ok = data_table["QC"]# * (data_table["R"] > 10)

latex_labels = {
    "TEFF_2": r"$T_{\rm eff}$ $[{\rm K}]$ $({\rm Literature})$",
    "TEFF_1": r"$T_{\rm eff}$ $[{\rm K}]$ $({\rm \it{RAVE}}{\rm -on})$",
    "LOGG_1": r"$\log{g}$ $({\rm \it{RAVE}}{\rm -on})$",
    "LOGG_2": r"$\log{g}$ $({\rm Literature})$",
    "FE_H": r"$[{\rm Fe/H}]$ $({\rm \it{RAVE}}{\rm -on})$",
    "FEH": r"$[{\rm Fe/H}]$ $({\rm Literature})$"
}


cannon_labels = ("TEFF_1", "LOGG_1", "FE_H")
literature_labels = ("TEFF_2", "LOGG_2", "FEH")

limits = {
    "TEFF_1": [3500, 7500],
    "LOGG_1": [0, 5.5],
    "FE_H": [-3.5, 0.75]
}

kwds = dict(cmap="plasma", vmin=np.nanmin(data_table["snr"]), vmax=np.nanmax(data_table["snr"]))


K = len(cannon_labels)
factor = 3.5
lbdim = 0.25 * factor
tdim = 0.1 * factor
rdim = 0.2 * factor
whspace = 0.05
yspace = factor
xspace = factor * K + factor * (K - 1) * whspace + lbdim * (K - 1)
xdim = lbdim + xspace + rdim
ydim = lbdim + yspace + tdim

fig, axes = plt.subplots(1, K, figsize=(xdim, ydim))
fig.subplots_adjust(
    left=lbdim/xdim, bottom=lbdim/ydim, right=(xspace + lbdim)/xdim,
    top=(yspace + lbdim)/ydim, wspace=whspace, hspace=whspace)


for i, (ax, cannon_label, literature_label) \
in enumerate(zip(axes, cannon_labels, literature_labels)):

    x = data_table[literature_label]
    y = data_table[cannon_label]
    c = data_table["snr"]

    #xerr = data_table["e_{}".format(literature_label)]
    yerr = data_table["E_{}".format(cannon_label).strip("_1")]

    ax.errorbar(x[ok], y[ok], yerr=yerr[ok], fmt=None, ecolor="#666666",
        zorder=-1)
    scat = ax.scatter(x[ok], y[ok], c=c[ok], s=50, **kwds)


_ = ax.scatter([-999], [-9999], c=[0], **kwds)



for ax, cannon_label, literature_label in zip(axes, cannon_labels, literature_labels):

    lims = limits[cannon_label]
    ax.plot(lims, lims, c="#666666", zorder=-1, linestyle=":")

    ax.set_xlim(lims)
    ax.set_ylim(lims)

    ax.xaxis.set_major_locator(MaxNLocator(6))
    ax.yaxis.set_major_locator(MaxNLocator(6))

    ax.set_xlabel(latex_labels[literature_label])
    ax.set_ylabel(latex_labels[cannon_label])


axes[0].set_xticks([4000, 5000, 6000, 7000])
axes[0].set_yticks([4000, 5000, 6000, 7000])

axes[-1].set_xticks([-3.5, -2.5, -1.5, -0.5, 0.5])
axes[-1].set_yticks([-3.5, -2.5, -1.5, -0.5, 0.5])

fig.tight_layout()


cbar = plt.colorbar(_, 
    cax=fig.add_axes([0.93, fig.subplotpars.bottom, 0.02, fig.subplotpars.top - fig.subplotpars.bottom]))
cbar.set_label(r"${\rm S/N}$ ${\rm RAVE}$ $[{\rm pixel}^{-1}]$")


fig.subplots_adjust(right=0.90)

fig.savefig("kordopatis-calibration.pdf", dpi=300)
fig.savefig("kordopatis-calibration.png")

for ref in set(data_table["REF"]):
    for cannon_label, literature_label in zip(cannon_labels, literature_labels):
        match = data_table["REF"] == ref
        x = data_table[cannon_label][match]
        y = data_table[literature_label][match]
        diff = y - x
        print(ref, np.isfinite(diff).sum(), cannon_label, np.nanmean(diff), np.nanstd(diff))
