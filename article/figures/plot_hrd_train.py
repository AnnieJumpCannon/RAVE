

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import MaxNLocator

"""
try:
    ms_training_set, giant_training_set, joint_training_set

except NameError:
    from rave_io import get_training_sets

    ms_training_set, giant_training_set, joint_training_set = get_training_sets()    

else:
    print("Using pre-loaded data")

"""
names = (
    r"${\rm Simple}$ ${\rm model}$",
    r"${\rm Main{-}sequence}$ ${\rm model}$", 
    r"${\rm Giant}$ ${\rm model}$", 
)

import AnniesLasso as tc
labelled_sets = (
    tc.load_model("rave-tgas-v46.model").labelled_set, # joint
    tc.load_model("rave-tgas-v43.model").labelled_set, # ms result
    tc.load_model("rave-tgas-v42.model").labelled_set, # giant result
)

labelled_sets[1]["TEFF"] = labelled_sets[1]["EPIC_TEFF"]
labelled_sets[1]["LOGG"] = labelled_sets[1]["EPIC_LOGG"]
labelled_sets[1]["FE_H"] = labelled_sets[1]["EPIC_FEH"]

M = len(labelled_sets)

factor = 3.5
cbar = 0.15 * factor
lbdim = 0.2 * factor
trdim = 0.1 * factor
whspace = 0.05
yspace = factor
xspace = factor * M + factor * (M - 1) * whspace + lbdim * (M - 1) + cbar
xdim = lbdim + xspace + trdim
ydim = lbdim + yspace + trdim

fig, axes = plt.subplots(1, M, figsize=(xdim, ydim))
fig.subplots_adjust(
    left=lbdim/xdim, bottom=lbdim/ydim, right=(xspace + lbdim)/xdim,
    top=(yspace + lbdim)/ydim, wspace=whspace, hspace=whspace)


xlim = (7500, 3500)
ylim = (5.5, 0)

vmin, vmax = (-3, 0.5)

# TODO: Add uncertainties.
# TODO: Show other training sets in the background of each axes as grey points.

for ax, labelled_set, name in zip(axes, labelled_sets, names):
    
    ax.scatter(labelled_set["TEFF"], labelled_set["LOGG"],
        c=labelled_set["FE_H"], cmap="plasma",
        vmin=vmin, vmax=vmax, s=50, edgecolor="k", linewidth=0.1, alpha=0.7)

    ax.text(0.05, 0.9, name,
        horizontalalignment="left", verticalalignment="bottom",
        transform=ax.transAxes, fontsize=14)
    ax.text(0.05, 0.82, r"${}$".format(len(labelled_set)) + r" ${\rm stars}$",
        horizontalalignment="left", verticalalignment="bottom", 
        transform=ax.transAxes, fontsize=14)
    
    ax.set_xlim(ax.get_xlim()[::-1])
    ax.set_ylim(ax.get_ylim()[::-1])

    ax.xaxis.set_major_locator(MaxNLocator(6))
    ax.yaxis.set_major_locator(MaxNLocator(6))

    ax.set_xlabel(r"$T_{\rm eff}$ $[{\rm K}]$")

    if ax.is_first_col():
        ax.set_ylabel(r"$\log{g}$")

    else:
        ax.set_yticklabels([])


scat = axes[0].scatter([0], [0], c=[0], cmap="plasma",
    vmin=vmin, vmax=vmax, alpha=1, s=50, edgecolor="none")

for ax in axes:
    ax.set_xticks([7000, 6000, 5000, 4000])
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)


cax, kw = mpl.colorbar.make_axes(list(axes), fraction=0.075, pad=0.025, aspect=10)
cbar = plt.colorbar(scat, cax=cax, ticks=[-3, -2, -1, 0])

cbar.set_label(r"$[{\rm Fe/H}]$")
#cbar.ax.yaxis.set_ticks([-3, -2, -1, 0])


fig.savefig("hrd-train-set.png")
fig.savefig("hrd-train-set.pdf", dpi=300)
