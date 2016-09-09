

import locale
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable

locale.setlocale(locale.LC_ALL, 'en_US')

try:
    cannon_dr1

except NameError:
    from rave_io import get_cannon_dr1
    cannon_dr1 = get_cannon_dr1()


else:
    print("Using pre-loaded data")




t = cannon_dr1

snrs = (100, 50, 10)
M = len(snrs)


factor = 3.5
lbdim = 0.2 * factor
trdim = 0.1 * factor
whspace = 0.05
yspace = factor * 2 + whspace
xspace = factor * M + factor * (M - 1) * whspace + lbdim * (M - 1)
xdim = lbdim + xspace + trdim
ydim = lbdim + yspace + trdim

fig, axes = plt.subplots(2, M, figsize=(xdim, ydim))
fig.subplots_adjust(
    left=lbdim/xdim, bottom=lbdim/ydim, right=(xspace + lbdim)/xdim,
    top=(yspace + lbdim)/ydim, wspace=whspace, hspace=whspace)

kwds = dict(rasterized=True, gridsize=65, extent=(3500, 7500, 0, 5.5), edgecolor="none", linewidths=0)

for i, snr in enumerate(snrs):
    # On top figure show density.
    # On bottom axes show metallicity.
    
    ok = cannon_dr1["QC"] * (cannon_dr1["SNR"] > snr)

    ax = axes[0, i]

    top_hexbin = ax.hexbin(t["TEFF"][ok], t["LOGG"][ok], 
        cmap="Blues", norm=LogNorm(), **kwds)

    K = locale.format("%d", sum(ok), grouping=True).replace(",", "$,$")

    ax.text(0.05, 0.9, r"$S/N > {:.0f}$".format(snr),
        horizontalalignment="left", verticalalignment="bottom",
        transform=ax.transAxes, fontsize=14)
    ax.text(0.05, 0.82, r"${}$".format(K) + r" ${\rm stars}$",
        horizontalalignment="left", verticalalignment="bottom", 
        transform=ax.transAxes, fontsize=14)

    ax.set_xticks([4000, 5000, 6000, 7000])
    
    ax.set_xlim(ax.get_xlim()[::-1])
    ax.set_ylim(ax.get_ylim()[::-1])

    ax.yaxis.set_major_locator(MaxNLocator(6))

    ax.set_xticklabels([])

    if ax.is_first_col():
        ax.set_ylabel(r"$\log{g}$")

    else:
        ax.set_yticklabels([])
    
    ax = axes[1, i]

    bottom_hexbin = ax.hexbin(t["TEFF"][ok], t["LOGG"][ok], C=t["FE_H"][ok],
        reduce_C_function=np.median, vmin=-2.5, vmax=0.25, cmap="plasma",
        **kwds)


    ax.set_xlabel(r"$T_{\rm eff}$ $[{\rm K}]$")

    if ax.is_first_col():
        ax.set_ylabel(r"$\log{g}$")

    else:
        ax.set_yticklabels([])

    ax.set_xticks([4000, 5000, 6000, 7000])

    ax.set_xlim(ax.get_xlim()[::-1])
    ax.set_ylim(ax.get_ylim()[::-1])

    ax.yaxis.set_major_locator(MaxNLocator(6))


#cax, kw = mpl.colorbar.make_axes(list(axes.flatten()), fraction=0.075, pad=0.025, aspect=10)
divider = make_axes_locatable(axes[0, 2])
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar_top = plt.colorbar(top_hexbin, cax=cax)
cbar_top.set_label(r"${\rm Counts}$")


divider = make_axes_locatable(axes[1, 2])
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar_bottom = plt.colorbar(bottom_hexbin, cax=cax, ticks=[-2.5, -2, -1.5, -1, -0.5, 0])
cbar_bottom.set_label(r"$[{\rm Fe/H}]$")

fig.tight_layout()

fig.savefig("hrd-test-set.png")
fig.savefig("hrd-test-set.pdf", dpi=300)