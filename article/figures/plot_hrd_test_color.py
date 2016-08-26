

import locale
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import MaxNLocator

locale.setlocale(locale.LC_ALL, 'en_US')

try:
    rave_cannon_dr1

except NameError:
    from rave_io import rave_cannon_dr1

else:
    print("Using pre-loaded data")




t = rave_cannon_dr1

#ok = (t["SNRK"] > 50) * (t["r_chi_sq_ms"] < 3) * (t["r_chi_sq_giant"] < 3) #* (t["WEIGHTED_VSINI"] < 1)

xlim = (7500, 3500)
ylim = (5.5, 0)
vmin, vmax = (-3, 0.5)

snrs = (100, 50, 25, 10)
M, N = (len(snrs), 50)


factor = 3.5
lbdim = 0.2 * factor
trdim = 0.1 * factor
whspace = 0.05
yspace = factor
xspace = factor * M + factor * (M - 1) * whspace + lbdim * (M - 1)
xdim = lbdim + xspace + trdim
ydim = lbdim + yspace + trdim

fig, axes = plt.subplots(1, M, figsize=(xdim, ydim))
fig.subplots_adjust(
    left=lbdim/xdim, bottom=lbdim/ydim, right=(xspace + lbdim)/xdim,
    top=(yspace + lbdim)/ydim, wspace=whspace, hspace=whspace)


for ax, snr in zip(axes, snrs):
    
    ok = (t["SNRK"] > snr) * (t["R_CHI_SQ"] < 3) * (t["R"] > 25) #* (t["VSINI"] < 1)

    ax.scatter(t["TEFF"][ok], t["LOGG"][ok], c=t["FE_H"][ok],
        vmin=vmin, vmax=vmax, alpha=0.25, s=50, edgecolor="none", cmap="plasma",
        rasterized=True)

    if ax.is_last_col():
        scat = ax.scatter([0], [0], c=[0], vmin=vmin, vmax=vmax, cmap="plasma")

    K = locale.format("%d", sum(ok), grouping=True).replace(",", "$,$")

    ax.text(0.05, 0.9, r"$S/N > {:.0f}$".format(snr),
        horizontalalignment="left", verticalalignment="bottom",
        transform=ax.transAxes, fontsize=14)
    ax.text(0.05, 0.82, r"${}$".format(K) + r" ${\rm stars}$",
        horizontalalignment="left", verticalalignment="bottom", 
        transform=ax.transAxes, fontsize=14)

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    ax.xaxis.set_major_locator(MaxNLocator(6))
    ax.xaxis.set_ticks([7000, 6000, 5000, 4000])
    ax.yaxis.set_major_locator(MaxNLocator(6))

    ax.set_xlabel(r"$T_{\rm eff}$ $({\rm K})$")

    if ax.is_first_col():
        ax.set_ylabel(r"$\log{g}$")

    else:
        ax.set_yticklabels([])

cax, kw = mpl.colorbar.make_axes(list(axes), fraction=0.075, pad=0.025, aspect=10)
cbar = plt.colorbar(scat, cax=cax, ticks=[-3, -2, -1, 0])

cbar.set_label(r"$[{\rm Fe/H}]$")


fig.savefig("hrd-test-set-color.png")
fig.savefig("hrd-test-set-color.pdf", dpi=300)