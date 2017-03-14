
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
    rave_cannon_dr1

except NameError: # Do you know who I am? That's Jeff Vader!

    from rave_io import get_cannon_dr1, get_ges_idr4, get_rave_kordopatis_dr4

    rave_cannon_dr1 = get_cannon_dr1()
    assert "Name" in rave_cannon_dr1.dtype.names
    #rave_cannon_dr1["Name"] = [each.split("/")[-2] + "_" + each.split("/")[-1].split(".rvsun.")[0] + "_" + each.split(".rvsun.")[1].split("-")[0] for each in rave_cannon_dr1["FILENAME"]]

    from astropy.table import join
    rave_cannon_dr1 = join(rave_cannon_dr1, get_rave_kordopatis_dr4(), keys=("Name", ))
    rave_cannon_dr1 = join(rave_cannon_dr1, get_ges_idr4(), keys=("Name", ))



asplund_2009_solar_abundances = {
    "O": 8.69, 
    "Al": 6.45, 
    "Mg": 7.6, 
    "Ca": 6.34, 
    "Si": 7.51, 
    "Fe": 7.5, 
    "Ni": 6.22,    
}

ges_species = ("O_1", "AL_1", "MG_1", "CA_1", "SI_1", "NI_1")

ok  = rave_cannon_dr1["QC"]


N, M = (3, 2)
assert N * M == len(ges_species)
factor = 3.5
bdim = 0.2 * factor
ldim = 0.25 * factor

tdim = 0.1 * factor
rdim = 0.05 + tdim
whspace = 0.05
xspace = factor * N + factor * (N - 1) * whspace + ldim * (N - 1)
yspace = factor * M + factor * (M - 1) * whspace + bdim * (M - 1)
xdim = ldim + xspace + rdim
ydim = bdim + yspace + tdim

fig, axes = plt.subplots(M, N, figsize=(xdim, ydim))
fig.subplots_adjust(
    left=ldim/xdim, bottom=bdim/ydim, right=(xspace + ldim)/xdim,
    top=(yspace + bdim)/ydim, wspace=whspace, hspace=whspace)

lims = (-1.75, 0.5)

for i, (ax, species) in enumerate(zip(axes.flatten(), ges_species)):

    element = species.split("_")[0]
    x   = rave_cannon_dr1[species.replace("_", "")] \
        - asplund_2009_solar_abundances[element.title()]
    y = rave_cannon_dr1["{}_H".format(element)]

    xerr = rave_cannon_dr1["E_{}".format(species.replace("_", ""))]
    yerr = rave_cannon_dr1["E_{}_H".format(element)]
    ax.errorbar(x[ok], y[ok], xerr[ok], yerr[ok], 
        fmt=None, ecolor="#666666", zorder=-1)
    scat = ax.scatter(x[ok], y[ok], c=rave_cannon_dr1["snr"][ok], s=50,
        vmin=np.nanmin(rave_cannon_dr1["snr"]),
        vmax=np.nanmax(rave_cannon_dr1["snr"]),
        cmap="plasma")

    Nstars = np.isfinite(x*y).sum()
    ax.text(0.05, 0.9, r"$[{{\rm {0}/H}}]$".format(element.title()),
        transform=ax.transAxes, fontsize=14,
        horizontalalignment="left")
    ax.text(0.05, 0.82, r"${0}$ ${{\rm stars}}$".format(Nstars),
        transform=ax.transAxes, fontsize=14)
    ax.text(0.05, 0.74, r"${{\rm Bias:}}$ ${0:.2f}$ ${{\rm dex}}$".format(
        np.nanmean(y - x)),
        horizontalalignment="left", transform=ax.transAxes, fontsize=14)
    ax.text(0.05, 0.66, r"${{\rm RMS:}}$ ${0:.2f}$ ${{\rm dex}}$".format(
        np.nanstd(y - x)),
        horizontalalignment="left", transform=ax.transAxes, fontsize=14)
    
    ax.plot(lims, lims, c="#666666", linestyle=":", zorder=-1)
    ax.set_xlim(lims)
    ax.set_ylim(lims)

    ax.xaxis.set_major_locator(MaxNLocator(6))
    ax.yaxis.set_major_locator(MaxNLocator(6))

    if not ax.is_first_col():
        ax.set_yticklabels([])
    else:
        ax.set_ylabel(r"$[{\rm X/H}]$ $({\rm \it{RAVE}}{\rm -on})$")

    if not ax.is_last_row():
        ax.set_xticklabels([])
    else:
        ax.set_xlabel(r"$[{\rm X/H}]$ $({\rm GES})$")


cbar = plt.colorbar(scat, 
    cax=fig.add_axes([0.93, fig.subplotpars.bottom, 0.02, fig.subplotpars.top - fig.subplotpars.bottom]))
cbar.set_label(r"${\rm S/N}$ ${\rm RAVE}$ $[{\rm pixel}^{-1}]$")

fig.subplots_adjust(right=0.90)

fig.savefig("ges-abundances.pdf", dpi=300)
fig.savefig("ges-abundances.png")

