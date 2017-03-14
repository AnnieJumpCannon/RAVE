

"""
Plot clusters (properly).
"""

import os
from astropy.table import Table, join
import matplotlib.pyplot as plt

try:
    unrave, dr4

except NameError:
    from rave_io import get_cannon_dr1, get_rave_kordopatis_dr4

    clusters = Table.read("../../RAVEDR4_OC.fits")

    unrave = join(get_cannon_dr1(), clusters, keys=("Name",))
    dr4 = join(get_rave_kordopatis_dr4(), clusters, keys=("Name",))

else:
    print("Warning: Using pre-loaded data")



# QC and labels for DR4
dr4_ok = np.ones(len(dr4), dtype=bool)
dr4_teff, dr4_logg, dr4_feh, dr4_e_teff, dr4_e_logg, dr4_e_feh \
    = ("TeffK_1", "loggK_1", "c_M_H_K_1", "e_TeffK_1", "e_loggK_1", "e__M_H_K_1")


unrave_ok = unrave["QC"]
unrave_teff, unrave_logg, unrave_feh, unrave_e_teff, unrave_e_logg, unrave_e_feh \
    = ("TEFF", "LOGG", "FE_H", "E_TEFF", "E_LOGG", "E_FE_H")

iso_teff, iso_logg = ("Teff", "logG")






kwds = dict(cmap="plasma", vmin=-0.5, vmax=0.5, s=80)

cluster_names = sorted(set(clusters["Cluster"]))
K = len(cluster_names)


factor = 3.5
lbdim = 0.2 * factor
rdim = 0.1 * factor
tdim = -0.20 * factor
whspace = 0.05
xspace = factor * 2 + whspace
yspace = factor * K + factor * (K - 1) * whspace + lbdim * (K - 1)
xdim = lbdim + xspace + rdim
ydim = lbdim + yspace + tdim

fig, axes = plt.subplots(K, 2, figsize=(xdim, ydim))
fig.subplots_adjust(
    left=lbdim/xdim, bottom=lbdim/ydim, right=(xspace + lbdim)/xdim,
    top=(yspace + lbdim)/ydim, wspace=whspace, hspace=whspace)


for i, (row, cluster_name) in enumerate(zip(axes, cluster_names)):

    isochrone, isochrone_path = (None, "{}_padova_iso.dat".format(cluster_name))
    if os.path.exists(isochrone_path):
        isochrone = Table.read(isochrone_path, format="ascii")
        isochrone[iso_teff] = 10**isochrone["logTe"]

    match = (dr4["Cluster"] == cluster_name) * dr4_ok
    assert sum(match) > 0

    # Show LHS RAVE DR4
    ax = row[0]
    ax.text(0.15, 0.82, r"${{\rm {0}}}$".format(cluster_name.strip().replace("NGC", "NGC\,")),
        horizontalalignment="left", transform=ax.transAxes, fontsize=16)
    Nstars = np.isfinite(dr4[dr4_teff][match] * dr4[dr4_logg][match]).sum()
    ax.text(0.15, 0.74, r"${0}$ ${{\rm stars}}$".format(Nstars),
        horizontalalignment="left", transform=ax.transAxes, fontsize=16)
    
    ax.errorbar(dr4[dr4_teff][match], dr4[dr4_logg][match],
        xerr=dr4[dr4_e_teff][match], yerr=dr4[dr4_e_logg][match], 
        fmt=None, ecolor="#666666", zorder=-1)
    ax.scatter(dr4[dr4_teff][match], dr4[dr4_logg][match], 
        c=dr4[dr4_feh][match], **kwds)

    if isochrone is not None:
        ax.plot(isochrone[iso_teff], isochrone[iso_logg], c="#666666", lw=2,
            zorder=-1)

    match = (unrave["Cluster"] == cluster_name) * unrave_ok
    assert sum(match) > 0

    # Show RHS unRAVE
    ax = row[1]
    ax.text(0.15, 0.82, r"${{\rm {0}}}$".format(cluster_name.strip().replace("NGC", "NGC\,")),
        horizontalalignment="left", transform=ax.transAxes, fontsize=16)
    Nstars = np.isfinite(unrave[unrave_teff][match] * unrave[unrave_logg][match]).sum()
    ax.text(0.15, 0.74, r"${0}$ ${{\rm stars}}$".format(Nstars),
        horizontalalignment="left", transform=ax.transAxes, fontsize=16)

    ax.errorbar(unrave[unrave_teff][match], unrave[unrave_logg][match],
        xerr=unrave[unrave_e_teff][match], yerr=unrave[unrave_e_logg][match], 
        fmt=None, ecolor="#666666", zorder=-1)
    scat = ax.scatter(unrave[unrave_teff][match], unrave[unrave_logg][match], 
        c=unrave[unrave_feh][match], **kwds)

    if isochrone is not None:
        ax.plot(isochrone[iso_teff], isochrone[iso_logg], c="k", lw=1.5,
            zorder=-1)


for ax in axes.flatten():
    ax.set_xlim(7500, 3500)
    ax.set_ylim(5.5, -0.5)

    ax.set_xticks([7000, 6000, 5000, 4000])
    ax.set_yticks([5, 4, 3, 2, 1, 0])

    if not ax.is_last_row():
        ax.set_xticklabels([])

    else:
        if ax.is_first_col():
            ax.set_xlabel(r"$T_{\rm eff}$ $[{\rm K}]$")

        else:
            ax.set_xlabel(r"$T_{\rm eff}$ $[{\rm K}]$")


    if ax.is_first_col():
        ax.set_ylabel(r"$\log{g}$")

    else:
        ax.set_yticklabels([])

    if ax.is_first_row():
        if ax.is_first_col():
            ax.set_title(r"${\rm RAVE}$ ${\rm DR4}$")
        else:
            ax.set_title(r"${\rm \it{RAVE}}{\rm -on}$")

fig.subplots_adjust(top=0.9)
cbar = plt.colorbar(scat, 
    cax=fig.add_axes([fig.subplotpars.left, 0.93, fig.subplotpars.right - fig.subplotpars.left, 0.02]),
    orientation='horizontal', ticks=[-2.5, -2, -1.5, -1, -0.5, 0, 0.5])
cbar.ax.xaxis.set_ticks_position('top')
cbar.ax.xaxis.set_label_position('top')

cbar.set_label(r"$[{\rm Fe/H}]$")



fig.savefig("open-clusters.pdf", dpi=300)
fig.savefig("open-clusters.png")


