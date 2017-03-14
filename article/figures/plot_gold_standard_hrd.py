
"""
Compare our results to the gold standard studies of Bensby et al. (2014), 
Reddy et al. (2003, 2006), and Valenti & Fischer (2005).
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable

try:
    bensby_2014, reddy_2003, reddy_2006, valenti_2005, rave_cannon_dr1

except NameError: # Do you know who I am? That's Jeff Vader!

    from rave_io import \
        (get_cannon_dr1, get_literature_bensby, get_literature_reddy_2003,
            get_literature_reddy_2006, get_literature_valenti_2005)

    from astropy.table import join

    rave_cannon_dr1 = get_cannon_dr1()

    """
    #rave_cannon_dr1["Name"] = rave_cannon_dr1["Name_1"]
    if "EPIC_TEFF" in rave_cannon_dr1.dtype.names:

        rave_cannon_dr1["TEFF"] = rave_cannon_dr1["EPIC_TEFF"]
        rave_cannon_dr1["LOGG"] = rave_cannon_dr1["EPIC_LOGG"]
        rave_cannon_dr1["FE_H"] = rave_cannon_dr1["EPIC_FEH"]

    if "SNRK" not in rave_cannon_dr1.dtype.names:
        from rave_io import get_rave_kordopatis_dr4
        rave_cannon_dr1 = join(rave_cannon_dr1, get_rave_kordopatis_dr4(), keys=("Name", ))
    
    if "R_CHI_SQ" not in rave_cannon_dr1.dtype.names:
        rave_cannon_dr1["R_CHI_SQ"] = rave_cannon_dr1["r_chi_sq"]
    
    OK = (rave_cannon_dr1["SNRK"] > 10) * (rave_cannon_dr1["R"] > 25) * (rave_cannon_dr1["R_CHI_SQ"] < 3)
    rave_cannon_dr1 = rave_cannon_dr1[OK].filled()
    """
    
    #ok = (rave_cannon_dr1["snr"] > 10) * (rave_cannon_dr1["r_chi_sq"] < 3)
    ok = rave_cannon_dr1["QC"]
    rave_cannon_dr1 = rave_cannon_dr1[ok].filled()

    bensby_2014 = join(rave_cannon_dr1, get_literature_bensby(), keys=("Name", ))
    reddy_2003 = join(rave_cannon_dr1, get_literature_reddy_2003(), keys=("Name", ))
    reddy_2006 = join(rave_cannon_dr1, get_literature_reddy_2006(), keys=("Name", ))
    valenti_2005 = join(rave_cannon_dr1, get_literature_valenti_2005(), keys=("Name", ))

else:
    print("Using pre-loaded data")


# Create convenience functions 

def get_data(table, cannon_label, translated_label):
    x = table[translated_label]
    #xerr #TODO
    y = table[cannon_label]
    # yerr #TODO
    return (x, y)


def bensby_2014_data(cannon_label):

    label_name = {
        "TEFF": "TeffH",
        "LOGG": "loggH",
        "FE_H": "Fe_H"
    }.get(cannon_label, cannon_label)

    return get_data(bensby_2014, cannon_label, label_name)


def reddy_2003_data(cannon_label):

    label_name = {
        "TEFF": "Teff",
        "LOGG": "logg",
        "FE_H": "__Fe_H_"
    }.get(cannon_label, cannon_label)

    return get_data(reddy_2003, cannon_label, label_name)


def reddy_2006_data(cannon_label):

    label_name = {
        "TEFF": "Teff",
        "LOGG": "logg",
        "FE_H": "__Fe_H_"
    }.get(cannon_label, cannon_label)

    return get_data(reddy_2006, cannon_label, label_name)


def valenti_2005_data(cannon_label):

    label_name = {
        "TEFF": "Teff",
        "LOGG": "log_g_",
        "FE_H": "__Fe_H_"
    }.get(cannon_label, cannon_label)

    return get_data(valenti_2005, cannon_label, label_name)






from collections import OrderedDict
comparisons = OrderedDict([
    [bensby_2014_data, dict(s=50, marker="s", )], #facecolor="#EEEFF7", marker="s", zorder=0)],
    [reddy_2003_data, dict(s=50, marker="o", )], #facecolor="#31353D", marker="o", )],
    [reddy_2006_data, dict(s=75, marker="*", )], # facecolor="#92CDCF", marker="o", )],
    [valenti_2005_data, dict(s=50, marker="^", )], #facecolor="#445878", marker="^", )],
])
comparison_labels = (
    r"${\rm Bensby}$ ${\rm et}$ ${\rm al.}$ $(2014)$",
    r"${\rm Reddy}$ ${\rm et}$ ${\rm al.}$ $(2003)$", 
    r"${\rm Reddy}$ ${\rm et}$ ${\rm al.}$ $(2006)$", 
    r"${\rm Valenti}$ ${\rm +}$ ${\rm Fischer}$ $(2005)$") 

label_names = ("TEFF", "LOGG", "FE_H")
N = len(label_names)

def latexify(label_name, axis):
    
    latex = {
        "TEFF": r"$T_{\rm eff}$",
        "LOGG": r"$\log{g}$",
        "FE_H": r"$[{\rm Fe/H}]$"
    }[label_name]

    suffix = r"$({\rm Literature})$" if axis == "x" else r"$({\rm \it{RAVE}}{\rm -on})$"

    return " ".join([latex, suffix])







# Make another 4-panel plot showing the HRD for:
# Bensby, Reddy, Valenti & Fishcer, .


comparison_labels = (
    r"${\rm Bensby}$ ${\rm et}$ ${\rm al.}$ $(2014)$",
    r"${\rm Reddy}$ ${\rm et}$ ${\rm al.}$ $(2003,2006)$", 
    r"${\rm Valenti}$ ${\rm +}$ ${\rm Fischer}$ $(2005)$",
    #r"${\rm Kordopatis}$ ${\rm et}$ ${\rm al.}$ $(2013{\rm ;}$ ${\rm RAVE}$ ${\rm DR4})$",
    r"${\rm \it{RAVE}}{\rm -on}$"
    )
M = len(comparison_labels)

link = True

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

kwds = dict(s=50, cmap="plasma", vmin=-1.5, vmax=0.5)

x, ux = bensby_2014_data("TEFF")
y, uy = bensby_2014_data("LOGG")
c0, uc = bensby_2014_data("FE_H")
axes[0].scatter(x, y, c=c0, **kwds)

if link:
    axes[0].scatter(ux, uy, c=uc, alpha=0.5, marker="^", zorder=-1, **kwds)
    for xi, uxi, uy, uyi in zip(x, ux, y, uy):
        axes[0].plot([xi, uxi], [uy, uyi], c="#666666", zorder=-2, linewidth=0.5,
            alpha=0.5)


x, ux = reddy_2003_data("TEFF")
y, uy = reddy_2003_data("LOGG")
c0, uc = reddy_2003_data("FE_H")
axes[1].scatter(x, y, c=c0, **kwds)

if link:
    axes[1].scatter(ux, uy, c=uc, alpha=0.5, marker="^", zorder=-1, **kwds)
    for xi, uxi, uy, uyi in zip(x, ux, y, uy):
        axes[1].plot([xi, uxi], [uy, uyi], c="#666666", zorder=-2, linewidth=0.5,
            alpha=0.5)



x, ux = reddy_2006_data("TEFF")
y, uy = reddy_2006_data("LOGG")
c0, uc = reddy_2006_data("FE_H")
axes[1].scatter(x, y, c=c0, **kwds)

if link:
    axes[1].scatter(ux, uy, c=uc, alpha=0.5, marker="^", zorder=-1, **kwds)
    for xi, uxi, uy, uyi in zip(x, ux, y, uy):
        axes[1].plot([xi, uxi], [uy, uyi], c="#666666", zorder=-2, linewidth=0.5,
            alpha=0.5)



x, ux = valenti_2005_data("TEFF")
y, uy = valenti_2005_data("LOGG")
c0, uc = valenti_2005_data("FE_H")
scat = axes[2].scatter(x, y, c=c0, **kwds)

if link:
    axes[2].scatter(ux, uy, c=uc, alpha=0.5, marker="^", zorder=-1, **kwds)
    for xi, uxi, uy, uyi in zip(x, ux, y, uy):
        axes[2].plot([xi, uxi], [uy, uyi], c="#666666", zorder=-2, linewidth=0.5,
            alpha=0.5)


"""
for comp in (bensby_2014_data, reddy_2003_data, reddy_2006_data, valenti_2005_data):

    _, x = comp("TeffK")
    _, y = comp("loggK")
    _, c = comp("c_M_H_K")

    axes[3].scatter(x, y, c=c, **kwds)
"""


for comp in (bensby_2014_data, reddy_2003_data, reddy_2006_data, valenti_2005_data):

    _, x = comp("TEFF")
    _, y = comp("LOGG")
    _, c = comp("FE_H")

    _, xerr = comp("E_TEFF")
    _, yerr = comp("E_LOGG")
    print(comp, np.nanmax(yerr))
    axes[-1].errorbar(x, y, xerr=xerr, yerr=yerr, fmt=None, ecolor="#666666", zorder=-1)
    axes[-1].scatter(x, y, c=c, **kwds)


xlim = np.array([ax.get_xlim() for ax in axes])
ylim = np.array([ax.get_ylim() for ax in axes])
for ax, label in zip(axes, comparison_labels):
    #ax.set_xlim(xlim.max(), xlim.min())
    #ax.set_ylim(ylim.max(), ylim.min())
    
    ax.set_xlim(7250, 4250)
    ax.set_ylim(5, 3)
    ax.xaxis.set_major_locator(MaxNLocator(6))
    ax.yaxis.set_major_locator(MaxNLocator(6))

    if ax.is_first_col():
        ax.set_ylabel(r"$\log{g}$")

    else:
        ax.set_yticklabels([])

    ax.set_xlabel(r"$T_{\rm eff}$ $[{\rm K}]$")

    ax.set_title(label)



cbar = plt.colorbar(scat, 
    cax=fig.add_axes([0.93,fig.subplotpars.bottom,0.02,0.9 - fig.subplotpars.bottom]),
    ticks=[-1.5, -1, -0.5, 0, 0.5])
cbar.set_label(r"$[{\rm Fe/H}]$")

fig.tight_layout()
fig.subplots_adjust(right=0.92)



fig.savefig("gold-standard-hrd.png")
fig.savefig("gold-standard-hrd.pdf", dpi=300)



