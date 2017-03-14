
"""
Compare our results to the gold standard studies of Bensby et al. (2014), 
Reddy et al. (2003, 2006), and Valenti & Fischer (2005).
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable

try:
    bensby_2014, reddy_2003, reddy_2006, valenti_2005

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
        "FE_H": "Fe_H",
        "E_TEFF": "e_Teff",
        "E_LOGG": "e_logg",
        "E_FE_H": "e_Fe_H"
    }.get(cannon_label, cannon_label)

    return get_data(bensby_2014, cannon_label, label_name)


def reddy_2003_data(cannon_label):

    label_name = {
        "TEFF": "Teff",
        "LOGG": "logg",
        "FE_H": "__Fe_H_",
    }.get(cannon_label, cannon_label)
    # No errors
    a, b =  get_data(reddy_2003, cannon_label, label_name)

    if cannon_label.startswith("E_"):
        a = np.nan * np.ones(len(a))
    return (a, b)


def reddy_2006_data(cannon_label):

    label_name = {
        "TEFF": "Teff",
        "LOGG": "logg",
        "FE_H": "__Fe_H_"
    }.get(cannon_label, cannon_label)

    # No errors..
    a, b = get_data(reddy_2006, cannon_label, label_name)
    if cannon_label.startswith("E_"):
        a = np.nan * np.ones(len(a))
    return a, b


def valenti_2005_data(cannon_label):

    label_name = {
        "TEFF": "Teff",
        "LOGG": "log_g_",
        "FE_H": "__Fe_H_"
    }.get(cannon_label, cannon_label)
    # No errors...

    a, b =  get_data(valenti_2005, cannon_label, label_name)
    if cannon_label.startswith("E_"):
        a = np.nan * np.ones(len(a))
    return a, b




from collections import OrderedDict
comparisons = OrderedDict([
    [bensby_2014_data, dict(s=75, marker="s", )], #facecolor="#EEEFF7", marker="s", zorder=0)],
    [reddy_2003_data, dict(s=75, marker="o", )], #facecolor="#31353D", marker="o", )],
    [reddy_2006_data, dict(s=100, marker="*", )], # facecolor="#92CDCF", marker="o", )],
    [valenti_2005_data, dict(s=75, marker="^", )], #facecolor="#445878", marker="^", )],
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






#fig, axes = plt.subplots(1, N)
from matplotlib import gridspec

fig = plt.figure(figsize=(14.15, 5.4))
gs = gridspec.GridSpec(2, N, width_ratios=[1] * N, height_ratios=[1, 3])

axes = [plt.subplot(gs[i]) for i in range(2 * N)]

handles = []

ylims = {
    "TEFF": 1200,
    "LOGG": 2,
    "FE_H": 1
}
limits = {
    "TEFF": (3500, 7500),
    "LOGG": (2.5, 5.5),
    "FE_H": (-2.5, 0.5)
}
for ax, label_name in zip(axes, list(label_names) + list(label_names)):

    allx = []
    ally = []
    diffs = []
    for comparison, kwds in comparisons.items():
        x, y = comparison(label_name) 
        xerr, yerr = comparison("E_{}".format(label_name))
        c0, c1 = comparison("snr")
        allx.extend(x)
        ally.extend(y)
        if ax.is_first_row():
            ax.errorbar(x, y - x, yerr=yerr, fmt=None, ecolor="#666666", zorder=-1)
            _ = ax.scatter(x, y - x, c=c0, cmap="plasma", **kwds)
            diffs.extend(y-x)

            if ax == axes[0]:
                handles.append(ax.scatter([-10000], [-10000], facecolor="#cccccc", **kwds))

        else:
            ax.errorbar(x, y, yerr=yerr, fmt=None, ecolor="#666666", zorder=-1)
            scat = ax.scatter(x, y, c=c0, cmap="plasma", **kwds)
    print(label_name, np.nanmedian(diffs), np.nanstd(diffs), np.isfinite(diffs).sum())

    if ax.is_first_row():
        ax.set_xticklabels([])

        ylim = ylims[label_name]
        ax.set_ylim(-ylim, +ylim)
        ax.axhline(0, c="#666666", zorder=-1, linestyle=":")

        ax.yaxis.set_major_locator(MaxNLocator(5))

    else:

        # calculate line coeffs

        x = np.array(allx)
        y = np.array(ally)
        yerr = np.ones(x.size)
 
        A = np.vstack((np.ones_like(x), x)).T
        C = np.diag(yerr * yerr)
        cov = np.linalg.inv(np.dot(A.T, np.linalg.solve(C, A)))
        b, m = np.dot(cov, np.dot(A.T, np.linalg.solve(C, y)))


        A = np.vstack((np.ones_like(y), y)).T
        C = np.diag(yerr*yerr)
        cov = np.linalg.inv(np.dot(A.T, np.linalg.solve(C, A)))
        b2, m2 = np.dot(cov, np.dot(A.T, np.linalg.solve(C, x)))

        print(label_name, m2, b2)
        
        
        
        #ax.plot(limits, m*limits + b, c='r', zorder=0)
        limit = limits[label_name]
        ax.plot(limit, limit, c="#666666", zorder=-1, linestyle=":")
        ax.set_xlim(limit)
        ax.set_ylim(limit)

        axes[axes.index(ax) - N].set_xlim(limit)

        ax.set_xlabel(latexify(label_name, "x"))
        ax.set_ylabel(latexify(label_name, "y"))
        

        ax.xaxis.set_major_locator(MaxNLocator(6))
        ax.yaxis.set_major_locator(MaxNLocator(6))


fig.legend(handles, comparison_labels, loc="upper center", 
    frameon=False, ncol=len(comparison_labels), numpoints=1, scatterpoints=1)

#divider = make_axes_locatable(axes[-1])
#cax = divider.append_axes("right", size="5%", pad=0.05)

fig.tight_layout()


cbar = plt.colorbar(scat, cax=fig.add_axes([0.93,fig.subplotpars.bottom,0.02,0.9 - fig.subplotpars.bottom]))
cbar.set_label(r"${\rm RAVE}$ $S/N$ $[{\rm pixel}^{-1}]$")

fig.subplots_adjust(top=0.90, right=0.91)

fig.savefig("gold-standard-comparison.png")
fig.savefig("gold-standard-comparison.pdf", dpi=300)



