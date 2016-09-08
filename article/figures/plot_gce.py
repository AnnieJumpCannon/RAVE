
"""
Plot [X/Fe] with respect to [Fe/H]
"""


import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import LogNorm, PowerNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable


try:
    rave_cannon_dr1

except NameError: # Do you know who I am? That's Jeff Vader!

    from rave_io import get_cannon_dr1, get_rave_kordopatis_dr4

    rave_cannon_dr1 = get_cannon_dr1()

    #from astropy.table import join
    #rave_cannon_dr1 = join(rave_cannon_dr1, get_rave_kordopatis_dr4(), keys=("Name", ))



x = rave_cannon_dr1["FE_H"] 


#ok  = (rave_cannon_dr1["snr"] > 10) * (rave_cannon_dr1["r_chi_sq"] < 3) \
#    * (rave_cannon_dr1["R"] > 10) * (rave_cannon_dr1["loggK"] < 3)

ok = rave_cannon_dr1["QC"]

default_kwds = dict(gridsize=(31, 9), cmap="Blues", norm=LogNorm())
elements_and_kwds = OrderedDict([
    ("O", dict(extent=(-1.2, 0.5, -0.05, 0.4))),
    ("AL", dict(extent=(-1.2, 0.5, -0.3, 0.5))),
    ("MG", dict(extent=(-1.2, 0.5, -0.1, 0.4))),
    ("CA", dict(extent=(-1.2, 0.5, -0.3, 0.4))),
    ("SI", dict(extent=(-1.2, 0.5, -0.1, 0.6))),
    ("NI", dict(extent=(-1.2, 0.5, -0.3, 0.2)))
])



fig, axes = plt.subplots(len(elements_and_kwds), 1, figsize=(7, 12))

for i, (element, updated_kwds) in enumerate(elements_and_kwds.items()):

    kwds = default_kwds.copy()
    kwds.update(updated_kwds)

    ax = axes[i]

    y = rave_cannon_dr1["{}_H".format(element)] - x 

    #x.hexbin(x[ok], y[ok], **kwds)
    #extent = kwds.pop("extent")
    #ax.hexbin(x[ok], y[ok], bins=[
    #    np.linspace(extent[0], extent[1], 31),
    #    np.linspace(extent[2], extent[3], 13)
    #    ],
    #    **kwds)
    cmap = plt.cm.Blues
    extent = kwds.get("extent")
    xbins = np.linspace(extent[0], extent[1], 50)#np.ptp(extent[:2])/0.04)
    ybins = np.linspace(extent[2], extent[3], 20)#np.ptp(extent[2:])/0.04)

    show = np.isfinite(x * y)
    ax.hist2d(x[ok * show], y[ok * show], bins=(xbins, ybins), norm=LogNorm(), cmap=cmap)
    #ax.hexbin(x[ok], y[ok], **kwds)
    ax.axhline(0, c="k", linewidth=2, linestyle="--")
    ax.axvline(0, c="k", linewidth=2, linestyle="--")

    ax.xaxis.set_major_locator(MaxNLocator(6))
    if ax.is_last_row():
        ax.set_xlabel(r"$[{\rm Fe/H}]$")

    else:
        ax.set_xticklabels([])

    ax.fill_between(extent[:2], [extent[2], extent[2]], [extent[3], extent[3]],
        facecolor=cmap(0), zorder=-1, edgecolor="none")

    ax.set_xlim(extent[:2])
    ax.set_ylim(extent[2:])

    ax.yaxis.set_major_locator(MaxNLocator(4))
    ax.set_ylabel(r"$[{{\rm {0}/Fe}}]$".format(element.title()))

fig.tight_layout()

for ax in axes:
    values = np.arange(-0.2, 0.8, 0.2)
    show_ticks = (ax.get_ylim()[1] >= values - 0.001) * (values >= ax.get_ylim()[0] + 0.001)
    ax.set_yticks(values[show_ticks])


"""
axes[-1]
y = np.mean([
        rave_cannon_dr1["MG_H"],
        rave_cannon_dr1["O_H"],
        rave_cannon_dr1["CA_H"],
    ], axis=0) - x

axes[-1].hexbin(x[ok], y[ok], extent=(-1.5, 0.5, -0.1, 0.3), **default_kwds)
"""
fig.savefig("gce.pdf", dpi=300)
fig.savefig("gce.png")