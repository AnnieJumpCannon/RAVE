
import os
import matplotlib.pyplot as plt
from astropy.table import Table, vstack
from matplotlib.colors import LogNorm
from matplotlib.ticker import MaxNLocator
from scipy.spatial import Delaunay


import AnniesLasso as tc


RESULTS_PATH = "/data/gaia-eso/arc/rave/results/"
RESULTS_PATH = "../../"

ms_results = Table.read(os.path.join(RESULTS_PATH, "rave-tgas-v43.fits.gz"))
giant_results = Table.read(os.path.join(RESULTS_PATH, "rave-tgas-v42.fits.gz"))
joint_results = Table.read(os.path.join(RESULTS_PATH, "rave-tgas-v46.fits.gz"))

for t in (ms_results, giant_results, joint_results):
    if "Name" not in t.dtype.names:
        t["Name"] = [each.split("/")[-2] + "_" + each.split("/")[-1].split(".rvsun.")[0] + "_" + each.split(".rvsun.")[1].split("-result")[0].replace("_result.pkl", "") for each in t["FILENAME"]]

    t.sort("Name")

assert np.all(ms_results["Name"] == joint_results["Name"])
assert np.all(giant_results["Name"] == joint_results["Name"])

ms_model = tc.load_model(os.path.join(RESULTS_PATH, "rave-tgas-v37.model"))

# HACK MAGIC BEGINS #
ms_model._labelled_set = Table.read(os.path.join(RESULTS_PATH, "rave-tgas-v16b-labelled-set-cut.fits"))
ms_model._labelled_set["TEFF"] = ms_model._labelled_set["EPIC_TEFF"]
ms_model._labelled_set["LOGG"] = ms_model._labelled_set["EPIC_LOGG"]
ms_model._labelled_set["FE_H"] = ms_model._labelled_set["EPIC_FEH"]

# HACK MAGIC ENDS


ms_results["TEFF"] = ms_results["EPIC_TEFF"]
ms_results["LOGG"] = ms_results["EPIC_LOGG"]
ms_results["FE_H"] = ms_results["EPIC_FEH"]




# Plot log(density) of the three models.
K = 3
factor = 3.5
lbdim = 0.2 * factor
trdim = 0.1 * factor
whspace = 0.05
yspace = factor
xspace = factor * K + factor * (K - 1) * whspace + lbdim * (K - 1)
xdim = lbdim + xspace + trdim
ydim = lbdim + yspace + trdim

fig, axes = plt.subplots(1, K, figsize=(xdim, ydim))
fig.subplots_adjust(
    left=lbdim/xdim, bottom=lbdim/ydim, right=(xspace + lbdim)/xdim,
    top=(yspace + lbdim)/ydim, wspace=whspace, hspace=whspace)


extent = (3000, 8000, 0.5, 5.5)

titles = (
    r"${\rm Simple}$ ${\rm model}$",
    r"${\rm Main{-}sequence}$ ${\rm model}$", 
    r"${\rm Giant}$ ${\rm branch}$ ${\rm model}$",
)
all_results = (joint_results, ms_results, giant_results)
for i, (ax, results, title) in enumerate(zip(axes, all_results, titles)):

    ax.hexbin(results["TEFF"], results["LOGG"], extent=extent,
        norm=LogNorm(), cmap="Blues", gridsize=35, linewidths=0.1,
        rasterized=True, edgecolor="#ffffff")

    ax.set_xlim(extent[:2][::-1])
    ax.set_ylim(extent[2:][::-1])

    ax.xaxis.set_major_locator(MaxNLocator(6))
    ax.yaxis.set_major_locator(MaxNLocator(6))

    ax.set_title(title)

    ax.set_xlabel(r"$T_{\rm eff}$ $[{\rm K}]$")
    if ax.is_first_col():
        ax.set_ylabel(r"$\log{g}$")

    else:
        ax.set_yticklabels([])


fig.tight_layout()

fig.savefig("test-set-density.pdf", dpi=300)
fig.savefig("test-set-density.png")








# OK, calculate how many sigma the stars are away from the joint model.



# Plot log(density) of the three models.
K = 2
factor = 3.5
lbdim = 0.2 * factor
trdim = 0.1 * factor
whspace = 0.05
xspace = factor
yspace = factor * K + factor * (K - 1) * whspace + lbdim * (K - 1)
xdim = lbdim + xspace + trdim
ydim = lbdim + yspace + trdim

fig, axes = plt.subplots(K, 1, figsize=(xdim, ydim))
fig.subplots_adjust(
    left=lbdim/xdim, bottom=lbdim/ydim, right=(xspace + lbdim)/xdim,
    top=(yspace + lbdim)/ydim, wspace=whspace, hspace=whspace)


x_mu, y_mu = (0, 0)
x_sigma, y_sigma = (90, 0.15)
x = ((ms_results["TEFF"] - joint_results["TEFF"]) - x_mu)/x_sigma
y = ((ms_results["LOGG"] - joint_results["LOGG"]) - y_mu)/y_sigma

ticks = (-10, -5, 0, 5, 10)

kwds = dict(extent=(-10, 10, -10, 10), gridsize=25, linewidths=0.1, 
    edgecolor="#000000", cmap="Blues", norm=LogNorm())

axes[0].hexbin(x, y, **kwds)
#axes[0].axhline(0, c="#FFFFFF", linewidth=0.5, linestyle="-")
#axes[0].axvline(0, c="#FFFFFF", linewidth=0.5, linestyle="-")

axes[0].set_xticks(ticks)
axes[0].set_yticks(ticks)

axes[0].set_xlabel(r"$(T_{{\rm eff},ms} - T_{{\rm eff},simple})/\delta_{T_{{\rm eff},ms}}$")
axes[0].set_ylabel(r"$(\log{g}_{ms} - \log{g}_{simple})/\delta_{\log{g},ms}$")


x_mu, y_mu = (0, 0)
x_sigma, y_sigma = (50, 0.15)

x2 = ((giant_results["TEFF"] - joint_results["TEFF"]) - x_mu)/x_sigma
y2 = ((giant_results["LOGG"] - joint_results["LOGG"]) - y_mu)/y_sigma

hexbin = axes[1].hexbin(x2, y2, **kwds)


axes[1].set_xticks(ticks)
axes[1].set_yticks(ticks)

axes[1].set_xlabel(r"$(T_{{\rm eff},giant} - T_{{\rm eff},simple})/\delta_{T_{{\rm eff},giant}}$")
axes[1].set_ylabel(r"$(\log{g}_{giant} - \log{g}_{simple})/\delta_{\log{g},giant}$")

fig.tight_layout()

fig.savefig("joint-model-differences.pdf", dpi=300)
fig.savefig("joint-model-differences.png")
