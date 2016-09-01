

"""
Join a MS and giant model.... in a bad way.
"""

import numpy as np
from astropy.table import Table
from matplotlib.colors import LogNorm


import AnniesLasso as tc

#try:
#    ms, giant, complete
#except NameError:

# Load the result files.
"""
ms = Table.read("results/rave-tgas-v16b.fits.gz")
giant = Table.read("results/rave-tgas-v23giant.fits.gz")
complete = Table.read("results/rave-tgas-v23all.fits.gz")
"""

"""
ms = Table.read("results/rave-tgas-v27.fits.gz")
giant = Table.read("results/rave-tgas-v23giant.fits.gz")
complete = Table.read("results/rave-tgas-v26.fits.gz")
"""

"""
ms = Table.read("results/rave-tgas-v31.fits.gz")
giant = Table.read("results/rave-tgas-v23giant.fits.gz")
complete = Table.read("results/rave-tgas-v26.fits.gz")
"""

ms = Table.read("results/rave-tgas-v37.fits.gz")
giant = Table.read("results/rave-tgas-v36.fits.gz")
complete = Table.read("results/rave-tgas-v26.fits.gz")


# Common columns
if "EPIC_TEFF" in ms.dtype.names:

    ms["TEFF"] = ms["EPIC_TEFF"]
    ms["LOGG"] = ms["EPIC_LOGG"]
    ms["FE_H"] = ms["EPIC_FEH"]

# Load the models
"""
ms_model = tc.load_model("rave-tgas-v16b.model")
giant_model = tc.load_model("rave-tgas-v23giants.model")
complete_model = tc.load_model("rave-tgas-v23all.model")
"""

"""
ms_model = tc.load_model("rave-tgas-v27.model")
giant_model = tc.load_model("rave-tgas-v23giants.model")
complete_model = tc.load_model("rave-tgas-v26.model")
"""

"""
ms_model = tc.load_model("rave-tgas-v31.model")
giant_model = tc.load_model("rave-tgas-v23giants.model")
complete_model = tc.load_model("rave-tgas-v26.model")
"""

ms_model = tc.load_model("rave-tgas-v37.model")
giant_model = tc.load_model("rave-tgas-v36.model")
complete_model = tc.load_model("rave-tgas-v26.model")



ms_model._labelled_set = Table.read("rave-tgas-v16b-labelled-set-cut.fits")
ms_model._labelled_set["TEFF"] = ms_model._labelled_set["EPIC_TEFF"]
ms_model._labelled_set["LOGG"] = ms_model._labelled_set["EPIC_LOGG"]
ms_model._labelled_set["FE_H"] = ms_model._labelled_set["EPIC_FEH"]

# Add flags for convex hull.
"""
giant["IN_CONVEX_HULL"] = giant_model.in_convex_hull(
    np.array([giant["TEFF"], giant["LOGG"], giant["FE_H"], giant["VSINI"]]).T)

ms["IN_CONVEX_HULL"] = ms_model.in_convex_hull(
    np.array([ms["TEFF"], ms["LOGG"], ms["FE_H"]]).T)#ms["VSINI"]]).T)

complete["IN_CONVEX_HULL"] = complete_model.in_convex_hull(
    np.array([complete["TEFF"], complete["LOGG"], complete["FE_H"], complete["VSINI"]]).T)
"""


complete.sort("Name")
ms.sort("Name")
giant.sort("Name")

"""
failed = [
    "20130404_0950m37_003",
    "20130404_1450m44_003"
]
ok = np.array([each not in failed for each in complete["Name"]])
complete = complete[ok]

ok = np.array([each not in failed for each in giant["Name"]])
giant = giant[ok]
"""

assert np.all(complete["Name"] == ms["Name"])
assert np.all(complete["Name"] == giant["Name"])

#else:
#    print("Using pre-loaded data!")



from scipy.spatial import Delaunay



subg = ms_model.labelled_set["LOGG"] < 4

bad1 = subg * (ms_model.labelled_set["TEFF"] < 4870)
bad2 = subg * ((ms_model.labelled_set["TEFF"] < 4965) * (ms_model.labelled_set["LOGG"] > 3.66))
keep = ~bad1 * ~bad2


ms_hull = Delaunay(ms_model.labels_array[:, :2])

ms_subg_convex_hull = Delaunay(ms_model.labels_array[:, :2][subg * keep])

giant_hull = Delaunay(giant_model.labels_array[:, :2])


in_ms_hull = ms_hull.find_simplex(np.array([ms["TEFF"], ms["LOGG"]]).T) >= 0
in_ms_subg_hull = ms_subg_convex_hull.find_simplex(np.array([ms["TEFF"], ms["LOGG"]]).T) >= 0

in_giant_hull = giant_hull.find_simplex(np.array([giant["TEFF"], giant["LOGG"]]).T) >= 0





ms_dist_to_model = np.abs(ms["LOGG"] - complete["LOGG"])
giant_dist_to_model = np.abs(giant["LOGG"] - complete["LOGG"])


# For stars in the subgiant area, if a star is outside the convex hull then we
# should set the distance to the parent model as being large.
#ms_subgiant_area = (5375 > ms["TEFF"]) * (ms["TEFF"] > 4830) \
#                 * (3.8 > ms["LOGG"])  * (ms["LOGG"] > 3.4)


ms_subgiant_area = (5500 > ms["TEFF"]) * (ms["LOGG"] < 4)

#giant_subgiant_area = (5375 > giant["TEFF"]) * (giant["TEFF"] > 4830) \
#                    * (3.8 > giant["LOGG"])  * (giant["LOGG"] > 3.4)          

giant_subgiant_area = (giant["LOGG"] > 3.5)# * (giant["TEFF"] > 4900)#5375)


"""
# Huber definition
for results in (ms, giant):
    results["is_giant_Huber"] = ((results["TEFF"] > 5000) * (results["LOGG"] < (13.463 - 0.00191*results["TEFF"]))) + (results["LOGG"] < 3.9)
    results["is_dwarf_Huber"] = results["LOGG"] > (1/4.671) * np.arctan((results["TEFF"] - 6300)/-67.172) + 3.876

    che = results["is_giant_Huber"] * results["is_dwarf_Huber"]
    results["is_dwarf_Huber"][che] = True
    results["is_giant_Huber"][che] = False

    results["is_subg_Huber"] = ((~results["is_giant_Huber"]) * (~results["is_dwarf_Huber"]))
"""



#ms_subgiant_area = ms["is_subg_Huber"] + ms_subgiant_area
#giant_subgiant_area = giant["is_subg_Huber"] * giant_subgiant_area


giant_dist_to_model[giant_subgiant_area * ~in_giant_hull] = np.inf
ms_dist_to_model[ms_subgiant_area * ~in_ms_subg_hull] = np.inf

# Both bad?
really_bad = (ms_subgiant_area * ~in_ms_subg_hull) * (giant_subgiant_area * ~in_giant_hull)
ms["TEFF"][really_bad] = np.nan
giant["TEFF"][really_bad] = np.nan


ms["USE"] = (giant_dist_to_model > ms_dist_to_model) * ~ms_subgiant_area
giant["USE"] = (giant_dist_to_model < ms_dist_to_model) * ~giant_subgiant_area

from astropy.table import vstack

joined = vstack([ms[ms["USE"]], giant[giant["USE"]]])

#ms["USE"] = ~ms_subgiant_area
#giant["USE"] = ~giant_subgiant_area

#giant["TEFF"][~np.isfinite(giant_dist_to_model)] = np.nan
#ms["TEFF"][~np.isfinite(ms_dist_to_model)] = np.nan


# Penalize by distance in logg.
beta = 1000.
giant_alpha = 200.
ms_alpha = 1000.0 

#beta = 500.
#giant_alpha =  beta/13.5
#ms_alpha = 12.5 * giant_alpha


#beta, giant_alpha, ms_alpha = 0, 0, 0

# And heavily penalize things in the subgiant area that are outside of the convex hull in Teff/logg

from scipy.spatial import Delaunay

ms_hull = Delaunay(ms_model.labels_array[:, :2])
giant_hull = Delaunay(giant_model.labels_array[:, :2])


in_ms_hull = ms_hull.find_simplex(np.array([ms["TEFF"], ms["LOGG"]]).T) >= 0
in_giant_hull = giant_hull.find_simplex(np.array([giant["TEFF"], giant["LOGG"]]).T) >= 0

print("Not adding 1000 in chisq")
#ms["chi_sq"][ms_subgiant_area * ~in_ms_hull] += 1000.
#giant["chi_sq"][giant_subgiant_area * ~in_giant_hull] += 1000.



# Join and calculate weights.
from astropy.table import join

t = join(ms, giant, keys=("Name", ), table_names=["ms", "giant"])

# Calculate weights.
#t["w_ms"] = np.exp(-0.5 * (t["chi_sq_ms"] + ms_alpha * ((ms["LOGG"] - ms_fiducial)**2 + (ms["LOGG"] - ms_model.vectorizer.fiducials[1])**2)))
#t["w_giant"] = np.exp(-0.5 * (t["chi_sq_giant"] + giant_alpha * ((giant["LOGG"] - giant_fiducial)**2 + (giant["LOGG"] - giant_model.vectorizer.fiducials[1])**2)))

t["w_ms"] = np.exp(
    -0.5 * ( \
            +     beta * (ms["LOGG"] - complete["LOGG"])**2 \
            + ms_alpha * (ms["LOGG"] - ms_model.vectorizer.fiducials[1])**2))

t["w_giant"] = np.exp(
    -0.5 * (  \
            +        beta * (giant["LOGG"] - complete["LOGG"])**2 \
            + giant_alpha * (giant["LOGG"] - giant_model.vectorizer.fiducials[1])**2))


print(np.sum(t["w_ms"] == 0), np.sum(t["w_giant"] == 0))


for label_name in ("TEFF", "LOGG", "FE_H"):#:"VSINI"):

    t["WEIGHTED_{}".format(label_name)] \
        = (t["w_ms"] * t["{}_ms".format(label_name)] + t["w_giant"] * t["{}_giant".format(label_name)]) \
        / (t["w_ms"] + t["w_giant"])

    t["WEIGHTED_{}".format(label_name)][ms["USE"]] = ms[label_name][ms["USE"]]
    t["WEIGHTED_{}".format(label_name)][giant["USE"]] = giant[label_name][giant["USE"]]

"""

N = len(complete)
columns = [complete["Name"]]
label_names = ("TEFF", "LOGG", "FE_H", "VSINI", "r_chi_sq", "IN_CONVEX_HULL")
for label_name in label_names:

    data = np.nan * np.ones(N, dtype=float)
    data[ms["USE"]] = ms[label_name][ms["USE"]]
    data[giant["USE"]] = giant[label_name][giant["USE"]]
    columns.append(data)

data = 0.5 * np.ones(N, dtype=float)
data[ms["USE"]] = 0
data[giant["USE"]] = 1
columns.append(data)

names = ["Name"] + map(str.upper, label_names) + ["SOURCE"]
t = Table(data=columns, names=names)
"""



rave_dr4 = Table.read("/data/gaia-eso/arc/rave-data-files/RAVE-DR4.fits")

t = join(t, rave_dr4, keys=("Name", ))


for label_name in ("TEFF", "LOGG", "FE_H"):#, "VSINI"):
    t[label_name] = t["WEIGHTED_{}".format(label_name)]

t["R_CHI_SQ"] = np.nanmax(np.array([t["r_chi_sq_ms"], t["r_chi_sq_giant"]]), axis=0)
#t.write("/data/gaia-eso/arc/rave-data-files/unrave-v0.2-16b_23giant.fits.gz")
#t.write("/data/gaia-eso/arc/rave-data-files/unrave-v0.3-27_23giant.fits.gz")

#t.write("/data/gaia-eso/arc/rave-data-files/unrave-v0.6-16b_36.fits.gz", overwrite=True)

t.write("/data/gaia-eso/arc/rave-data-files/unrave-v0.7-37_36.fits.gz", overwrite=True)


ok = (t["SNRK"] > 50) * (t["r_chi_sq_ms"] < 3) * (t["r_chi_sq_giant"] < 3) #* (t["WEIGHTED_VSINI"] < 1)

ok = (t["SNRK"] > 50) * (t["r_chi_sq_ms"] < 3) * (t["r_chi_sq_giant"] < 3) #* (t["WEIGHTED_VSINI"] < 1)




N = 70
fig, ax = plt.subplots()

ax.hexbin(t["WEIGHTED_TEFF"][ok], t["WEIGHTED_LOGG"][ok], gridsize=N,
    extent=(3000, 7500, 0, 5.5),
    cmap="Blues", norm=LogNorm(), edgecolor="#ffffff", linewidths=0.0)

#ax.hist2d(t["WEIGHTED_TEFF"][ok], t["WEIGHTED_LOGG"][ok],
#    bins=(np.linspace(3500, 7000, N), np.linspace(0.5, 5, N)))

#    norm=LogNorm())
#    norm=LogNorm(),
#    cmap="plasma")

raise a

ax.text(0.05, 0.90, r"${:.0f}$ ${{\rm stars}}$".format(sum(ok)), horizontalalignment="left",
    verticalalignment="bottom", transform=ax.transAxes)

ax.set_xlim(ax.get_xlim()[::-1])
ax.set_ylim(ax.get_ylim()[::-1])

ax.set_xlabel(r"$T_{\rm eff}$ $({\rm K})$")
ax.set_ylabel(r"$\log{g}$")
assert sum(ok) == 208522



