
"""
Join models together based on a main-sequence model, a giant model, and a joint
model.
"""

import os
import matplotlib.pyplot as plt
from astropy.table import Table, vstack
from matplotlib.colors import LogNorm
from scipy.spatial import Delaunay


import AnniesLasso as tc

output_filename, overwrite = \
    ("/data/gaia-eso/arc/rave-data-files/unrave-v0.8-36-37-39.fits.gz", True)

RESULTS_PATH = "/data/gaia-eso/arc/rave/results/"
RESULTS_PATH = ""

ms_results = Table.read(os.path.join(RESULTS_PATH, "rave-tgas-v37.fits.gz"))
giant_results = Table.read(os.path.join(RESULTS_PATH, "rave-tgas-v36.fits.gz"))
joint_results = Table.read(os.path.join(RESULTS_PATH, "rave-tgas-v41.fits.gz"))

for t in (ms_results, giant_results, joint_results):
    if "Name" not in t.dtype.names:
        t["Name"] = [each.split("/")[-2] + "_" + each.split("/")[-1].split(".rvsun.")[0] + "_" + each.split(".rvsun.")[1].split("-result")[0].replace("_result.pkl", "") for each in t["FILENAME"]]

    t.sort("Name")

assert np.all(ms_results["Name"] == joint_results["Name"])
assert np.all(giant_results["Name"] == joint_results["Name"])

ms_model = tc.load_model("rave-tgas-v37.model")

# HACK MAGIC BEGINS #
ms_model._labelled_set = Table.read("rave-tgas-v16b-labelled-set-cut.fits")
ms_model._labelled_set["TEFF"] = ms_model._labelled_set["EPIC_TEFF"]
ms_model._labelled_set["LOGG"] = ms_model._labelled_set["EPIC_LOGG"]
ms_model._labelled_set["FE_H"] = ms_model._labelled_set["EPIC_FEH"]

# HACK MAGIC ENDS


ms_results["TEFF"] = ms_results["EPIC_TEFF"]
ms_results["LOGG"] = ms_results["EPIC_LOGG"]
ms_results["FE_H"] = ms_results["EPIC_FEH"]


# In the giant and main-sequence results, ignore anything with \chi_{red} > 3
for table in (ms_results, giant_results):
    bad = table["r_chi_sq"] > 3
    for label in ("TEFF", "LOGG", "FE_H"):
        table[label][bad] = np.nan



bad = joint_results["LOGG"] > 3.5
for label in ("TEFF", "LOGG", "FE_H"):
    giant_results[label][bad] = np.nan

label_names = ("TEFF", "LOGG", "FE_H")


subg = ms_model.labelled_set["LOGG"] < 4

# Exclude out a few bad eggs from the convex hull
bad1 = subg * (ms_model.labelled_set["TEFF"] < 4870)
bad2 = subg * ((ms_model.labelled_set["TEFF"] < 4965) * (ms_model.labelled_set["LOGG"] > 3.66))
keep = ~bad1 * ~bad2

ms_subg_convex_hull = Delaunay(ms_model.labels_array[:, :2][subg * keep])
in_ms_subg_hull = ms_subg_convex_hull.find_simplex(
    np.array([ms_results["TEFF"], ms_results["LOGG"]]).T) >= 0

bad = (ms_results["LOGG"] < 4) * (ms_results["TEFF"] < 5000) * ~in_ms_subg_hull
for label in ("TEFF", "LOGG", "FE_H"):
    ms_results[label][bad] = np.nan


"""


subg = ms_model.labelled_set["LOGG"] < 4

# Exclude out a few bad eggs from the convex hull
bad1 = subg * (ms_model.labelled_set["TEFF"] < 4870)
bad2 = subg * ((ms_model.labelled_set["TEFF"] < 4965) * (ms_model.labelled_set["LOGG"] > 3.66))
keep = ~bad1 * ~bad2

ms_subg_convex_hull = Delaunay(ms_model.labels_array[:, :2][subg * keep])
in_ms_subg_hull = ms_subg_convex_hull.find_simplex(
    np.array([ms_results["TEFF"], ms_results["LOGG"]]).T) >= 0

bad = (ms_results["LOGG"] < 4) * (ms_results["TEFF"] < 5000) * ~in_ms_subg_hull
for label in ("TEFF", "LOGG", "FE_H"):
    ms_results[label][bad] = np.nan


label_names = ("TEFF", "LOGG", "FE_H")
"""


"""
fig, axes = plt.subplots(2, 3)
common_kwds = dict(cmap="jet", gridsize=50, linewidths=0.5, edgecolor="#ffffff")
specific_kwds = {
    "TEFF": dict(extent=(3500, 7000, 3500, 7000), C=joint_results["LOGG"], vmin=0, vmax=5),
    "LOGG": dict(extent=(0, 5, 0, 5),  C=joint_results["TEFF"], vmin=3500, vmax=7000),
    "FE_H": dict(extent=(-4, 0.5, -4, 0.5))
}


specific_kwds = {
    "TEFF": dict(extent=(3500, 7000, 3500, 7000), norm=LogNorm()),
    "LOGG": dict(extent=(0, 5, 0, 5), norm=LogNorm()),
    "FE_H": dict(extent=(-4, 0.5, -4, 0.5), norm=LogNorm())
}


ok = joint_results["snr"] > 1
valid = (
    True,
    joint_results["TEFF"] < 5200,
    True
)

for i, (col, label_name) in enumerate(zip(axes.T, label_names)):

    x = joint_results[label_name]

    kwds = common_kwds.copy()
    kwds.update(specific_kwds[label_name])

    for j, (ax, model_results) in enumerate(zip(col, (ms_results, giant_results))):

        use = ok * valid[j]

        _ = ax.hexbin(x/141., model_results[label_name]/141, **kwds)

    raise a

"""


# OK, calculate how many sigma the stars are away from the joint model.



fig, axes = plt.subplots(1, 2, figsize=(11.55, 4.55))


x_mu, y_mu = (0, 0)
x_sigma, y_sigma = (90, 0.15)
x = ((ms_results["TEFF"] - joint_results["TEFF"]) - x_mu)/x_sigma
y = ((ms_results["LOGG"] - joint_results["LOGG"]) - y_mu)/y_sigma


axes[0].hexbin(x, y, gridsize=100, extent=(-3, +3, -3, +3), norm=LogNorm(), linewidths=0.1)


# TODO WITHOUT THIS PENALISATION TERM THE GOLD STANDARD COMPARISON IS BAD.
x_mu, y_mu = (0, 1.0)
x_sigma, y_sigma = (50, 0.15)

x2 = ((giant_results["TEFF"] - joint_results["TEFF"]) - x_mu)/x_sigma
y2 = ((giant_results["LOGG"] - joint_results["LOGG"]) - y_mu)/y_sigma


axes[1].hexbin(x2, y2, gridsize=100, extent=(-3, +3, -3, +3), norm=LogNorm(), linewidths=0.1)



ms_distance = np.sqrt(x**2 + y**2)
giant_distance = np.sqrt(x2**2 + y2**2)


# Weight by relative distance, or just take the closest of the two?
combined_data = {
    "Name": giant_results["Name"],
    "snr": joint_results["snr"],
    "r_chi_sq": joint_results["r_chi_sq"]
}

#w1 = np.exp(-0.5 * ms_distance**2)
#w2 = np.exp(-0.5 * (giant_distance**2))

w1 = 1.0/((ms_distance)**2)
w2 = 1.0/((giant_distance)**2)# + 10**2)

#assert np.isfinite(ms_distance).sum() == np.isfinite(w1).sum()
#assert np.isfinite(giant_distance).sum() == np.isfinite(w2).sum()


weights = np.array([w1, w2])
weights[~np.isfinite(weights)] = 0
weights = weights/np.sum(weights, axis=0)


ms_bad = ms_distance > 10
giant_bad = giant_distance > 10


for label_name in label_names:

    values = np.array([ms_results[label_name], giant_results[label_name]])

    values[0, ms_bad] = np.nan
    values[1, giant_bad] = np.nan

    foo = values * weights
    #assert np.isfinite(foo).sum() == np.isfinite(values).sum()

    weights2 = weights.copy()
    weights2[~np.isfinite(foo)] = 0

    combined_data[label_name] = np.nansum(foo, axis=0)/np.sum(weights2, axis=0)



combined_table = Table(data=combined_data)





N = 70
fig, ax = plt.subplots()

ok = combined_table["snr"] > 10
ax.hexbin(combined_table["TEFF"][ok], combined_table["LOGG"][ok], gridsize=N,
    extent=(3000, 7500, 0, 5.5),
    cmap="Blues", norm=LogNorm(), edgecolor="#ffffff", linewidths=0.0)




N = 70
fig, ax = plt.subplots()

ok = combined_table["snr"] > 10
ax.hexbin(combined_table["TEFF"][ok], combined_table["LOGG"][ok], gridsize=N,
    extent=(3000, 7500, 0, 5.5), C=weights2[0][ok], edgecolor="#ffffff", linewidths=0.0)




raise a

is_giant = giant_distance < ms_distance
giant_subset = giant_results[is_giant]
ms_subset = ms_results[~is_giant]
del ms_subset["COV"]
del giant_subset["COV"]

combined_table = vstack([giant_subset, ms_subset])

raise a






# In the main-sequence model, then anything with LOGG < 4 and teff < 5000 *must*
# be within the convex hull of the training set.

subg = ms_model.labelled_set["LOGG"] < 4

# Exclude out a few bad eggs from the convex hull
bad1 = subg * (ms_model.labelled_set["TEFF"] < 4870)
bad2 = subg * ((ms_model.labelled_set["TEFF"] < 4965) * (ms_model.labelled_set["LOGG"] > 3.66))
keep = ~bad1 * ~bad2

ms_subg_convex_hull = Delaunay(ms_model.labels_array[:, :2][subg * keep])
in_ms_subg_hull = ms_subg_convex_hull.find_simplex(
    np.array([ms_results["TEFF"], ms_results["LOGG"]]).T) >= 0

bad = (ms_results["LOGG"] < 4) * (ms_results["TEFF"] < 5000) * ~in_ms_subg_hull
for label in ("TEFF", "LOGG", "FE_H"):
    ms_results[label][bad] = np.nan


# Sort all in the same way.
for table in (ms_results, giant_results, joint_results):
    if "Name" not in table.dtype.names:
        table["Name"] = [each.split("/")[-2] + "_" + each.split("/")[-1].split(".rvsun.")[0] + "_" + each.split(".rvsun.")[1].split("-")[0] for each in table["FILENAME"]]

ms_results.sort("Name")
giant_results.sort("Name")
joint_results.sort("Name")

assert np.all(joint_results["Name"] == ms_results["Name"])
assert np.all(joint_results["Name"] == giant_results["Name"])

# Identify stars in the joint model that are giants by the Huber definition
is_giant = (joint_results["LOGG"] < 3.8) \
    + ((joint_results["TEFF"] > 5000) \
        * (joint_results["LOGG"] < (13.363 - 0.00191*joint_results["TEFF"])))

ms_subset = ms_results[~is_giant * np.isfinite(ms_results["TEFF"])]
giant_subset = giant_results[is_giant * np.isfinite(giant_results["TEFF"])]

# Fix columns and covariance matrices.
ms_subset["E_VSINI"] = ms_subset["COV"][:, -1, -1]**0.5
for column in ("EPIC_TEFF", "EPIC_LOGG", "EPIC_FEH"):
    del ms_subset[column]


K, L = (giant_subset["COV"].shape[-1], 3)
new_cov = np.nan * np.ones((len(ms_subset), K, K))
new_cov[:, :L, :L] = ms_subset["COV"][:, :L, :L]

del ms_subset["COV"]
ms_subset["COV"] = new_cov

print("Corrected covariance matrices based on 4-label to K-label model")

combined_table = vstack([giant_subset, ms_subset])
combined_table.write(output_filename, overwrite=overwrite)
print("Written to {}".format(output_filename))

N = 70
fig, ax = plt.subplots()

ok = combined_table["snr"] > 25
ax.hexbin(combined_table["TEFF"][ok], combined_table["LOGG"][ok], gridsize=N,
    extent=(3000, 7500, 0, 5.5),
    cmap="Blues", norm=LogNorm(), edgecolor="#ffffff", linewidths=0.0)

