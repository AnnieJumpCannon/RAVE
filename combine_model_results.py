
"""
Join models together based on a main-sequence model, a giant model, and a joint
model.
"""

import matplotlib.pyplot as plt
from astropy.table import Table, vstack
from matplotlib.colors import LogNorm
from scipy.spatial import Delaunay


import AnniesLasso as tc

output_filename, overwrite = \
    ("/data/gaia-eso/arc/rave-data-files/unrave-v0.8-36-37-39.fits.gz", True)

ms_results = Table.read("/data/gaia-eso/arc/rave/results/rave-tgas-v37.fits.gz")
giant_results = Table.read("/data/gaia-eso/arc/rave/results/rave-tgas-v36.fits.gz")
joint_results = Table.read("/data/gaia-eso/arc/rave/results/rave-tgas-v39.fits.gz")

ms_model = tc.load_model("/data/gaia-eso/arc/rave/rave-tgas-v37.model")

# HACK MAGIC BEGINS #
ms_model._labelled_set = Table.read("/data/gaia-eso/arc/rave/rave-tgas-v16b-labelled-set-cut.fits")
ms_model._labelled_set["TEFF"] = ms_model._labelled_set["EPIC_TEFF"]
ms_model._labelled_set["LOGG"] = ms_model._labelled_set["EPIC_LOGG"]
ms_model._labelled_set["FE_H"] = ms_model._labelled_set["EPIC_FEH"]

ms_results["TEFF"] = ms_results["EPIC_TEFF"]
ms_results["LOGG"] = ms_results["EPIC_LOGG"]
ms_results["FE_H"] = ms_results["EPIC_FEH"]
# HACK MAGIC ENDS

# In the giant and main-sequence results, ignore anything with \chi_{red} > 3
for table in (ms_results, giant_results):
    bad = table["r_chi_sq"] > 3
    for label in ("TEFF", "LOGG", "FE_H"):
        table[label][bad] = np.nan

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

