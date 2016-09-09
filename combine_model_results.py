
"""
Join models together based on a main-sequence model, a giant model, and a joint
model.
"""

import os
import matplotlib.pyplot as plt
from astropy.table import Table, vstack
from matplotlib.colors import LogNorm
from scipy.spatial import Delaunay
from collections import OrderedDict

import AnniesLasso as tc


RESULTS_PATH = "/data/gaia-eso/arc/rave/results/"
RESULTS_PATH = ""

ms_results = Table.read(os.path.join(RESULTS_PATH, "rave-tgas-v43.fits.gz"))
giant_results = Table.read(os.path.join(RESULTS_PATH, "rave-tgas-v42.fits.gz"))
joint_results = Table.read(os.path.join(RESULTS_PATH, "rave-tgas-v46.fits.gz"))


for t in (ms_results, giant_results, joint_results):
    if "RAVE_OBS_ID" not in t.dtype.names:
        t["RAVE_OBS_ID"] = [each.split("/")[-2] + "_" + each.split("/")[-1].split(".rvsun.")[0] + "_" + each.split(".rvsun.")[1].split("-result")[0].replace("_result.pkl", "") for each in t["FILENAME"]]

    t.sort("RAVE_OBS_ID")


assert np.all(ms_results["RAVE_OBS_ID"] == joint_results["RAVE_OBS_ID"])
assert np.all(giant_results["RAVE_OBS_ID"] == joint_results["RAVE_OBS_ID"])

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
ms_results["E_TEFF"] = ms_results["E_EPIC_TEFF"]
ms_results["E_LOGG"] = ms_results["E_EPIC_LOGG"]
ms_results["E_FE_H"] = ms_results["E_EPIC_FEH"]


# In the giant and main-sequence results, ignore anything with \chi_{red} > 3
for table in (ms_results, giant_results):
    bad = table["r_chi_sq"] > 3
    for label in ("TEFF", "LOGG", "FE_H"):
        table[label][bad] = np.nan
        table["E_{}".format(label)][bad] = np.nan


bad = joint_results["LOGG"] > 3.5
for label in ("TEFF", "LOGG", "FE_H"):
    giant_results[label][bad] = np.nan
    giant_results["E_{}".format(label)][bad] = np.nan

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
    ms_results["E_{}".format(label)][bad] = np.nan


misclassified = Table.read("dr5-misclassified.fits")
indices = []
for rave_obs_id in misclassified["RAVE_OBS_ID"]:
    match = ms_results["RAVE_OBS_ID"] == rave_obs_id.strip()
    indices.append(np.where(match)[0][0])

indices = np.array(indices)
for label in ("TEFF", "LOGG", "FE_H"):
    ms_results[label][indices] = np.nan
    ms_results["E_{}".format(label)][indices] = np.nan


# Remove stars outside the convex hull of the lower main-sequence, since they are
# actually giants 

"""
ms_lower_ms_convex_hull = Delaunay(ms_model.labels_array[:, :2][ms_model.labelled_set["TEFF"] < 4000])
in_lower_ms_convex_hull = ms_lower_ms_convex_hull.find_simplex(
    np.array([ms_results["TEFF"], ms_results["LOGG"]]).T) >= 0

bad = (ms_results["TEFF"] < 4000) * ~in_lower_ms_convex_hull
for label in ("TEFF", "LOGG", "FE_H"):
    ms_results[label][bad] = np.nan
"""



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
#z = joint_results["FE_H"]**2

axes[0].hexbin(x, y, gridsize=100, extent=(-3, +3, -3, +3), norm=LogNorm(), linewidths=0.1)


# TODO WITHOUT THIS PENALISATION TERM THE GOLD STANDARD COMPARISON IS BAD.
x_mu, y_mu = (0, 0)
x_sigma, y_sigma = (50, 0.15)

x2 = ((giant_results["TEFF"] - joint_results["TEFF"]) - x_mu)/x_sigma
y2 = ((giant_results["LOGG"] - joint_results["LOGG"]) - y_mu)/y_sigma
#z2 = ((giant_results["FE_H"] - joint_results["FE_H"]) - 0)/0.15
z2 = np.abs(joint_results["FE_H"])**2

axes[1].hexbin(x2, y2, gridsize=100, extent=(-3, +3, -3, +3), norm=LogNorm(), linewidths=0.1)



ms_distance = np.sqrt(x**2 + y**2)
giant_distance = np.sqrt(x2**2 + y2**2)


# Weight by relative distance, or just take the closest of the two?
combined_data = OrderedDict([
    #("RAVE_OBS_ID", giant_results["RAVE_OBS_ID"]),
    ("RAVE_OBS_ID", giant_results["RAVE_OBS_ID"]),
    ("TEFF", None),
    ("LOGG", None),
    ("FE_H", None),
    ("O_H",  None),
    ("MG_H", None),
    ("AL_H", None),
    ("SI_H", None),
    ("CA_H", None),
    ("NI_H", None),
    ("E_TEFF", None),
    ("E_LOGG", None),
    ("E_FE_H", None),
    ("E_O_H",  None),
    ("E_MG_H", None),
    ("E_AL_H", None),
    ("E_SI_H", None),
    ("E_CA_H", None),
    ("E_NI_H", None),
    ("SNR", joint_results["snr"]),
    ("R_CHI_SQ", joint_results["r_chi_sq"]),
    ("QC", np.ones(len(joint_results), dtype=bool))
])

    #("COV", None),
    #("snr", joint_results["snr"]),
    #("r_chi_sq", joint_results["r_chi_sq"]),
    

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


for i, label_name in enumerate(label_names):

    values = np.array([ms_results[label_name], giant_results[label_name]])

    values[0, ms_bad] = np.nan
    values[1, giant_bad] = np.nan

    foo = values * weights
    #assert np.isfinite(foo).sum() == np.isfinite(values).sum()

    weights2 = weights.copy()
    weights2[~np.isfinite(foo)] = 0

    combined_data[label_name] = np.nansum(foo, axis=0)/np.sum(weights2, axis=0)

    # And the errors.
    #errors = np.array([ms_results["COV"][:, i, i]**0.5, giant_results["COV"][:, i, i]**0.5])
    errors = np.array([ms_results["E_{}".format(label_name)], giant_results["E_{}".format(label_name)]])
    errors[0, ms_bad] = np.nan
    errors[1, giant_bad] = np.nan

    bar = errors * weights

    combined_data["E_{}".format(label_name)] = np.nansum(bar, axis=0)/np.sum(weights2, axis=0)


# Need abundances for things that are *probably* giants.
is_giant = (weights2[1] > 0.95)

abundance_label_names = ["MG_H", "AL_H", "O_H", "NI_H", "SI_H", "CA_H"]
for i, label_name in enumerate(abundance_label_names):
    combined_data[label_name] = np.nan * np.ones(len(joint_results))
    combined_data[label_name][is_giant] = giant_results[label_name][is_giant]

    # And the errors.
    combined_data["E_{}".format(label_name)] = np.nan * np.ones(len(joint_results))
    #combined_data["E_{}".format(label_name)][is_giant] = giant_results["COV"][:, i+3, i+3][is_giant]**0.5
    combined_data["E_{}".format(label_name)][is_giant] = giant_results["E_{}".format(label_name)][is_giant]


"""
# Construct a common-shape covariance matrix for each table.
K, L = (giant_results["COV"].shape[-1], 3)
new_cov = np.nan * np.ones((len(ms_results), K, K))
new_cov[:, :L, :L] = ms_results["COV"][:, :L, :L]
ms_results["COV"] = new_cov

# Produce a weighted covariance matrix.
giant_weights = weights2[1]
for label_name in ["TEFF", "LOGG", "FE_H"] + abundance_label_names:

"""


combined_table = Table(data=combined_data)


dr5_catalog_with_correct_raveids = Table.read("RAVEDR5_PublicCut_20160905.csv.gz", format="csv")
keep_columns = (
    # Positional information
    "RAVE_OBS_ID", "RAdeg", "DEdeg",

    # Velocity information.
    "HRV", "eHRV", "StdDev_HRV", "MAD_HRV", "CorrelationCoeff", "PeakHeight",
    "PeakWidth", "CorrectionRV",

    # Skyline information.
    "SkyRV",
    "eSkyRV",
    "SkyCorrelationCoeff",

    # Some flags from SPARV. Others removed if we already provide them (e.g SNR)
    "Teff_SPARV",
    "Vrot_SPARV",
    "ZeroPointFLAG",

    # Morphological information.
    "c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8", "c9", "c10", "c11", "c12", 
    "c13", "c14", "c15", "c16", "c17", "c18", "c19", "c20")



# Fix string columns.
from astropy.table import Column
for column in ("StdDev_HRV", "MAD_HRV", "Vrot_SPARV", "Teff_SPARV"):
    data = np.array([[e, np.nan][e == "NULL"] for e in dr5_catalog_with_correct_raveids[column].copy()], dtype=float)

    index = list(dr5_catalog_with_correct_raveids.dtype.names).index(column)
    del dr5_catalog_with_correct_raveids[column]


    dr5_catalog_with_correct_raveids.add_column(
        Column(data, name=column),
        index=index)


for column in dr5_catalog_with_correct_raveids.dtype.names:
    if column not in keep_columns:
        del dr5_catalog_with_correct_raveids[column]



# Rename columns, since units will go into the metadata of the FITS table.
rename_columns = [
    ("RAdeg", "RA"),
    ("DEdeg", "DEC"),
]
for before, after in rename_columns:
    dr5_catalog_with_correct_raveids.rename_column(before, after)


from astropy.table import join
combined_table = join(dr5_catalog_with_correct_raveids, combined_table, keys=("RAVE_OBS_ID", ))

# QC_FLAG
# 0 = No problems
# 1 = Don't trust the stellar parameters.
# 2 = Don't trust the stellar parameters or the abundances.
# It's a bitmask flag
# Remove hot stars based on SPARV Teff values


combined_table["QC"] = \
      (combined_table["SNR"] >= 10) \
    * (combined_table["R_CHI_SQ"] < 3) \
    * (combined_table["Teff_SPARV"] <= 7000)

metal_poor = combined_table["FE_H"] < -1.5
for label_name in abundance_label_names:
    combined_table[label_name][metal_poor] = np.nan
    combined_table["E_{}".format(label_name)][metal_poor] = np.nan
    


"""
# Make a figure
ok = (combined_table["SNR"] >= 10) * (combined_table["R_CHI_SQ"] < 3)
fig, ax = plt.subplots()
hexbin = ax.hexbin(combined_table["TEFF"][ok], combined_table["LOGG"][ok], 
    C=combined_tpwable["Teff_SPARV"][ok], reduce_C_function=np.nanmax,
    gridsize=50, vmin=4000, vmax=10000, extent=(3500, 7500, 0, 5.5),
    cmap="plasma",)

ax.set_xlim(ax.get_xlim()[::-1])
ax.set_ylim(ax.get_ylim()[::-1])

ax.set_xlabel(r"$T_{\rm eff}$ $[{\rm K}]$")
ax.set_ylabel(r"$\log{g}$")

cbar = plt.colorbar(hexbin, ticks=[4000, 5000, 6000, 7000, 8000, 9000, 10000])
cbar.set_label(r"$T_{{\rm eff},SPARV}$ $[{\rm K}]$")

fig.tight_layout()
fig.savefig("article/figures/hot-stars.pdf", dpi=300)
fig.savefig("article/figures/hot-stars.png")
"""

# Remove the SPARV temperature.
del combined_table["Teff_SPARV"]



#for column in ("TEFF", "LOGG", "FE_H", "O_H", "MG_H", "AL_H", "SI_H", "CA_H", "NI_H"):
#    combined_table[column][hot] = np.nan
#    combined_table["E_{}".format(column)] = np.nan



error_floor = {
    "TEFF": 70, #100,
    "LOGG": 0.12, #0.15,
    "FE_H": 0.06, #0.08,
    "MG_H": 0.07, #0.08,
    "O_H": 0.07, #0.08,
    "AL_H": 0.08, #0.08,
    "CA_H": 0.07, #0.08,
    "SI_H": 0.07, #0.08,
    "NI_H": 0.06, #0.08,
}

for label_name, floor in error_floor.items():
    combined_table["E_{}".format(label_name)] = np.sqrt(combined_table["E_{}".format(label_name)]**2 + floor**2)

combined_table.write("unrave-v0.97.fits.gz", overwrite=True)



N = 70
fig, ax = plt.subplots()

ok = combined_table["QC"]
ax.hexbin(combined_table["TEFF"][ok], combined_table["LOGG"][ok], gridsize=N,
    extent=(3000, 7500, 0, 5.5),
    cmap="Blues", norm=LogNorm(), edgecolor="#ffffff", linewidths=0.0)

raise a

"""
# Plot log(density) of the three models.
K = 1
factor = 3.5
lbdim = 0.2 * factor
tdim = 0.1 * factor
rdim = 0.25 * factor
whspace = 0.05
xspace = factor
yspace = factor * K + factor * (K - 1) * whspace + lbdim * (K - 1)
xdim = lbdim + xspace + rdim
ydim = lbdim + yspace + tdim

fig, ax = plt.subplots(K, 1, figsize=(xdim, ydim))
fig.subplots_adjust(
    left=lbdim/xdim, bottom=lbdim/ydim, right=(xspace + lbdim)/xdim,
    top=(yspace + lbdim)/ydim, wspace=whspace, hspace=whspace)


N = 35

ok = combined_table["snr"] > 10
_ = ax.hexbin(combined_table["TEFF"][ok], combined_table["LOGG"][ok], gridsize=N,
    extent=(3000, 8000, 0, 5.5), C=weights2[0][ok], edgecolor="#ffffff", linewidths=0.1,
    cmap="Spectral")

ax.set_xlim(ax.get_xlim()[::-1])
ax.set_ylim(ax.get_ylim()[::-1])
ax.set_xlabel(r"$T_{\rm eff}$ $[{\rm K}]$")
ax.set_ylabel(r"$\log{g}$")

cbar = plt.colorbar(_, ticks=[0, 0.25, 0.5, 0.75, 1.0])
cbar.set_label(r"$w_{ms}/(w_{ms} + w_{giant})$")

fig.tight_layout()
fig.savefig("article/figures/model-weights.pdf", dpi=300)
fig.savefig("article/figures/model-weights.png", dpi=300)

"""





raise a



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

