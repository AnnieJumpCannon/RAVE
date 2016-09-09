
"""
Compare pair-wise differences in multiple visits.

"""

import os
from astropy.table import Table, join
from matplotlib.ticker import MaxNLocator

DATA_PATH = "../../"

individual_results = Table.read(os.path.join(DATA_PATH, "unrave-v0.97.fits.gz"))

OK = individual_results["QC"]
individual_results = individual_results[OK]

if "RAVE_OBS_ID" not in individual_results.dtype.names:
    individual_results["RAVE_OBS_ID"] \
        = [each.split("/")[-2] + "_" + each.split("/")[-1].split(".rvsun.")[0] + "_" + each.split(".rvsun.")[1].split("-")[0] for each in individual_results["FILENAME"]]


# Match the individual results to RAVE DR5.
dr5_catalog = Table.read(
    os.path.join(DATA_PATH, "RAVEDR5_PublicCut_20160905.csv.gz"), format="csv")
individual_results = join(individual_results, dr5_catalog, keys=("RAVE_OBS_ID", ))



label_names = ('TEFF', 'LOGG', 'FE_H', 'O_H', 'MG_H', 'AL_H', 'SI_H', 'CA_H', 'NI_H')

# Group by RAVE_ID 
individual_results = individual_results.group_by("RAVEID")


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

pairwise_metrics = {}
for label_name in label_names:
    pairwise_metrics[label_name] = []

indices = individual_results.groups.indices
for j, si in enumerate(indices[:-1]):
    ei = indices[j + 1]
    N = ei - si

    if N == 1: continue
    
    for label_name in label_names:

        value = np.repeat(individual_results[label_name].data[si:ei], N).reshape(-1, N)
        error = np.repeat(individual_results["E_{}".format(label_name)].data[si:ei], N).reshape(-1, N)

        diff = value - value.T
        summed_error = (error**2 + error.T**2 + 2*error_floor[label_name]**2)**0.5

        # Just take the upper diagonal of the matrix so we get all pair-wise
        # comparisons, just once.
        pairwise_metric = (diff/summed_error)[np.triu_indices(N, 1)]
        pairwise_metrics[label_name].extend(pairwise_metric)



for label_name in label_names:
    finite = np.isfinite(pairwise_metrics[label_name])
    pairwise_metrics[label_name] = np.array(pairwise_metrics[label_name])[finite]



latex_labels = { 
    "TEFF": r"$T_{\rm eff}$",
    "LOGG": r"$\log{g}$",
    "FE_H": r"$[{\rm Fe/H}]$",
    "MG_H": r"$[{\rm Mg/H}]$",
    "O_H": r"$[{\rm O/H}]$",
    "AL_H": r"$[{\rm Al/H}]$",
    "CA_H": r"$[{\rm Ca/H}]$",
    "SI_H": r"$[{\rm Si/H}]$",
    "NI_H": r"$[{\rm Ni/H}]$",

}
K = len(label_names)

x = np.linspace(-6, 6, 30)
xi = np.linspace(-6, 6, 100)



fig, axes = plt.subplots(3, 3, figsize=(7.66, 7.06))

for ax, label_name in zip(axes.flatten(), label_names):

    not_outlier = np.abs(pairwise_metrics[label_name]) < 10

    ax.hist(pairwise_metrics[label_name][not_outlier], bins=x, facecolor="#CCCCCC",
        edgecolor="#666666", normed=True)
    ax.plot(xi,  np.exp(-0.5 * xi**2)/np.sqrt(2*np.pi), c='r')
    
    print(label_name, np.std(pairwise_metrics[label_name][not_outlier]))

    ax.text(0.10, 0.80, latex_labels[label_name], fontsize=14,
        transform=ax.transAxes)
    if label_name == "TEFF":
        ax.text(0.9, 0.80, r"${:.0f}$ ${{\rm K}}$".format(error_floor[label_name]),
            transform=ax.transAxes, horizontalalignment="right", fontsize=14)
    else:
        ax.text(0.9, 0.80, r"${:.2f}$ ${{\rm dex}}$".format(error_floor[label_name]),
            transform=ax.transAxes, horizontalalignment="right", fontsize=14)

    #if label_name == "TEFF":
    #    ax.text(0.15, 0.73, r"$\sigma_{{T_{{\rm eff}},floor}} = {0:.0f}$ ${{\rm K}}$".format(error_floor[label_name]),
    #        fontsize=16, transform=ax.transAxes)
    #else:
    #    ax.text(0.15, 0.73, r"$\sigma_{{{0},floor}} = {1:.2f}$".format(error_floor[label_name]),
    #        fontsize=16, transform=ax.transAxes)

    ax.set_xlim(-5, 5)
    ax.xaxis.set_major_locator(MaxNLocator(6))

    ax.set_ylim(0, 0.6)
    ax.set_yticks([0, 0.25, 0.5])
    if ax.is_last_row():
        ax.set_xlabel(r"$\eta$")
    else:   
        ax.set_xticklabels([])

    
    if not ax.is_first_col():
        ax.set_yticklabels([])

fig.tight_layout()

fig.savefig("pairwise-metrics.pdf", dpi=300)
fig.savefig("pairwise-metrics.png", dpi=300)

#40, 0.08, 0.08
# 1.67, 1.44, 1.51

# 100
# 1.44, 1.33, 1.23

# Show the distributions of the pairwise metrics.

