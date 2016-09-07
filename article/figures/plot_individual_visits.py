
"""
Compare the results from stacked spectra to individual visits.

"""

import os
from astropy.table import Table, join

DATA_PATH = "/data/gaia-eso/arc/rave-data-files/"

stacked_results = Table.read(os.path.join(DATA_PATH, "stacked-spectra-v44.fits.gz"))
individual_results = Table.read(os.path.join(DATA_PATH, "rave-tgas-v44.fits.gz"))
stacked_results["RAVEID"] \
    = [os.path.basename(p).split("-result-")[0] for p in stacked_results["FILENAME"]]

individual_results["RAVE_OBS_ID"] \
    = [each.split("/")[-2] + "_" + each.split("/")[-1].split(".rvsun.")[0] + "_" + each.split(".rvsun.")[1].split("-")[0] for each in individual_results["FILENAME"]]


# Match the individual results to RAVE DR5.
dr5_catalog = Table.read(
    os.path.join(DATA_PATH, "RAVEDR5_PublicCut_20160905.csv.gz"), format="csv")
individual_results = join(individual_results, dr5_catalog, keys=("RAVE_OBS_ID", ))



label_names = ("TEFF", "LOGG", "FE_H")


# Get the list of results for the stacked spectra.
snrs = []
differences = {}
for label_name in label_names:
    differences[label_name] = []

# Metadata
matches = []
stacked_snrs = []

for i, rave_id in enumerate(stacked_results["RAVEID"]):

    match = (individual_results["RAVEID"] == rave_id)
    N = sum(match)
    print(i, rave_id, N)
    stacked_snrs.extend([stacked_results["snr"][i]] * N)
    assert N > 0

    # Store a list of matches, so each entry in snrs 
    snrs.extend(individual_results["snr"][match])
    for label_name in label_names:
        # Stacked result - individual result
        differences[label_name].extend(
            stacked_results[label_name][i] - individual_results[label_name][match])

    matches.extend([N] * N)
    stacked_snrs.extend([stacked_results["snr"][i]] * N)


snrs = np.array(snrs)
for label_name in label_names:
    differences[label_name] = np.array(differences[label_name])

matches = np.array(matches)
stacked_snrs = np.array(stacked_snrs)

assert snrs.size == matches.size == stacked_snrs.size



# Plot the average RMS in the difference as a function of SNR.
K = len(label_names)

# Chose the metric to apply in each SNR bin.
metric = lambda _: np.nanmean(np.abs(_))

snr_bins = np.linspace(0, 100, 11)
snr_limits = [min(snr_bins), max(snr_bins)]


ok = np.isfinite(snrs) * (snr_limits[1] >= snrs) * (snrs >= snr_limits[0]) \
    * (matches > 2) * (stacked_snrs >= 100)


common_kwds = dict(gridsize=20, cmap="Blues", extent=snr_limits + [0, 1])

specific_kwds = {
    "TEFF": dict(extent=snr_limits + [0, 300]),
    "LOGG": dict(extent=snr_limits + [0, 2.5]),
    "FE_H": dict(extent=snr_limits + [0, 2.5])
}


fig, axes = plt.subplots(K)
for j, (ax, label_name) in enumerate(zip(axes, label_names)):

    kwds = common_kwds.copy()
    kwds.update(specific_kwds[label_name])

    x = snrs[ok]
    y = differences[label_name][ok]

    indices = np.digitize(x, list(snr_bins) + [np.inf])
    y_aggregate = np.array([metric(y[indices == k]) for k in range(0, 1 + max(indices))])

    ax.hexbin(x, np.abs(y), **kwds)

    ax.scatter(snr_bins, y_aggregate, facecolor="r")
    ax.plot(snr_bins, y_aggregate, drawstyle="steps-mid")

    ax.set_xlim(snr_limits)
    ax.set_ylim(kwds["extent"][2:])

