
"""
Compare the results from stacked spectra to individual visits.

"""

import os
import cPickle as pickle
from astropy.table import Table, join
from matplotlib.ticker import MaxNLocator

DATA_PATH = "../../"

stacked_results = Table.read(os.path.join(DATA_PATH, "stacked-unrave-v0.95.fits.gz"))
individual_results = Table.read(os.path.join(DATA_PATH, "unrave-v0.95.fits.gz"))


pickled_filename = "individual-visits.pkl"
label_names = ("TEFF", "LOGG", "FE_H", "O_H", "MG_H", "AL_H", "SI_H", "CA_H", "NI_H")

if not os.path.exists(pickled_filename):


    # Match the individual results to RAVE DR5.
    dr5_catalog = Table.read(
        os.path.join(DATA_PATH, "RAVEDR5_PublicCut_20160905.csv.gz"), format="csv")
    individual_results = join(individual_results, dr5_catalog, keys=("RAVE_OBS_ID", ))

    individual_results = individual_results[individual_results["QC"]]


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
        #assert N > 0

        if 1 > N:
            continue

        # Store a list of matches, so each entry in snrs 
        snrs.extend(individual_results["SNR"][match])
        for label_name in label_names:
            # Stacked result - individual result
            differences[label_name].extend(
                stacked_results[label_name][i] - individual_results[label_name][match])

        matches.extend([N] * N)
        stacked_snrs.extend([stacked_results["SNR"][i]] * N)
        

    snrs = np.array(snrs)
    for label_name in label_names:
        differences[label_name] = np.array(differences[label_name])

    matches = np.array(matches)
    stacked_snrs = np.array(stacked_snrs)

    assert snrs.size == matches.size == stacked_snrs.size

    with open(pickled_filename, "wb") as fp:
        pickle.dump((snrs, differences, stacked_snrs, matches), fp, -1)

else:
    with open(pickled_filename, "rb") as fp:
        snrs, differences, stacked_snrs, matches = pickle.load(fp)


# Plot the average RMS in the difference as a function of SNR.
K = len(label_names)

# Chose the metric to apply in each SNR bin.
metric = lambda _: np.nanstd(_)

snr_bins = np.arange(0, 101, 10)


latex_labels = {
    "TEFF": r"$\sigma_{T_{\rm eff}}$ $[{\rm K}]$",
    "LOGG": r"$\sigma_{\log{g}}$ $[{\rm dex}]$",
    "FE_H": r"$\sigma_{[{\rm Fe/H}]}$ $[{\rm dex}]$",
    "O_H": r"$\sigma_{[{\rm O/H}]}$ $[{\rm dex}]$",
    "MG_H": r"$\sigma_{[{\rm Mg/H}]}$ $[{\rm dex}]$",
    "AL_H": r"$\sigma_{[{\rm Al/H}]}$ $[{\rm dex}]$",
    "SI_H": r"$\sigma_{[{\rm Si/H}]}$ $[{\rm dex}]$",
    "CA_H": r"$\sigma_{[{\rm Ca/H}]}$ $[{\rm dex}]$",
    "NI_H": r"$\sigma_{[{\rm Ni/H}]}$ $[{\rm dex}]$",
}

ok = (matches > 1) * (stacked_snrs >= 100)


fig, axes = plt.subplots(3, 3)
for j, label_name in enumerate(label_names):
    #if j >= len(axes) - 1:
    #    ax = axes[-1]
    #    label = r"${\rm " + label_name.split("_")[0].title() + r"}$"

    #else:
    #    ax = axes[j]
    #    label = None
    ax = axes.flatten()[j]

    x = snrs[ok]
    y = differences[label_name][ok]

    # Aggregate into each bin.
    y_aggregate = []
    for k, s in enumerate(snr_bins[:-1]):
        e = snr_bins[k + 1]

        in_bin = (e >= x) * (x > s)
        y_aggregate.append(metric(y[in_bin]))


    y_aggregate = np.array(y_aggregate)


    x2 = np.repeat(snr_bins[1:], 2)
    y2 = np.repeat(y_aggregate[1:], 2)

    ax.fill_between(
        x2[1:-1], np.zeros_like(y2),
        y2, facecolor="#CCCCCC", edgecolor="#CCCCCC")
    ax.plot(snr_bins[:-1] + np.diff(snr_bins)[0], y_aggregate, drawstyle="steps",
        linewidth=2, c="#666666")



    ax.set_xlim(min(snr_bins), max(snr_bins))
    ax.xaxis.set_major_locator(MaxNLocator(6))

    ax.set_ylim(0, ax.get_ylim()[1])
    ax.yaxis.set_major_locator(MaxNLocator(4))

    if ax.is_last_row():
        ax.set_xlabel(r"$S/N$ $[{\rm pixel}^{-1}]$")
    else:
        ax.set_xticklabels([])

    ax.set_ylabel(latex_labels[label_name])

    if j > 1:
        ax.set_ylim(0, 0.35)
        ax.set_yticks([0, 0.1, 0.2, 0.3])


fig.tight_layout()
fig.savefig("repeat-visits.pdf", dpi=300)
fig.savefig("repeat-visits.png")

