
"""
Plot histograms of parameters (as requested by referee #1).
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from matplotlib.ticker import MaxNLocator


catalog = Table.read("RAVE-on-v1.0.fits.gz")

passed_qc = catalog["QC"]
catalog = catalog[passed_qc]


label_names = ('TEFF', 'LOGG', 'FE_H', 'O_H', 'MG_H', 'AL_H', 'SI_H', 'CA_H', 'NI_H')

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

fig, axes = plt.subplots(3, 3, figsize=(7.66, 7.06))

B = 30

for ax, label_name in zip(axes.flatten(), label_names):

    finite = np.isfinite(catalog[label_name])

    if label_name == "TEFF":
        bins = np.linspace(3500, 7000, B)

    elif label_name == "LOGG":
        bins = np.linspace(0, 5, B)

    else:
        bins = np.linspace(-1.2, 0.6, B)

    ax.hist(catalog[label_name][finite], bins=bins, facecolor="#CCCCCC",
        edgecolor="#666666", normed=True)
    
    #ax.text(0.10, 0.80, latex_labels[label_name], fontsize=14,
    #    transform=ax.transAxes)

    ax.xaxis.set_major_locator(MaxNLocator(4))

    ax.set_xlabel(latex_labels[label_name])

    ax.set_yticks([])    
    ax.set_yticklabels([])

    if ax.is_first_col():
        ax.set_ylabel(r"$N$ $({\rm normalized})$")

fig.tight_layout()

fig.savefig("label-histograms.pdf", dpi=300)
fig.savefig("label-histograms.png", dpi=300)
