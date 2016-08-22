

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.ticker import MaxNLocator

try:
    data_table

except NameError:
    from rave_io import (rave_cannon_dr1, literature_pastel)

    from astropy.table import join

    data_table = join(rave_cannon_dr1, literature_pastel, keys=("Name",))

else:
    print("Warning: Using pre-loaded data!")

"""
# Kordopatis et al compared against:


Ruchti et al (2011); 229 spectra            --> Acquired
CFLIB (Pastel database; 224+169 spectra)    --> Acquired
Fulbright et al (2010) (163 spectra)        --> Data not in Vizier
M67 (Pancino et al 2010; 16 spectra)        --> Stellar parameters not in Vizier
IC4651 (Pasquini et al 2004; 6 spectra)     --> Only 6 stars; uninformative
"""


# Plot labels against literature.

label_names = ("TEFF", "LOGG", "FE_H")
pastel_label_names = {
    "TEFF": "PASTEL_TEFF",
    "LOGG": "PASTEL_LOGG",
    "FE_H": "PASTEL_FE_H",
}

label_limits = {
    "TEFF": (3000, 8000),
    "LOGG": (0, 5),
    "FE_H": (-3, 0.75)
}

latex_labels = {
    "TEFF": r"$T_{\rm eff}$",
    "LOGG": r"$\log{g}$",
    "FE_H": r"$[{\rm Fe/H}]$"
}

"""
exclude_comparisons = ("2013MNRAS.429..126R", )
ok = np.ones(len(data_table), dtype=bool)
for bibcode in exclude_comparisons:
    ok *= (data_table["bibcode"] != bibcode)
"""

# Exclude ones we think are bad.
ok = data_table["OK"]

data_table = data_table[ok]

# Assign markers random colours and orientations based on the author's name,
# just to keep things common between figures.
# This idea and execution was inspired by David W Hogg's Platypus repository:
# https://github.com/davidwhogg/Platypus
data_table.sort(["bibcode", "Author"])

unique_authors = np.unique(data_table["Author"])
author_ids = np.array([
    np.where(author == unique_authors)[0][0] for author in data_table["Author"]])
M = len(unique_authors)

cmap = plt.get_cmap("gist_earth")
angles = (120. / M) * author_ids
colors = cmap(author_ids/float(M))

K = len(label_names)
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


for i, (ax, label_name) in enumerate(zip(axes, label_names)):


    x = data_table[label_name]
    #xerr
    y = data_table[pastel_label_names[label_name]]
    # yerr

    for j in range(len(x)):
        ax.scatter([x[j]], [y[j]],
            c=colors[j], alpha=0.75, marker=(3, 1, angles[j]), s=50, 
            linewidths=0.5, edgecolors="#000000")

    ax.set_xlim(label_limits[label_name])
    ax.set_ylim(label_limits[label_name])
    ax.xaxis.set_major_locator(MaxNLocator(6))
    ax.yaxis.set_major_locator(MaxNLocator(6))
    ax.set_xlabel(" ".join([latex_labels[label_name], r"$({\rm UNRAVE})$"]))
    ax.set_ylabel(" ".join([latex_labels[label_name], r"$({\rm Literature})$"]))

    [_.set_rotation(30) for _ in ax.get_xticklabels()]
    [_.set_rotation(30) for _ in ax.get_yticklabels()]

fig.tight_layout()

fig.savefig("literature-comparison.pdf", dpi=300)
fig.savefig("literature-comparison.png")
