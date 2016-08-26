
"""
Make label comparisons with Reddy et al. (2003, 2006).
"""

import numpy as np
import matplotlib.pyplot as plt


try:
    vf_data

except NameError: # Do you know who I am?

    from rave_io import get_cannon_dr1, get_literature_valenti_2005


    rave_cannon_dr1 = get_cannon_dr1()
    #OK = (data["SNRK"] > 10) * (data["R_CHI_SQ"] < 3) * (data["R"] > 25)
    #rave_cannon_dr1 = rave_cannon_dr1[OK]

    valenti_2005 = get_literature_valenti_2005()
    
    from astropy.table import join

    vf_data = join(rave_cannon_dr1, valenti_2005, keys=("Name", ))
    
    ok = (vf_data["SNRK"] > 10) * (vf_data["R_CHI_SQ"] < 3) * (vf_data["R"] > 25)
    vf_data = vf_data[ok].filled()

else:
    print("Using pre-loaded data!")

latex_labels = {
    "TEFF": r"$T_{\rm eff}$",
    "LOGG": r"$\log{g}$",
    "FE_H": r"$[{\rm Fe/H}]$"
}

def scatter_comparison(axis, literature_label, label_name, c=None):
    """
    Show a scatter plot on the given axis.
    """

    x = vf_data[literature_label]
    y = vf_data[label_name]
    c = vf_data["__Fe_H__2"]

    axis.scatter(x, y, c=c)

    limits = np.array([ax.get_xlim(), ax.get_ylim()]).flatten()
    limits = [np.min(limits), np.max(limits)]    
    axis.plot(limits, limits, c="#666666", zorder=-1, linestyle=":")
    axis.set_xlim(limits)
    axis.set_ylim(limits)

    diff = y - x
    print(label_name, np.nanmean(diff), np.nanstd(diff))

    axis.set_xlabel(" ".join([latex_labels[label_name], r"$({\rm Valenti}+$ $2005)$"]))
    axis.set_ylabel(" ".join([latex_labels[label_name], r"$({\rm unRAVE})$"]))


# Compare teff, logg, [Fe/H]
fig, axes = plt.subplots(1, 3)

labels = [
    ("TEFF", "Teff"),
    ("LOGG", "log_g_"),
    ("FE_H", "__Fe_H__2")
]

for ax, (cannon_label, literature_label) in zip(axes, labels):
    scatter_comparison(ax, literature_label, cannon_label)



# Compare abundances.


