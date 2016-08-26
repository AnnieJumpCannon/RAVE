
"""
Compare our results to the gold standard studies of Bensby et al. (2014), 
Reddy et al. (2003, 2006), and Valenti & Fischer (2005).
"""

import numpy as np
import matplotlib.pyplot as plt

try:
    bensby_2014, reddy_2003, reddy_2006, valenti_2005

except NameError: # Do you know who I am? That's Jeff Vader!

    from rave_io import \
        (get_cannon_dr1, get_literature_bensby, get_literature_reddy_2003,
            get_literature_reddy_2006, get_literature_valenti_2005)

    from astropy.table import join

    rave_cannon_dr1 = get_cannon_dr1()
    
    OK = (rave_cannon_dr1["SNRK"] > 10) * (rave_cannon_dr1["R_CHI_SQ"] < 3) * (rave_cannon_dr1["R"] > 25)
    rave_cannon_dr1 = rave_cannon_dr1[OK].filled()

    bensby_2014 = join(rave_cannon_dr1, get_literature_bensby(), keys=("Name", ))
    reddy_2003 = join(rave_cannon_dr1, get_literature_reddy_2003(), keys=("Name", ))
    reddy_2006 = join(rave_cannon_dr1, get_literature_reddy_2006(), keys=("Name", ))
    valenti_2005 = join(rave_cannon_dr1, get_literature_valenti_2005(), keys=("Name", ))

else:
    print("Using pre-loaded data")


# Create convenience functions 

def get_data(table, cannon_label, translated_label):
    x = table[translated_label]
    #xerr #TODO
    y = table[cannon_label]
    # yerr #TODO
    return (x, y)


def bensby_2014_data(cannon_label):

    label_name = {
        "TEFF": "Teff",
        "LOGG": "logg",
        "FE_H": "Fe_H"
    }[cannon_label]

    return get_data(bensby_2014, cannon_label, label_name)


def reddy_2003_data(cannon_label):

    label_name = {
        "TEFF": "Teff",
        "LOGG": "logg",
        "FE_H": "__Fe_H__2"
    }[cannon_label]

    return get_data(reddy_2003, cannon_label, label_name)


def reddy_2006_data(cannon_label):

    label_name = {
        "TEFF": "Teff",
        "LOGG": "logg",
        "FE_H": "__Fe_H__2"
    }[cannon_label]

    return get_data(reddy_2006, cannon_label, label_name)


def valenti_2005_data(cannon_label):

    label_name = {
        "TEFF": "Teff",
        "LOGG": "log_g_",
        "FE_H": "__Fe_H__2"
    }[cannon_label]

    return get_data(valenti_2005, cannon_label, label_name)



comparisons = {
    bensby_2014_data: dict(facecolor="b", marker="s"),
    reddy_2006_data: dict(facecolor="r", marker="o"),
    reddy_2003_data: dict(facecolor="w", edgecolor="k", linewidth=2, marker="o"),
    valenti_2005_data: dict(facecolor="g", marker="^")
}

label_names = ("TEFF", "LOGG", "FE_H")
N = len(label_names)

fig, axes = plt.subplots(1, N)

for ax, label_name in zip(axes, label_names):


    for comparison, kwds in comparisons.items():
        x, y = comparison(label_name)

        _, c = comparison("FE_H")

        scat = ax.scatter(x, y, c=c, vmin=-2.5, vmax=0.5, **kwds)

        comp = c > -0.7
        diff = y[comp] - x[comp]
        print(comparison, label_name, np.nanmean(diff), np.nanstd(diff))

    limits = np.array([ax.get_xlim(), ax.get_ylim()])
    limits = [limits.min(), limits.max()]
    
    ax.plot(limits, limits, c="#666666", zorder=-1, linestyle=":")
    ax.set_xlim(limits)
    ax.set_ylim(limits)


    ax.set_xlabel(label_name)

plt.colorbar(scat)