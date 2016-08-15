

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.ticker import MaxNLocator

import AnniesLasso as tc

try:
    cannon_model_path

except NameError:
    from rave_io import (cannon_loocv_results, cannon_model_path, )

else:
    print("Warning: Using pre-loaded data!")


label_limits = {
    "TEFF": (3000, 8000),
    "LOGG": (0, 5),
    "FE_H": (-2, 0.75)
}

latex_labels = {
    "TEFF": r"$T_{\rm eff}$",
    "LOGG": r"$\log{g}$",
    "FE_H": r"$[{\rm Fe/H}]$"
}



model = tc.load_model(cannon_model_path)

K = len(model.vectorizer.label_names)

N = (1, 5)


assert K <= (N[0] * N[1])

fig, axes = plt.subplots(N[0], N[1])

c = model.labelled_set["SNRK"]

for i, (ax, label_name) in enumerate(zip(axes, model.vectorizer.label_names)):


    x = model.labelled_set[label_name]
    #xerr = model.labelled_set["E_{}".format(label_name)]
    y = cannon_loocv_results[label_name]
    #yerr = cannon_loocv_results["E_{}".format(label_name)]

    ax.scatter(x, y, c=c, cmap="viridis", s=50)

    _ = limits[label_name]
    ax.plot(_, _, c="#666666", zorder=-1)
    ax.set_xlim(_)
    ax.set_ylim(_)
    ax.xaxis.set_major_locator(MaxNLocator(6))
    ax.yaxis.set_major_locator(MaxNLocator(6))
    ax.set_xlabel(" ".join([latex_labels[label_name], r"$({\rm Labelled})$"]))
    ax.set_ylabel(" ".join([latex_labels[label_name], r"$({\rm LOOCV})$"]))

    [_.set_rotation(30) for _ in ax.get_xticklabels()]
    [_.set_rotation(30) for _ in ax.get_yticklabels()]



fig.tight_layout()

fig.savefig("loocv-results.pdf", dpi=300)
fig.savefig("loocv-results.png", dpi=300)
