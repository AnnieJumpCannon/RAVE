import os

import numpy as np
import cPickle as pickle
import matplotlib.pyplot as plt

from astropy.table import Table
from glob import glob
from matplotlib.ticker import MaxNLocator

scale_factors = [0.5, 1, 2, 5, 10]
Lambdas = [0, 0.001, 0.01, 0.1, 1.0, 10.0, 25.0, 50.0, 100.0]
labels = ("TEFF", "LOGG", "FE_H", "MG_H", "AL_H", "O_H", "NI_H", "SI_H", "CA_H")



# y axis should be some metric(s)
# x axis should be Lambda




#labelled_set = Table.read("{}-labelled-set.fits".format(model_name))
labelled_set = Table.read("/data/gaia-eso/arc/rave/preprocessing/giant_trainingsets_APOGEERAVE_V2.fits")

ok = np.ones(len(labelled_set), dtype=bool)
for label_name in labels:
    ok *= np.isfinite(labelled_set[label_name])
    ok *= (labelled_set[label_name] > -5)

labelled_set = labelled_set[ok]
labelled_set = np.array([labelled_set[label_name] for label_name in labels]).T



def metric(labelled_set, test_set):
    return np.nanstd(labelled_set - test_set, axis=0)


color_choices = ["r", "b", "g", "m", "y"]
color_choices = ["#3498DB"]

files = glob("/data/gaia-eso/arc/rave/rave-tgas-v34-test-labelled-set-?.????????e*-1.00.pkl")

metric_index = 0

offset = 5
x = []
colors = []
y = []

for filename in files:

    with open(filename, "rb") as fp:
        test_set = pickle.load(fp)

    basename = os.path.basename(filename)
    Lambda = basename[32:41 + offset]
    scale_factor = basename[42 + offset:46 + offset]

    x.append(float(Lambda))
    colors.append(float(scale_factor))

    y.append(metric(labelled_set, test_set))

x = np.array(x)
y = np.array(y)
colors = np.array(colors)

fig, ax = plt.subplots()
N_metrics = len(y[0])

marker_colors = dict(zip(set(colors), color_choices))

for color in set(colors):

    base_match = (colors == color) * (x == 0)

    match = colors == color
    x_ = x[match]
    y_ = 100 * (y[match] - y[base_match])/(y[base_match])

    indices = np.argsort(x_)
    x_ = x_[indices]
    y_ = y_[indices]

    

    """
    for index in range(9):

        ax.scatter(x_, y_[:, index], facecolor=marker_colors[color])
        ax.plot(x_, y_[:, index], c=marker_colors[color], label=color)
    """ 
    ax.scatter(x_, np.mean(y_, axis=1), facecolor=marker_colors[color], s=50, zorder=100)
    ax.plot(x_, np.mean(y_, axis=1), c=marker_colors[color], label=color, linewidth=2, zorder=50)
    ax.fill_between(x_, 
        np.min(y_, axis=1), #np.mean(y_, axis=1) - np.std(y_, axis=1), #np.min(y_, axis=1),
        np.max(y_, axis=1), #np.mean(y_, axis=1) + np.std(y_, axis=1), #np.max(y_, axis=1),
        facecolor=marker_colors[color], edgecolor=marker_colors[color],
        alpha=0.5, linewidth=2, zorder=10)



ax.semilogx()
ax.axhline(0, c="k", linestyle=":", zorder=-100)

ax.set_xlim(0, 101)
ax.set_ylim(-50, 50)

ax.yaxis.set_major_locator(MaxNLocator(6))

ax.set_xlabel(r"${\rm Regularization}$ ${\rm hyperparameter,}$ $\Lambda$")
ax.set_ylabel(r"${\rm Change}$ ${\rm in}$ ${\rm label}$ ${\rm RMS}$ ${\rm w.r.t.}$ $\Lambda = 0{\rm ,}$ $\Delta\sigma$ $[{\rm \%}]$")

fig.tight_layout()


fig.savefig("set-giant-regularization.pdf", dpi=300)
fig.savefig("set-giant-regularization.png")

raise a


#plt.legend()

# Show model sparsity.
tolerance = 1e-3

offset = -18

import AnniesLasso as tc

Lambda = []
scale_factor = []
sparsity = []
for filename in glob("rave-tgas-v34-?.????????e???-*.model"):
    
    Lambda.append(float(filename[32 + offset:46 + offset]))
    scale_factor.append(float(filename[47 + offset:51 + offset]))
    print(filename, Lambda[-1], scale_factor[-1])

    model = tc.load_model(filename)
    N = np.sum(np.abs(model.theta[:, 11:]) > tolerance)
    sparsity.append(float(N)/model.theta[:, 11:].size)

Lambda = np.array(Lambda)
sparsity = np.array(sparsity)
scale_factor = np.array(scale_factor)
ok = Lambda > 0

fig, ax = plt.subplots()
scat = ax.scatter(Lambda[ok], sparsity[ok], c=scale_factor[ok])

ax.semilogx()
#ax.semilogy()
ax.set_xlim(10.0**-3, 10.**2.0)
ax.set_ylim(0, 1)

cbar = plt.colorbar(scat)


