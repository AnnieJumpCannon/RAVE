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


import AnniesLasso as tc

# y axis should be some metric(s)
# x axis should be Lambda






#def metric(labelled_set, test_set):
#    return np.nanstd(labelled_set - test_set, axis=0)


comparisons = [
    (
        tc.load_model("rave-tgas-v34.model"), 
        glob("rave-tgas-v34-test-labelled-set-?.????????e*-1.00.pkl"),
        "rave-tgas-v34-test-labelled-set-0.00000000e+00-1.00.pkl",
        "#3498DB",
        r"${\rm Giant}$ ${\rm branch}$ ${\rm model}$"

    ),
    (
        tc.load_model("rave-tgas-v16b.model"),
        glob("rave-tgas-v16b-test-labelled-set-*+vsini.pkl"),
        'rave-tgas-v16b-test-labelled-set-1.00000000e-03+vsini.pkl',
        "r",
        r"${\rm Main{-}sequence}$ ${\rm model}$"
    )
]



# Plot log(density) of the three models.
K = 1
factor = 3.5
lbdim = 0.2 * factor
trdim = 0.1 * factor
whspace = 0.05
yspace = factor
xspace = factor * K + factor * (K - 1) * whspace + lbdim * (K - 1)
xdim = lbdim + xspace + trdim
ydim = lbdim + yspace + trdim

fig, ax = plt.subplots(1, K, figsize=(xdim, ydim))
fig.subplots_adjust(
    left=lbdim/xdim, bottom=lbdim/ydim, right=(xspace + lbdim)/xdim,
    top=(yspace + lbdim)/ydim, wspace=whspace, hspace=whspace)


offset = 0
for i, (model, filenames, base_filename, color, title) in enumerate(comparisons):


    x = []
    y = []

    for filename in filenames:

        with open(filename, "rb") as fp:
            test_set = pickle.load(fp)

        basename = os.path.basename(filename)
        
        _ = basename.split("-test-labelled-set-")[1]
        Lambda = float(_[:_.index("e") + 4])


        
        x.append(Lambda)
        y.append(np.nanstd(model.labels_array - test_set, axis=0))

        if i == 0 and "-5.736" in filename:
            print(i, filename, np.nanmean(model.labels_array - test_set, axis=0), np.nanstd(model.labels_array - test_set, axis=0))

    x = np.array(x)
    y = np.array(y)

    with open(base_filename, "rb") as fp:
        base_test_set = pickle.load(fp)

    y_comparison = np.nanstd(model.labels_array - base_test_set, axis=0)

    y = 100 * (y - y_comparison)/y_comparison

    indices = np.argsort(x)
    x = x[indices]
    y = y[indices]

    

    """
    for index in range(9):

        ax.scatter(x_, y_[:, index], facecolor=marker_colors[color])
        ax.plot(x_, y_[:, index], c=marker_colors[color], label=color)
    """ 
    ax.scatter(x, np.mean(y, axis=1), facecolor=color, s=50, zorder=100 - i)
    ax.plot(x, np.mean(y, axis=1), c=color, linewidth=2, zorder=-i*50,
        label=title)
    ax.fill_between(x, 
        np.min(y, axis=1), #np.mean(y_, axis=1) - np.std(y_, axis=1), #np.min(y_, axis=1),
        np.max(y, axis=1), #np.mean(y_, axis=1) + np.std(y_, axis=1), #np.max(y_, axis=1),
        facecolor=color, edgecolor=color,
        alpha=0.5, linewidth=2, zorder=-i*10)


    ax.semilogx()
    ax.axhline(0, c="k", linestyle=":", zorder=-100)


    ax.set_xlim(0, 1001)
    ax.set_ylim(-50, 50)

    #ax.text(0.05, 0.92, title, fontsize=14, transform=ax.transAxes)

    ax.yaxis.set_major_locator(MaxNLocator(6))

    if ax.is_first_col():
        ax.set_ylabel(r"${\rm Change}$ ${\rm in}$ ${\rm label}$ ${\rm RMS}$ ${\rm w.r.t.}$ $\Lambda = 0{\rm ,}$ $\Delta\sigma$ $[{\rm \%}]$")
    
    else:
        ax.set_yticklabels([])

    ax.set_xlabel(r"${\rm Regularization}$ ${\rm hyperparameter,}$ $\Lambda$")
    

plt.legend(loc="lower right", frameon=False)

fig.tight_layout()

fig.savefig("set-hyperparameters.pdf", dpi=300)
fig.savefig("set-hyperparameters.png", dpi=300)

raise a

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
    y.append(metric(labelled_set, test_set))


x = np.array(x)
y = np.array(y)
colors = np.array(colors)

N_metrics = len(y[0])





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


