
"""
Plot model coefficients.
"""


import AnniesLasso as tc


giant_model = tc.load_model("/data/gaia-eso/arc/rave/rave-tgas-v23giants.model") # giants
ms_model = tc.load_model("/data/gaia-eso/arc/rave/rave-tgas-v16b.model") # dwarfs

kwds = [
    dict(color="r", alpha=0.8, label="Giant model - Combined model"),
    dict(color="g", alpha=0.8, label="Dwarf model - Combined model"),
    #dict(color="b", alpha=0.8),
    dict(color="m", alpha=0.8)
]

# For comparison sake
ms_model.vectorizer._label_names = ["TEFF", "LOGG", "FE_H", "VSINI"] 

models = (giant_model, ms_model)

comparison_model = tc.load_model("/data/gaia-eso/arc/rave/rave-tgas-v26.model")   # both

N = comparison_model.theta.shape[1]


fig, axes = plt.subplots(4, 1, sharex=True, figsize=(9, 13))

offset = 0
for i in range(N):


    term = comparison_model.vectorizer.get_human_readable_label_vector().split(" + ")[i]

    if "FE_H" not in term:
        offset += 1
        continue

    ax = axes[i - offset]

    for k, (model, kwd) in enumerate(zip(models, kwds)):

        # Get the index.
        

        try:
            index = model.vectorizer.get_human_readable_label_vector().split(" + ").index(term)
        
        except ValueError:
            print("Term {} not in model".format(term))
            continue



        y1 = model.theta[:, index].copy()

        # Scale the coefficients back by their isotropic scales.s
        
        if i > 0:
            for j, z in model.vectorizer.terms[index - 1]:
                y1 = y1 * model.vectorizer.scales[j]
                print(term, k, model, j, model.vectorizer.scales[j])

        y2 = comparison_model.theta[:, i].copy()
        if i > 0:
            for j, k in comparison_model.vectorizer.terms[i - 1]:
                y2 = y2 * comparison_model.vectorizer.scales[j]
        

        """
        if "EPIC_TEFF" in term:
            y1 = y1 / (model.theta[:, 1] * model.vectorizer.scales[0])
            y2 = y2 / (comparison_model.theta[:, 1] * comparison_model.vectorizer.scales[0])

        if "EPIC_LOGG" in term:
            print("OK LOGG", term)
            y1 = y1 / (model.theta[:, 2] * model.vectorizer.scales[1])
            y2 = y2 / (comparison_model.theta[:, 2] * comparison_model.vectorizer.scales[1])
        """
        #y2 = 0
        ax.plot(model.dispersion, y1 - y2, **kwd)
        kwd2 = kwd.copy()
        del kwd2["label"]
        ax.plot(model.dispersion, -(y1 - y2), linestyle=":", **kwd2)

    #ax.plot(comparison_model.dispersion, y2, c='b')

    ax.axhline(0, c="#666666", zorder=-10, linestyle=":")
    #ax.set_ylabel(term, rotation=0, labelpad=50)
    ax.set_title(term.replace("EPIC_", ""))

    ax.set_ylabel(r"$\Delta\theta_{" + str(int(i)) + r"}$")

    #print(" // ".join([model.vectorizer.get_human_readable_label_vector().split(" + ")[i] for model in models]))



fig.tight_layout()
fig.subplots_adjust(wspace=0, hspace=0)


raise a

def objective_function(xdata, a, b):
    y = xdata * a + b
    return y
    diff = (ms_model.theta[:, 3] - y)
    print(a,b, np.sum(diff**2))
    return diff

import scipy.optimize as op

foo = op.curve_fit(objective_function)

    

