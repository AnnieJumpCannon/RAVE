

"""
Stack RAVE spectra from repeat visits.
"""

import cPickle as pickle
import os
import numpy as np
from astropy.table import Table

parent_spectrum_dir = "/data/gaia-eso/arc/rave/pre-normalized-spectra-with-correct-errors"
stacked_spectrum_dir = os.path.join(parent_spectrum_dir, "stacked-spectra")

if not os.path.exists(stacked_spectrum_dir):
    os.mkdir(stacked_spectrum_dir)

dr5 = Table.read("/data/gaia-eso/arc/rave-data-files/rave-dr5-positions.fits")
dr5 = dr5.filled()


def get_spectrum_path(rave_obs_id):

    date, field, fibre = rave_obs_id.split("_")
    year = date[:4]

    return os.path.join(parent_spectrum_dir, year, date,
        "{0}.rvsun.{1}.pkl".format(field, fibre.strip()))



for group in dr5.group_by("GroupID").groups:
    if group["GroupID"][0] < 0 or group["GroupSize"][0] < 2: continue

    group_id = group["GroupID"][0]

    flux = []
    ivar = []
    subset = np.ones(len(group), dtype=bool)

    for i, visit in enumerate(group):

        spectrum_path = get_spectrum_path(visit["RAVE_OBS_ID"])
        if not os.path.exists(spectrum_path):
            print("Could not find {} in group {}".format(spectrum_path, group_id))
            subset[i] = False
            raise WTFError
            continue

        with open(spectrum_path, "rb") as fp:
            visit_flux, visit_ivar = pickle.load(fp)

        flux.append(visit_flux)
        ivar.append(visit_ivar)

    flux = np.array(flux)
    ivar = np.array(ivar)

    if flux.shape[0] < 2:
        print("Skipping group {} because only not enough spectra found".format(
            group_id))
        continue

    # Produce a stacked spectrum.
    stacked_ivar = np.sum(ivar, axis=0)
    stacked_flux = np.sum(flux * ivar, axis=0)/stacked_ivar

    assert np.any(np.isfinite(stacked_flux))
    if not np.all(np.isfinite(stacked_ivar)):
        print("Warning: {} pixels in {} had non-finite inverse variance".format(
            np.sum(~np.isfinite(stacked_ivar)), group["RAVEID"][0]))
        stacked_ivar[~np.isfinite(stacked_ivar)] = 0

    assert np.all(np.isfinite(stacked_ivar))

    stacked_spectrum_path = os.path.join(
        stacked_spectrum_dir, "{}.pkl".format(group["RAVEID"][0].strip()))

    with open(stacked_spectrum_path, "wb") as fp:
        pickle.dump((stacked_flux, stacked_ivar), fp, -1)

    print("Created {}".format(stacked_spectrum_path))
