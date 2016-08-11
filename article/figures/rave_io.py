
import numpy as np
from astropy.table import Table

rave_kordopatis_dr4 = Table.read("../../RAVE-DR4.fits")
rave_cannon_dr1 = Table.read("../../results/rave-tgas-v3a+b-xm-dr4.fits")

w1 = np.exp(-0.5 * rave_cannon_dr1["chi-sq_1"])
w2 = np.exp(-0.5 * rave_cannon_dr1["chi-sq_2"])
for label in ("TEFF", "LOGG", "FE_H"):
    rave_cannon_dr1[label] \
        = (rave_cannon_dr1["{}_1".format(label)] * w1 \
            + rave_cannon_dr1["{}_2".format(label)] * w2)/(w1 + w2)
