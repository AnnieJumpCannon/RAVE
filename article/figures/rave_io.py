
import os
import numpy as np
from astropy.table import Table

DATA_PATH = "/data/gaia-eso/arc/rave-data-files/"

rave_kordopatis_dr4 = Table.read(os.path.join(DATA_PATH, "RAVE-DR4.fits"))
rave_cannon_dr1 = Table.read(os.path.join(DATA_PATH, "rave-tgas-v1.fits"))
ges_idr4 = Table.read(os.path.join(DATA_PATH, "literature-ges-idr4-xm.fits"))

literature_pastel = Table.read(os.path.join(DATA_PATH, "literature-pastel-database-xm.fits"))