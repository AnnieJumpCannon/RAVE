
import os
import numpy as np
from astropy.table import Table

DATA_PATH = "/data/gaia-eso/arc/rave-data-files/"

rave_kordopatis_dr4 = Table.read(os.path.join(DATA_PATH, "RAVE-DR4.fits"))
rave_cannon_dr1 = Table.read(os.path.join(DATA_PATH, "rave-tgas-v1.fits"))


if "OK" not in rave_cannon_dr1.dtype.names:
    try:
        rave_cannon_dr1["OK"] = (rave_cannon_dr1["SNRK"] > 50) * (rave_cannon_dr1["R_CHI_SQ"] < 3)
    except KeyError:
        print("No 'OK' subset")

    else:
        print("Setting OK as stars with S/N > 50 and R_CHI_SQ < 3")

else:
    print("Using pre-defined subset of 'OK'")


ges_idr4 = Table.read(os.path.join(DATA_PATH, "literature-ges-idr4-xm.fits"))
ges_idr4["Name"] = [each.strip() for each in ges_idr4["Name"]]

literature_pastel = Table.read(os.path.join(DATA_PATH, "literature-pastel-database-xm.fits"))


def get_training_sets():

    ms_training_set = Table.read(os.path.join(DATA_PATH, "rave-ms-labelled-set.fits"))
    ms_training_set["TEFF"] = ms_training_set["EPIC_TEFF"]
    ms_training_set["LOGG"] = ms_training_set["EPIC_LOGG"]
    ms_training_set["FE_H"] = ms_training_set["EPIC_FEH"]

    giant_training_set = Table.read(os.path.join(DATA_PATH, "rave-giant-labelled-set.fits"))

    joint_training_set = Table.read(os.path.join(DATA_PATH, "rave-joint-labelled-set.fits"))

    return (ms_training_set, giant_training_set, joint_training_set)
