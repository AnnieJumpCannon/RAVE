
import os
import numpy as np
from astropy.table import Table

DATA_PATH = "/data/gaia-eso/arc/rave-data-files/"


def get_rave_kordopatis_dr4():
    return Table.read(os.path.join(DATA_PATH, "RAVE-DR4.fits"))



def get_cannon_dr1(filename=None):

    #filename = "unrave-v0.2-16b_23giant.fits.gz"
    #filename = "unrave-v0.3-27_23giant.fits.gz"
    #filename = "unrave-v0.4-29_23giant.fits.gz"
    #filename = "rave-tgas-v31.fits.gz"
    #filename = "rave-tgas-v29.fits.gz"
    #filename = "rave-tgas-v37.fits.gz"
    filename = filename or "unrave-v0.7-37_36.fits.gz"

    rave_cannon_dr1 = Table.read(os.path.join(DATA_PATH, filename))
    
    #rave_cannon_dr1["TEFF"] = rave_cannon_dr1["EPIC_TEFF"]
    #rave_cannon_dr1["LOGG"] = rave_cannon_dr1["EPIC_LOGG"]
    #rave_cannon_dr1["FE_H"] = rave_cannon_dr1["EPIC_FEH"]
    
    print("LOADED FROM FILENAME {}".format(filename))
    if "OK" not in rave_cannon_dr1.dtype.names:
        try:
            rave_cannon_dr1["OK"] = (rave_cannon_dr1["SNRK"] > 10) * (rave_cannon_dr1["R_CHI_SQ"] < 3)
        except KeyError:
            print("No 'OK' subset")

        else:
            print("Setting OK as stars with S/N > 50 and R_CHI_SQ < 3")

    else:
        print("Using pre-defined subset of 'OK'")

    return rave_cannon_dr1


def get_ges_idr4():

    ges_idr4 = Table.read(os.path.join(DATA_PATH, "literature-ges-idr4-xm.fits"))
    ges_idr4["Name"] = [each.strip() for each in ges_idr4["Name"]]
    return ges_idr4


def get_literature_pastel():
    return Table.read(os.path.join(DATA_PATH, "literature-pastel-database-xm.fits"))


def get_literature_bensby():
    data = Table.read(os.path.join(DATA_PATH, "literature-bensby-2014.fits"))
    data["Name"] = [each.strip() for each in data["Name"]]
    return data



def get_literature_valenti_2005():
    data = Table.read(os.path.join(DATA_PATH, "literature-valenti-2005.fits"))
    data["Name"] = [each.strip() for each in data["Name"]]
    return data



def get_literature_reddy_2003():
    data = Table.read(os.path.join(DATA_PATH, "literature-reddy-2003.fits"))
    data["Name"] = [each.strip() for each in data["Name"]]
    return data

def get_literature_reddy_2006():
    data = Table.read(os.path.join(DATA_PATH, "literature-reddy-2006.fits"))
    data["Name"] = [each.strip() for each in data["Name"]]
    return data


def get_training_sets():

    ms_training_set = Table.read(os.path.join(DATA_PATH, "rave-ms-labelled-set.fits"))
    ms_training_set["TEFF"] = ms_training_set["EPIC_TEFF"]
    ms_training_set["LOGG"] = ms_training_set["EPIC_LOGG"]
    ms_training_set["FE_H"] = ms_training_set["EPIC_FEH"]

    giant_training_set = Table.read(os.path.join(DATA_PATH, "rave-giant-labelled-set.fits"))

    joint_training_set = Table.read(os.path.join(DATA_PATH, "rave-joint-labelled-set.fits"))

    return (ms_training_set, giant_training_set, joint_training_set)
