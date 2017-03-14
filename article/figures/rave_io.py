
import os
import numpy as np
from astropy.table import Table
from scipy.io import readsav

DATA_PATH = "/data/gaia-eso/arc/rave-data-files/"
DATA_PATH = "/Users/arc/research/rave/"

def get_rave_kordopatis_dr4():
    return Table.read(os.path.join(DATA_PATH, "RAVE-DR4.fits"))

def get_kordopatis_comparisons():
    data = readsav(os.path.join(DATA_PATH, "RAVE_DR5_calibration_data.save"))
    
    return Table(data={
        "TEFF": data["calibration_data"]["TEFF"][0],
        "LOGG": data["calibration_data"]["LOGG"][0],
        "FEH": data["calibration_data"]["FEH"][0],
        "REF": data["calibration_data"]["REF"][0],
        "Name": [each.strip() for each in data["calibration_data"]["DR5_OBS_ID"][0]]
        })




def get_cannon_dr1(filename=None):

    filename = filename or os.path.join(DATA_PATH, "RAVE-on-v1.0.fits")

    rave_cannon_dr1 = Table.read(os.path.join(DATA_PATH, filename))

    # Just so I don't have to rewrite a bunch of code.
    rave_cannon_dr1['Name'] = rave_cannon_dr1["RAVE_OBS_ID"]
    rave_cannon_dr1['snr'] = rave_cannon_dr1["SNR"]
    rave_cannon_dr1['r_chi_sq'] = rave_cannon_dr1["R_CHI_SQ"]

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
