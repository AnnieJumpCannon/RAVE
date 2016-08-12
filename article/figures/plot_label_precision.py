
"""
Plot label precision as a function of S/N for repeated stars in the sample.
"""

import numpy as np
from collections import Counter

try:
    combined_table

except NameError:
    from rave_io import (rave_cannon_dr1, rave_kordopatis_dr4)

    from astropy.table import join

    combined_table = join(rave_cannon_dr1, rave_kordopatis_dr4, keys=("Name", ))

else:
    print("Warning! Using pre-loaded data")

# TODO: Stack all the spectra and do an analysis on the joint spectrum so that
# we have some reasonable high S/N spectra to compare each group against.

# Identify the groups.
combined_table = combined_table.group_by("RAVE")
M = len(combined_table.groups)

snr_column = "SNRK"
labels = ["TEFF", "LOGG", "FE_H", "ALPHA"]

deltas = {}
snr_values = []
snr_reference_value = []
for i, group in enumerate(combined_table.groups):

    print(i, M)

    N = len(group)
    if N == 1: continue

    # Get the highest S/N star as a reference.
    index = np.argmax(group[snr_column])
    ok = np.ones(N, dtype=bool)
    ok[index] = False

    # Record the SNR values for prosperity
    snr_values.append(group[snr_column][ok])
    snr_reference_value.append(group[snr_column][index])
    

    for label in labels:
        deltas.setdefault(label, [])   
        deltas[label].append(group[label][ok] - group[label][index])

        if N > 2:
            print(N, label, np.std(deltas[label][-1]))

