""" script to analyze the PSD cut for IBD events and NC events, that pass all cuts (except of PSD cut):

    not all events are analyzed like in analyze_PSD_cut.py, but only the events that pass all cuts
    (analyzed with analyze_spectrum_v2.py)

    For each time window, the TTR values of events that pass all cuts (e.g. TTR_beforePSD_IBDevents_100ns_to_600ns.txt
    and TTR_IBDlike_NCevents_100ns_to_600ns.txt) are read and then the IBD suppression is calculated for different NC
    suppressions (very similar to analyze_PSD_cut_v2).

"""
import datetime
import numpy as np
import os
import sys
from matplotlib import pyplot as plt
from NC_background_functions import pulse_shape
from NC_background_functions import tot_efficiency

""" parameters for tail to total method: """
# INFO-me: parameters should agree with the bin-width of the time window!
# start of the tail in ns:
start_tail = [50.0, 50.0, 50.0, 100.0, 100.0, 100.0, 125.0, 125.0, 125.0, 150.0, 150.0, 150.0, 175.0, 175.0, 175.0,
              200.0, 200.0, 200.0, 225.0, 225.0, 225.0, 250.0, 250.0, 250.0, 275.0, 275.0, 275.0, 300.0, 300.0, 300.0,
              325.0, 325.0, 325.0, 350.0, 350.0]

# end of the tail in ns:
stop_tail = [600.0, 800.0, 1000.0, 600.0, 800.0, 1000.0, 600.0, 800.0, 1000.0, 600.0, 800.0, 1000.0,
             600.0, 800.0, 1000.0, 600.0, 800.0, 1000.0, 600.0, 800.0, 1000.0, 600.0, 800.0, 1000.0,
             600.0, 800.0, 1000.0, 600.0, 800.0, 1000.0, 600.0, 800.0, 1000.0, 600.0, 800.0]

# set input path, where TTR values are saved:
input_path = "/home/astro/blum/juno/atmoNC/data_NC/output_detsim_v2/" \
             "DCR_results_16000mm_10MeVto100MeV_500nsto1ms_mult1_2400PEto3400PE_dist500mm_R17700mm_PSD90/"

# loop over the different time windows:
for index in range(len(start_tail)):

    print("tail start = {0:.1f} ns".format(start_tail[index]))
    print("tail end = {0:.1f} ns".format(stop_tail[index]))

    # load TTR values of IBD events, that pass all cuts:
    ttr_IBD = np.loadtxt(input_path + "TTR_beforePSD_IBDevents_{0:.0f}ns_to_{1:.0f}ns.txt".format(start_tail[index],
                                                                                                  stop_tail[index]))

    # load TTR values of NC events, that pass all cuts (IBD-like events):
    ttr_NC = np.loadtxt(input_path + "TTR_IBDlike_NCevents_{0:.0f}ns_to_{1:.0f}ns.txt".format(start_tail[index],
                                                                                              stop_tail[index]))

    # check the efficiency of PSD for different cut-efficiencies of NC events:
    supp_NC = [93.0, 94.0, 95.0, 96.0, 97.0, 98.0, 99.0]

    # calculate the IBD suppression and the TTR cut value for all NC suppressions:
    for index1 in range(len(supp_NC)):

        # calculate IBD suppression and corresponding TTR cut value:
        supp_IBD, ttr_cut_value = tot_efficiency(ttr_IBD, ttr_NC, supp_NC[index1])

        print(supp_NC[index1])
        print(supp_IBD)
        print(ttr_cut_value)











