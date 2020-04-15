""" script to analyze the PSD cut for IBD events and NC events, that pass all cuts (except of PSD cut):

    not all events are analyzed like in analyze_PSD_cut.py, but only the events that pass all cuts
    (analyzed with analyze_spectrum_v2.py)

    For each time window, the TTR values of events that pass all cuts (e.g. TTR_beforePSD_IBDevents_100ns_to_600ns.txt
    and TTR_IBDlike_NCevents_100ns_to_600ns.txt) are read and then the IBD suppression is calculated for different NC
    suppressions (very similar to analyze_PSD_cut_v2).

"""
import numpy as np
from NC_background_functions import tot_efficiency

# flag if PSD efficiency is independent of energy or not:
PSD_energy_independent = False

""" parameters for tail to total method: """
# INFO-me: parameters should agree with the bin-width of the time window!
# start of the tail in ns:
start_tail = [225.0, 225.0, 225.0, 250.0, 250.0, 250.0, 275.0, 275.0, 275.0, 300.0, 300.0, 300.0,
              325.0, 325.0, 325.0, 350.0, 350.0]

# end of the tail in ns:
stop_tail = [600.0, 800.0, 1000.0, 600.0, 800.0, 1000.0, 600.0, 800.0, 1000.0, 600.0, 800.0, 1000.0,
             600.0, 800.0, 1000.0, 600.0, 800.0]

# set input path, where TTR values are saved:
input_path = "/home/astro/blum/juno/atmoNC/data_NC/output_detsim_v2/" \
             "DCR_results_16000mm_10MeVto100MeV_1000nsto1ms_mult1_1800keVto2550keV_dist500mm_R17700mm_PSD90/"

if PSD_energy_independent:
    # if PSD_energy_independent = True: PSD efficiency calculation does not depend on energy of prompt signal!

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

else:
    # if PSD_energy_independent = False: PSD efficiency calculation does depend on energy of prompt signal!

    # loop over the different time windows:
    for index in range(len(start_tail)):

        print("")
        print("tail start = {0:.1f} ns".format(start_tail[index]))
        print("tail end = {0:.1f} ns".format(stop_tail[index]))

        # load TTR values of IBD events, that pass all cuts:
        ttr_IBD = np.loadtxt(input_path + "TTR_beforePSD_IBDevents_{0:.0f}ns_to_{1:.0f}ns.txt".format(start_tail[index],
                                                                                                      stop_tail[index]))

        # load TTR values of NC events, that pass all cuts (IBD-like events):
        ttr_NC = np.loadtxt(input_path + "TTR_IBDlike_NCevents_{0:.0f}ns_to_{1:.0f}ns.txt".format(start_tail[index],
                                                                                                  stop_tail[index]))

        # load filenumber, evtID and Evis of IBD events that pass all cuts:
        IBD_array = np.loadtxt(input_path + "IBD_filenumber_evtID_Evis_pass_all_cuts_wo_PSD.txt")

        # load filenumber, evtID and Evis of NC events that pass all cuts (IBD-like events):
        NC_array = np.loadtxt(input_path + "atmoNC_filenumber_evtID_Evis_pass_all_cuts_wo_PSD.txt")

        # get visible energy of IBD events that pass all cuts:
        Evis_IBD = IBD_array[:, 2]

        # get visible energy of NC events that pass all cuts:
        Evis_NC = NC_array[:, 2]

        # preallocate arrays, where TTR values of save depending on their energies:
        ttr_IBD_10_20 = []
        ttr_IBD_20_30 = []
        ttr_IBD_30_40 = []
        ttr_IBD_40_100 = []
        # loop over ttr_IBD and fill arrays depending on their energy:
        for index2 in range(len(ttr_IBD)):
            if Evis_IBD[index2] <= 20.0:
                ttr_IBD_10_20.append(ttr_IBD[index2])
            elif 20.0 < Evis_IBD[index2] <= 30.0:
                ttr_IBD_20_30.append(ttr_IBD[index2])
            elif 30.0 < Evis_IBD[index2] <= 40.0:
                ttr_IBD_30_40.append(ttr_IBD[index2])
            else:
                ttr_IBD_40_100.append(ttr_IBD[index2])

        # preallocate arrays, where TTR values of save depending on their energies:
        ttr_NC_10_20 = []
        ttr_NC_20_30 = []
        ttr_NC_30_40 = []
        ttr_NC_40_100 = []
        # loop over ttr_NC and fill arrays depending on their energy:
        for index2 in range(len(ttr_NC)):
            if Evis_NC[index2] <= 20.0:
                ttr_NC_10_20.append(ttr_NC[index2])
            elif 20.0 < Evis_NC[index2] <= 30.0:
                ttr_NC_20_30.append(ttr_NC[index2])
            elif 30.0 < Evis_NC[index2] <= 40.0:
                ttr_NC_30_40.append(ttr_NC[index2])
            else:
                ttr_NC_40_100.append(ttr_NC[index2])

        # check the efficiency of PSD for different cut-efficiencies of NC events:
        # supp_NC = [95.0, 97.0, 98.0, 99.0, 99.9]
        supp_NC = [99.99]

        # calculate the IBD suppression and the TTR cut value for all NC suppressions:
        for index1 in range(len(supp_NC)):

            # calculate IBD suppression and corresponding TTR cut value:
            # supp_IBD_10_20, ttr_cut_value_10_20 = tot_efficiency(ttr_IBD_10_20, ttr_NC_10_20, supp_NC[index1])

            # print("10 MeV to 20 MeV:")
            # print(supp_NC[index1])
            # print(supp_IBD_10_20)
            # print(ttr_cut_value_10_20)

            # calculate IBD suppression and corresponding TTR cut value:
            # supp_IBD_20_30, ttr_cut_value_20_30 = tot_efficiency(ttr_IBD_20_30, ttr_NC_20_30, supp_NC[index1])

            # print("20 MeV to 30 MeV:")
            # print(supp_NC[index1])
            # print(supp_IBD_20_30)
            # print(ttr_cut_value_20_30)

            # calculate IBD suppression and corresponding TTR cut value:
            # supp_IBD_30_40, ttr_cut_value_30_40 = tot_efficiency(ttr_IBD_30_40, ttr_NC_30_40, supp_NC[index1])

            # print("30 MeV to 40 MeV:")
            # print(supp_NC[index1])
            # print(supp_IBD_30_40)
            # print(ttr_cut_value_30_40)

            # calculate IBD suppression and corresponding TTR cut value:
            supp_IBD_40_100, ttr_cut_value_40_100 = tot_efficiency(ttr_IBD_40_100, ttr_NC_40_100, supp_NC[index1])

            # print("40 MeV to 100 MeV:")
            print(supp_NC[index1])
            print(supp_IBD_40_100)
            print(ttr_cut_value_40_100)
            print("")








