""" Script to do a pulse shape analysis of atmospheric NC neutrino events that could mimic an IBD signal in JUNO
    detector.

    The hittime distribution of the prompt signal of the NC event, which passes all the cuts (and therefore mimics an
    IBD signal), are save in folder output_detsim (the hittime distributions are analyzed with
    prompt_signal_preselected_evts.py).

    As reference the hittime distribution of positron events are also analyzed for kinetic energy of positrons
    of 10 MeV and 100 MeV.

    The hittime distributions of positrons are equal to the hittime distribution of the prompt signal of real IBD events
    and therefore act as reference.

"""
import datetime
import os
import numpy as np
from matplotlib import pyplot as plt
from NC_background_functions import pulse_shape
from NC_background_functions import tot_efficiency

# get the date and time, when the script was run:
date = datetime.datetime.now()
now = date.strftime("%Y-%m-%d %H:%M")

# path, where output is saved:
output_path = "/home/astro/blum/juno/atmoNC/data_NC/output_PSD/"

# Set flag, if plots should be saved:
DISPLAY_PLOTS = True

""" parameters for tail to total method: """
# INFO-me: parameters should agree with the bin-width of the time window!
# start of the tail in ns:
# start_tail = np.arange(320, 365, 5)
start_tail = np.array([335, 340, 350])

# end of the tail in ns:
# stop_tail = np.arange(400, 720, 20)
stop_tail = np.array([540, 600])

print("start_tail = {0}".format(start_tail))
print("stop_tail = {0}".format(stop_tail))

""" parameters that define the time window of prompt signal: """
# start of the time window in ns:
start_time = 0.0
# end of the time window in ns:
end_time = 2000.0

""" preallocate values: """
# best positron efficiency for NC efficiency of 90 % and the corresponding tail-to-total value and values of start and
# stop of tail:
best_pos_eff_for_NC_90 = 100.0
tot_value_90 = 0.0
tail_start_90 = 0.0
tail_stop_90 = 0.0
# best positron efficiency for NC efficiency of 95 % and the corresponding tail-to-total value and values of start and
# stop of tail:
best_pos_eff_for_NC_95 = 100.0
tot_value_95 = 0.0
tail_start_95 = 0.0
tail_stop_95 = 0.0
# best positron efficiency for NC efficiency of 96 % and the corresponding tail-to-total value and values of start and
# stop of tail:
best_pos_eff_for_NC_96 = 100.0
tot_value_96 = 0.0
tail_start_96 = 0.0
tail_stop_96 = 0.0
# best positron efficiency for NC efficiency of 97 % and the corresponding tail-to-total value and values of start and
# stop of tail:
best_pos_eff_for_NC_97 = 100.0
tot_value_97 = 0.0
tail_start_97 = 0.0
tail_stop_97 = 0.0
# best positron efficiency for NC efficiency of 98 % and the corresponding tail-to-total value and values of start and
# stop of tail:
best_pos_eff_for_NC_98 = 100.0
tot_value_98 = 0.0
tail_start_98 = 0.0
tail_stop_98 = 0.0
# best positron efficiency for NC efficiency of 99 % and the corresponding tail-to-total value and values of start and
# stop of tail:
best_pos_eff_for_NC_99 = 100.0
tot_value_99 = 0.0
tail_start_99 = 0.0
tail_stop_99 = 0.0

# """ analyze the hittime distribution of the 10 MeV positrons: """
# print("analyze 10 MeV positrons...")
# # path, where hittime distributions of 10 MeV positrons are saved:
# input_path_positron10 = "/home/astro/blum/juno/atmoNC/data_NC/output_PSD/positron_hittime/"
#
# # kinetic energy of positron:
# kinetic_energy_10 = 10
#
# # number of events that are analyzed:
# num_events_analyzed_10 = 0
#
# # preallocate array, where tail-to-total ratios are stored:
# array_tot_ratio_10 = []
#
# # loop over all files in folder input_path_positron10, that start with 'file' and end with '10_MeV.txt'
# # (files where hittime distribution is saved, each file is equal to one event):
# for file_10 in os.listdir(input_path_positron10):
#     if file_10.startswith("file") and file_10.endswith("10_MeV.txt"):
#
#         # increment num_event_analyzed_10:
#         num_events_analyzed_10 += 1
#
#         # get the file name:
#         file_name_10 = input_path_positron10 + file_10
#
#         # read txt file:
#         file_data_10 = np.loadtxt(file_name_10)
#         # 0th entry in file_data_10 is minimum of time window in ns:
#         min_time_10 = file_data_10[0]
#         # 1st entry in file_data_10 is maximum of time window in ns:
#         max_time_10 = file_data_10[1]
#         # 2nd entry in file_data_10 is bin-width in ns:
#         bin_width = file_data_10[2]
#
#         # the rest of file_data_10 is the hittime distribution histogram in nPE per bin:
#         number_pe_per_bin_10 = file_data_10[3:]
#
#         # time window corresponding to number_pe_per_bin_10:
#         time_window_10 = np.arange(min_time_10, max_time_10 + bin_width, bin_width)
#
#         # analyze the hittime distribution of these event:
#         tot_ratio_10, npe_norm_10 = pulse_shape(time_window_10, number_pe_per_bin_10, start_tail, stop_tail)
#
#         # append tail-to-total ratio to array:
#         array_tot_ratio_10.append(tot_ratio_10)
#
#     else:
#         continue

# """ analyze the hittime distribution of the 100 MeV positrons: """
# print("analyze 100 MeV positrons...")
# # path, where hittime distributions of 100 MeV positrons are saved:
# input_path_positron100 = "/home/astro/blum/juno/atmoNC/data_NC/output_PSD/positron_hittime/"
#
# # kinetic energy of positron:
# kinetic_energy_100 = 100
#
# # number of events that are analyzed:
# num_events_analyzed_100 = 0
#
# # preallocate array, where tail-to-total ratios are stored:
# array_tot_ratio_100 = []
#
# # loop over all files in folder input_path_positron100, that start with 'file' and end with '100_MeV.txt'
# # (files where hittime distribution is saved, each file is equal to one event):
# for file_100 in os.listdir(input_path_positron100):
#     if file_100.startswith("file") and file_100.endswith("100_MeV.txt"):
#
#         # increment num_event_analyzed_100:
#         num_events_analyzed_100 += 1
#
#         # get the file name:
#         file_name_100 = input_path_positron100 + file_100
#
#         # read txt file:
#         file_data_100 = np.loadtxt(file_name_100)
#         # 0th entry in file_data_100 is minimum of time window in ns:
#         min_time_100 = file_data_100[0]
#         # 1st entry in file_data_100 is maximum of time window in ns:
#         max_time_100 = file_data_100[1]
#         # 2nd entry in file_data_100 is bin-width in ns:
#         bin_width = file_data_100[2]
#
#         # the rest of file_data_100 is the hittime distribution histogram in nPE per bin:
#         number_pe_per_bin_100 = file_data_100[3:]
#
#         # time window corresponding to number_pe_per_bin_100:
#         time_window_100 = np.arange(min_time_100, max_time_100 + bin_width, bin_width)
#
#         # analyze the hittime distribution of these event:
#         tot_ratio_100, npe_norm_100 = pulse_shape(time_window_100, number_pe_per_bin_100, start_tail, stop_tail)
#
#         # append tail-to-total ratio to array:
#         array_tot_ratio_100.append(tot_ratio_100)
#
#     else:
#         continue

# loop over different start values of the tail:
for index in range(len(start_tail)):

    # loop over different stop values of the tail:
    for index1 in range(len(stop_tail)):

        """ analyze the hittime distribution of positrons with kinetic energy uniformly distributed from 10 MeV to 
        100 MeV: """
        print("analyze positrons...")

        # path, where hittime distributions of 100 MeV positrons are saved:
        input_path_positron = "/home/astro/blum/juno/atmoNC/data_NC/output_PSD/positron_hittime/"

        # number of events that are analyzed:
        num_events_analyzed_positron = 0

        # preallocate array, where tail-to-total ratios are stored:
        array_tot_ratio_positron = []

        # preallocate array, where a average hittime-distribution of positrons are stored (number of pe per bin):
        # length of average hittime distribution (number of bins):
        length_average_hittime = 500
        hittime_average_positron = np.zeros(length_average_hittime)

        # loop over all files in folder input_path_positron, that start with 'file' and end with 'positron.txt'
        # (files where hittime distribution is saved, each file is equal to one event):
        for file_positron in os.listdir(input_path_positron):
            if file_positron.startswith("file") and file_positron.endswith("positron.txt"):

                # increment num_event_analyzed_positron:
                num_events_analyzed_positron += 1

                # get the file name:
                file_name_positron = input_path_positron + file_positron

                # read txt file:
                file_data_positron = np.loadtxt(file_name_positron)
                # 0th entry in file_data_positron is minimum of time window in ns:
                min_time_positron = file_data_positron[0]
                # 1st entry in file_data_positron is maximum of time window in ns:
                max_time_positron = file_data_positron[1]
                # 2nd entry in file_data_positron is bin-width in ns:
                bin_width = file_data_positron[2]

                # the rest of file_data_positron is the hittime distribution histogram in nPE per bin:
                number_pe_per_bin_positron = file_data_positron[3:]

                # check if max_time_positron is greater than end_time:
                if max_time_positron > end_time:
                    # prompt signal is longer than time window.
                    print("max_time_positron {0:.2f} ns > end_time {1:.1f} ns in file {2}".format(max_time_positron,
                                                                                                  end_time,
                                                                                                  file_name_positron))

                # time window:
                time_window_positron = np.arange(min_time_positron, end_time + bin_width, bin_width)

                # compare len(time_window_positron) with len(number_pe_per_bin_positron):
                missing_zeros = len(time_window_positron) - len(number_pe_per_bin_positron)

                # append the missing_zeros to number_pe_per_bin_positron:
                number_pe_per_bin_positron = np.pad(number_pe_per_bin_positron, (0, missing_zeros), 'constant',
                                                    constant_values=(0.0, 0.0))

                # analyze the hittime distribution of these event:
                tot_ratio_positron, npe_norm_positron = pulse_shape(time_window_positron, number_pe_per_bin_positron,
                                                                    start_tail[index], stop_tail[index1])

                # append tail-to-total ratio to array:
                array_tot_ratio_positron.append(tot_ratio_positron)

                # append zeros to npe_norm_positron to get a average length of the hittimes:
                npe_norm_positron = np.pad(npe_norm_positron, (0, length_average_hittime - len(npe_norm_positron)),
                                           'constant', constant_values=(0.0, 0.0))

                # add the normalized hittime distribution (npe_norm_positron) to the average hittime distribution
                # (hittime_average_positron):
                hittime_average_positron = hittime_average_positron + npe_norm_positron

            else:
                continue

        # to get the average hittime distribution with a maximum of 1, normalize hittime_average_positron with
        # max(hittime_average_positron):
        hittime_average_positron = hittime_average_positron / max(hittime_average_positron)

        """ analyze the hittime distribution of NC events that can mimic IBD signal: """
        print("analyze NC events...")
        # path, where hittime distributions of preselected NC events are saved:
        input_path_NCevents = "/home/astro/blum/juno/atmoNC/data_NC/output_detsim/"

        # number of events that are analyzed:
        num_events_analyzed_NC = 0

        # preallocate array, where tail-to-total ratios are stored:
        array_tot_ratio_NC = []

        # preallocate array, where a average hittime-distribution of NC events are stored (number of pe per bin):
        hittime_average_NC = np.zeros(length_average_hittime)

        # loop over all files in folder input_path_NCevents, that start with "file" and end with ".txt"
        # (file where hittime distribution is saved, each file is equal to one event that mimics IBD signal):
        for file_NC in os.listdir(input_path_NCevents):
            if file_NC.startswith("file") and file_NC.endswith(".txt"):

                # increment num_event_analyzed_NC:
                num_events_analyzed_NC += 1

                # get the file name:
                file_name_NC = input_path_NCevents + file_NC

                # read txt file:
                file_data_NC = np.loadtxt(file_name_NC)
                # 0th entry in file_data_NC is minimum of time window in ns:
                min_time_NC = file_data_NC[0]
                # 1st entry in file_data_NC is maximum of time window in ns:
                max_time_NC = file_data_NC[1]
                # 2nd entry in file_data_NC is bin-width in ns:
                bin_width = file_data_NC[2]

                # the rest of file_data_NC is the hittime distribution histogram in nPE per bin:
                number_pe_per_bin_NC = file_data_NC[3:]

                # check if max_time_NC is greater than end_time:
                if max_time_NC > end_time:
                    # prompt signal is longer than time window. -> Set max_time_NC = end_time:
                    print("max_time_NC {0:.2f} ns > end_time {1:.1f} ns in file {2}".format(max_time_NC, end_time,
                                                                                            file_name_NC))
                    max_time_NC = end_time

                # time window corresponding to number_pe_per_bin_NC:
                time_window_NC = np.arange(min_time_NC, end_time + bin_width, bin_width)

                # compare len(time_window_NC) with len(number_pe_per_bin_NC):
                missing_zeros = len(time_window_NC) - len(number_pe_per_bin_NC)

                # append the missing_zeros to number_pe_per_bin_positron:
                number_pe_per_bin_NC = np.pad(number_pe_per_bin_NC, (0, missing_zeros), 'constant',
                                              constant_values=(0.0, 0.0))

                # analyze the hittime distribution of these event:
                tot_ratio_NC, npe_norm_NC = pulse_shape(time_window_NC, number_pe_per_bin_NC, start_tail[index],
                                                        stop_tail[index1])

                # append tail-to-total ratio to array:
                array_tot_ratio_NC.append(tot_ratio_NC)

                # append zeros to npe_norm_NC to get a average length of the hittimes:
                npe_norm_NC = np.pad(npe_norm_NC, (0, length_average_hittime - len(npe_norm_NC)), 'constant',
                                     constant_values=(0.0, 0.0))

                # add the normalized hittime distribution (npe_norm_NC) to the average hittime distribution
                # (hittime_average_NC):
                hittime_average_NC = hittime_average_NC + npe_norm_NC

            else:
                continue

        # to get the average hittime distribution with a maximum of 1, normalize hittime_average_NC with
        # max(hittime_average_NC):
        hittime_average_NC = hittime_average_NC / max(hittime_average_NC)

        """ calculate the efficiency of the tail-to-total pulse shape analysis (for the values of start_tail[index] and 
        stop_tail[index1]): """
        if max(array_tot_ratio_positron) >= max(array_tot_ratio_NC):
            maximum_tot_value = max(array_tot_ratio_positron)
        else:
            maximum_tot_value = max(array_tot_ratio_NC)

        # check the efficiency of PSD for different cut-efficiencies of NC events:
        efficiency_NC_99 = 99.0
        efficiency_NC_98 = 98.0
        efficiency_NC_97 = 97.0
        efficiency_NC_96 = 96.0
        efficiency_NC_95 = 95.0
        efficiency_NC_90 = 90.0

        # calculate the cut-efficiencies for positrons depending of efficiency_NC:
        efficiency_pos_99, tot_cut_value_99 = tot_efficiency(array_tot_ratio_positron, array_tot_ratio_NC,
                                                             efficiency_NC_99)
        efficiency_pos_98, tot_cut_value_98 = tot_efficiency(array_tot_ratio_positron, array_tot_ratio_NC,
                                                             efficiency_NC_98)
        efficiency_pos_97, tot_cut_value_97 = tot_efficiency(array_tot_ratio_positron, array_tot_ratio_NC,
                                                             efficiency_NC_97)
        efficiency_pos_96, tot_cut_value_96 = tot_efficiency(array_tot_ratio_positron, array_tot_ratio_NC,
                                                             efficiency_NC_96)
        efficiency_pos_95, tot_cut_value_95 = tot_efficiency(array_tot_ratio_positron, array_tot_ratio_NC,
                                                             efficiency_NC_95)
        efficiency_pos_90, tot_cut_value_90 = tot_efficiency(array_tot_ratio_positron, array_tot_ratio_NC,
                                                             efficiency_NC_90)

        # check efficiency_pos to get the "best" (in this case smallest) value:
        if efficiency_pos_99 < best_pos_eff_for_NC_99:
            best_pos_eff_for_NC_99 = efficiency_pos_99
            # also store the corresponding tot-value, start and stop time of the tail:
            tot_value_99 = tot_cut_value_99
            tail_start_99 = start_tail[index]
            tail_stop_99 = stop_tail[index1]

        # check efficiency_pos to get the "best" (in this case smallest) value:
        if efficiency_pos_98 < best_pos_eff_for_NC_98:
            best_pos_eff_for_NC_98 = efficiency_pos_98
            # also store the corresponding tot-value, start and stop time of the tail:
            tot_value_98 = tot_cut_value_98
            tail_start_98 = start_tail[index]
            tail_stop_98 = stop_tail[index1]

        # check efficiency_pos to get the "best" (in this case smallest) value:
        if efficiency_pos_97 < best_pos_eff_for_NC_97:
            best_pos_eff_for_NC_97 = efficiency_pos_97
            # also store the corresponding tot-value, start and stop time of the tail:
            tot_value_97 = tot_cut_value_97
            tail_start_97 = start_tail[index]
            tail_stop_97 = stop_tail[index1]

        # check efficiency_pos to get the "best" (in this case smallest) value:
        if efficiency_pos_96 < best_pos_eff_for_NC_96:
            best_pos_eff_for_NC_96 = efficiency_pos_96
            # also store the corresponding tot-value, start and stop time of the tail:
            tot_value_96 = tot_cut_value_96
            tail_start_96 = start_tail[index]
            tail_stop_96 = stop_tail[index1]

        # check efficiency_pos to get the "best" (in this case smallest) value:
        if efficiency_pos_95 < best_pos_eff_for_NC_95:
            best_pos_eff_for_NC_95 = efficiency_pos_95
            # also store the corresponding tot-value, start and stop time of the tail:
            tot_value_95 = tot_cut_value_95
            tail_start_95 = start_tail[index]
            tail_stop_95 = stop_tail[index1]

        # check efficiency_pos to get the "best" (in this case smallest) value:
        if efficiency_pos_90 < best_pos_eff_for_NC_90:
            best_pos_eff_for_NC_90 = efficiency_pos_90
            # also store the corresponding tot-value, start and stop time of the tail:
            tot_value_90 = tot_cut_value_90
            tail_start_90 = start_tail[index]
            tail_stop_90 = stop_tail[index1]

        print("tail start = {0:.1f} ns".format(start_tail[index]))
        print("tail end = {0:.1f} ns".format(stop_tail[index1]))
        print("NC efficiency = {0:.1f} %, positron efficiency = {1:.2f} %, tot-value = {2:.5f}"
              .format(efficiency_NC_90, efficiency_pos_90, tot_cut_value_90))
        # print("$\\epsilon_{pos}$ / $\\epsilon_{NC}$" + " = {0:.4f}\n".format(eff_ratio_90))
        print("NC efficiency = {0:.1f} %, positron efficiency = {1:.2f} %, tot-value = {2:.5f}"
              .format(efficiency_NC_95, efficiency_pos_95, tot_cut_value_95))
        # print("$\\epsilon_{pos}$ / $\\epsilon_{NC}$" + " = {0:.4f}\n".format(eff_ratio_95))
        print("NC efficiency = {0:.1f} %, positron efficiency = {1:.2f} %, tot-value = {2:.5f}"
              .format(efficiency_NC_96, efficiency_pos_96, tot_cut_value_96))
        # print("$\\epsilon_{pos}$ / $\\epsilon_{NC}$" + " = {0:.4f}\n".format(eff_ratio_96))
        print("NC efficiency = {0:.1f} %, positron efficiency = {1:.2f} %, tot-value = {2:.5f}"
              .format(efficiency_NC_97, efficiency_pos_97, tot_cut_value_97))
        # print("$\\epsilon_{pos}$ / $\\epsilon_{NC}$" + " = {0:.4f}\n".format(eff_ratio_97))
        print("NC efficiency = {0:.1f} %, positron efficiency = {1:.2f} %, tot-value = {2:.5f}"
              .format(efficiency_NC_98, efficiency_pos_98, tot_cut_value_98))
        # print("$\\epsilon_{pos}$ / $\\epsilon_{NC}$" + " = {0:.4f}\n".format(eff_ratio_98))
        print("NC efficiency = {0:.1f} %, positron efficiency = {1:.2f} %, tot-value = {2:.5f}"
              .format(efficiency_NC_99, efficiency_pos_99, tot_cut_value_99))
        # print("$\\epsilon_{pos}$ / $\\epsilon_{NC}$" + " = {0:.4f}\n".format(eff_ratio_99))
        #
        # print(efficiency_pos_90)
        # print(efficiency_pos_95)
        # print(efficiency_pos_96)
        # print(efficiency_pos_97)
        # print(efficiency_pos_98)
        # print(efficiency_pos_99)

        if DISPLAY_PLOTS:
            """ Display tot ratios in histograms: """
            h1 = plt.figure(1, figsize=(15, 8))
            First_bin = 0.0
            Last_bin = 0.05
            Bin_width = (Last_bin-First_bin) / 200
            Bins = np.arange(First_bin, Last_bin+Bin_width, Bin_width)

            plt.hist(array_tot_ratio_positron, bins=Bins, histtype="step", align='mid', color="r", linewidth=1.5,
                     label="positrons with kinetic energy between 10 MeV and 100 MeV (entries = {0:d})"
                     .format(num_events_analyzed_positron))

            plt.hist(array_tot_ratio_NC, bins=Bins, histtype="step", align='mid', color="b", linewidth=1.5,
                     label="prompt signal of NC events that mimic IBD signal (entries = {0:d})"
                     .format(num_events_analyzed_NC))

            plt.xlabel("tail-to-total ratio")
            plt.ylabel("events")
            plt.title("Tail-to-total value for prompt signals of positrons and NC events" +
                      "\n(tail window {0:0.1f} ns to {1:0.1f} ns)".format(start_tail[index], stop_tail[index1]))
            plt.legend()
            plt.grid()
            plt.savefig(output_path + "tot_ratio_tail_{0:0.0f}ns_to_{1:0.0f}ns.png".format(start_tail[index],
                                                                                           stop_tail[index1]))
            plt.close()

            """ Display tot ratios in histograms: """
            h2 = plt.figure(2, figsize=(15, 8))
            First_bin = 0.0
            Last_bin = 0.05
            Bin_width = (Last_bin-First_bin) / 200
            Bins = np.arange(First_bin, Last_bin+Bin_width, Bin_width)
            n_pos_1, bins_pos_1, patches_pos_1 = plt.hist(array_tot_ratio_positron, bins=Bins, histtype="step",
                                                          align='mid',
                                                          color="r", linewidth=1.5,
                                                          label="positrons with kinetic energy between 10 MeV and "
                                                                "100 MeV "
                                                                "(entries = {0:d})"
                                                          .format(num_events_analyzed_positron))
            n_nc_1, bins_nc_1, patches_nc_1 = plt.hist(array_tot_ratio_NC, bins=Bins, histtype="step", align='mid',
                                                       color="b",
                                                       linewidth=1.5,
                                                       label="prompt signal of NC events that mimic IBD signal "
                                                             "(entries = {0:d})"
                                                       .format(num_events_analyzed_NC))
            plt.vlines(tot_cut_value_99, 0, max(n_pos_1)+max(n_pos_1)/10, colors="k", linestyles="-",
                       label="$\\epsilon_{NC}$ = "+"{0:0.2f} %\n".format(efficiency_NC_99)+"$\\epsilon_{IBD}$ = " +
                             "{0:0.2f} %\n".format(efficiency_pos_99)+"tot value = {0:.5f}".format(tot_value_99))
            plt.vlines(tot_cut_value_95, 0, max(n_pos_1)+max(n_pos_1)/10, colors="k", linestyles="--",
                       label="$\\epsilon_{NC}$ = "+"{0:0.2f} %\n".format(efficiency_NC_95)+"$\\epsilon_{IBD}$ = " +
                             "{0:0.2f} %\n".format(efficiency_pos_95)+"tot value = {0:.5f}".format(tot_value_95))
            plt.vlines(tot_cut_value_90, 0, max(n_pos_1)+max(n_pos_1)/10, colors="k", linestyles=":",
                       label="$\\epsilon_{NC}$ = "+"{0:0.2f} %\n".format(efficiency_NC_90)+"$\\epsilon_{IBD}$ = " +
                             "{0:0.2f} %\n".format(efficiency_pos_90)+"tot value = {0:.5f}".format(tot_value_90))
            plt.xlabel("tail-to-total ratio")
            plt.ylabel("events")
            plt.title("Tail-to-total value for prompt signals of positrons and NC events" +
                      "\n(tail window {0:0.1f} ns to {1:0.1f} ns)".format(start_tail[index], stop_tail[index1]))
            plt.legend()
            plt.grid()
            plt.savefig(output_path + "tot_ratio_tail_{0:0.0f}ns_to_{1:0.0f}ns_efficiency.png"
                        .format(start_tail[index], stop_tail[index1]))
            plt.close()

            """ Display the average hittime distributions of positrons and IBD-like NC events: """
            h3 = plt.figure(3, figsize=(15, 8))
            bin_edges = np.arange(0.0, bin_width*length_average_hittime, bin_width)
            plt.semilogy(bin_edges, hittime_average_positron, linestyle="steps", color="r",
                         label="average positron hittime distribution")
            plt.semilogy(bin_edges, hittime_average_NC, linestyle="steps", color="b",
                         label="average hittime distribution of IBD-like NC events")
            plt.xlabel("hittime in ns")
            plt.ylabel("probability per bin (bin-width = {0:0.1f} ns)".format(bin_width))
            plt.xlim(xmin=0.0, xmax=end_time)
            plt.ylim(ymin=1e-4, ymax=2.0)
            plt.title("Average hittime distribution of prompt signals")
            plt.legend()
            plt.grid()
            plt.savefig(output_path + "average_hittimes.png")
            plt.close()

""" print the best positron efficiencies with the corresponding start and stop time of the tail for each NC efficiency:
"""
print("NC efficiency = {0:.1f} %, best positron efficiency = {1:.2f} %, tail start = {2:.1f} ns, tail end = {3:.1f} ns"
      ", tail-to-total value = {4:.5f}"
      .format(efficiency_NC_90, best_pos_eff_for_NC_90, tail_start_90, tail_stop_90, tot_value_90))
print("NC efficiency = {0:.1f} %, best positron efficiency = {1:.2f} %, tail start = {2:.1f} ns, tail end = {3:.1f} ns"
      ", tail-to-total value = {4:.5f}"
      .format(efficiency_NC_95, best_pos_eff_for_NC_95, tail_start_95, tail_stop_95, tot_value_95))
print("NC efficiency = {0:.1f} %, best positron efficiency = {1:.2f} %, tail start = {2:.1f} ns, tail end = {3:.1f} ns"
      ", tail-to-total value = {4:.5f}"
      .format(efficiency_NC_96, best_pos_eff_for_NC_96, tail_start_96, tail_stop_96, tot_value_96))
print("NC efficiency = {0:.1f} %, best positron efficiency = {1:.2f} %, tail start = {2:.1f} ns, tail end = {3:.1f} ns"
      ", tail-to-total value = {4:.5f}"
      .format(efficiency_NC_97, best_pos_eff_for_NC_97, tail_start_97, tail_stop_97, tot_value_97))
print("NC efficiency = {0:.1f} %, best positron efficiency = {1:.2f} %, tail start = {2:.1f} ns, tail end = {3:.1f} ns"
      ", tail-to-total value = {4:.5f}"
      .format(efficiency_NC_98, best_pos_eff_for_NC_98, tail_start_98, tail_stop_98, tot_value_98))
print("NC efficiency = {0:.1f} %, best positron efficiency = {1:.2f} %, tail start = {2:.1f} ns, tail end = {3:.1f} ns"
      ", tail-to-total value = {4:.5f}"
      .format(efficiency_NC_99, best_pos_eff_for_NC_99, tail_start_99, tail_stop_99, tot_value_99))

