""" script to do pulse shape analysis of the fast neutron background:

    Neutron hittimes from script hittime_distribution_fastneutron.py, that are saved in folder
    /home/astro/blum/PhD/work/MeVDM_JUNO/fast_neutrons/hittimes/, are analyzed with this script.

    As input for the different values defining the pulse shape analysis (tail start, tail end, tot value and efficiency
    of positron hittimes) the results from script pulse_shape_analysis_v1.py (summarized in file PSD_results.ods) are
    taken.

    Procedure to get the efficiency of how many fast neutron background events can be cut away:
    1.  calculate the tail-to-total values for each neutron hittime distribution with function pulse_shape() also used
        in script pulse_shape_analysis_v1.py for start and stop value of tail given by PSD_results.ods.
    2.  Take the tail-to-total value corresponding to the IBD efficiency from PSD_results.ods and calculate the
        fast neutron efficiency due to this tail-to-total value.

    To compare the results, also the positron hittimes and NC hittimes are analyzed in the same way like in
    pulse_shape_analysis_v1.py

"""
import datetime
import os
import sys
import numpy as np
from matplotlib import pyplot as plt
from NC_background_functions import pulse_shape


def fast_n_efficiency(array_tot_fn, tot_value_pos):
    """
    calculate the number of values in array_tot_fn, that are greater than tot_values_pos
    :param array_tot_fn: array, where the tail-to-total values of fast neutron hittimes are stored
    :param tot_value_pos: tail-to-total value from positron and NC PSD
    :return:
    """
    # number of analyzed hittime distributions:
    number_events = len(array_tot_fn)

    # sort array_tot_fn in ascending order (from small to large):
    array_tot_fn.sort()

    # loop over array_tot_fn until you reach tot_value_pos:
    for index9 in range(number_events):
        if array_tot_fn[index9] >= tot_value_pos:
            # you reach tot_value_pos -> break from the loop
            break

    # index9 events are smaller than tot_value_pos -> How many events are greater than tot_value_pos?
    number_cut_away = float(number_events - index9)

    # calculate the fast neutron efficiency (How many events are cut away) in percent:
    fn_eff = number_cut_away / float(number_events) * 100.0

    return fn_eff


# get the date and time, when the script was run:
date = datetime.datetime.now()
now = date.strftime("%Y-%m-%d %H:%M")

# path, where output is saved:
output_path = "/home/astro/blum/PhD/work/MeVDM_JUNO/fast_neutrons/"

# Set flag, if plots should be saved:
DISPLAY_PLOTS = True

""" parameters for tail to total method from pulse shape analysis of IBD events and NC events (PSD_results.ods): """
# INFO-me: parameters should agree with the bin-width of the time window!
# start of the tail in ns:
start_tail = np.array([335, 350, 340])

# end of the tail in ns:
stop_tail = np.array([600, 540, 540])

# tail-to-total value corresponding to tail window:
tot_value_positron = np.array([0.00661, 0.00497, 0.00525])

# best positron (IBD) efficiencies in %:
best_positron_efficiencies = np.array([2.3, 5.77, 14.16])

# corresponding NC (IBD-like) efficiencies in %:
NC_efficiencies = np.array([95, 96, 97])

# check if array have same length:
if (len(start_tail) != len(stop_tail) != len(tot_value_positron) != len(best_positron_efficiencies)
        != len(NC_efficiencies)):
    sys.exit("ERROR: input parameters have not the same length!!!")

print("start_tail = {0}".format(start_tail))
print("stop_tail = {0}".format(stop_tail))

""" parameters that define the time window of prompt signal: """
# start of the time window in ns:
start_time = 0.0
# end of the time window in ns:
end_time = 2000.0

# loop over different start values of the tail:
for index in range(len(start_tail)):

    """ analyze the hittime distribution of neutrons: """
    print("analyze neutrons...")

    # path, where hittime distributions neutrons are saved:
    input_path_neutron = "/home/astro/blum/PhD/work/MeVDM_JUNO/fast_neutrons/hittimes/"

    # number of events that are analyzed:
    num_events_analyzed_neutron = 0

    # preallocate array, where tail-to-total ratios are stored:
    array_tot_ratio_neutron = []

    # preallocate array, where a average hittime-distribution of neutrons are stored (number of pe per bin):
    # length of average hittime distribution (number of bins):
    length_average_hittime = 500
    hittime_average_neutron = np.zeros(length_average_hittime)

    # loop over all files in folder input_path_neutron, that start with 'file' and end with 'neutron.txt'
    # (files where hittime distribution is saved, each file is equal to one event):
    for file_neutron in os.listdir(input_path_neutron):
        if file_neutron.startswith("file") and file_neutron.endswith("neutron.txt"):

            # get the file name:
            file_name_neutron = input_path_neutron + file_neutron

            # read txt file:
            file_data_neutron = np.loadtxt(file_name_neutron)
            # 0th entry in file_data_neutron is minimum of time window in ns:
            min_time_neutron = file_data_neutron[0]
            # 1st entry in file_data_neutron is maximum of time window in ns:
            max_time_neutron = file_data_neutron[1]
            # 2nd entry in file_data_neutron is bin-width in ns:
            bin_width = file_data_neutron[2]

            # the rest of file_data_neutron is the hittime distribution histogram in nPE per bin:
            number_pe_per_bin_neutron = file_data_neutron[3:]

            # check if max_time_neutron is greater than end_time:
            if max_time_neutron > end_time:
                # prompt signal is longer than time window.
                print("max_time_neutron {0:.2f} ns > end_time {1:.1f} ns in file {2}".format(max_time_neutron,
                                                                                             end_time,
                                                                                             file_name_neutron))

            # time window:
            time_window_neutron = np.arange(min_time_neutron, end_time + bin_width, bin_width)

            # compare len(time_window_neutron) with len(number_pe_per_bin_neutron):
            missing_zeros = len(time_window_neutron) - len(number_pe_per_bin_neutron)

            # append the missing_zeros to number_pe_per_bin_neutron:
            number_pe_per_bin_neutron = np.pad(number_pe_per_bin_neutron, (0, missing_zeros), 'constant',
                                               constant_values=(0.0, 0.0))

            # analyze the hittime distribution of these event:
            tot_ratio_neutron, npe_norm_neutron = pulse_shape(time_window_neutron, number_pe_per_bin_neutron,
                                                              start_tail[index], stop_tail[index])

            # check if tot-value is not 0:
            if tot_ratio_neutron == 0:
                continue

            # increment number of analyzed events:
            num_events_analyzed_neutron += 1

            # append tail-to-total ratio to array:
            array_tot_ratio_neutron.append(tot_ratio_neutron)

            # append zeros to npe_norm_neutron to get a average length of the hittimes:
            npe_norm_neutron = np.pad(npe_norm_neutron, (0, length_average_hittime - len(npe_norm_neutron)),
                                      'constant', constant_values=(0.0, 0.0))

            # add the normalized hittime distribution (npe_norm_neutron) to the average hittime distribution
            # (hittime_average_neutron):
            hittime_average_neutron = hittime_average_neutron + npe_norm_neutron

        else:
            continue

    # array_tot_ratio_neutron contains the tot-values of each neutron hittime distribution!

    # to get the average hittime distribution with a maximum of 1, normalize hittime_average_neutron with
    # max(hittime_average_neutron):
    hittime_average_neutron = hittime_average_neutron / max(hittime_average_neutron)

    """ calculate the fast neutron efficiency due to the tot-value from positron and NC pulse shape analysis: """
    fast_n_eff = fast_n_efficiency(array_tot_ratio_neutron, tot_value_positron[index])

    print("\ntail start = {0:.1f} ns, tail end = {1:.1f} ns, tot value = {2:.5f}"
          .format(start_tail[index], stop_tail[index], tot_value_positron[index]))
    print("positron (IBD) efficiency = {0:.3f} %".format(best_positron_efficiencies[index]))
    print("NC efficiency = {0:.3f} %".format(NC_efficiencies[index]))
    print("Fast Neutron efficiency = {0:.3f} %\n".format(fast_n_eff))

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
                                                                start_tail[index], stop_tail[index])

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
                                                    stop_tail[index])

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

    if DISPLAY_PLOTS:
        # display tot-values for positrons, NC and fast neutron events for the given configuration:
        h1 = plt.figure(1, figsize=(15, 8))
        First_bin = 0.0
        Last_bin = 0.05
        Bin_width = (Last_bin-First_bin) / 200
        Bins = np.arange(First_bin, Last_bin+Bin_width, Bin_width)

        plt.hist(array_tot_ratio_positron, bins=Bins, histtype="step", align='mid', color="r",
                 label="positrons with kinetic energy between 10 MeV and 100 MeV (entries = {0:d})"
                 .format(num_events_analyzed_positron))

        plt.hist(array_tot_ratio_NC, bins=Bins, histtype="step", align='mid', color="b",
                 label="prompt signal of NC events that mimic IBD signal (entries = {0:d})"
                 .format(num_events_analyzed_NC))

        plt.hist(array_tot_ratio_neutron, bins=Bins, histtype="step", align='mid', color="g",
                 label="prompt signal of neutrons representing fast neutron events (entries = {0:d})"
                 .format(num_events_analyzed_neutron))

        plt.xlabel("tail-to-total ratio")
        plt.ylabel("events")
        plt.title("Tail-to-total value for prompt signals of positron, NC and fast neutron events" +
                  "\n(tail window {0:0.1f} ns to {1:0.1f} ns)".format(start_tail[index], stop_tail[index]))
        plt.legend()
        plt.grid()
        plt.savefig(output_path + "tot_ratio_tail_{0:.0f}_PosNCfastN.png".format(NC_efficiencies[index]))
        plt.close()

        # display tot-values for positrons, NC and fast neutron events for the given configuration with efficiencies:
        h2 = plt.figure(2, figsize=(15, 8))
        First_bin = 0.0
        Last_bin = 0.05
        Bin_width = (Last_bin-First_bin) / 200
        Bins = np.arange(First_bin, Last_bin+Bin_width, Bin_width)

        n_pos_1, bins_pos_1, patches_pos_1 = plt.hist(array_tot_ratio_positron, bins=Bins, histtype="step", align='mid',
                                                      color="r", linewidth=1.5,
                                                      label="positrons with kinetic energy between 10 MeV and 100 MeV "
                                                            "(entries = {0:d})"
                                                      .format(num_events_analyzed_positron))

        n_NC_1, bins_NC_1, patches_NC_1 = plt.hist(array_tot_ratio_NC, bins=Bins, histtype="step", align='mid',
                                                   color="b", linewidth=1.5,
                                                   label="prompt signal of NC events that mimic IBD signal "
                                                         "(entries = {0:d})"
                                                   .format(num_events_analyzed_NC))

        n_n_1, bins_n_1, patches_n_1 = plt.hist(array_tot_ratio_neutron, bins=Bins, histtype="step", align='mid',
                                                color="g", linewidth=1.5,
                                                label="prompt signal of neutrons representing fast neutron "
                                                      "events (entries = {0:d})"
                                                .format(num_events_analyzed_neutron))

        plt.vlines(tot_value_positron[index], 0, max(n_pos_1)+max(n_pos_1)/10, colors="k", linestyles="--",
                       label="$\\epsilon_{IBD}$ = "+"{0:0.2f} %\n".format(best_positron_efficiencies[index])+
                             "$\\epsilon_{NC}$ = "+"{0:0.2f} %\n".format(NC_efficiencies[index])+
                             "$\\epsilon_{fastN}$ = "+"{0:0.2f} %\n".format(fast_n_eff)+
                             "tot value = {0:.5f}".format(tot_value_positron[index]))

        plt.xlabel("tail-to-total ratio")
        plt.ylabel("events")
        plt.title("Tail-to-total value for prompt signals of positron, NC and fast neutron events" +
                  "\n(tail window {0:0.1f} ns to {1:0.1f} ns)".format(start_tail[index], stop_tail[index]))
        plt.legend()
        plt.grid()
        plt.savefig(output_path + "tot_ratio_tail_{0:.0f}_PosNCfastN_efficiencies.png".format(NC_efficiencies[index]))
        plt.close()

        """ Display the average hittime distributions of positrons and IBD-like NC events: """
        h3 = plt.figure(3, figsize=(15, 8))
        bin_edges = np.arange(0.0, bin_width*length_average_hittime, bin_width)
        plt.semilogy(bin_edges, hittime_average_positron, linestyle="steps", color="r",
                     label="average positron hittime distribution")
        plt.semilogy(bin_edges, hittime_average_NC, linestyle="steps", color="b",
                     label="average hittime distribution of IBD-like NC events")
        plt.semilogy(bin_edges, hittime_average_neutron, linestyle="steps", color="g",
                     label="average hittime distribution of fast neutron events")
        plt.xlabel("hittime in ns")
        plt.ylabel("probability per bin (bin-width = {0:0.1f} ns)".format(bin_width))
        plt.xlim(xmin=0.0, xmax=end_time)
        plt.ylim(ymin=1e-4, ymax=2.0)
        plt.title("Average hittime distribution of prompt signals")
        plt.legend()
        plt.grid()
        plt.savefig(output_path + "average_hittimes_PosNCfastN.png")
        plt.close()











