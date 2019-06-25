""" Script to do a pulse shape analysis of atmospheric NC neutrino events that could mimic an IBD signal in JUNO
    detector.
    The hittime distribution of the prompt signal of the NC event, which passes all the cuts (and therefore mimics an
    IBD signal), are save in folder output_detsim (the hittime distributions are analyzed with
    prompt_signal_preselected_evts.py).
    As reference the hittime distribution of prompt positron events are also analyzed for kinetic energy of positrons
    of 10 MeV and 100 MeV.

    The hittime distributions of positrons are equal to the hittime distribution of the prompt signal of real IBD events
    and therefore act as reference.

"""
import datetime
import os
import numpy as np
from matplotlib import pyplot as plt


def pulse_shape(hittime, npe, tail_start, tail_end):
    """
    function to analyzed the hittime distribution of one event with the tail to total method (charge integration method)

    1. time of flight correction: shift the distribution that is starts at 0 ns
    2. normalize the hittime distribution to 1
    3. calculate the 'charge' of the whole distribution
    4. calculate the 'charge' of the tail of the distribution (define by tail_start and tail_end)
    5. calculate the ration between charge of tail and charge of total distribution


    :param hittime: hittime array, which defines the time window (array in ns)
    :param npe: number of pe per bin (corresponds to hittime array) (array)
    :param tail_start: defines the start value of the tail in ns
    :param tail_end: defines the stop value of the tail in ns
    :return:
    """
    """ time of flight correction: """
    # get maximum value of npe:
    maximum_npe = np.max(npe)
    # calculate 10 % of maximum_pe:
    start_condition = 0.1 * maximum_npe

    # loop over npe until one entry is higher than start_condition:
    for index in range(len(npe)):
        if npe[index] >= start_condition:
            break

    # do time of flight correction for hittime and npe:
    hittime = hittime[index:]
    npe = npe[index:]

    """ normalize hittime distribution to 1: """
    # calculate the integral of the whole hittime distribution:
    integral_npe = np.trapz(npe, hittime)
    # normalize npe distribution to 1:
    npe_norm = npe/integral_npe

    """ integral (charge) of total distribution: """
    # should be 1 because of normalization:
    integral_total = np.trapz(npe_norm, hittime)

    """ integral (charge) of the tail of the distribution: """
    # get the index of hittime, which correspond to tail_start:
    for index1 in range(len(hittime)):
        if hittime[index1] == tail_start:
            index_tail_start = index1
        elif hittime[index1] == tail_end:
            index_tail_end = index1
        else:
            continue

    # define the time window of the tail of the distribution:
    hittime_tail = hittime[index_tail_start:index_tail_end+1]
    # get the corresponding npe_norm array:
    npe_tail = npe_norm[index_tail_start:index_tail_end+1]

    # integrate the tail of the distribution:
    integral_tail = np.trapz(npe_tail, hittime_tail)

    """ tail to total ratio: """
    # calculate the ratio between integral od tail of distribution and of total distribution:
    tot_ratio = integral_tail / integral_total

    return tot_ratio, npe_norm


def tot_efficiency(tot_values_positron, tot_values_nc, eff_nc):
    """
    function to get the efficiency, how many positron events are cut away, depending on eff_nc
    :param tot_values_positron: list/array of all tail-to-total values of positron events
    :param tot_values_nc: list/array of all tail-to-total values of NC events
    :param eff_nc: efficiency, how many NC events should be cut away (in percent). For example: eff_nc=99% -> 99% of all
    NC events are cut away
    :return:
    """
    # calculate the percentile of tot_values_nc defined by (100-eff_nc) (value of tot, where (100-eff_nc) % of all
    # tot-values are smaller)
    tot_cut_value = np.percentile(tot_values_nc, 100 - eff_nc)

    # use tot_cut_value to calculate, how many positron events are cut away:
    num_positron_cut = 0
    # loop over tot_values_positron:
    for index in range(len(tot_values_positron)):
        if tot_values_positron[index] >= tot_cut_value:
            # positron events is cut away by PSD:
            num_positron_cut += 1
        else:
            continue

    # get total number of positron events:
    num_pos_total = len(tot_values_positron)

    # calculate efficiency, how many positron events are cut away by PSD based on eff_nc (in percent):
    eff_positron = float(num_positron_cut) / float(num_pos_total) * 100

    return eff_positron, tot_cut_value


# get the date and time, when the script was run:
date = datetime.datetime.now()
now = date.strftime("%Y-%m-%d %H:%M")

# path, where output is saved:
output_path = "/home/astro/blum/juno/atmoNC/data_NC/output_PSD/"

""" parameters for tail to total method: """
# INFO-me: parameters should agree with the bin-width of the time window!
# start of the tail in ns:
start_tail = 450.0
# end of the tail in ns:
stop_tail = 500.0

""" analyze the hittime distribution of the 10 MeV positrons: """
print("analyze 10 MeV positrons...")
# path, where hittime distributions of 10 MeV positrons are saved:
input_path_positron10 = "/home/astro/blum/juno/atmoNC/data_NC/output_PSD/positron_hittime/"

# kinetic energy of positron:
kinetic_energy_10 = 10

# number of events that are analyzed:
num_events_analyzed_10 = 0

# preallocate array, where tail-to-total ratios are stored:
array_tot_ratio_10 = []

# loop over all files in folder input_path_positron10, that start with 'file' and end with '10_MeV.txt'
# (files where hittime distribution is saved, each file is equal to one event):
for file_10 in os.listdir(input_path_positron10):
    if file_10.startswith("file") and file_10.endswith("10_MeV.txt"):

        # increment num_event_analyzed_10:
        num_events_analyzed_10 += 1

        # get the file name:
        file_name_10 = input_path_positron10 + file_10

        # read txt file:
        file_data_10 = np.loadtxt(file_name_10)
        # 0th entry in file_data_10 is minimum of time window in ns:
        min_time_10 = file_data_10[0]
        # 1st entry in file_data_10 is maximum of time window in ns:
        max_time_10 = file_data_10[1]
        # 2nd entry in file_data_10 is bin-width in ns:
        bin_width = file_data_10[2]

        # the rest of file_data_10 is the hittime distribution histogram in nPE per bin:
        number_pe_per_bin_10 = file_data_10[3:]

        # time window corresponding to number_pe_per_bin_10:
        time_window_10 = np.arange(min_time_10, max_time_10 + bin_width, bin_width)

        # analyze the hittime distribution of these event:
        tot_ratio_10, npe_norm_10 = pulse_shape(time_window_10, number_pe_per_bin_10, start_tail, stop_tail)

        # append tail-to-total ratio to array:
        array_tot_ratio_10.append(tot_ratio_10)

    else:
        continue

""" analyze the hittime distribution of the 100 MeV positrons: """
print("analyze 100 MeV positrons...")
# path, where hittime distributions of 100 MeV positrons are saved:
input_path_positron100 = "/home/astro/blum/juno/atmoNC/data_NC/output_PSD/positron_hittime/"

# kinetic energy of positron:
kinetic_energy_100 = 100

# number of events that are analyzed:
num_events_analyzed_100 = 0

# preallocate array, where tail-to-total ratios are stored:
array_tot_ratio_100 = []

# loop over all files in folder input_path_positron100, that start with 'file' and end with '100_MeV.txt'
# (files where hittime distribution is saved, each file is equal to one event):
for file_100 in os.listdir(input_path_positron100):
    if file_100.startswith("file") and file_100.endswith("100_MeV.txt"):

        # increment num_event_analyzed_100:
        num_events_analyzed_100 += 1

        # get the file name:
        file_name_100 = input_path_positron100 + file_100

        # read txt file:
        file_data_100 = np.loadtxt(file_name_100)
        # 0th entry in file_data_100 is minimum of time window in ns:
        min_time_100 = file_data_100[0]
        # 1st entry in file_data_100 is maximum of time window in ns:
        max_time_100 = file_data_100[1]
        # 2nd entry in file_data_100 is bin-width in ns:
        bin_width = file_data_100[2]

        # the rest of file_data_100 is the hittime distribution histogram in nPE per bin:
        number_pe_per_bin_100 = file_data_100[3:]

        # time window corresponding to number_pe_per_bin_100:
        time_window_100 = np.arange(min_time_100, max_time_100 + bin_width, bin_width)

        # analyze the hittime distribution of these event:
        tot_ratio_100, npe_norm_100 = pulse_shape(time_window_100, number_pe_per_bin_100, start_tail, stop_tail)

        # append tail-to-total ratio to array:
        array_tot_ratio_100.append(tot_ratio_100)

    else:
        continue

""" analyze the hittime distribution of positrons with kinetic energy uniformly distributed from 10 MeV to 100 MeV: """
print("analyze positrons...")
# path, where hittime distributions of 100 MeV positrons are saved:
input_path_positron = "/home/astro/blum/juno/atmoNC/data_NC/output_PSD/positron_hittime/"

# number of events that are analyzed:
num_events_analyzed_positron = 0

# preallocate array, where tail-to-total ratios are stored:
array_tot_ratio_positron = []

# preallocate array, where a average hittime-distribution of positrons are stored (number of pe per bin):
# length of average hittime distribution (number of bins):
length_average_hittime = 300
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

        # time window corresponding to number_pe_per_bin_positron:
        time_window_positron = np.arange(min_time_positron, max_time_positron + bin_width, bin_width)

        # analyze the hittime distribution of these event:
        tot_ratio_positron, npe_norm_positron = pulse_shape(time_window_positron, number_pe_per_bin_positron,
                                                            start_tail, stop_tail)

        # append tail-to-total ratio to array:
        array_tot_ratio_positron.append(tot_ratio_positron)

        # append zeros to the end of npe_norm_positron to get the same length of the array like
        # hittime_average_positron:
        npe_norm_positron = np.pad(npe_norm_positron, (0, length_average_hittime-len(npe_norm_positron)), 'constant',
                                   constant_values=(0.0, 0.0))

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

        # time window corresponding to number_pe_per_bin_NC:
        time_window_NC = np.arange(min_time_NC, max_time_NC + bin_width, bin_width)

        # analyze the hittime distribution of these event:
        tot_ratio_NC, npe_norm_NC = pulse_shape(time_window_NC, number_pe_per_bin_NC, start_tail, stop_tail)

        # append tail-to-total ratio to array:
        array_tot_ratio_NC.append(tot_ratio_NC)

        # append zeros to the end of npe_norm_NC to get the same length of the array like hittime_average_NC:
        npe_norm_NC = np.pad(npe_norm_NC, (0, length_average_hittime-len(npe_norm_NC)), 'constant',
                             constant_values=(0.0, 0.0))

        # add the normalized hittime distribution (npe_norm_NC) to the average hittime distribution
        # (hittime_average_NC):
        hittime_average_NC = hittime_average_NC + npe_norm_NC

    else:
        continue

# to get the average hittime distribution with a maximum of 1, normalize hittime_average_NC with
# max(hittime_average_NC):
hittime_average_NC = hittime_average_NC / max(hittime_average_NC)

""" calculate the efficiency of the tail-to-total pulse shape analysis (for the values of start_tail and stop_tail): """
if max(array_tot_ratio_positron) >= max(array_tot_ratio_NC):
    maximum_tot_value = max(array_tot_ratio_positron)
else:
    maximum_tot_value = max(array_tot_ratio_NC)

# check the efficiency of PSD for different cut-efficiencies of NC events:
efficiency_NC_99 = 99
efficiency_NC_95 = 95
efficiency_NC_90 = 90

# calculate the cut-efficiencies for positrons depending of efficiency_NC:
efficiency_pos_99, tot_cut_value_99 = tot_efficiency(array_tot_ratio_positron, array_tot_ratio_NC, efficiency_NC_99)
efficiency_pos_95, tot_cut_value_95 = tot_efficiency(array_tot_ratio_positron, array_tot_ratio_NC, efficiency_NC_95)
efficiency_pos_90, tot_cut_value_90 = tot_efficiency(array_tot_ratio_positron, array_tot_ratio_NC, efficiency_NC_90)

print(efficiency_pos_99)
print(tot_cut_value_99)
print(efficiency_pos_95)
print(tot_cut_value_95)
print(efficiency_pos_90)
print(tot_cut_value_90)


""" Display tot ratios in histograms: """
h1 = plt.figure(1, figsize=(15, 8))
First_bin = 0.0
Last_bin = maximum_tot_value
Bin_width = (Last_bin-First_bin) / 200
Bins = np.arange(First_bin, Last_bin+Bin_width, Bin_width)
# plt.hist(array_tot_ratio_10, bins=Bins, histtype="step", align='mid', color="r", label="positrons with {1:0.0f} MeV "
#                                                                                        "kinetic energy "
#                                                                                        "(entries = {0:d})"
#          .format(num_events_analyzed_10, kinetic_energy_10))

# plt.hist(array_tot_ratio_100, bins=Bins, histtype="step", align='mid', color="g", label="positrons with {1:0.0f} "
#                                                                                              "MeV kinetic energy "
#                                                                                              "(entries = {0:d})"
#          .format(num_events_analyzed_100, kinetic_energy_100))

plt.hist(array_tot_ratio_positron, bins=Bins, histtype="step", align='mid', color="r",
         label="positrons with kinetic energy between 10 MeV and 100 MeV (entries = {0:d})"
         .format(num_events_analyzed_positron))

plt.hist(array_tot_ratio_NC, bins=Bins, histtype="step", align='mid', color="b",
         label="prompt signal of NC events that mimic IBD signal (entries = {0:d})"
         .format(num_events_analyzed_NC))

plt.xlabel("tail-to-total ratio")
plt.ylabel("events")
plt.title("Tail-to-total value for prompt signals of positrons and NC events" +
          "\n(tail window {0:0.1f} ns to {1:0.1f} ns)".format(start_tail, stop_tail))
# plt.title("Tail-to-total value for positrons of different kinetic energy"+
#           "\n(tail window {0:0.1f} ns to {1:0.1f} ns)".format(start_tail, stop_tail))
plt.legend()
plt.grid()
plt.savefig(output_path + "tot_ratio_tail_{0:0.0f}ns_to_{1:0.0f}ns.png".format(start_tail, stop_tail))
plt.close()

""" Display tot ratios in histograms: """
h2 = plt.figure(2, figsize=(15, 8))
First_bin = 0.0
Last_bin = maximum_tot_value
Bin_width = (Last_bin-First_bin) / 200
Bins = np.arange(First_bin, Last_bin+Bin_width, Bin_width)
n_pos_1, bins_pos_1, patches_pos_1 = plt.hist(array_tot_ratio_positron, bins=Bins, histtype="step", align='mid',
                                              color="r", linewidth=1.5,
                                              label="positrons with kinetic energy between 10 MeV and 100 MeV "
                                                    "(entries = {0:d})"
                                              .format(num_events_analyzed_positron))
n_nc_1, bins_nc_1, patches_nc_1 = plt.hist(array_tot_ratio_NC, bins=Bins, histtype="step", align='mid', color="b",
                                           linewidth=1.5,
                                           label="prompt signal of NC events that mimic IBD signal (entries = {0:d})"
                                           .format(num_events_analyzed_NC))
plt.vlines(tot_cut_value_99, 0, max(n_pos_1)+max(n_pos_1)/10, colors="k", linestyles="-",
           label="$\\epsilon_{NC}$ = "+"{0:0.2f} %\n".format(efficiency_NC_99)+"$\\epsilon_{IBD}$ = "+"{0:0.2f} %"
           .format(efficiency_pos_99))
plt.vlines(tot_cut_value_95, 0, max(n_pos_1)+max(n_pos_1)/10, colors="k", linestyles="--",
           label="$\\epsilon_{NC}$ = "+"{0:0.2f} %\n".format(efficiency_NC_95)+"$\\epsilon_{IBD}$ = "+"{0:0.2f} %"
           .format(efficiency_pos_95))
plt.vlines(tot_cut_value_90, 0, max(n_pos_1)+max(n_pos_1)/10, colors="k", linestyles=":",
           label="$\\epsilon_{NC}$ = "+"{0:0.2f} %\n".format(efficiency_NC_90)+"$\\epsilon_{IBD}$ = "+"{0:0.2f} %"
           .format(efficiency_pos_90))
plt.xlabel("tail-to-total ratio")
plt.ylabel("events")
plt.title("Tail-to-total value for prompt signals of positrons and NC events" +
          "\n(tail window {0:0.1f} ns to {1:0.1f} ns)".format(start_tail, stop_tail))
# plt.title("Tail-to-total value for positrons of different kinetic energy"+
#           "\n(tail window {0:0.1f} ns to {1:0.1f} ns)".format(start_tail, stop_tail))
plt.legend()
plt.grid()
plt.savefig(output_path + "tot_ratio_tail_{0:0.0f}ns_to_{1:0.0f}ns_efficiency.png".format(start_tail, stop_tail))
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
plt.xlim(xmin=0.0, xmax=500.0)
plt.ylim(ymin=1e-3, ymax=2.0)
plt.title("Average hittime distribution of prompt signals")
plt.legend()
plt.grid()
plt.savefig(output_path + "average_hittimes.png")
plt.close()
