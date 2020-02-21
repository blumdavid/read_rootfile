""" Script to do a pulse shape analysis of atmospheric NC neutrino events and of real IBD signals.

    Read already calculated pulse shapes of prompt signals (hittime distribution) of NC and IBD events
    (calculated with analyze_prompt_delayed_cut_v2.py).

    Calculate the TTR (tail-to-total ratio) for a specific tail window (tail start and tail end) for each pulse shape.

    Create histogram of all TTR values of NC and IBD events.

    Store the filenumber and evtID of each event, that passes the PSD cut which is specified by the TTR cut value.

    Also save the number of events that pass and the number of events that are rejected to txt file.

"""
import datetime
import numpy as np
import os
import sys
from matplotlib import pyplot as plt
from NC_background_functions import pulse_shape
from NC_background_functions import tot_efficiency

# get the date and time, when the script was run:
date = datetime.datetime.now()
now = date.strftime("%Y-%m-%d %H:%M")

# flag if plots should be create and saved:
CREATE_PLOTS = True
# flag if arrays with TTR values should be saved:
SAVE_ARRAY_TTR = True

""" parameters for tail to total method: """
# INFO-me: parameters should agree with the bin-width of the time window!
# start of the tail in ns:
# start_tail = np.arange(50, 650, 100)
start_tail = np.array([225, 275, 325])

# end of the tail in ns:
stop_tail = np.arange(600, 1200, 200)
# stop_tail = np.array([1000])

print("start_tail = {0}".format(start_tail))
print("stop_tail = {0}".format(stop_tail))

""" parameters: """
# integer to split the analysis:
analysis_part = 0
# start of the time window in ns:
start_time = 0.0
# end of the time window in ns:
end_time = 1000.0
# bin-width of the time window in ns:
binwidth = 5.0
# first file of NC events to be read:
start_filenumber_NC = 0
# last file of NC events to be read:
stop_filenumber_NC = 999
# first file of IBD events to be read:
start_filenumber_IBD = 0
# last file of IBD events to be read:
stop_filenumber_IBD = 199
# number of event per file:
num_evts_per_file = 100

# path, where output is saved:
output_path = "/home/astro/blum/juno/atmoNC/data_NC/output_PSD_v2/DCR/"

""" preallocate values: """
# best IBD efficiency for NC efficiency of 90 % and the corresponding tail-to-total value and values of start and
# stop of tail:
best_IBD_eff_for_NC_90 = 100.0
ttr_value_90 = 0.0
tail_start_90 = 0.0
tail_stop_90 = 0.0
best_array_TTR_NC_90 = []
best_array_TTR_IBD_90 = []
# best IBD efficiency for NC efficiency of 95 % and the corresponding tail-to-total value and values of start and
# stop of tail:
best_IBD_eff_for_NC_95 = 100.0
ttr_value_95 = 0.0
tail_start_95 = 0.0
tail_stop_95 = 0.0
best_array_TTR_NC_95 = []
best_array_TTR_IBD_95 = []
# best IBD efficiency for NC efficiency of 96 % and the corresponding tail-to-total value and values of start and
# stop of tail:
best_IBD_eff_for_NC_96 = 100.0
ttr_value_96 = 0.0
tail_start_96 = 0.0
tail_stop_96 = 0.0
best_array_TTR_NC_96 = []
best_array_TTR_IBD_96 = []
# best IBD efficiency for NC efficiency of 97 % and the corresponding tail-to-total value and values of start and
# stop of tail:
best_IBD_eff_for_NC_97 = 100.0
ttr_value_97 = 0.0
tail_start_97 = 0.0
tail_stop_97 = 0.0
best_array_TTR_NC_97 = []
best_array_TTR_IBD_97 = []
# best IBD efficiency for NC efficiency of 98 % and the corresponding tail-to-total value and values of start and
# stop of tail:
best_IBD_eff_for_NC_98 = 100.0
ttr_value_98 = 0.0
tail_start_98 = 0.0
tail_stop_98 = 0.0
best_array_TTR_NC_98 = []
best_array_TTR_IBD_98 = []
# best IBD efficiency for NC efficiency of 99 % and the corresponding tail-to-total value and values of start and
# stop of tail:
best_IBD_eff_for_NC_99 = 100.0
ttr_value_99 = 0.0
tail_start_99 = 0.0
tail_stop_99 = 0.0
best_array_TTR_NC_99 = []
best_array_TTR_IBD_99 = []

# loop over different start values:
for index in range(len(start_tail)):

    # loop over different stop values of the tail:
    for index1 in range(len(stop_tail)):

        """ analyze the pulse shape of prompt signal of real IBD events: """
        print("analyze IBD events...")

        # path, where IBD pulse shapes are stored:
        input_path_IBD = "/home/astro/blum/juno/IBD_events/hittimes/"

        # number of events that are analyzed:
        num_events_analyzed_IBD = 0

        # preallocate array, where tail-to-total ratios are stored:
        array_TTR_IBD = []
        filenumber_IBD = []
        evtID_IBD = []

        # preallocate array, where a average hittime-distribution of IBD events are stored (number of pe per bin):
        # length of average hittime distribution (number of bins):
        length_average_hittime = int((end_time - start_time + binwidth) / binwidth) + 10
        hittime_average_IBD = np.zeros(length_average_hittime)

        # loop over all IBD pulse shapes:
        for filenumber in range(start_filenumber_IBD, stop_filenumber_IBD + 1, 1):
            # print(filenumber)
            for event in range(num_evts_per_file):

                # increment num_events_analyzed_IBD:
                num_events_analyzed_IBD += 1

                # read pulse shape:
                pulse_shape_data_IBD = np.loadtxt(input_path_IBD + "file{0:d}_evt{1:d}_prompt_signal_DCR.txt"
                                                  .format(filenumber, event))
                # pulse_shape_data_IBD = np.loadtxt(input_path_IBD + "file{0:d}_evt{1:d}_prompt_signal.txt"
                #                                   .format(filenumber, event))

                # get reconstructed positions x, y, z in mm:
                x_reco_IBD = pulse_shape_data_IBD[0]
                y_reco_IBD = pulse_shape_data_IBD[1]
                z_reco_IBD = pulse_shape_data_IBD[2]

                # 0th entry in pulse_shape_data_IBD is minimum of time window in ns:
                min_time_IBD = pulse_shape_data_IBD[3]
                # 1st entry in pulse_shape_data_IBD is maximum of time window in ns:
                max_time_IBD = pulse_shape_data_IBD[4]
                # 2nd entry in pulse_shape_data_IBD is bin-width in ns:
                bin_width_IBD = pulse_shape_data_IBD[5]

                # check if bin_width_IBD is equal to binwidth:
                if bin_width_IBD != binwidth:
                    sys.exit("ERROR: set binwidth ({0:.1f} ns) != binwidth from pulse shape file ({1:.1f} ns)"
                             .format(binwidth, bin_width_IBD))

                # the rest of pulse_shape_data_IBD is the hittime distribution histogram in nPE per bin. Take only the
                # prompt signal defined by start_time and end_time:
                nPE_per_bin_IBD = pulse_shape_data_IBD[6:
                                                       (int((end_time + binwidth + np.abs(min_time_IBD)) / binwidth)+3)]

                # set the time window corresponding to nPE_per_bin_IBD:
                time_window_IBD = np.arange(min_time_IBD, end_time+binwidth, binwidth)

                # analyze pulse shape of this event:
                ttr_IBD, npe_norm_IBD = pulse_shape(time_window_IBD, nPE_per_bin_IBD, start_tail[index],
                                                    stop_tail[index1])

                # append tail-to-total ratio to array:
                array_TTR_IBD.append(ttr_IBD)
                # append filenumber and eventID to arrays:
                filenumber_IBD.append(filenumber)
                evtID_IBD.append(event)

                # append zeros to npe_norm_IBD to get a average length of the hittimes:
                npe_norm_IBD = np.pad(npe_norm_IBD, (0, length_average_hittime - len(npe_norm_IBD)), 'constant',
                                      constant_values=(0.0, 0.0))

                # add the normalized pulse_shape (npe_norm_IBD) to the average hittime distribution
                # (hittime_average_positron):
                hittime_average_IBD = hittime_average_IBD + npe_norm_IBD

        # to get the average hittime distribution with a maximum of 1, normalize hittime_average_IBD with
        # max(hittime_average_IBD):
        hittime_average_IBD = hittime_average_IBD / max(hittime_average_IBD)

        """ analyze the pulse shape of prompt signal of NC events: """
        print("analyze NC events...")

        # path, where NC pulse shapes are stored:
        input_path_NC = "/home/astro/blum/juno/atmoNC/data_NC/output_detsim_v2/hittimes/"

        # number of events that are analyzed:
        num_events_analyzed_NC = 0

        # preallocate array, where tail-to-total ratios are stored:
        array_TTR_NC = []
        filenumber_NC = []
        evtID_NC = []

        # preallocate array, where a average hittime-distribution of NC events are stored (number of pe per bin):
        # length of average hittime distribution (number of bins):
        length_average_hittime = int((end_time - start_time + binwidth) / binwidth) + 10
        hittime_average_NC = np.zeros(length_average_hittime)

        # loop over all NC pulse shapes:
        for filenumber in range(start_filenumber_NC, stop_filenumber_NC + 1, 1):
            # print(filenumber)
            for event in range(num_evts_per_file):

                # check if the file exists:
                flag_file_exists = os.path.isfile(input_path_NC + "file{0:d}_evt{1:d}_prompt_signal_DCR.txt"
                                                  .format(filenumber, event))
                # flag_file_exists = os.path.isfile(input_path_NC + "file{0:d}_evt{1:d}_prompt_signal.txt"
                #                                   .format(filenumber, event))

                if flag_file_exists:
                    # read pulse shape:
                    pulse_shape_data_NC = np.loadtxt(input_path_NC + "file{0:d}_evt{1:d}_prompt_signal_DCR.txt"
                                                     .format(filenumber, event))
                    # pulse_shape_data_NC = np.loadtxt(input_path_NC + "file{0:d}_evt{1:d}_prompt_signal.txt"
                    #                                  .format(filenumber, event))
                else:
                    # file does not exist (event without initial particles)
                    # go to next pulse shape
                    continue

                # get reconstructed positions x, y, z in mm:
                x_reco_NC = pulse_shape_data_NC[0]
                y_reco_NC = pulse_shape_data_NC[1]
                z_reco_NC = pulse_shape_data_NC[2]

                # 0th entry in pulse_shape_data_NC is minimum of time window in ns:
                min_time_NC = pulse_shape_data_NC[3]
                # 1st entry in pulse_shape_data_NC is maximum of time window in ns:
                max_time_NC = pulse_shape_data_NC[4]
                # 2nd entry in pulse_shape_data_NC is bin-width in ns:
                bin_width_NC = pulse_shape_data_NC[5]

                # check if bin_width_NC is equal to binwidth:
                if bin_width_NC != binwidth:
                    sys.exit("ERROR: set binwidth ({0:.1f} ns) != binwidth from pulse shape file ({1:.1f} ns)"
                             .format(binwidth, bin_width_NC))

                # the rest of pulse_shape_data_NC is the hittime distribution histogram in nPE per bin. Take only the
                # prompt signal defined by start_time and end_time:
                nPE_per_bin_NC = pulse_shape_data_NC[6:(int((end_time + binwidth + np.abs(min_time_NC)) / binwidth)+3)]

                # check if nPE_per_bin_NC has entries above 0:
                if max(nPE_per_bin_NC) < 500:
                    # only analyze event with max(nPE_per_bin) >= 500 PE
                    continue

                # set the time window corresponding to nPE_per_bin_NC:
                time_window_NC = np.arange(min_time_NC, end_time+binwidth, binwidth)

                # increment num_events_analyzed_NC:
                num_events_analyzed_NC += 1

                # analyze pulse shape of this event:
                ttr_NC, npe_norm_NC = pulse_shape(time_window_NC, nPE_per_bin_NC, start_tail[index], stop_tail[index1])

                # append tail-to-total ratio to array:
                array_TTR_NC.append(ttr_NC)
                # append filenumber and evtID to arrays:
                filenumber_NC.append(filenumber)
                evtID_NC.append(event)

                # append zeros to npe_norm_IBD to get a average length of the hittimes:
                npe_norm_NC = np.pad(npe_norm_NC, (0, length_average_hittime - len(npe_norm_NC)), 'constant',
                                     constant_values=(0.0, 0.0))

                # add the normalized pulse_shape (npe_norm_NC) to the average hittime distribution
                # (hittime_average_positron):
                hittime_average_NC = hittime_average_NC + npe_norm_NC

        # to get the average hittime distribution with a maximum of 1, normalize hittime_average_NC with
        # max(hittime_average_NC):
        hittime_average_NC = hittime_average_NC / max(hittime_average_NC)

        """ calculate the efficiency of the tail-to-total pulse shape analysis (for the values of start_tail[index] and 
        stop_tail[index1]): """
        if max(array_TTR_IBD) >= max(array_TTR_NC):
            maximum_tot_value = max(array_TTR_IBD)
        else:
            maximum_tot_value = max(array_TTR_NC)

        # check the efficiency of PSD for different cut-efficiencies of NC events:
        efficiency_NC_99 = 99.0
        efficiency_NC_98 = 98.0
        efficiency_NC_97 = 97.0
        efficiency_NC_96 = 96.0
        efficiency_NC_95 = 95.0
        efficiency_NC_90 = 90.0

        # calculate the cut-efficiencies for positrons depending of efficiency_NC:
        efficiency_IBD_99, ttr_cut_value_99 = tot_efficiency(array_TTR_IBD, array_TTR_NC, efficiency_NC_99)
        efficiency_IBD_98, ttr_cut_value_98 = tot_efficiency(array_TTR_IBD, array_TTR_NC, efficiency_NC_98)
        efficiency_IBD_97, ttr_cut_value_97 = tot_efficiency(array_TTR_IBD, array_TTR_NC, efficiency_NC_97)
        efficiency_IBD_96, ttr_cut_value_96 = tot_efficiency(array_TTR_IBD, array_TTR_NC, efficiency_NC_96)
        efficiency_IBD_95, ttr_cut_value_95 = tot_efficiency(array_TTR_IBD, array_TTR_NC, efficiency_NC_95)
        efficiency_IBD_90, ttr_cut_value_90 = tot_efficiency(array_TTR_IBD, array_TTR_NC, efficiency_NC_90)

        print("tail start = {0:.1f} ns".format(start_tail[index]))
        print("tail end = {0:.1f} ns".format(stop_tail[index1]))
        print("NC efficiency = {0:.1f} %".format(efficiency_NC_90))
        print(efficiency_IBD_90)
        print(ttr_cut_value_90)
        print("NC efficiency = {0:.1f} %".format(efficiency_NC_95))
        print(efficiency_IBD_95)
        print(ttr_cut_value_95)
        print("NC efficiency = {0:.1f} %".format(efficiency_NC_96))
        print(efficiency_IBD_96)
        print(ttr_cut_value_96)
        print("NC efficiency = {0:.1f} %".format(efficiency_NC_97))
        print(efficiency_IBD_97)
        print(ttr_cut_value_97)
        print("NC efficiency = {0:.1f} %".format(efficiency_NC_98))
        print(efficiency_IBD_98)
        print(ttr_cut_value_98)
        print("NC efficiency = {0:.1f} %".format(efficiency_NC_99))
        print(efficiency_IBD_99)
        print(ttr_cut_value_99)

        if CREATE_PLOTS:

            if SAVE_ARRAY_TTR:
                np.savetxt(output_path + 'TTR_IBD_{0:.0f}ns_{1:.0f}ns_{2:d}.txt'
                           .format(start_tail[index], stop_tail[index1], analysis_part),
                           np.c_[filenumber_IBD, evtID_IBD, array_TTR_IBD], fmt='%.5f',
                           header='filenumber | evtID | TTR values of {0:.0f} IBD events'.format(len(array_TTR_IBD)))

                np.savetxt(output_path + 'TTR_atmoNC_{0:.0f}ns_{1:.0f}ns_{2:d}.txt'
                           .format(start_tail[index], stop_tail[index1], analysis_part),
                           np.c_[filenumber_NC, evtID_NC, array_TTR_NC], fmt='%.5f',
                           header='filenumber | evtID | TTR values of {0:.0f} atmo. NC events'
                           .format(len(array_TTR_NC)))

            """ Display tot ratios in histograms: """
            h1 = plt.figure(1, figsize=(15, 8))
            First_bin = 0.0
            Last_bin = maximum_tot_value
            Bin_width = (Last_bin-First_bin) / 200
            Bins = np.arange(First_bin, Last_bin+Bin_width, Bin_width)

            plt.hist(array_TTR_IBD, bins=Bins, histtype="step", align='mid', color="r", linewidth=1.5,
                     label="prompt signal of IBD events (entries = {0:d})".format(num_events_analyzed_IBD))

            plt.hist(array_TTR_NC, bins=Bins, histtype="step", align='mid', color="b", linewidth=1.5,
                     label="prompt signal of NC events (entries = {0:d})".format(num_events_analyzed_NC))

            plt.xlabel("tail-to-total ratio")
            plt.ylabel("events")
            plt.title("Tail-to-total ratio of prompt signals of IBD and NC events" +
                      "\n(tail window {0:0.1f} ns to {1:0.1f} ns)".format(start_tail[index], stop_tail[index1]))
            plt.legend()
            plt.grid()
            plt.savefig(output_path + "TTR_tail_{0:0.0f}ns_to_{1:0.0f}ns.png".format(start_tail[index],
                                                                                     stop_tail[index1]))
            plt.close()

            """ Display tot ratios in histograms: """
            h2 = plt.figure(2, figsize=(15, 8))
            First_bin = 0.0
            Last_bin = maximum_tot_value
            Bin_width = (Last_bin-First_bin) / 200
            Bins = np.arange(First_bin, Last_bin+Bin_width, Bin_width)
            n_pos_1, bins_pos_1, patches_pos_1 = plt.hist(array_TTR_IBD, bins=Bins, histtype="step",
                                                          align='mid',
                                                          color="r", linewidth=1.5,
                                                          label="prompt signal of IBD events "
                                                                "(entries = {0:d})"
                                                          .format(num_events_analyzed_IBD))
            n_nc_1, bins_nc_1, patches_nc_1 = plt.hist(array_TTR_NC, bins=Bins, histtype="step", align='mid',
                                                       color="b",
                                                       linewidth=1.5,
                                                       label="prompt signal of NC events (entries = {0:d})"
                                                       .format(num_events_analyzed_NC))
            plt.vlines(ttr_cut_value_99, 0, max(n_pos_1)+max(n_pos_1)/10, colors="k", linestyles="-",
                       label="$\\epsilon_{NC}$ = "+"{0:0.2f} %\n".format(efficiency_NC_99)+"$\\epsilon_{IBD}$ = " +
                             "{0:0.2f} %\n".format(efficiency_IBD_99)+"tot value = {0:.5f}".format(ttr_cut_value_99))
            plt.vlines(ttr_cut_value_95, 0, max(n_pos_1)+max(n_pos_1)/10, colors="k", linestyles="--",
                       label="$\\epsilon_{NC}$ = "+"{0:0.2f} %\n".format(efficiency_NC_95)+"$\\epsilon_{IBD}$ = " +
                             "{0:0.2f} %\n".format(efficiency_IBD_95)+"tot value = {0:.5f}".format(ttr_cut_value_95))
            plt.xlabel("tail-to-total ratio")
            plt.ylabel("events")
            plt.title("Tail-to-total value for prompt signals of IBD and NC events" +
                      "\n(tail window {0:0.1f} ns to {1:0.1f} ns)".format(start_tail[index], stop_tail[index1]))
            plt.legend()
            plt.grid()
            plt.savefig(output_path + "TTR_tail_{0:0.0f}ns_to_{1:0.0f}ns_efficiency.png"
                        .format(start_tail[index], stop_tail[index1]))
            plt.close()

""" Display the average hittime distributions of positrons and IBD-like NC events: """
h3 = plt.figure(3, figsize=(15, 8))
bin_edges = np.arange(start_time, end_time+binwidth+10*binwidth, binwidth)
plt.semilogy(bin_edges, hittime_average_IBD, linestyle="steps", color="r",
             label="average pulse shape of IBD events")
plt.semilogy(bin_edges, hittime_average_NC, linestyle="steps", color="b",
             label="average pulse shape of NC events")
plt.xlabel("hittime in ns")
plt.ylabel("probability per bin (bin-width = {0:0.1f} ns)".format(binwidth))
plt.xlim(xmin=start_time, xmax=end_time)
plt.ylim(ymin=1e-4, ymax=2.0)
plt.title("Average pulse shape of prompt signals")
plt.legend()
plt.grid()
plt.savefig(output_path + "average_pulse_shape.png")
plt.close()

""" print the best IBD efficiencies with the corresponding start and stop time of the tail for each NC efficiency: """
print("NC efficiency = {0:.1f} %, best IBD efficiency = {1:.2f} %, tail start = {2:.1f} ns, tail end = {3:.1f} ns"
      ", tail-to-total value = {4:.5f}"
      .format(efficiency_NC_90, best_IBD_eff_for_NC_90, tail_start_90, tail_stop_90, ttr_value_90))
print("NC efficiency = {0:.1f} %, best IBD efficiency = {1:.2f} %, tail start = {2:.1f} ns, tail end = {3:.1f} ns"
      ", tail-to-total value = {4:.5f}"
      .format(efficiency_NC_95, best_IBD_eff_for_NC_95, tail_start_95, tail_stop_95, ttr_value_95))
print("NC efficiency = {0:.1f} %, best IBD efficiency = {1:.2f} %, tail start = {2:.1f} ns, tail end = {3:.1f} ns"
      ", tail-to-total value = {4:.5f}"
      .format(efficiency_NC_96, best_IBD_eff_for_NC_96, tail_start_96, tail_stop_96, ttr_value_96))
print("NC efficiency = {0:.1f} %, best IBD efficiency = {1:.2f} %, tail start = {2:.1f} ns, tail end = {3:.1f} ns"
      ", tail-to-total value = {4:.5f}"
      .format(efficiency_NC_97, best_IBD_eff_for_NC_97, tail_start_97, tail_stop_97, ttr_value_97))
print("NC efficiency = {0:.1f} %, best IBD efficiency = {1:.2f} %, tail start = {2:.1f} ns, tail end = {3:.1f} ns"
      ", tail-to-total value = {4:.5f}"
      .format(efficiency_NC_98, best_IBD_eff_for_NC_98, tail_start_98, tail_stop_98, ttr_value_98))
print("NC efficiency = {0:.1f} %, best IBD efficiency = {1:.2f} %, tail start = {2:.1f} ns, tail end = {3:.1f} ns"
      ", tail-to-total value = {4:.5f}"
      .format(efficiency_NC_99, best_IBD_eff_for_NC_99, tail_start_99, tail_stop_99, ttr_value_99))

""" save different number of events of PSD cut to txt file: """
np.savetxt(output_path + "numbers_PSD.txt",
           np.array([(stop_filenumber_NC-start_filenumber_NC+1)*num_evts_per_file, num_events_analyzed_NC,
                     (stop_filenumber_IBD-start_filenumber_IBD+1)*num_evts_per_file, num_events_analyzed_IBD]),
           fmt='%i',
           header='number of events from analyze_PSD_cut_v2.py ({0}):'
                  '\nanalyzed atmoNC files: file{1:d}_evt0_pulse_shape.txt to file{2:d}_evt99_pulse_shape.txt,'
                  '\nanalyzed IBD files: file{3:d}_evt0_pulse_shape.txt to file{4:d}_evt99_pulse_shape.txt,'
                  '\nValues below:'
                  '\nnumber of total NC events,'
                  '\nnumber of analyzed NC events,'
                  '\nnumber of total IBD events,'
                  '\nnumber of analyzed IBD events:'
                  .format(now, start_filenumber_NC, stop_filenumber_NC, start_filenumber_IBD, stop_filenumber_IBD))

