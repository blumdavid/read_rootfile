""" Script to get the energy of the prompt signal of preselected events of atmospheric NC neutrino background that are
    simulated with JUNO detector simulation.

    1.  read only the preselected events (preselection done with script preselection_detsim_user.py and saved in folder
        /home/astro/blum/juno/atmoNC/data_NC/output_preselection/preselection_detsim/ in files evtID_preselected_{}.txt)

    2.  create histogram with hit-times for these events and analyze two parts of this hittime distribution:
        2.1:    check, if there is really only one delayed signal from neutron capture
        2.2:    calculate the number of pe of the prompt signal

    3.  do a cut on the prompt signal: use function conversion_npe_to_evis() of atmoNC_spectrum.py and convert the
        number of pe of the prompt signal to visible energy in MeV and do a cut on the visible energy (10 MeV to 100
        MeV)

    3.  Save number of pe of the prompt signal and number of pe of delayed signal of each event in txt file.


    This txt file can than be analyzed further with script atmoNC_spectrum.py:

    4.  Convert the energy of the prompt signal from number of pe to visible energy in the detector in MeV

    5.  Put all these calculated visible energies into histogram to get the spectrum of atmospheric NC neutrino
        background as function of the visible energy

    6.  Consider the event rate of NC interactions on C12 inside the detector (calculated with cross-sections and
        neutrino fluxes) and calculate the 'real' spectrum of atmospheric NC background, JUNO will measure after 10
        years of data taking
"""
import datetime
import ROOT
import sys
import numpy as np
from matplotlib import pyplot as plt


def conversion_npe_to_evis(number_photo_electron):
    """
    Function to calculate the visible energy in MeV for a given number of p.e. for the prompt signal.
    This function is the result of linear fit from script check_conversion_npe_mev.py.

    :param number_photo_electron: number of photo-electrons of the prompt signal
    :return: quenched deposited energy (visible energy in MeV)
    """
    # TODO-me: if this function is changed, also change function in atmoNC_spectrum!!!!

    # TODO-me: parameters of the fit has to be checked!!!!!!!!!
    # first fit parameter (slope) in MeV/nPE:
    parameter_a = 0.0007872

    energy = parameter_a * number_photo_electron

    return energy


def analyze_delayed_signal(npe_per_time, bins_time, first_index, threshold, threshold2, min_pe_delayed, max_pe_delayed,
                           evt):
    """
    function to analyze the time window, where a delayed signal could be. This time window is check for a possible
    delayed signal (signal has to be greater than 'threshold' and nPE of signal peak must be between 'min_pe_delayed'
    and 'max_pe_delayed')

    :param npe_per_time: number of pe as function of hittime for time window of delayed signal (array of integer)
    :param bins_time: array of bins, which contains the information about the hittime in ns (array of float)
    :param first_index: first index of npe_per_time, that should be analyzed (integer)
    :param threshold: threshold in nPE per bin, that a delayed signal must have (integer)
    :param threshold2: threshold in nPE per bin, that specifies the minimal nPE corresponding to a peak
    :param min_pe_delayed:  minimum number of PE for delayed energy cut of neutron capture on Hydrogen
                            (values from check_delayed_energy.py)
    :param max_pe_delayed:  maximum number of PE for delayed energy cut of neutron capture on Hydrogen
                            (values from check_delayed_energy.py)
    :param evt: event ID of the event

    :return:
    """
    # preallocate number of delayed signal and index_after_peak1 and number of pe in the delayed signal:
    num_del_signal = 0
    index_after_peak1 = len(npe_per_time)
    num_pe_del = 0
    # preallocate flag, if event is rejected by delayed energy cut or not:
    flag_rejected_delayed_energy_cut = 0

    # loop over values of the histogram bins of delayed window:
    for index4 in range(first_index, len(npe_per_time)):
        # check if number of PE in this bin is above the threshold:
        if npe_per_time[index4] > threshold:
            # possible delayed signal (signal in delayed window):

            # check if index4 = first_value -> first value of npe_per_time above threshold -> break from for loop and
            # return index_after_peak1 = index4 + 1:
            if index4 == first_index:
                index_after_peak1 = index4 + 1
                print("WARNING: npe_per_time[0] > threshold in evtID = {0:d}".format(evt))
                break

            # calculate number of PEs in this signal peak:
            # preallocate sum of PE per bin in the signal peak:
            sum_pe_peak = 0

            # add nPE of npe_per_time[index4] to sum_pe_peak:
            sum_pe_peak = sum_pe_peak + npe_per_time[index4]

            # check previous bins:
            for num in range(1, 100, 1):
                # check if value in previous bins_time[index4 - num] is above threshold2:
                if npe_per_time[index4 - num] > threshold2:
                    # add nPE of this bin to sum_pe_peak:
                    sum_pe_peak = sum_pe_peak + npe_per_time[index4 - num]

                else:
                    # hittime, when signal peak begins, in ns:
                    begin_peak = bins_time[index4 - num]
                    break

            # check following bins:
            for num in range(1, 1000, 1):
                # check if index4 + num is in the range of npe_per_time array:
                if (index4 + num) >= len(npe_per_time):
                    print("Warning: iteration reaches last index of npe_per_time in evtID = {0:d}".format(evt))
                    break

                # check if value in following bins_time[index4 + num] is above threshold2:
                if npe_per_time[index4 + num] > threshold2:
                    # add nPE of this bin to sum_pe_peak:
                    sum_pe_peak = sum_pe_peak + npe_per_time[index4 + num]

                else:
                    # hittime, when signal peak ends, in ns:
                    end_peak = bins_hittime[index4 + num]

                    # get index, where npe_hittime is below threshold2:
                    index_after_peak1 = index4 + num
                    break

            # print("number of pe in delayed signal = {0:d}".format(sum_pe_peak))

            # set number of pe of delayed signal:
            num_pe_del = sum_pe_peak

            # first peak is analyzed:
            # check, if number of PE in this signal peak agree with delayed energy cut:
            if min_pe_delayed < sum_pe_peak < max_pe_delayed:
                # PE of signal peak agree with delayed energy cut
                # set delayed flag:
                num_del_signal = 1
            else:
                flag_rejected_delayed_energy_cut = 1
                print("number of pe in delayed signal = {0:d}".format(sum_pe_peak))

            # after analyzing the first peak, break out of for-loop (possible second delayed signal will be checked
            # below)
            break
        else:
            # go to next bin:
            continue

    return num_del_signal, index_after_peak1, num_pe_del, flag_rejected_delayed_energy_cut


# get the date and time, when the script was run:
date = datetime.datetime.now()
now = date.strftime("%Y-%m-%d %H:%M")

""" set the number of the first file and number of the last file that should be read: """
start_number = 0
stop_number = 999
# number of entries in the input files:
Number_entries_input = 100
# set the path of the input root files:
input_path_root = "/local/scratch1/pipc51/astro/blum/detsim_output_data/"
# set the path of evtID_preselected_{}.txt files:
input_path_preselect = "/home/astro/blum/juno/atmoNC/data_NC/output_preselection/preselection_detsim/"
# set the path of the output, where the txt file with the number of pe of each preselected events is saved:
output_path = "/home/astro/blum/juno/atmoNC/data_NC/output_detsim/"

""" define time window and bin width: """
# set time window of whole signal in ns:
min_time = 0
max_time = 1000000
# set time in ns, where the prompt signal should be 0:
time_limit_prompt = 500
# Set bin-width of hittime histogram in ns:
binwidth = 5.0

""" parameters for prompt energy cut: """
# minimal visible energy of prompt signal in MeV:
min_energy = 10.0
# maximal visible energy of prompt signal in MeV:
max_energy = 100.0
# preallocate number of events that are rejected by prompt energy cut:
number_rejected_prompt_cut = 0
# preallocate number of events where nPE of prompt signal is below min_energy:
number_rejected_prompt_cut_min = 0
# preallocate number of events where nPE of prompt signal is above max_energy:
number_rejected_prompt_cut_max = 0

""" thresholds and cuts for delayed signal: """
# Set threshold of number of PE per bin for possible delayed signal (bin-width = 5 ns):
threshold1_del = 50
# set threshold2 of number of PEs per bin (signal peak is summed as long as nPE is above threshold2):
threshold2_del = 0
# min and max number of PE for delayed energy cut (from check_delayed_energy.py):
min_PE_delayed = 2300.0
max_PE_delayed = 3600.0
# preallocate number of events that are rejected by delayed energy cut:
number_rejected_delayed_energy_cut = 0

# loop over the files:
for index in range(start_number, stop_number+1, 1):
    # read evtID_preselected_index.txt file:
    evtID_pre_arr = np.loadtxt(input_path_preselect + "evtID_preselected_{0:d}.txt".format(index))

    # load user_atmoNC_index.root file:
    rfile = ROOT.TFile(input_path_root + "user_atmoNC_{0:d}.root".format(index))
    print("... read {0}...".format(rfile))

    # get the "evt"-TTree from the TFile:
    rtree_evt = rfile.Get("evt")
    # get the number of events in the 'evt' Tree:
    number_events_evt = rtree_evt.GetEntries()
    # check number of events:
    if number_events_evt != Number_entries_input:
        sys.exit("ERROR: number of events in root file ({0:d}) != {1:d}"
                 .format(number_events_evt, Number_entries_input))

    # preallocate array, where total nPE of prompt signal per event for one file is saved:
    number_pe_total = np.array([])
    # preallocate array, where total nPE of delayed signal per event for one file is saved:
    number_pe_total_del = np.array([])

    # loop over the length of evtID_pre_arr and read only the preselected events from the root file:
    for index1 in range(len(evtID_pre_arr)):
        # get evtID of preselected event:
        event_id = int(evtID_pre_arr[index1])
        # get event of 'evt'-tree:
        rtree_evt.GetEntry(event_id)
        # get evtID of the tree and compare with evtID of evtID_pre_arr:
        evt_id = int(rtree_evt.GetBranch('evtID').GetLeaf('evtID').GetValue())
        if evt_id != event_id:
            sys.exit("ERROR: evtID of tree ({0:d}) != evtID of evtID_preselected.txt ({1:d})".format(evt_id, event_id))

        print("\nanalyze event {0:d}".format(evt_id))

        # get number of photons of this event:
        n_photons = int(rtree_evt.GetBranch('nPhotons').GetLeaf('nPhotons').GetValue())

        # preallocate empty array to build default hittime-histogram:
        hittime_array = []

        # loop over every photon in the event:
        for index2 in range(n_photons):

            # get PMT ID, where photon is absorbed:
            pmt_id = int(rtree_evt.GetBranch('pmtID').GetLeaf('pmtID').GetValue(index2))

            # only 20 inch PMTs (PMT ID of 20 inch PMTs are below 21000, PMT ID of 3 inch PMTs start at 290000):
            if pmt_id < 25000:
                # get nPE for this photon:
                n_pe = int(rtree_evt.GetBranch('nPE').GetLeaf('nPE').GetValue(index2))
                # check, if photon produces only 1 PE:
                if n_pe != 1:
                    print("{1:d} PE for 1 photon in event {0:d} in file user_atmoNC_{2:d}.root"
                          .format(event_id, n_pe, index))

                # get hittime of this photon:
                hit_time = float(rtree_evt.GetBranch('hitTime').GetLeaf('hitTime').GetValue(index2))

                # append hittime to array:
                hittime_array.append(hit_time)

            else:
                continue

        """ analyze prompt signal: """
        # build histogram, where hittimes are saved:
        # set bin-edges of hittime histogram in ns:
        bins_hittime = np.arange(min_time, max_time + 2 * binwidth, binwidth)
        # build hittime histogram:
        npe_per_hittime, bin_edges_hittime = np.histogram(hittime_array, bins_hittime)

        # get index of bins_hittime corresponding to min_time (should be index = 0):
        index_min_hittime_prompt = int(min_time / binwidth)

        # Where does prompt signal end?
        # get index of bins_hittime corresponding to time_limit_prompt
        index_time_limit_prompt = int(time_limit_prompt / binwidth)
        # check if npe_per_hittime is 0 for this index:
        if npe_per_hittime[index_time_limit_prompt] == 0:
            # prompt signal already 0:
            index_max_hittime_prompt = index_time_limit_prompt
        else:
            # prompt signal not yet 0.
            # loop over npe_per_hittime from index_time_limit_prompt until npe_per_hittime is 0:
            for index3 in range(index_time_limit_prompt, index_time_limit_prompt+200):
                if npe_per_hittime[index3] == 0:
                    index_max_hittime_prompt = index3
                    break

        # calculate nPE as function of hittime only for prompt time window (from min_hittime_prompt to
        # max_hittime_prompt+1, last index should be included):
        npe_per_hittime_prompt = npe_per_hittime[index_min_hittime_prompt:index_max_hittime_prompt+1]
        # bin edges of hittime histogram only for prompt time window:
        bins_hittime_prompt = bin_edges_hittime[index_min_hittime_prompt:index_max_hittime_prompt+1]

        # get the minimum and maximum time of the prompt signal time window in ns:
        min_time_prompt = bins_hittime_prompt[0]
        max_time_prompt = bins_hittime_prompt[-1]

        # sum up the values of npe_per_hittime_prompt to get the total number of pe of the prompt signal:
        number_pe_prompt = np.sum(npe_per_hittime_prompt)

        # convert the total number of pe to quenched deposited energy in MeV:
        quenched_deposit_energy = conversion_npe_to_evis(number_pe_prompt)

        # check, if energy is in the correct time window:
        if quenched_deposit_energy < min_energy or quenched_deposit_energy > max_energy:
            # event is rejected:
            number_rejected_prompt_cut += 1
            # check further:
            if quenched_deposit_energy < min_energy:
                # e_vis to small:
                number_rejected_prompt_cut_min += 1
            elif quenched_deposit_energy > max_energy:
                # e_vis too large:
                number_rejected_prompt_cut_max += 1

            # go to next event:
            continue

        # append number of pe to array:
        number_pe_total = np.append(number_pe_total, number_pe_prompt)

        """ analyze delayed signal: """
        # INFO-me: nPE of delayed signal is only saved as info -> NO cut is applied to delayed signal in this script!
        print("analyze delayed signal")
        # get index of bins_hittime corresponding to the end of the prompt signal window:
        index_min_hittime_del = index_max_hittime_prompt
        # get index of bins_hittime corresponding to max_time:
        index_max_hittime_del = int(max_time / binwidth)

        # calculate nPE as function of hittime only for delayed time window (from min_hittime_del to max_hittime_del+1):
        npe_per_hittime_del = npe_per_hittime[index_min_hittime_del:index_max_hittime_del+1]
        # bin edges of hittime histogram only for delayed time window:
        bins_hittime_del = bin_edges_hittime[index_min_hittime_del:index_max_hittime_del+1]

        # get the minimum and maximum time of the delayed signal time window in ns:
        min_time_delayed = bins_hittime_del[0]
        max_time_delayed = bins_hittime_del[-1] + binwidth

        # preallocate number of possible delayed signals:
        number_delayed_signal = 0
        # preallocate first index of npe_per_hittime_del:
        index_first_del = 0
        # preallocate number of pe of delayed signal:
        number_pe_delayed = 0
        # preallocate flag, if event is rejected by delayed energy cut (0 means it passes the cut):
        number_rejected_delayed = 0

        # analyze npe_per_hittime_del for possible delayed signals. As long as number_delayed_signal<2 and as long as
        # index has not reached the end of npe_per_hittime_del, check event for possible delayed signals
        while number_delayed_signal < 2 and index_first_del < len(npe_per_hittime_del):
            number_delayed, index_first_del, num_pe_delayed, number_rejected_delayed = \
                analyze_delayed_signal(npe_per_hittime_del, bins_hittime_del, index_first_del, threshold1_del,
                                       threshold2_del, min_PE_delayed, max_PE_delayed, event_id)

            number_delayed_signal = number_delayed_signal + number_delayed
            number_pe_delayed = number_pe_delayed + num_pe_delayed

        # add number_rejected_delayed (0 if event passes delayed energy cut, 1 if event is rejected by delayed energy
        # cut) to number_rejected_delayed_energy_cut:
        number_rejected_delayed_energy_cut += number_rejected_delayed

        # append number of pe of delayed signal to array:
        number_pe_total_del = np.append(number_pe_total_del, number_pe_delayed)

        # check number_delayed_signal:
        if number_delayed_signal == 0:
            print("---------------ERROR: no delayed signal in event {0:d}".format(event_id))
            h1 = plt.figure(1)
            plt.step(bins_hittime_del, npe_per_hittime_del, label="nPE of peak = {0:d}".format(number_pe_delayed))
            plt.xlabel("hit-time in ns")
            plt.ylabel("number of p.e. per bin (bin-width = {0:0.2f} ns)".format(binwidth))
            plt.title("Hit-time distribution of delayed time window of event {0:d}\nNO delayed signal".format(event_id))
            plt.xlim(xmin=min_time_delayed, xmax=max_time_delayed)
            plt.legend()
            plt.grid()
            plt.savefig(output_path + "file{1:d}_evt{0:d}_no_delayed_signal.png".format(event_id, index))
            plt.close()
            # plt.show()

        elif number_delayed_signal > 1:
            print("+++++++++++++++ERROR: more than one delayed signal in event {0:d}".format(event_id))
            h1 = plt.figure(1)
            plt.step(bins_hittime_del, npe_per_hittime_del, label="nPE of peak = {0:d}".format(number_pe_delayed))
            plt.xlabel("hit-time in ns")
            plt.ylabel("number of p.e. per bin (bin-width = {0:0.2f} ns)".format(binwidth))
            plt.title("Hit-time distribution of delayed time window of event {0:d}\n"
                      "More than 1 delayed signals".format(event_id))
            plt.xlim(xmin=min_time_delayed, xmax=max_time_delayed)
            plt.legend()
            plt.grid()
            plt.savefig(output_path + "file{1:d}_evt{0:d}_more_than_1_delayed_signal.png".format(event_id, index))
            plt.close()
            # plt.show()

        # save hittime distribution of the IBD-like events to png and txt file:
        h2 = plt.figure(2)
        plt.step(bins_hittime_prompt, npe_per_hittime_prompt, label="number of pe = {0:d}\n"
                                                                    "visible energy = {1:.3f} MeV"
                 .format(number_pe_prompt, quenched_deposit_energy))
        plt.xlabel("hit-time in ns")
        plt.ylabel("number of p.e. per bin (bin-width = {0:0.2f} ns)".format(binwidth))
        plt.title("Hit-time distribution of prompt time window of event {0:d}".format(event_id))
        plt.xlim(xmin=min_time_prompt, xmax=max_time_prompt)
        plt.legend()
        plt.grid()
        plt.savefig(output_path + "file{1:d}_evt{0:d}_prompt_signal.png".format(event_id, index))
        plt.close()
        # plt.show()

        # save npe_per_hittime_prompt to txt file:
        # build list, where 0th entry is start-hittime in ns, 1st entry is last-hittime in ns, 2nd entry is binwidth in
        # ns and the following entries are nPE of each hittime-bin of prompt signal:
        npe_per_hittime_prompt_save = [min_time_prompt, max_time_prompt, binwidth]
        npe_per_hittime_prompt_save.extend(npe_per_hittime_prompt)
        np.savetxt(output_path + "file{0:d}_evt{1:d}_prompt_signal.txt".format(index, event_id),
                   npe_per_hittime_prompt_save, fmt='%1.2f',
                   header="Number of pe as function of the hittime of the prompt signal of file user_atmoNC_{0:d}.root,"
                          "\npreselected event {1:d} (analyzed with script prompt_signal_preselected_evts.py, {2}):"
                          "\ntime window of hittime: from {3:.3f} ns to {4:.3f} ns with bin-width = {5:0.3f} ns,"
                          "\nEnergy cut on prompt signal is applied: {6:0.1f} MeV <= E_vis <= {7:0.1f} MeV,"
                          "\nConversion function E_vis = 0.0007872 * nPE:"
                   .format(index, event_id, now, min_time_prompt, max_time_prompt, binwidth, min_energy, max_energy))

    # save array number_pe_total to txt file:
    np.savetxt(output_path + "number_pe_file{0:d}.txt".format(index), number_pe_total, fmt="%i",
               header="Total number of pe of prompt signal for every preselected event (that pass prompt energy cut) in"
                      " file user_atmoNC_{0:d}.root\n"
                      "({1}). Number of pe is analyzed with script prompt_signal_preselected_evts.py.\n"
                      "Time window of prompt signal is defined from {2:0.2f} ns to around {3:0.2f} ns (number of events"
                      " = {4:d}).\n"
                      "Energy cut on prompt signal is applied: {5:0.1f} MeV <= E_vis <= {6:0.1f} MeV,\n"
                      "Conversion function E_vis = 0.0007872 * nPE:"
               .format(index, now, min_time_prompt, time_limit_prompt, len(number_pe_total), min_energy, max_energy))

    # save array number_pe_total_del to txt file:
    np.savetxt(output_path + "number_pe_delayed_file{0:d}.txt".format(index), number_pe_total_del, fmt="%i",
               header="Total number of pe of delayed signal for every preselected event (that pass prompt energy cut) "
                      "in file user_atmoNC_{0:d}.root\n"
                      "({1}). Number of pe is analyzed with script prompt_signal_preselected_evts.py.\n"
                      "Time window of delayed signal is defined from around {2:0.2f} ns to {3:0.2f} ns.\n"
                      "Threshold of delayed signal is set to {4:0.2f} PE and {5:0.2f} PE (number of events = {6:d}).\n"
                      "A not very strict energy cut on the delayed signal ({7:0.0f} PE <= nPE delayed <= {8:0.0f} PE) "
                      "is also applied:"
               .format(index, now, time_limit_prompt, max_time_delayed, threshold1_del, threshold2_del,
                       len(number_pe_total_del), min_PE_delayed, max_PE_delayed))

print("number of events (of preselect. evts) rejected by prompt energy cut = {0:d}".format(number_rejected_prompt_cut))
print("number of events (of preselect. evts) with nPE of prompt signal < min_energy = {0:d}"
      .format(number_rejected_prompt_cut_min))
print("number of events (of preselect. evts) with nPE of prompt signal > max_energy = {0:d}"
      .format(number_rejected_prompt_cut_max))
print("number of events (of preselect. evts) rejected by delayed energy cut = {0:d}"
      .format(number_rejected_delayed_energy_cut))
