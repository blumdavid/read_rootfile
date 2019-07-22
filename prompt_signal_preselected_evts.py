""" Script to get the energy and hittime distribution of the prompt signal of preselected events of atmospheric
    NC neutrino background that are simulated with JUNO detector simulation.

    1.  read only the preselected events (preselection done with script preselection_detsim_user.py and saved in folder
        /home/astro/blum/juno/atmoNC/data_NC/output_preselection/preselection_detsim/ in files evtID_preselected_{}.txt)

    2.  Calculate hittime distribution (with time-of-flight correction and PMT time resolution) for each event:
        Procedure to get the hittime distribution with vertex reconstruction and time smearing of PMTs (same procedure
        like in script hittime_distribution_positron.py):

        2.1.    calculate time of flight:
            2.1.1   for every photon, that hits a PMT (20inch and 3inch), take the PMT position (via PMT ID from file
                    PMT_position.root) and calculate the time-of-flight distance with the reconstructed position from
                    file evtID_preselected_{}.txt
            2.1.2   with the time-of-flight distance, calculate the time-of-flight of this photon from production to
                    PMT by considering an effective speed of light in the LS.

        2.2.    consider TTS of PMTs:
            2.2.1   for every photon, that hits a PMT (20inch and 3inch), take the time resolution (sigma) of the PMT
                    (via PMT ID either from file PmtData.root for the 20inch PMTs or set TTS = 5 ns for 3inch PMTs.)
            2.2.2   the TTS of the PMT is FWHM. Therefore calculate sigma from TTS (FWHM = 2*sqrt(2*ln(2)) * sigma).
            2.2.3   smear hittime of detsim with gaussian of sigma (time resolution) around the value of detsim hittime
                    to get the smeared hittime

        2.3.    for every photon, calculate the 'real' hittime (= smeared hittime - time_of_flight) and store it in
                array

        2.4.    Do points 2.1 to 2.3 for every photon. Then you get the correct hittime of this event. Build histogram
                with correct hittimes.

    3.  Take the prompt signal of the corrected hittime histogram and do a cut on the prompt signal:
        use function conversion_npe_to_evis() and convert the number of pe of the prompt signal to visible energy in
        MeV and do a cut on the visible energy (10 MeV to 100 MeV).
        Only analyze events further that pass prompt energy cut.

    4.  Analyze delayed signal:


    5.  Save the hittime histogram of the prompt signal to txt file and png file for further analysis with script
        pulse_shape_analysis.py

    6.  Save number of pe of the prompt signal and number of pe of delayed signal of each event in txt file.
"""
import datetime
import ROOT
import sys
import NC_background_functions
import numpy as np
from matplotlib import pyplot as plt

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
min_time = -50
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
min_PE_delayed = 2805.53
max_PE_delayed = 3731.04
# preallocate number of events that are rejected by delayed energy cut:
number_rejected_delayed_energy_cut = 0
# preallocate array, where npe of delayed signal, that wouldn't pass the delayed energy cut are saved:
number_pe_delayed_rejected_array = np.array([])

""" load position of the PMTs and corresponding PMT ID from file PMT_position.root: """
file_PMT_position = "/home/astro/blum/juno/atmoNC/PMT_information/PMT_position.root"
# array with PMT ID and corresponding x, y, z position in mm:
pmtID_pos_file, x_pos_pmt, y_pos_pmt, z_pos_pmt = NC_background_functions.get_pmt_position(file_PMT_position)

""" load 'time resolution' in ns of the 20 inch PMTs and corresponding PMT ID from file PmtData.root: """
file_PMT_time = "/home/astro/blum/juno/atmoNC/PMT_information/PmtData.root"
# array with PMT ID and corresponding sigma in ns:
pmtID_time_file, sigma_time_20inch = NC_background_functions.get_20inchpmt_tts(file_PMT_time)
# set TTS (FWHM) of the 3inch PMTs in ns:
tts_3inch = 5.0
# calculate time resolution (sigma) for the 3inch PMTs in ns:
sigma_time_3inch = tts_3inch / (2 * np.sqrt(2 * np.log(2)))
# set effective speed of light in the liquid scintillator in mm/ns (see page 7 of c_effective_JUNO-doc-3144-v2.pdf in
# folder /home/astro/blum/PhD/paper/Pulse_Shape_Discrimination/). Effective refraction index in LS n_eff = 1.54.
# c/n_eff = 299792458 m / 1.54 s ~ 194670427 m/s = 194670427 * 10**(-6) mm/ns ~ 194.67 mm/ns:
c_effective = 194.67

# loop over the files:
for index in range(start_number, stop_number+1, 1):
    # read evtID_preselected_{}.txt file:
    evtID_pre_arr, x_reco_arr, y_reco_arr, z_reco_arr = np.loadtxt(input_path_preselect + "evtID_preselected_{0:d}.txt"
                                                                   .format(index), unpack=True)

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

        # print("\nanalyze event {0:d}".format(evt_id))

        """ calculate the real hittime distribution (time of flight correction with reconstructed position and time 
        smearing with TTS for each hit): """
        # get number of photons of this event:
        n_photons = int(rtree_evt.GetBranch('nPhotons').GetLeaf('nPhotons').GetValue())

        # preallocate empty array to build default hittime-histogram:
        hittime_array = []

        # loop over every photon in the event:
        for index2 in range(n_photons):

            # get nPE for this photon:
            n_pe = int(rtree_evt.GetBranch('nPE').GetLeaf('nPE').GetValue(index2))
            # check, if photon produces only 1 PE:
            if n_pe != 1:
                print("{1:d} PE for 1 photon in event {0:d} in file user_atmoNC_{2:d}.root"
                      .format(event_id, n_pe, index))

            # get the pmtID of the hit PMT:
            pmtID = int(rtree_evt.GetBranch('pmtID').GetLeaf('pmtID').GetValue(index2))

            # get hittime of this photon:
            hit_time = float(rtree_evt.GetBranch('hitTime').GetLeaf('hitTime').GetValue(index2))

            # get position of the PMT with specific pmtID (pmtID is ascending number from 0 to 17738 (17739 large PMTs)
            # and from 300000 to 336571 (36572 small PMTs)).
            # For large PMTs -> For 20inch PMTs, the pmtID is equal to index of x,y,z_pos_pmt array.
            # For small PMTs -> For 3inch PMTs, the pmtID - (300000 - 17739) is equal to index of x,y,z_pos_pmt array.
            # check if PMT is 20 inch or 3inch (pmtID < 50000 means 20inch PMT):
            if pmtID < 50000:
                # 20inch PMT:
                # get PMT position in mm from arrays:
                x_pmt = x_pos_pmt[pmtID]
                y_pmt = y_pos_pmt[pmtID]
                z_pmt = z_pos_pmt[pmtID]
            else:
                # 3inch PMT:
                # calculate index of pos_pmt array that correspond to pmtID of 3inch PMTs (for example:
                # first small PMT: 300000-282261 = 17739, last small PMT: 336571-282261 = 54310)
                index_3inch = pmtID - 282261
                # get PMT position in mm from arrays:
                x_pmt = x_pos_pmt[index_3inch]
                y_pmt = y_pos_pmt[index_3inch]
                z_pmt = z_pos_pmt[index_3inch]

            # calculate distance between reconstructed position of event and position of PMT (in mm):
            distance_tof = np.sqrt((x_reco_arr[index1] - x_pmt)**2 + (y_reco_arr[index1] - y_pmt)**2 +
                                   (z_reco_arr[index1] - z_pmt)**2)

            # calculate time of flight in ns:
            time_of_flight = distance_tof / c_effective

            """ time resolution of PMT: """
            # get time resolution of PMT with specific pmtID (pmtID is ascending number from 0 to 17738 (17739 large
            # PMTs)) -> For 20inch PMTs, the pmtID is equal to index of sigma_time_20inch array.
            # check if PMT is 20 inch or 3inch (pmtID < 50000 means 20inch PMT):
            if pmtID < 50000:
                # 20inch PMT:
                # get time resolution (sigma) of PMT in ns from array:
                sigma_pmt = sigma_time_20inch[pmtID]

            else:
                # 3inch PMT:
                sigma_pmt = sigma_time_3inch

            # consider time resolution of PMT by generating normal distributed random number with mu = hit_time and
            # sigma = sigma_pmt (only the hit_time at the PMT must be smeared, not the time-of-flight):
            hittime_tts = np.random.normal(hit_time, sigma_pmt)

            """ calculate the 'real' hittime of the photon in ns: """
            hittime_real = hittime_tts - time_of_flight
            if hittime_real < min_time:
                print("------")
                print(hittime_real)
                print(pmtID)
                print(sigma_pmt)

            # append hittime to array:
            hittime_array.append(hittime_real)

        """ analyze prompt signal: """
        # build histogram, where hittimes are saved:
        # set bin-edges of hittime histogram in ns:
        bins_hittime = np.arange(min_time, max_time + 2 * binwidth, binwidth)
        # build hittime histogram:
        npe_per_hittime, bin_edges_hittime = np.histogram(hittime_array, bins_hittime)

        # get index of bins_hittime corresponding to min_time (should be index = 0):
        index_min_hittime_prompt = 0

        # Where does prompt signal end?
        # get index of bins_hittime corresponding to time_limit_prompt
        index_time_limit_prompt = int((time_limit_prompt + np.abs(min_time)) / binwidth)
        # check if npe_per_hittime (and the following two bins) are 0 for this index:
        if (npe_per_hittime[index_time_limit_prompt] == npe_per_hittime[index_time_limit_prompt+1]
                == npe_per_hittime[index_time_limit_prompt+2] == 0):
            # prompt signal already 0:
            index_max_hittime_prompt = index_time_limit_prompt
        else:
            # prompt signal not yet 0.
            # loop over npe_per_hittime from index_time_limit_prompt until npe_per_hittime (and the following two bins)
            # are 0:
            for index3 in range(index_time_limit_prompt, index_time_limit_prompt+200):
                if npe_per_hittime[index3] == npe_per_hittime[index3+1] == npe_per_hittime[index3+2] == 0:
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
        quenched_deposit_energy = NC_background_functions.conversion_npe_to_evis(number_pe_prompt)

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
        # print("analyze delayed signal")
        # get index of bin_edges_hittime corresponding to the end of the prompt signal window:
        index_min_hittime_del = index_max_hittime_prompt

        # calculate nPE as function of hittime only for delayed time window (from min_hittime_del to end):
        npe_per_hittime_del = npe_per_hittime[index_min_hittime_del:]
        # bin edges of hittime histogram only for delayed time window:
        bins_hittime_del = bin_edges_hittime[index_min_hittime_del:-1]

        # get the minimum and maximum time of the delayed signal time window in ns:
        min_time_delayed = bins_hittime_del[0]
        max_time_delayed = bins_hittime_del[-1] + binwidth

        # preallocate number of possible delayed signals:
        number_delayed_signal = 0
        # preallocate first index of npe_per_hittime_del:
        index_first_del = 0
        # preallocate number of pe of delayed signal:
        number_pe_delayed = 0

        # analyze npe_per_hittime_del for possible delayed signals. As long as number_delayed_signal<2 and as long as
        # index has not reached the end of npe_per_hittime_del, check event for possible delayed signals
        while number_delayed_signal < 2 and index_first_del < len(npe_per_hittime_del):

            number_delayed, index_first_del, num_pe_delayed, begin_pulse, end_pulse = \
                NC_background_functions.analyze_delayed_signal(npe_per_hittime_del, bins_hittime_del, index_first_del,
                                                               threshold1_del, threshold2_del, min_PE_delayed,
                                                               max_PE_delayed, event_id)

            number_delayed_signal += number_delayed
            number_pe_delayed += num_pe_delayed

        if number_delayed_signal != 1:
            # 0 or more than one delayed signals that pass the delayed energy cut (-> no or more than one neutron
            # capture on hydrogen):
            number_rejected_delayed_energy_cut += 1
            # append npe of delayed signal, that wouldn't pass delayed energy cut to array:
            number_pe_delayed_rejected_array = np.append(number_pe_delayed_rejected_array, number_pe_delayed)

        # append number of pe of delayed signal to array:
        number_pe_total_del = np.append(number_pe_total_del, number_pe_delayed)

        # check number_delayed_signal:
        if number_delayed_signal == 0:
            print("---------------ERROR: no delayed signal in event {0:d}".format(event_id))
            h1 = plt.figure(1)
            plt.step(bins_hittime_del, npe_per_hittime_del, label="nPE of peak = {0:.0f}".format(number_pe_delayed))
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
                   header="Number of pe as function of the hittime of the prompt signal (time-of-flight correction "
                          "and TTS smearing) of file user_atmoNC_{0:d}.root,"
                          "\npreselected event {1:d} (analyzed with script prompt_signal_preselected_evts.py, {2}):"
                          "\ntime window of hittime: from {3:.3f} ns to {4:.3f} ns with bin-width = {5:0.3f} ns,"
                          "\nEnergy cut on prompt signal is applied: {6:0.1f} MeV <= E_vis <= {7:0.1f} MeV,"
                          "\nConversion function E_vis = 0.0007475 * nPE:"
                   .format(index, event_id, now, min_time_prompt, max_time_prompt, binwidth, min_energy, max_energy))

    # save array number_pe_total to txt file:
    np.savetxt(output_path + "number_pe_file{0:d}.txt".format(index), number_pe_total, fmt="%i",
               header="Total number of pe of prompt signal for every preselected event (that pass prompt energy cut) in"
                      " file user_atmoNC_{0:d}.root\n"
                      "({1}). Number of pe is analyzed with script prompt_signal_preselected_evts.py.\n"
                      "Time window of prompt signal is defined from {2:0.2f} ns to around {3:0.2f} ns (number of events"
                      " = {4:d}).\n"
                      "Energy cut on prompt signal is applied: {5:0.1f} MeV <= E_vis <= {6:0.1f} MeV,\n"
                      "Conversion function E_vis = 0.0007475 * nPE:"
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
print("number of events (of preselect. evts), that would be rejected by delayed energy cut = {0:d}"
      .format(number_rejected_delayed_energy_cut))
