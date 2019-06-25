""" Script to read the root files from positron simulation and saved their hittime distributions to txt file.
    These hittime distributions can then be analyzed further with pulse_shape_analysis.py as reference to the prompt
    signal of IBD-like NC events (to compare hittime distributions of positrons and NC events).

"""
import datetime
import os
import ROOT
import sys
from array import array
import NC_background_functions
import numpy as np
from matplotlib import pyplot as plt


def get_hittime_from_rootfile_fixenergy(input_path, output_path, first_file, last_file, num_evts, kin_energy, min_t,
                                        max_t, t_limit, bin_width, radius, now):
    """
    function to read events of the root files of positron simulation and save hittime distribution to png and txt file.

    function is used for positron with fixed kinetic energy.

    Prompt signal is analyzed like in prompt_signal_preselected_evts.py

    :param input_path: path, where root-files are saved
    :param output_path: path, where png and txt files of hittime distribution are saved
    :param first_file: number of the first file to read
    :param last_file: number of the last file to read
    :param num_evts: number of events per file
    :param kin_energy: kinetic energy of positron in MeV
    :param min_t: minimum of time window of whole signal in ns
    :param max_t: maximum of time window of whole signal in ns
    :param t_limit: time in ns, where prompt signal should be 0
    :param bin_width: bin-width of hittime-distribution in ns
    :param radius: cut radius in mm
    :param now: actual time
    :return:
    """
    # preallocate array, where total nPE of prompt signal per event of all files is stored:
    num_pe_total = np.array([])
    # number of events that are analyzed (pass volume cut):
    number_analyzed = 0

    # loop over root files with positron simulation:
    for index in range(first_file, last_file + 1, 1):
        # load user_positron_{}.root file:
        rfile = ROOT.TFile(input_path + "user_positron_{0:d}_MeV_{1:d}.root".format(kin_energy, index))
        print("... read {0}...".format(rfile))

        # get the "evt"-TTree from the TFile:
        rtree_evt = rfile.Get("evt")
        # get geninfo tree from TFile:
        rtree_geninfo = rfile.Get("geninfo")
        # get the number of events in the 'evt' Tree:
        num_events_evt = rtree_evt.GetEntries()
        # check number of events:
        if num_events_evt != num_evts:
            sys.exit("ERROR: number of events in root file ({0:d}) != {1:d}"
                     .format(num_events_evt, num_evts))

        # loop over the events:
        for event in range(num_evts):
            # check volume cut:
            rtree_geninfo.GetEntry(event)
            # get number of initial particles:
            n_init_particles = int(rtree_geninfo.GetBranch('nInitParticles').GetLeaf('nInitParticles').GetValue())
            if n_init_particles != 1:
                sys.exit("ERROR: more than 1 initial particles in event {0:d}".format(event))
            # get initial x, y, z position:
            x_init = float(rtree_geninfo.GetBranch('InitX').GetLeaf('InitX').GetValue())
            y_init = float(rtree_geninfo.GetBranch('InitY').GetLeaf('InitY').GetValue())
            z_init = float(rtree_geninfo.GetBranch('InitZ').GetLeaf('InitZ').GetValue())
            # calculate distance to center:
            r_init = np.sqrt(x_init**2 + y_init**2 + z_init**2)

            if r_init >= radius:
                print("file {0:d}, event = {1:d}: r_init = {2:0.2f} mm".format(index, event, r_init))
                continue

            # get event of 'evt'-tree:
            rtree_evt.GetEntry(event)
            # get evtID of the tree and compare with event:
            evt_id = int(rtree_evt.GetBranch('evtID').GetLeaf('evtID').GetValue())
            if evt_id != event:
                sys.exit("ERROR: evtID of tree ({0:d}) != {1:d}".format(evt_id, event))

            print("\nanalyze event {0:d}".format(evt_id))

            # increment number_analyzed:
            number_analyzed += 1

            # get number of photons of this event:
            n_photons = int(rtree_evt.GetBranch('nPhotons').GetLeaf('nPhotons').GetValue())

            # preallocate empty array to build default hittime-histogram:
            hittime_array = []

            # loop over every photon in the event:
            for index1 in range(n_photons):

                # get PMT ID, where photon is absorbed:
                pmt_id = int(rtree_evt.GetBranch('pmtID').GetLeaf('pmtID').GetValue(index1))

                # only 20 inch PMTs (PMT ID of 20 inch PMTs are below 21000, PMT ID of 3 inch PMTs start at 290000):
                if pmt_id < 25000:
                    # get nPE for this photon:
                    n_pe = int(rtree_evt.GetBranch('nPE').GetLeaf('nPE').GetValue(index1))
                    # check, if photon produces only 1 PE:
                    if n_pe != 1:
                        print("{1:d} PE for 1 photon in event {0:d} in file user_positron_{3:d}_MeV_{2:d}.root"
                              .format(evt_id, n_pe, index, kin_energy))

                    # get hittime of this photon:
                    hit_time = float(rtree_evt.GetBranch('hitTime').GetLeaf('hitTime').GetValue(index1))

                    # append hittime to array:
                    hittime_array.append(hit_time)

                else:
                    continue

            """ analyze prompt signal: """
            # build histogram, where hittimes are saved:
            # set bin-edges of hittime histogram in ns:
            bins_hittime = np.arange(min_t, max_t + 2 * bin_width, bin_width)
            # build hittime histogram:
            npe_per_hittime, bin_edges_hittime = np.histogram(hittime_array, bins_hittime)

            # get index of bins_hittime corresponding to min_time (should be index = 0):
            index_min_hittime_prompt = int(min_t / bin_width)

            # Where does prompt signal end?
            # get index of bins_hittime corresponding to t_limit:
            index_time_limit_prompt = int(t_limit / bin_width)
            # check if npe_per_hittime is 0 for this index:
            if npe_per_hittime[index_time_limit_prompt] == 0:
                # prompt signal already 0:
                index_max_hittime_prompt = index_time_limit_prompt
            else:
                # prompt signal not yet 0.
                # loop over npe_per_hittime from index_time_limit_prompt until npe_per_hittime is 0:
                for index2 in range(index_time_limit_prompt, index_time_limit_prompt+200):
                    if npe_per_hittime[index2] == 0:
                        index_max_hittime_prompt = index2
                        break

            # calculate nPE as function of hittime only for prompt time window (from min_hittime_prompt to
            # max_hittime_prompt+1 (last index should be included)):
            npe_per_hittime_prompt = npe_per_hittime[index_min_hittime_prompt:index_max_hittime_prompt+1]
            # bin edges of hittime histogram only for prompt time window:
            bins_hittime_prompt = bin_edges_hittime[index_min_hittime_prompt:index_max_hittime_prompt+1]

            # get the minimum and maximum time of the prompt signal time window in ns:
            min_time_prompt = bins_hittime_prompt[0]
            max_time_prompt = bins_hittime_prompt[-1]

            # sum up the values of npe_per_hittime_prompt to get the total number of pe of the prompt signal:
            number_pe_prompt = np.sum(npe_per_hittime_prompt)

            # append number of pe to array:
            num_pe_total = np.append(num_pe_total, number_pe_prompt)

            h1 = plt.figure(1)
            plt.step(bins_hittime_prompt, npe_per_hittime_prompt, label="number of pe = {0:d}".format(number_pe_prompt))
            plt.xlabel("hit-time in ns")
            plt.ylabel("number of p.e. per bin (bin-width = {0:0.2f} ns)".format(bin_width))
            plt.title("Hit-time distribution of prompt time window of event {0:d}".format(evt_id))
            plt.xlim(xmin=min_time_prompt, xmax=max_time_prompt)
            plt.legend()
            plt.grid()
            plt.savefig(output_path + "file{1:d}_evt{0:d}_positron_{2:d}_MeV.png".format(evt_id, index, kin_energy))
            plt.close()
            # plt.show()

            # save npe_per_hittime_prompt to txt file:
            # build list, where 0th entry is start-hittime in ns, 1st entry is last-hittime in ns, 2nd entry is binwidth
            # in ns and the following entries are nPE of each hittime-bin of prompt signal:
            npe_per_hittime_prompt_save = [min_time_prompt, max_time_prompt, bin_width]
            npe_per_hittime_prompt_save.extend(npe_per_hittime_prompt)
            np.savetxt(output_path + "file{0:d}_evt{1:d}_positron_{2:d}_MeV.txt".format(index, evt_id, kin_energy),
                       npe_per_hittime_prompt_save, fmt='%1.2f',
                       header="Number of pe as function of the hittime of the prompt positron signal of file "
                              "user_positron_{6:d}_MeV_{0:d}.root,"
                              "\nevent = {1:d}, (analyzed with hittime_distribution_positron.py, {2}):"
                              "\ntime window of hittime: from {3:.3f} ns to {4:.3f} ns with bin-width = {5:0.3f} ns:"
                       .format(index, evt_id, now, min_time_prompt, max_time_prompt, bin_width, kin_energy))

    return num_pe_total, number_analyzed


def get_hittime_from_rootfile(input_path, output_path, first_file, last_file, num_evts, min_t, max_t, t_limit,
                              bin_width, radius, now):
    """
    function to read events of the root files of positron simulation and save hittime distribution to png and txt file

    function is used for positron with uniformly distributed kinetic energy (not fixed energy).

    Prompt signal is analyzed like in prompt_signal_preselected_evts.py

    :param input_path: path, where root-files are saved
    :param output_path: path, where png and txt files of hittime distribution are saved
    :param first_file: number of the first file to read
    :param last_file: number of the last file to read
    :param num_evts: number of events per file
    :param min_t: minimum of time window of whole signal in ns
    :param max_t: maximum of time window of whole signal in ns
    :param t_limit: time in ns, where prompt signal should be 0
    :param bin_width: bin-width of hittime-distribution in ns
    :param radius: cut radius in mm
    :param now: actual time
    :return:
    """
    # preallocate array, where total nPE of prompt signal per event of all files is stored:
    num_pe_total = np.array([])
    # number of events that are analyzed (pass volume cut):
    number_analyzed = 0

    # loop over root files with positron simulation:
    for index in range(first_file, last_file + 1, 1):
        # load user_positron_{}.root file:
        rfile = ROOT.TFile(input_path + "user_positron_{0:d}.root".format(index))
        print("... read {0}...".format(rfile))

        # get the "evt"-TTree from the TFile:
        rtree_evt = rfile.Get("evt")
        # get geninfo tree from TFile:
        rtree_geninfo = rfile.Get("geninfo")
        # get the number of events in the 'evt' Tree:
        num_events_evt = rtree_evt.GetEntries()
        # check number of events:
        if num_events_evt != num_evts:
            sys.exit("ERROR: number of events in root file ({0:d}) != {1:d}"
                     .format(num_events_evt, num_evts))

        # loop over the events:
        for event in range(num_evts):
            # check volume cut:
            rtree_geninfo.GetEntry(event)
            # get number of initial particles:
            n_init_particles = int(rtree_geninfo.GetBranch('nInitParticles').GetLeaf('nInitParticles').GetValue())
            if n_init_particles != 1:
                sys.exit("ERROR: more than 1 initial particles in event {0:d}".format(event))
            # get initial x, y, z position:
            x_init = float(rtree_geninfo.GetBranch('InitX').GetLeaf('InitX').GetValue())
            y_init = float(rtree_geninfo.GetBranch('InitY').GetLeaf('InitY').GetValue())
            z_init = float(rtree_geninfo.GetBranch('InitZ').GetLeaf('InitZ').GetValue())
            # calculate distance to center:
            r_init = np.sqrt(x_init**2 + y_init**2 + z_init**2)

            if r_init >= radius:
                print("file {0:d}, event = {1:d}: r_init = {2:0.2f} mm".format(index, event, r_init))
                continue

            # get event of 'evt'-tree:
            rtree_evt.GetEntry(event)
            # get evtID of the tree and compare with event:
            evt_id = int(rtree_evt.GetBranch('evtID').GetLeaf('evtID').GetValue())
            if evt_id != event:
                sys.exit("ERROR: evtID of tree ({0:d}) != {1:d}".format(evt_id, event))

            print("\nanalyze event {0:d}".format(evt_id))

            # increment number_analyzed:
            number_analyzed += 1

            # get number of photons of this event:
            n_photons = int(rtree_evt.GetBranch('nPhotons').GetLeaf('nPhotons').GetValue())

            # preallocate empty array to build default hittime-histogram:
            hittime_array = []

            # loop over every photon in the event:
            for index1 in range(n_photons):

                # get PMT ID, where photon is absorbed:
                pmt_id = int(rtree_evt.GetBranch('pmtID').GetLeaf('pmtID').GetValue(index1))

                # only 20 inch PMTs (PMT ID of 20 inch PMTs are below 21000, PMT ID of 3 inch PMTs start at 290000):
                if pmt_id < 25000:
                    # get nPE for this photon:
                    n_pe = int(rtree_evt.GetBranch('nPE').GetLeaf('nPE').GetValue(index1))
                    # check, if photon produces only 1 PE:
                    if n_pe != 1:
                        print("{1:d} PE for 1 photon in event {0:d} in file user_positron_{2:d}.root"
                              .format(evt_id, n_pe, index))

                    # get hittime of this photon:
                    hit_time = float(rtree_evt.GetBranch('hitTime').GetLeaf('hitTime').GetValue(index1))

                    # append hittime to array:
                    hittime_array.append(hit_time)

                else:
                    continue

            """ analyze prompt signal: """
            # build histogram, where hittimes are saved:
            # set bin-edges of hittime histogram in ns:
            bins_hittime = np.arange(min_t, max_t + 2 * bin_width, bin_width)
            # build hittime histogram:
            npe_per_hittime, bin_edges_hittime = np.histogram(hittime_array, bins_hittime)

            # get index of bins_hittime corresponding to min_time (should be index = 0):
            index_min_hittime_prompt = int(min_t / bin_width)

            # Where does prompt signal end?
            # get index of bins_hittime corresponding to t_limit:
            index_time_limit_prompt = int(t_limit / bin_width)
            # check if npe_per_hittime is 0 for this index:
            if npe_per_hittime[index_time_limit_prompt] == 0:
                # prompt signal already 0:
                index_max_hittime_prompt = index_time_limit_prompt
            else:
                # prompt signal not yet 0.
                # loop over npe_per_hittime from index_time_limit_prompt until npe_per_hittime is 0:
                for index2 in range(index_time_limit_prompt, index_time_limit_prompt+200):
                    if npe_per_hittime[index2] == 0:
                        index_max_hittime_prompt = index2
                        break

            # calculate nPE as function of hittime only for prompt time window (from min_hittime_prompt to
            # max_hittime_prompt+1 (last index should be included)):
            npe_per_hittime_prompt = npe_per_hittime[index_min_hittime_prompt:index_max_hittime_prompt+1]
            # bin edges of hittime histogram only for prompt time window:
            bins_hittime_prompt = bin_edges_hittime[index_min_hittime_prompt:index_max_hittime_prompt+1]

            # get the minimum and maximum time of the prompt signal time window in ns:
            min_time_prompt = bins_hittime_prompt[0]
            max_time_prompt = bins_hittime_prompt[-1]

            # sum up the values of npe_per_hittime_prompt to get the total number of pe of the prompt signal:
            number_pe_prompt = np.sum(npe_per_hittime_prompt)

            # append number of pe to array:
            num_pe_total = np.append(num_pe_total, number_pe_prompt)

            h1 = plt.figure(1)
            plt.step(bins_hittime_prompt, npe_per_hittime_prompt, label="number of pe = {0:d}".format(number_pe_prompt))
            plt.xlabel("hit-time in ns")
            plt.ylabel("number of p.e. per bin (bin-width = {0:0.2f} ns)".format(bin_width))
            plt.title("Hit-time distribution of prompt time window of event {0:d}".format(evt_id))
            plt.xlim(xmin=min_time_prompt, xmax=max_time_prompt)
            plt.legend()
            plt.grid()
            plt.savefig(output_path + "file{1:d}_evt{0:d}_positron.png".format(evt_id, index))
            plt.close()
            # plt.show()

            # save npe_per_hittime_prompt to txt file:
            # build list, where 0th entry is start-hittime in ns, 1st entry is last-hittime in ns, 2nd entry is binwidth
            # in ns and the following entries are nPE of each hittime-bin of prompt signal:
            npe_per_hittime_prompt_save = [min_time_prompt, max_time_prompt, bin_width]
            npe_per_hittime_prompt_save.extend(npe_per_hittime_prompt)
            np.savetxt(output_path + "file{0:d}_evt{1:d}_positron.txt".format(index, evt_id),
                       npe_per_hittime_prompt_save, fmt='%1.2f',
                       header="Number of pe as function of the hittime of the prompt positron signal of file "
                              "user_positron_{0:d}.root,"
                              "\nevent = {1:d}, (analyzed with hittime_distribution_positron.py, {2}):"
                              "\ntime window of hittime: from {3:.3f} ns to {4:.3f} ns with bin-width = {5:0.3f} ns:"
                       .format(index, evt_id, now, min_time_prompt, max_time_prompt, bin_width))

    return num_pe_total, number_analyzed


# get the date and time, when the script was run:
date = datetime.datetime.now()
NOW = date.strftime("%Y-%m-%d %H:%M")

""" define time window and bin width: """
# set time window of whole signal in ns:
min_time = 0
max_time = 1000000
# set time in ns, where the prompt signal should be 0:
time_limit_prompt = 500
# Set bin-width of hittime histogram in ns:
binwidth = 5.0

""" set parameter for volume cut (must be the same like for NC events): """
# cut-radius in mm:
radius_cut = 16000

# path, where root files of positron simulation are saved:
input_path_positron = "/local/scratch1/pipc51/astro/blum/positron_output/"
# path, where hittime distributions (png and txt) are saved:
output_path_positron = "/home/astro/blum/juno/atmoNC/data_NC/output_PSD/positron_hittime/"
# event per root file:
number_evts_positron = 10

""" analyze 10 MeV positron hittime distributions """
# first file of positron simulation for 10 MeV:
first_file_10 = 0
# last file of positron simulation for 10 MeV:
last_file_10 = 99
# kinetic energy of positrons in MeV:
energy_positron_10 = 10
# total number of positron events with 10 MeV:
number_evts_total_10 = (last_file_10 - first_file_10 + 1) * number_evts_positron

# # get hittime distribution (return values: number of pe of each event, number of events that are analyzed):
# number_pe_10, number_analyzed_10 = get_hittime_from_rootfile_fixenergy(input_path_positron, output_path_positron,
#                                                                        first_file_10, last_file_10,
#                                                                        number_evts_positron, energy_positron_10,
#                                                                        min_time, max_time, time_limit_prompt,
#                                                                        binwidth, radius_cut, NOW)
#
# # check the distribution of total nPE for all events with positrons of 10 MeV kinetic energy:
# # calculate mean of number of PE:
# mean_10_MeV = np.mean(number_pe_10)
# # display number_pe_10 in histogram:
# h1 = plt.figure(1, figsize=(15, 8))
# n_PE_10MeV, bins_10MeV, patches1 = plt.hist(number_pe_10, align='mid', bins=100,
#                                             label="{0:d} positrons with kinetic energy = {1:d} MeV"
#                                             .format(number_evts_total_10, energy_positron_10))
# plt.vlines(mean_10_MeV, ymin=0, ymax=max(n_PE_10MeV), label="mean = {0:0.2f} nPE".format(mean_10_MeV))
# plt.xlabel("number of PE (per positron)", fontsize=13)
# plt.ylabel("entries per bin", fontsize=13)
# plt.title("Number of PE of 10 MeV positrons", fontsize=18)
# plt.legend()
# plt.grid()
# plt.savefig(output_path_positron + "hist_nPE_10_MeV.png")
# plt.close()

""" analyze 100 MeV positron hittime distributions """
# first file of positron simulation for 100 MeV:
first_file_100 = 0
# last file of positron simulation for 100 MeV:
last_file_100 = 99
# kinetic energy of positrons in MeV:
energy_positron_100 = 100
# total number of positron events with 10 MeV:
number_evts_total_100 = (last_file_100 - first_file_100 + 1) * number_evts_positron

# number_pe_100, number_analyzed_100 = get_hittime_from_rootfile_fixenergy(input_path_positron, output_path_positron,
#                                                                          first_file_100, last_file_100,
#                                                                          number_evts_positron, energy_positron_100,
#                                                                          min_time, max_time, time_limit_prompt,
#                                                                          binwidth, radius_cut, NOW)
#
# # check the distribution of total nPE for all events with positrons of 100 MeV kinetic energy:
# # calculate mean of number of PE:
# mean_100_MeV = np.mean(number_pe_100)
# # display number_pe_100 in histogram:
# h2 = plt.figure(2, figsize=(15, 8))
# n_PE_100MeV, bins_100MeV, patches2 = plt.hist(number_pe_100, align='mid', bins=100,
#                                               label="{0:d} positrons with kinetic energy = {1:d} MeV"
#                                               .format(number_evts_total_100, energy_positron_100))
# plt.vlines(mean_100_MeV, ymin=0, ymax=max(n_PE_100MeV), label="mean = {0:0.2f} nPE".format(mean_100_MeV))
# plt.xlabel("number of PE (per positron)", fontsize=13)
# plt.ylabel("entries per bin", fontsize=13)
# plt.title("Number of PE of 100 MeV positrons", fontsize=18)
# plt.legend()
# plt.grid()
# plt.savefig(output_path_positron + "hist_nPE_100_MeV.png")
# plt.close()

""" analyze positron hittime distribution for kinetic energy from 10 MeV to 100 MeV (uniformly distributed): """
# first file of positron simulation:
first_file_positron = 0
# last file of positron simulation:
last_file_positron = 99
# number of events per file:
number_evts_per_file = 100
# total number of positron events:
number_evts_total = (last_file_positron - first_file_positron + 1) * number_evts_per_file

number_pe_positron, number_analyzed_positron = get_hittime_from_rootfile(input_path_positron, output_path_positron,
                                                                         first_file_positron, last_file_positron,
                                                                         number_evts_per_file, min_time, max_time,
                                                                         time_limit_prompt, binwidth, radius_cut, NOW)

# check the distribution of total nPE for all events with positrons of kinetic energy between 10 and 100 MeV:
# display number_pe_positron in histogram:
h3 = plt.figure(2, figsize=(15, 8))
n_PE_positron, bins_positron, patches3 = plt.hist(number_pe_positron, align='mid', bins=100,
                                                  label="{0:d} positrons with kinetic energy between {1:d} MeV and "
                                                        "{2:d} MeV"
                                                  .format(number_evts_total, energy_positron_10, energy_positron_100))
plt.xlabel("number of PE (per positron)", fontsize=13)
plt.ylabel("entries per bin", fontsize=13)
plt.title("Number of PE of positrons with kinetic energy from 10 MeV to 100 MeV", fontsize=18)
plt.legend()
plt.grid()
plt.savefig(output_path_positron + "hist_nPE_positron.png")
plt.close()

print("total number of events = {0:d}".format(number_evts_total))
print("number of analyzed events = {0:d}".format(number_analyzed_positron))

