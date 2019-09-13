""" Script to analyze the volume cut and its efficiency (either of positron or NC events):

    Version 2 (v2) means, that:
    - the volume cut efficiency is calculated independently of the other cuts.
    - filenumber and evtID of each event passing the volume cut is stored in file
    - all values of R_reco are stored in a histogram with the fiducial volume cut parameter R_FV

    1.  Set the cut parameter (R_FV)

    2.  read all NC events (user_atmoNC_0.root to user_atmoNC_999.root, each file contains 100 events.
        Therefore, 100000 events are analyzed.)
        or read all IBD events (user_IBD_hepevt_0.root to user_IBD_hepevt_199.root, each file contains 100
        events. Therefore, 20000 events are analyzed.)

    3.  analyze the events: do volume cut on the reconstructed position of the event

    4.  store filenumber and event number of the events, that pass the volume cut on recon. position into array and
        save it to file

    5.  store R_reco of all events in histogram and saved histogram with R_FV

    6.  analyze events also for the MC truth (initial) position

    7.  store also filenumber and event number of the events, that pass the volume cut on initial position (MC truth)
        into array and save it to file

    8.  calculate number of leak-in and leak-out events

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

""" file information: """
# first file to be read:
start_number = 0
# last file to be read:
stop_number = 199
# number of entries in the input files:
Number_entries_input = 100
# set string, that define, if positrons or NC events are analyzed:
event_type = "IBD"
# set the path of the input root files:
if event_type == "atmoNC":
    input_path = "/local/scratch1/pipc51/astro/blum/detsim_output_data/user_atmoNC_"
elif event_type == "IBD":
    input_path = "/local/scratch1/pipc51/astro/blum/IBD_hepevt/user_IBD_hepevt_"

# set the path of the output, where the txt file with visible energy of events that pass cuts is saved:
if event_type == "atmoNC":
    output_path = "/home/astro/blum/juno/atmoNC/data_NC/output_detsim_v2/"
elif event_type == "IBD":
    output_path = "/home/astro/blum/juno/IBD_events/"

""" define parameters depending on the cuts: """
# fiducial volume cut in mm:
R_cut_mm = 16000.0

# number of analyzed events (total number of events):
number_total_events = (stop_number + 1 - start_number) * Number_entries_input

""" preallocate variables: """
# number of events without initial particles:
number_without_particles = 0
# number of events before the volume cut (number of events, that were analyzed):
number_analyzed = 0
# number of events that pass volume cut on prompt recon. position:
number_volume_pass_prompt_reco = 0
# number of events that are rejected by volume cut on prompt recon. position:
number_volume_rejected_prompt_reco = 0
# number of events that pass volume cut on prompt initial (MC truth) position:
number_volume_pass_prompt_init = 0
# number of events that are rejected by volume cut on prompt initial (MC truth) position:
number_volume_rejected_prompt_init = 0
# number of leak-in events (initial position outside, but recon. position inside fiducial volume):
number_leak_in = 0
# number of leak-out events (initial position inside, but recon. position outside fiducial volume):
number_leak_out = 0
# array, where filenumber of events that pass the cut (reconstructed (real) position) are stored:
array_filenumber = []
# array, where corresponding evtID of events that pass the cut (resonctructed (real) position) are stored:
array_evtID = []
# array, where all recon. positions are stored to build a histogram:
array_R_reco = []
# array, where filenumber of events that pass the cut (MC truth/initial (ideal) position) are stored:
array_filenumber_MCtruth = []
# array, where corresponding evtID of events that pass the cut (MC truth/initial (ideal) position) are stored:
array_evtID_MCtruth = []

# loop over the files that are read:
for filenumber in range(start_number, stop_number+1):

    # file name of the input file:
    input_name = input_path + "{0:d}.root".format(filenumber)

    # load user_atmoNC_index.root file:
    rfile = ROOT.TFile(input_name)
    # get the "geninfo"-TTree from the TFile:
    rtree_geninfo = rfile.Get("geninfo")
    # get the "prmtrkdep"-TTree from the TFile:
    rtree_prmtrkdep = rfile.Get("prmtrkdep")

    # get the number of events in the geninfo Tree:
    number_events_geninfo = rtree_geninfo.GetEntries()
    # get the number of events in the prmtrkdep tree:
    number_events_prmtrkdep = rtree_prmtrkdep.GetEntries()
    if number_events_geninfo == number_events_prmtrkdep:
        number_events = number_events_geninfo
    else:
        sys.exit("ERROR: number of events in the Trees are NOT equal!!")

    # check if number_events is equal to number_entries_input (if not, the detector simulation was incorrect!!):
    if number_events != Number_entries_input:
        sys.exit("ERROR: number of events {0:d} are not equal to {1:d} -> Detector Simulation not correct!"
                 .format(number_events, Number_entries_input))

    # loop over every event in the file:
    for event in range(number_events):

        # sum of Qedep of all initial particles (in MeV):
        Qedep_sum = 0

        """ get reconstructed position of event """
        # get the current event in the TTree:
        rtree_geninfo.GetEntry(event)
        # get the value of the event ID:
        event_id = int(rtree_geninfo.GetBranch('evtID').GetLeaf('evtID').GetValue())

        # get the value of the number of initial particles:
        nInitParticles_geninfo = int(rtree_geninfo.GetBranch('nInitParticles').GetLeaf('nInitParticles').GetValue())

        # check if there are initial particles in this event:
        if nInitParticles_geninfo < 1:
            # there are NO particles in the event:
            number_without_particles += 1
            # go to next event
            continue

        # increment number_analyzed:
        number_analyzed += 1

        # get initial position of the first initial particle:
        init_x = float(rtree_geninfo.GetBranch('InitX').GetLeaf('InitX').GetValue(0))
        init_y = float(rtree_geninfo.GetBranch('InitY').GetLeaf('InitY').GetValue(0))
        init_z = float(rtree_geninfo.GetBranch('InitZ').GetLeaf('InitZ').GetValue(0))

        # check, if initial positions of the other initial particles are equal to position of first initial:
        for index in range(1, nInitParticles_geninfo):
            # get initial x position:
            init_x_test = float(rtree_geninfo.GetBranch('InitX').GetLeaf('InitX').GetValue(index))
            # get initial y position:
            init_y_test = float(rtree_geninfo.GetBranch('InitY').GetLeaf('InitY').GetValue(index))
            # get initial z position:
            init_z_test = float(rtree_geninfo.GetBranch('InitZ').GetLeaf('InitZ').GetValue(index))

            if init_x != init_x_test or init_y != init_y_test or init_z != init_z_test:
                sys.exit("ERROR: initial positions are NOT equal (event = {0:d}, filenumber = {1:d})"
                         .format(event, filenumber))

        # calculate the initial position of the event:
        r_initial = np.sqrt(init_x**2 + init_y**2 + init_z**2)
        # do volume cut on initial position:
        if r_initial < R_cut_mm:
            # initial position passes the cut:
            number_volume_pass_prompt_init += 1
            # append filenumber and evtID of this event to the arrays:
            array_filenumber_MCtruth.append(filenumber)
            array_evtID_MCtruth.append(event_id)
        else:
            # initial position is rejected by cut:
            number_volume_rejected_prompt_init += 1

        # get current event in TTree:
        rtree_prmtrkdep.GetEntry(event)
        # get number of initial particles in prmtrkdep:
        nInitParticles_prmtrkdep = int(rtree_geninfo.GetBranch('nInitParticles').GetLeaf('nInitParticles').GetValue())

        # get quenched deposit energy of all initial particles in MeV:
        for index in range(nInitParticles_prmtrkdep):
            Qedep = float(rtree_prmtrkdep.GetBranch('Qedep').GetLeaf('Qedep').GetValue(index))
            Qedep_sum += Qedep

        # Calculate vertex smeared position.
        if Qedep_sum != 0:
            # Smear initial x,y and z position with function position_smearing(). (returns reconstructed position
            # in mm):
            x_reconstructed = NC_background_functions.position_smearing(init_x, Qedep_sum)
            y_reconstructed = NC_background_functions.position_smearing(init_y, Qedep_sum)
            z_reconstructed = NC_background_functions.position_smearing(init_z, Qedep_sum)

        else:
            # Qedep_sum = 0, use initial position:
            x_reconstructed = init_x
            y_reconstructed = init_y
            z_reconstructed = init_z

        # calculate reconstructed distance to detector center in mm:
        r_reconstructed = np.sqrt(x_reconstructed**2 + y_reconstructed**2 + z_reconstructed**2)

        """ volume cut on prompt recon. signal """
        # apply volume cut:
        if r_reconstructed < R_cut_mm:
            # recon. position passes volume cut:
            number_volume_pass_prompt_reco += 1
            # append filenumber and evtID of this event to the arrays:
            array_filenumber.append(filenumber)
            array_evtID.append(event_id)

        else:
            # recon. position is rejected by cut:
            number_volume_rejected_prompt_reco += 1

        """ store r_reconstructed in array to build histogram: """
        array_R_reco.append(r_reconstructed)

        """ check leak-in/leak-out efficiency: """
        if r_reconstructed < R_cut_mm and r_initial >= R_cut_mm:
            # leak-in event:
            number_leak_in += 1

        if r_reconstructed >= R_cut_mm and r_initial < R_cut_mm:
            # leak-out event:
            number_leak_out += 1

""" build histogram from array_R_reco: """
h1 = plt.figure(1, figsize=(15, 8))
bin_width = 500.0
bins_radius = np.arange(0.0, 18000+bin_width, bin_width)
values, bins, patches = plt.hist(array_R_reco, bins_radius, align='mid', histtype='step', linewidth='1.5',
                                 label='entries = {0:d},\n'
                                       'events passing cut = {1:d},\n'
                                       'events rejected by cut = {2:d}'
                                 .format(number_total_events, number_volume_pass_prompt_reco,
                                         number_volume_rejected_prompt_reco))
plt.vlines(R_cut_mm, ymin=0.0, ymax=(max(values)+max(values)*0.1), colors='k', linestyles='solid',
           label='fiducial volume cut $R_{FV}$ = '+'{0:.0f} mm'.format(R_cut_mm))
plt.xlim(xmin=0.0, xmax=18000)
plt.xlabel("distance to detector center in mm")
plt.ylabel("number of events per bin (bin-width = {0:.0f} mm)".format(bin_width))
plt.title("Reconstructed position of {0} events".format(event_type))
plt.legend(loc='upper left')
plt.grid()
plt.savefig(output_path + "histo_volume_cut_{0}_{1:.0f}mm.png".format(event_type, R_cut_mm))
plt.close()

""" save array_filenumber and array_evtID of reconstructed position (real) to txt file: """
np.savetxt(output_path + 'filenumber_evtID_volume_cut_{0}_{1:.0f}mm.txt'.format(event_type, R_cut_mm),
           np.c_[array_filenumber, array_evtID], fmt='%i',
           header='filenumber | evtID of events that pass volume cut')

""" save array_filenumber_MCtruth and array_evtID_MCtruth of initial/MC truth position (ideal) to txt file: """
np.savetxt(output_path + 'filenumber_evtID_volume_cut_MCtruth_{0}_{1:.0f}mm.txt'.format(event_type, R_cut_mm),
           np.c_[array_filenumber_MCtruth, array_evtID_MCtruth], fmt='%i',
           header='filenumber | evtID of events that pass volume cut')

""" save different number of events to txt file: """
np.savetxt(output_path + "numbers_volume_cut_{0}_{1:.0f}mm.txt".format(event_type, R_cut_mm),
           np.array([number_total_events, number_without_particles, number_analyzed,
                     number_volume_pass_prompt_reco, number_volume_rejected_prompt_reco,
                     number_volume_pass_prompt_init, number_volume_rejected_prompt_init,
                     number_leak_in, number_leak_out]), fmt='%i',
           header='number of events from analyze_volume_cut_v2.py ({0}):'
                  '\n{1} events are analyzed for R_FV = {2:.0f} mm.'
                  '\nValues below:'
                  '\nnumber of total events,'
                  '\nnumber of events without initial particles,'
                  '\nnumber of events before cut (these events were analyzed),'
                  '\nnumber of events that pass volume cut on prompt recon. position,'
                  '\nnumber of events that are rejected by volume cut on prompt recon. position,'
                  '\nnumber of events that pass volume cut on prompt initial (MC truth) position,'
                  '\nnumber of events that are rejected by volume cut on prompt initial (MC truth) position,'
                  '\nnumber of leak-in events (counted too much) (initial position outside, but recon. position inside '
                  'fiducial volume),'
                  '\nnumber of leak-out events (counted too less) (initial position inside, but recon. position '
                  'outside fiducial volume),'
                  .format(now, event_type, R_cut_mm))



