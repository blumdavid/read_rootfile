""" Script to check the efficiency of the fiducial volume cut for atmoNC events.

    The radius of JUNO detector is 17.7m, but only a fiducial volume of the detector is used for detection to suppress
    background from the surroundings (normally 16 m).

    The total volume cut efficiency (how many events are inside fiducial volume after vertex reconstruction) is done
    with script preselection_detsim_user.py. (73.1 % of all events inside 16 m (03.07.2019))

    With this script the leak-in/leak-out efficiency can be calculated. This means, how many events with initial
    position < 16 m are reconstructed outside 16 m (leak-out), and, how many events with initial position >= 16 m, are
    reconstructed inside 16 m (leak-in).

    In the yellow book, this efficiency (leak-out/leak-in) is 91.8 % for r = 17 m (page 39)
    (-> this means the number of leak-out events is smaller than the number of leak-in events
    -> there are more events inside fiducial volume because of reconstruction than in reality)

"""
import datetime
import numpy as np
import ROOT
import sys
import NC_background_functions


def check_volume_cut(input_path, number_entries_input, radius_cut):
    """

    :param input_path: file name with path to input root files from tut_detsim.py: user_atmoNC_{}.root
    :param number_entries_input:  number of entries, that the input files should have (integer), normally = 100
    :param radius_cut: radius, which defines the fiducial volume, in mm
    :return:
    """
    # load the ROOT file:
    rfile = ROOT.TFile(input_path)
    # get the "geninfo"-TTree from the TFile:
    rtree_geninfo = rfile.Get("geninfo")
    # get the "prmtrkdep"-TTree from TFile:
    rtree_prmtrkdep = rfile.Get("prmtrkdep")

    # get the number of events in the geninfo Tree:
    number_events = rtree_geninfo.GetEntries()
    # check if number_events_geninfo is equal to number_entries_input (if not, the detector simulation was incorrect!!):
    if number_events != number_entries_input:
        sys.exit("ERROR: number of events are not equal to {0:d} -> Detector Simulation not correct!"
                 .format(number_entries_input))

    # number of events with reconstructed position inside fiducial volume (r_reco < Radius_cut):
    n_reco_inside = 0
    # number of leak-out events (r_init < Radius_cut, but r_reco >= Radius_cut):
    n_leak_out = 0
    # number of leak-in events (r_init >= Radius_cut, but r_reco < Radius_cut):
    n_leak_in = 0

    # loop over every event, i.e. every entry, in the TTree:
    for event in range(number_events):

        """ read prmtrkdep tree: """
        rtree_prmtrkdep.GetEntry(event)
        # get nInitParticles from prmtrkdep tree:
        n_par_prmtrkdep = int(rtree_prmtrkdep.GetBranch('nInitParticles').GetLeaf('nInitParticles').GetValue())
        # preallocate sum of Qedep of initial particles:
        qedep_sum = 0
        # to get total quenched deposited energy, sum over initial particles:
        for index1 in range(n_par_prmtrkdep):
            # get deposit energy of initial neutron in MeV:
            qedep = float(rtree_prmtrkdep.GetBranch('Qedep').GetLeaf('Qedep').GetValue(index1))
            # add qedep to qedep_sum:
            qedep_sum += qedep

        """ read 'geninfo' tree: """
        # get the current event in the Tree:
        rtree_geninfo.GetEntry(event)

        # get event ID:
        evt_id = int(rtree_geninfo.GetBranch('evtID').GetLeaf('evtID').GetValue())
        if evt_id != event:
            sys.exit("ERROR: event ID's are not equal in file {0}".format(rfile))

        # get number of particles in the event:
        n_par_geninfo = int(rtree_geninfo.GetBranch('nInitParticles').GetLeaf('nInitParticles').GetValue())

        # check if there are initial particles:
        if n_par_geninfo == 0:
            # no initial particle -> go to next event
            continue

        # preallocate array for initial position:
        initx = np.array([])
        inity = np.array([])
        initz = np.array([])

        # loop over number of particles in the event:
        for index1 in range(n_par_geninfo):
            # get initial position of the particle in mm :
            initx = np.append(initx, float(rtree_geninfo.GetBranch('InitX').GetLeaf('InitX').GetValue(index1)))
            inity = np.append(inity, float(rtree_geninfo.GetBranch('InitY').GetLeaf('InitY').GetValue(index1)))
            initz = np.append(initz, float(rtree_geninfo.GetBranch('InitZ').GetLeaf('InitZ').GetValue(index1)))

        # set 0th entry of array as initial position in mm:
        x_init = initx[0]
        y_init = inity[0]
        z_init = initz[0]

        # check if all initial position are equal:
        for index1 in range(n_par_geninfo):
            if x_init != initx[index1] or y_init != inity[index1] or z_init != initz[index1]:
                sys.exit("ERROR: initial positions are not equal for all initial particles (event {0:d}, file {1})"
                         .format(event, rfile))

        # smear the initial position with vertex reconstruction:
        x_reco = NC_background_functions.position_smearing(x_init, qedep_sum)
        y_reco = NC_background_functions.position_smearing(y_init, qedep_sum)
        z_reco = NC_background_functions.position_smearing(z_init, qedep_sum)

        # calculate r_init in mm:
        r_init = np.sqrt(x_init**2 + y_init**2 + z_init**2)
        # calculate r_reco in mm:
        r_reco = np.sqrt(x_reco**2 + y_reco**2 + z_reco**2)

        if r_reco < radius_cut:
            # reco. position inside fiducial volume:
            n_reco_inside += 1

        if r_init < radius_cut and r_reco >= radius_cut:
            n_leak_out += 1

        if r_init >= radius_cut and r_reco < radius_cut:
            n_leak_in += 1

    return number_events, n_reco_inside, n_leak_out, n_leak_in


# get the date and time, when the script was run:
date = datetime.datetime.now()
now = date.strftime("%Y-%m-%d %H:%M")

# Set cut radius, that defines the fiducial volume, in mm (normally volume cut R < 16 m = 16000 mm):
Radius_cut = 16000

""" analyze the user_atmoNC_{}.root file from atmospheric Neutral current simulation: """
# set the path of the input files (filename must be 'user_atmoNC_{}.root'):
Input_path = "/local/scratch1/pipc51/astro/blum/detsim_output_data/"

# set path, where results should be saved:
Output_path = "/home/astro/blum/juno/atmoNC/data_NC/output_volume_cut/"

# set the number of the first file and number of the last file that should be read:
start_number = 0
stop_number = 999
# number of entries in the input files:
Number_entries_input = 100

# number of total events, that are analyzed:
Number_events = 0
# number of events with reconstructed position inside fiducial volume (r_reco < Radius_cut):
number_reco_inside = 0
# number of leak-out events (r_init < Radius_cut, but r_reco >= Radius_cut):
number_leak_out = 0
# number of leak-in events (r_init >= Radius_cut, but r_reco < Radius_cut):
number_leak_in = 0

# loop over files:
for file_num in range(start_number, stop_number+1, 1):
    # path to file:
    input_file = Input_path + "user_atmoNC_{0:d}.root".format(file_num)
    # print("Start reading {0} ...".format(input_file))

    num_events, num_reco_inside, num_leak_out, num_leak_in = check_volume_cut(input_file,
                                                                              Number_entries_input, Radius_cut)

    # add values:
    Number_events += num_events
    number_reco_inside += num_reco_inside
    number_leak_out += num_leak_out
    number_leak_in += num_leak_in

# calculate total volume cut efficiency in percent (efficiency: reco. events inside fiducial volume / total events):
efficiency_total = float(number_reco_inside) / float(Number_events) * 100.0

# calculate leak efficiency in percent ():
efficiency_leak = 100.0 + float((number_leak_out - number_leak_in)) / float(number_reco_inside) * 100.0

print("NC events:")
print("total number of events = {0:d}".format(Number_events))
print("number of reco. events inside fiducial volume = {0:d}".format(number_reco_inside))
print("total volume cut efficiency = {0:0.4f} %".format(efficiency_total))
print("")
print("number of leak out events = {0:d}".format(number_leak_out))
print("number of leak in events = {0:d}".format(number_leak_in))
print("leak efficiency (1 - (number leak-out - number leak-in) / number_reco) = {0:.4f} %".format(efficiency_leak))
print("-------------------")

""" analyze the user_positron_volcut_{}.root file from simulation of positrons with momentum from 10 MeV to 100 MeV 
    -> represent the prompt signal of real IBD-like events: """
# set the path of the input files (filename must be 'user_positron_{}.root'):
Input_path_pos = "/local/scratch1/pipc51/astro/blum/positron_output_volcut/"

# set path, where results should be saved:
Output_path_pos = "/home/astro/blum/juno/atmoNC/data_NC/output_volume_cut/"

# set the number of the first file and number of the last file that should be read:
start_number_pos = 0
stop_number_pos = 99
# number of entries in the input files:
Number_entries_input_pos = 1000

# number of total events, that are analyzed:
Number_events_pos = 0
# number of events with reconstructed position inside fiducial volume (r_reco < Radius_cut):
number_reco_inside_pos = 0
# number of leak-out events (r_init < Radius_cut, but r_reco >= Radius_cut):
number_leak_out_pos = 0
# number of leak-in events (r_init >= Radius_cut, but r_reco < Radius_cut):
number_leak_in_pos = 0

# loop over files:
for file_num in range(start_number_pos, stop_number_pos+1, 1):
    # path to file:
    input_file_pos = Input_path_pos + "user_positron_volcut_{0:d}.root".format(file_num)
    # print("Start reading {0} ...".format(input_file))

    num_events_pos, num_reco_inside_pos, num_leak_out_pos, num_leak_in_pos = check_volume_cut(input_file_pos,
                                                                                              Number_entries_input_pos,
                                                                                              Radius_cut)

    # add values:
    Number_events_pos += num_events_pos
    number_reco_inside_pos += num_reco_inside_pos
    number_leak_out_pos += num_leak_out_pos
    number_leak_in_pos += num_leak_in_pos

# calculate total volume cut efficiency in percent (efficiency: reco. events inside fiducial volume / total events):
efficiency_total_pos = float(number_reco_inside_pos) / float(Number_events_pos) * 100.0

# calculate leak efficiency in percent ():
efficiency_leak_pos = 100.0 + float((number_leak_out_pos - number_leak_in_pos)) / float(number_reco_inside_pos) * 100.0

print("Positron events:")
print("total number of events = {0:d}".format(Number_events_pos))
print("number of reco. events inside fiducial volume = {0:d}".format(number_reco_inside_pos))
print("total volume cut efficiency = {0:0.4f} %".format(efficiency_total_pos))
print("")
print("number of leak out events = {0:d}".format(number_leak_out_pos))
print("number of leak in events = {0:d}".format(number_leak_in_pos))
print("leak efficiency (1 - (number leak-out - number leak-in) / number_reco) = {0:.4f} %".format(efficiency_leak_pos))



