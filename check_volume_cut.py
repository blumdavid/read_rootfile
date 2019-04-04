""" Script to check the efficiency of the fiducial volume cut for atmoNC events.

    The radius of JUNO detector is 17.7m, but only a fiducial volume of the detector is used for detection to suppress
    background from the surroundings (normally 17 m).

    The efficiency of this volume cut is calculated with 3 different parameters:
    1. with edepX, edepY, edepZ from 'evt'-tree of user_atmoNC_.root
    2. with InitX, InitY, InitZ from 'geninfo'-tree of user_atmoNC_.root
    3. with ExitX, ExitY, ExitZ from 'geninfo'-tree of user_atmoNC_.root

    Then these 3 efficiency can be compared
"""
import datetime
import numpy as np
import ROOT
import sys


def check_volume_cut(input_path, file_number, number_entries_input, radius_cut):
    """

    :param input_path: path to input root files from tut_detsim.py: user_atmoNC_{}.root
    :param file_number: number of the input root file: e.g. 2 -> user_atmoNC_2.root
    :param number_entries_input:  number of entries, that the input files should have (integer), normally = 100
    :param radius_cut: radius, which defines the fiducial volume, in mm
    :return:
    """
    # load the ROOT file:
    rfile = ROOT.TFile(input_path + "user_atmoNC_{0:d}.root".format(file_number))
    # get the "evt"-TTree from the TFile:
    rtree_evt = rfile.Get("evt")
    # get the "geninfo"-TTree from the TFile:
    rtree_geninfo = rfile.Get("geninfo")

    # get the number of events in the 'evt' Tree:
    number_events_evt = rtree_evt.GetEntries()
    # get the number of events in the geninfo Tree:
    number_events_geninfo = rtree_geninfo.GetEntries()
    if number_events_geninfo == number_events_evt:
        number_events = number_events_geninfo
    else:
        sys.exit("ERROR: number of events in Trees are NOT equal!!")

    # check if number_events is equal to number_entries_input (if not, the detector simulation was incorrect!!):
    if number_events != number_entries_input:
        sys.exit("ERROR: number of events are not equal to {0:d} -> Detector Simulation not correct!"
                 .format(number_entries_input))

    # preallocate number of events INSIDE fiducial volume for edepX, edepY, edepZ:
    number_cut_edep_xyz = 0
    # preallocate number of events INSIDE fiducial volume for InitX, InitY, InitZ:
    number_cut_init_xyz = 0
    # preallocate number of events INSIDE fiducial volume for ExitX, ExitYm ExitZ:
    number_cut_exit_xyz = 0

    # loop over every event, i.e. every entry, in the TTree:
    for event in range(number_events):

        # set flag, that init-position of all particles in event is inside fiducial volume:
        flag_init_inside = True
        # set flag, that exit-position of all particles in event is inside fiducial volume:
        flag_exit_inside = True

        """ first read 'evt' Tree: """
        # get the current event in the Tree:
        rtree_evt.GetEntry(event)

        # get the value of the event ID:
        evt_id_evt = int(rtree_evt.GetBranch('evtID').GetLeaf('evtID').GetValue())

        # get x-position, where energy is deposited, in mm:
        edepx = float(rtree_evt.GetBranch('edepX').GetLeaf('edepX').GetValue())

        # get y-position, where energy is deposited, in mm:
        edepy = float(rtree_evt.GetBranch('edepY').GetLeaf('edepY').GetValue())

        # get z-position, where energy is deposited, in mm:
        edepz = float(rtree_evt.GetBranch('edepZ').GetLeaf('edepZ').GetValue())

        # calculate distance to detector center in mm:
        r_edepxyz = np.sqrt(edepx**2 + edepy**2 + edepz**2)

        if r_edepxyz < radius_cut:
            # deposit energy of whole event inside fiducial volume:
            number_cut_edep_xyz = number_cut_edep_xyz + 1

        """ read 'geninfo' tree: """
        # get the current event in the Tree:
        rtree_geninfo.GetEntry(event)

        # get event ID:
        evt_id_geninfo = int(rtree_geninfo.GetBranch('evtID').GetLeaf('evtID').GetValue())

        if evt_id_evt == evt_id_geninfo:
            evt_id = evt_id_evt
        else:
            sys.exit("Event ID's of the three trees are NOT equal for 1 event!")

        # get number of particles in the event:
        n_par_geninfo = int(rtree_geninfo.GetBranch('nInitParticles').GetLeaf('nInitParticles').GetValue())

        # preallocate array of distance to center of initial and exit position in mm:
        r_init_array = np.array([])
        r_exit_array = np.array([])

        # loop over number of particles in the event:
        for index1 in range(n_par_geninfo):
            # get initial position of the particle in mm :
            initx = float(rtree_geninfo.GetBranch('InitX').GetLeaf('InitX').GetValue(index1))
            inity = float(rtree_geninfo.GetBranch('InitY').GetLeaf('InitY').GetValue(index1))
            initz = float(rtree_geninfo.GetBranch('InitZ').GetLeaf('InitZ').GetValue(index1))

            # calculate distance to center for initial position in mm:
            r_init = np.sqrt(initx**2 + inity**2 + initz**2)
            r_init_array = np.append(r_init_array, r_init)

            # get exit position of the particles in mm:
            exitx = float(rtree_geninfo.GetBranch('ExitX').GetLeaf('ExitX').GetValue(index1))
            exity = float(rtree_geninfo.GetBranch('ExitY').GetLeaf('ExitY').GetValue(index1))
            exitz = float(rtree_geninfo.GetBranch('ExitZ').GetLeaf('ExitZ').GetValue(index1))

            # calculate distance to center for exit position in mm:
            r_exit = np.sqrt(exitx**2 + exity**2 + exitz**2)
            r_exit_array = np.append(r_exit_array, r_exit)

        for index2 in range(len(r_init_array)):
            if r_init_array[index2] >= radius_cut:
                # particle outside of the fiducial volume:
                flag_init_inside = False
                break
            else:
                # particle inside of fiducial volume -> check next particle position
                continue

        for index3 in range(len(r_exit_array)):
            if r_exit_array[index3] >= radius_cut:
                # particle outside of the fiducial volume:
                flag_exit_inside = False
                break
            else:
                # particle inside of fiducial volume -> check next particle position
                continue

        # check flags:
        if flag_init_inside:
            # all particles of event are inside fiducial volume:
            number_cut_init_xyz = number_cut_init_xyz + 1

        if flag_exit_inside:
            # all particles of event are inside fiducial volume:
            number_cut_exit_xyz = number_cut_exit_xyz + 1

    return number_events, number_cut_edep_xyz, number_cut_init_xyz, number_cut_exit_xyz


# get the date and time, when the script was run:
date = datetime.datetime.now()
now = date.strftime("%Y-%m-%d %H:%M")

# set the path of the input files (filename must be 'user_atmoNC_{}.root'):
Input_path = "/local/scratch1/pipc51/astro/blum/detsim_output_data/"

# set path, where results should be saved:
Output_path = "/home/astro/blum/juno/atmoNC/data_NC/output_volume_cut/"

# set the number of the first file and number of the last file that should be read:
start_number = 0
stop_number = 99
# number of entries in the input files:
Number_entries_input = 100

# Set cut radius, that defines the fiducial volume, in mm (normally volume cut R < 17 m = 17000 mm):
Radius_cut = 17000

# number of total events, that are analyzed:
Number_events = 0
# preallocate number of events INSIDE fiducial volume for edepX, edepY, edepZ:
Number_cut_edep_xyz = 0
# preallocate number of events INSIDE fiducial volume for InitX, InitY, InitZ:
Number_cut_init_xyz = 0
# preallocate number of events INSIDE fiducial volume for ExitX, ExitYm ExitZ:
Number_cut_exit_xyz = 0

# loop over files:
for file_num in range(start_number, stop_number+1, 1):
    # path to file:
    input_file = "user_atmoNC_{0:d}.root".format(file_num)
    print("Start reading {0} ...".format(input_file))

    num_events, num_cut_edep_xyz, num_cut_init_xyz, num_cut_exit_xyz = check_volume_cut(Input_path, file_num,
                                                                                        Number_entries_input,
                                                                                        Radius_cut)

    # add values:
    Number_events = Number_events + num_events
    Number_cut_edep_xyz = Number_cut_edep_xyz + num_cut_edep_xyz
    Number_cut_init_xyz = Number_cut_init_xyz + num_cut_init_xyz
    Number_cut_exit_xyz = Number_cut_exit_xyz + num_cut_exit_xyz

# calculate the efficiencies in percent (efficiency: events inside fiducial volume / total events):
efficiency_edep_xyz = float(Number_cut_edep_xyz) / float(Number_events) * 100
efficiency_init_xyz = float(Number_cut_init_xyz) / float(Number_events) * 100
efficiency_exit_xyz = float(Number_cut_exit_xyz) / float(Number_events) * 100

# save results in txt file:
np.savetxt(Output_path + "vol_cut_atmoNC_{0:d}_to_{1:d}.txt".format(start_number, stop_number),
           np.array([Number_events, Number_cut_edep_xyz, Number_cut_init_xyz, Number_cut_exit_xyz, efficiency_edep_xyz,
                     efficiency_init_xyz, efficiency_exit_xyz]), fmt='%.3f',
           header="Results of script check_volume_cut.py ({0}):\n"
                  "Input-files: user_atmoNC_{1:d}.root to user_atmoNC_{2:d}.root;\n"
                  "Cut radius = {3:d} mm;\n"
                  "\n"
                  "Results:\n"
                  "Total number of events;\n"
                  "Number of events inside fid. volume from 'edepXYZ' of 'evt'-tree;\n"
                  "Number of events inside fid. volume from 'InitXYZ' of 'geninfo'-tree;\n"
                  "Number of events inside fid. volume from 'ExitXYZ' of 'geninfo'-tree;\n"
                  "Efficiency of cut on 'edepXYZ' in percent;\n"
                  "Efficiency of cut on 'InitXYZ' in percent;\n"
                  "Efficiency of cut on 'ExitXYZ' in percent:".format(now, start_number, stop_number, Radius_cut))





