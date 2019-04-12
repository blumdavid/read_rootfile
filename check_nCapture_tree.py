""" script to check the information saved in the nCapture Tree in detsim user output root file:

    What information is saved in the different parameters?


"""
import datetime
import ROOT
import sys
import numpy as np
from matplotlib import pyplot as plt


def check_n_capture(input_file, start_num, stop_num):
    """
    function to read the nCapture Tree of the input root file and return information about nCapture Tree
    :param input_file: detsim user output root file (e.g. user_p_atmoNC3_evt83.root)
    :param start_num: first event of the root file that is read
    :param stop_num: last event of the root file that is read
    :return:
    """
    # load the ROOT file:
    rfile = ROOT.TFile(input_file)
    # get the "nCapture"-TTree from the TFile:
    rtree_n_capture = rfile.Get("nCapture")
    # get the 'geninfo' TTree from the TFile:
    rtree_geninfo = rfile.Get("geninfo")

    # preallocate the arrays:
    """
    arr_neutron_trkid = np.array([])
    arr_neutron_kine = np.array([])
    arr_neutron_capture_t = np.array([])
    arr_nc_start_x = np.array([])
    arr_nc_start_y = np.array([])
    arr_nc_start_z = np.array([])
    arr_nc_stop_x = np.array([])
    arr_nc_stop_y = np.array([])
    arr_nc_stop_z = np.array([])
    arr_nc_track_length = np.array([])

    arr_trkid = np.array([])
    arr_pdgid = np.array([])
    arr_kine = np.array([])
    arr_px = np.array([])
    arr_py = np.array([])
    arr_pz = np.array([])
    arr_x = np.array([])
    arr_y = np.array([])
    arr_z = np.array([])
    arr_t = np.array([])
    """

    # loop over events, i.e. entries, in the TTree:
    for event in range(start_num, stop_num+1, 1):

        print("--------------read event {0:d}:--------------".format(event))

        # get the current event in the geninfo tree:
        rtree_geninfo.GetEntry(event)

        # get the value of the event ID:
        evt_id_geninfo = int(rtree_geninfo.GetBranch('evtID').GetLeaf('evtID').GetValue())
        print("\n--------read geninfo Tree (evtID={0:d}):--------".format(evt_id_geninfo))

        # get the number of initial particles:
        n_init_particles = int(rtree_geninfo.GetBranch('nInitParticles').GetLeaf('nInitParticles').GetValue())

        # get initial PDG IDs:
        for index in range(n_init_particles):
            init_pdgid = int(rtree_geninfo.GetBranch('InitPDGID').GetLeaf('InitPDGID').GetValue(index))
            print("Initial PDGID = {0:d}".format(init_pdgid))

        # set neutron_capture flag:
        flag_n_capture = False

        # get the current event in the nCapture Tree:
        rtree_n_capture.GetEntry(event)

        # get the value of the event ID:
        evt_id = int(rtree_n_capture.GetBranch('evtID').GetLeaf('evtID').GetValue())
        print("\n--------read nCapture Tree (evtID={0:d}):--------".format(evt_id))

        # get the number of entries of the variable NeutronTrkid:
        num_neutron_trkid = rtree_n_capture.GetBranch('NeutronTrkid').GetLeaf('NeutronTrkid').GetLen()

        # get the number of entries of the variable trkid:
        num_trkid = rtree_n_capture.GetBranch('trkid').GetLeaf('trkid').GetLen()

        if num_neutron_trkid == 0 and num_trkid == 0:

            print("no neutron and no secondaries in event\n")
            continue

        elif num_neutron_trkid == 0 or num_trkid == 0:
            # just to check if this is possible:
            print("WARNING: num_NeutronTrkid == 0 or num_trkid == 0!!!\n")
            continue

        else:
            # set neutron capture flag:
            flag_n_capture = True

        if flag_n_capture:
            # there is at least one neutron and one secondary in the event:
            print("\n------ Info about Neutron Capture: ------")
            # read variables from Neutron Capture:
            for num in range(num_neutron_trkid):
                # get NeutronTrkid:
                neutron_trkid = int(rtree_n_capture.GetBranch("NeutronTrkid").GetLeaf("NeutronTrkid").GetValue(num))
                print("NeutronTrkid = {0:d}".format(neutron_trkid))

                # get NeutronKine in MeV:
                neutron_kine = float(rtree_n_capture.GetBranch("NeutronKine").GetLeaf("NeutronKine").GetValue(num))
                print("NeutronKine = {0:.2f} MeV".format(neutron_kine))

                # get NeutronCaptureT in ns:
                neutron_capture_t = float(rtree_n_capture.GetBranch("NeutronCaptureT").GetLeaf("NeutronCaptureT")
                                          .GetValue(num))
                print("NeutronCaptureT = {0:.2f} ns".format(neutron_capture_t))

                # get NCStartX in mm:
                nc_start_x = float(rtree_n_capture.GetBranch("NCStartX").GetLeaf("NCStartX").GetValue(num))
                print("NCStartX = {0:.2f} mm".format(nc_start_x))

                # get NCStartY in mm:
                nc_start_y = float(rtree_n_capture.GetBranch("NCStartY").GetLeaf("NCStartY").GetValue(num))
                print("NCStartY = {0:.2f} mm".format(nc_start_y))

                # get NCStartZ in mm:
                nc_start_z = float(rtree_n_capture.GetBranch("NCStartZ").GetLeaf("NCStartZ").GetValue(num))
                print("NCStartZ = {0:.2f} mm".format(nc_start_z))

                # get NCStopX in mm:
                nc_stop_x = float(rtree_n_capture.GetBranch("NCStopX").GetLeaf("NCStopX").GetValue(num))
                print("NCStopX = {0:.2f} mm".format(nc_stop_x))

                # get NCStopY in mm:
                nc_stop_y = float(rtree_n_capture.GetBranch("NCStopY").GetLeaf("NCStopY").GetValue(num))
                print("NCStopY = {0:.2f} mm".format(nc_stop_y))

                # get NCStopZ in mm:
                nc_stop_z = float(rtree_n_capture.GetBranch("NCStopZ").GetLeaf("NCStopZ").GetValue(num))
                print("NCStopZ = {0:.2f} mm".format(nc_stop_z))

                # calculate distance between start and stop point in mm:
                nc_distance = np.sqrt((nc_stop_x - nc_start_x)**2 + (nc_stop_y - nc_start_y)**2 +
                                      (nc_stop_z - nc_start_z)**2)
                print("distance Stop - Start = {0:.2f} mm".format(nc_distance))

                # get NCTrackLength in mm:
                nc_track_length = float(rtree_n_capture.GetBranch("NCTrackLength").GetLeaf("NCTrackLength").GetValue(num))
                print("NCTrackLength = {0:.2f} mm\n".format(nc_track_length))

            print("------ Info about Secondaries: ------")
            # read variables from secondaries:
            for num in range(num_trkid):
                # get trkid:
                trkid = int(rtree_n_capture.GetBranch("trkid").GetLeaf("trkid").GetValue(num))
                print("\ntrkid = {0:d}".format(trkid))

                # get pdgid:
                pdgid = int(rtree_n_capture.GetBranch("pdgid").GetLeaf("pdgid").GetValue(num))
                print("pdgid = {0:d}".format(pdgid))

                # get kinetic energy in MeV:
                kine = float(rtree_n_capture.GetBranch("kine").GetLeaf("kine").GetValue(num))
                print("kine = {0:.2f} MeV".format(kine))

                # get momentum in x, y, z direction in MeV:
                px = float(rtree_n_capture.GetBranch("px").GetLeaf("px").GetValue(num))
                print("px = {0:.2f} MeV".format(px))
                py = float(rtree_n_capture.GetBranch("py").GetLeaf("py").GetValue(num))
                print("py = {0:.2f} MeV".format(py))
                pz = float(rtree_n_capture.GetBranch("pz").GetLeaf("pz").GetValue(num))
                print("pz = {0:.2f} MeV".format(pz))

                # calculate total momentum in MeV:
                p_total = np.sqrt(px**2 + py**2 + pz**2)
                print("total momentum = {0:.2f} MeV".format(p_total))

                # get x, y, z position in mm:
                x = float(rtree_n_capture.GetBranch("x").GetLeaf("x").GetValue(num))
                print("x = {0:.2f} mm".format(x))
                y = float(rtree_n_capture.GetBranch("y").GetLeaf("y").GetValue(num))
                print("y = {0:.2f} mm".format(y))
                z = float(rtree_n_capture.GetBranch("z").GetLeaf("z").GetValue(num))
                print("z = {0:.2f} mm".format(z))

                # get time in ns:
                t = float(rtree_n_capture.GetBranch("t").GetLeaf("t").GetValue(num))
                print("t = {0:.2f} ns".format(t))

    return


# get the date and time, when the script was run:
date = datetime.datetime.now()
now = date.strftime("%Y-%m-%d %H:%M")

# set the path of the input files (filename must be 'user_atmoNC_{}.root'):
# Input_path = "/local/scratch1/pipc51/astro/blum/detsim_output_data/"
Input_path = "/home/astro/blum/juno/atmoNC/data_NC/secondary_anamgr/"

# set path, where results should be saved:
Output_path = "/home/astro/blum/juno/atmoNC/data_NC/secondary_anamgr/"

# set the file number, the number of the first event and number of the last event that should be read:
file_number = 0
start_number = 0
stop_number = 2

# path to file:
# Input_file = Input_path + ".root"
Input_file = Input_path + "user_atmoNC_{0:d}_secAnaMgr.root".format(file_number)

print("Start reading {0} ...".format(Input_file))

# check_n_capture(Input_file, start_number, stop_number)
























