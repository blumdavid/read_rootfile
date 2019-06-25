""" script to check the information saved in the new "secondary" Tree in detsim user output root file:

    What information is saved in the different parameters?


"""
import datetime
import ROOT
import sys
import numpy as np
from matplotlib import pyplot as plt
from check_nCapture_tree import check_n_capture


def read_secondaries_tree(input_file, first_evt, last_evt):
    """ function to read the secondaries tree from user_atmoNC_{}.root to check the different secondary particles, which
    are produced in the detector.


    :param input_file: input root file with user output
    :param first_evt: first event in the file to be read
    :param last_evt: last event in the file to be read

    :return:
    """
    # load the ROOT file:
    rfile = ROOT.TFile(input_file)
    # get the "geninfo"-Tree from the TFile:
    rtree_geninfo = rfile.Get("geninfo")
    # get the "secondaries"-TTree from the TFile:
    rtree_sec = rfile.Get("secondaries")

    # loop over events, i.e. entries, in the TTree:
    for event in range(first_evt, last_evt+1, 1):

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

        """ read secondaries tree: """
        # set secondaries flag:
        flag_secondaries = False

        # get the current event in the geninfo tree:
        rtree_sec.GetEntry(event)

        # get the value of the event ID:
        evt_id = int(rtree_sec.GetBranch('evtID').GetLeaf('evtID').GetValue())
        print("\n--------read secondaries Tree (evtID={0:d}):--------".format(evt_id))

        # get the number of primaries:
        num_primaries = int(rtree_sec.GetBranch('numPrimaries').GetLeaf('numPrimaries').GetValue())

        # get the number of entries of the variable SecTrkID:
        num_secondaries = int(rtree_sec.GetBranch('numSec').GetLeaf('numSec').GetValue())

        if num_primaries == 0 and num_secondaries == 0:

            print("NO particles in event\n")
            continue

        elif num_secondaries == 0:

            print("NO secondaries in event\n")
            continue

        else:
            # set secondaries flag:
            flag_secondaries = True

        if flag_secondaries:
            # loop over primary particles:
            for num in range(num_primaries):
                # there are secondaries in the event:
                print("\n------------ Info about primaries: ------------")

                # get PrimaryTrkID:
                pri_trkid = int(rtree_sec.GetBranch("PrimaryTrkID").GetLeaf("PrimaryTrkID").GetValue(num))
                print("PrimaryTrkID = {0:d}".format(pri_trkid))

                # get PrimaryPDGID:
                pri_pdgid = int(rtree_sec.GetBranch("PrimaryPDGID").GetLeaf("PrimaryPDGID").GetValue(num))
                print("PrimaryPDGID = {0:d}".format(pri_pdgid))

                # get kinetic energy:
                pri_kine = float(rtree_sec.GetBranch("PrimaryKine").GetLeaf("PrimaryKine").GetValue(num))
                print("PrimaryKine = {0:.2f} MeV".format(pri_kine))

                # get stop position:
                pri_stopx = float(rtree_sec.GetBranch("PrimaryStopX").GetLeaf("PrimaryStopX").GetValue(num))
                pri_stopy = float(rtree_sec.GetBranch("PrimaryStopY").GetLeaf("PrimaryStopY").GetValue(num))
                pri_stopz = float(rtree_sec.GetBranch("PrimaryStopZ").GetLeaf("PrimaryStopZ").GetValue(num))
                print("PrimaryStopX = {0:.2f} mm".format(pri_stopx))
                print("PrimaryStopY = {0:.2f} mm".format(pri_stopy))
                print("PrimaryStopZ = {0:.2f} mm".format(pri_stopz))

                # get time:
                pri_time = float(rtree_sec.GetBranch("PrimaryTime").GetLeaf("PrimaryTime").GetValue(num))
                print("PrimaryTime = {0:.2f} ns".format(pri_time))

                # get track length:
                pri_trklength = float(rtree_sec.GetBranch("PrimaryTrkLen").GetLeaf("PrimaryTrkLen").GetValue(num))
                print("PrimaryTrkLen = {0:.2f} mm".format(pri_trklength))

                # loop over possible secondaries:
                for num1 in range(num_secondaries):

                    # get parent track ID ParentTrkID:
                    parent_trkid = int(rtree_sec.GetBranch("ParentTrkID").GetLeaf("ParentTrkID").GetValue(num1))

                    if pri_trkid != parent_trkid:
                        # secondary do NOT correspond to primary:
                        continue
                    else:
                        print("\n------ Info about secondaries: ------")
                        # secondary correspond to primary:
                        sec_trkid = int(rtree_sec.GetBranch("SecTrkID").GetLeaf("SecTrkID").GetValue(num1))
                        print("SecTrkID = {0:d}".format(sec_trkid))

                        # get SecPDGID:
                        sec_pdgid = int(rtree_sec.GetBranch("SecPDGID").GetLeaf("SecPDGID").GetValue(num1))
                        print("SecPDGID = {0:d}".format(sec_pdgid))

                        # get SecProcess (number of process):
                        sec_process = int(rtree_sec.GetBranch("SecProcess").GetLeaf("SecProcess").GetValue(num1))
                        # get name of process:
                        sec_process_name = get_process_name(sec_process)
                        print("SecProcess Name = {0}".format(sec_process_name))

                        # get SecKine:
                        sec_kine = float(rtree_sec.GetBranch("SecKine").GetLeaf("SecKine").GetValue(num1))
                        print("SecKine = {0:.2f} MeV".format(sec_kine))

                        # get position of secondary:
                        sec_posx = float(rtree_sec.GetBranch("SecPosX").GetLeaf("SecPosX").GetValue(num1))
                        sec_posy = float(rtree_sec.GetBranch("SecPosY").GetLeaf("SecPosY").GetValue(num1))
                        sec_posz = float(rtree_sec.GetBranch("SecPosZ").GetLeaf("SecPosZ").GetValue(num1))
                        print("SecPosX = {0:.2f} mm".format(sec_posx))
                        print("SecPosY = {0:.2f} mm".format(sec_posy))
                        print("SecPosZ = {0:.2f} mm".format(sec_posz))

                        # get SecTime:
                        sec_time = float(rtree_sec.GetBranch("SecTime").GetLeaf("SecTime").GetValue(num1))
                        print("SecTime = {0:.2f} ns".format(sec_time))

    return


def get_process_name(process_number):
    """ Function to get the process name from the process number set in SecondariesAnaMgr.cc

    :param process_number: process_number (integer)
    :return: process_name (string)
    """
    if process_number == 1:
        # Ionisation process of charged hadrons and ions, including low energy extensions:
        process_name = "hLowEIoni"
    elif process_number == 2:
        # Process for Proton Inelastic scattering:
        process_name = "ProtonInelastic"
    elif process_number == 3:
        # Final state production model for hadron nuclear elastic scattering:
        process_name = "LElastic"
    elif process_number == 10:
        # Electron/positron processes: Low Energy electromagnetic process, electron Ionisation:
        process_name = "LowEnergyIoni"
    elif process_number == 11:
        # Process of e+ annihilation into 2 gammas:
        process_name = "annihil"
    elif process_number == 12:
        # Ionisation process for e-/e+:
        process_name = "eIoni"
    elif process_number == 13:
        # Bremsstrahlung for e-/e+:
        process_name = "eBrem"
    elif process_number == 14:
        # Low Energy electromagnetic process, electron Bremsstrahlung:
        process_name = "LowEnBrem"
    elif process_number == 15:
        # Process for positron nuclear inelastic scattering
        process_name = "PositronNuclear"
    elif process_number == 20:
        # neutron: neutron capture:
        process_name = "nCapture"
    elif process_number == 21:
        # Process of neutron elastic scattering:
        process_name = "neutronElastic"
    elif process_number == 22:
        # Process for Neutron Inelastic scattering:
        process_name = "NeutronInelastic"
    elif process_number == 30:
        # Muon: Ionisation process for muons:
        process_name = "muIoni"
    elif process_number == 31:
        # MuonMinus captured at rest (interesting for michel electrons (see InterestingProcessAnaMgr.cc)):
        process_name = "muMinusCaptureAtRest"
    elif process_number == 40:
        # Decay:
        process_name = "Decay"
    elif process_number == 41:
        # nuclear radioactive decay process. It simulates the decays of radioactive nuclei:
        process_name = "RadioactiveDecay"
    elif process_number == 50:
        # Pion: Process for PionMinus Inelastic scattering:
        process_name = "PionMinusInelastic"
    elif process_number == 51:
        # Process for PionPlus Inelastic scattering:
        process_name = "PionPlusInelastic"
    elif process_number == 60:
        # Kaon: Process for KaonZeroLong Inelastic scattering:
        process_name = "KaonZeroLInelastic"
    elif process_number == 61:
        # Process for KaonPlus Inelastic scattering:
        process_name = "KaonPlusInelastic"
    elif process_number == 70:
        # deuteron: Process for Deuteron Inelastic scattering:
        process_name = "DeuteronInelastic"
    elif process_number == 80:
        # Lambda: Process for Lambda Inelastic scattering:
        process_name = "LambdaInelastic"
    elif process_number == 90:
        # Photon/gamma: Process of Compton scattering:
        process_name = "LowEnCompton"
    elif process_number == 91:
        # Photon/gamma: Process of photoelectric effect:
        process_name = "LowEnPhotoElec"
    elif process_number == 92:
        # Photon/gamma: Process of photon conversion (pair production):
        process_name = "LowEnConversion"
    elif process_number == 93:
        # Photon/gamma: Process for photon nuclear inelastic scattering:
        process_name = "PhotonInelastic"
    else:
        # process name and number not yet implemented (process_number = 100 means "new" process)
        process_name = "WARNING: not yet implemented"

    return process_name


def check_secondaries_tree(input_file, first_evt, last_evt, print_all, print_process, print_initials):
    """
    function to read the secondaries tree and get the processes, which lead to neutron capture. This is to understand,
    which processes can cause neutron capture and therefore a possible delayed signal for IBD-like signals.

    :param input_file: input root file with user output
    :param first_evt: first event in the file to be read
    :param last_evt: last event in the file to be read
    :param print_all: print all information into terminal (boolean)
    :param print_process: print only the information of the processes, where there is neutron capture, but no initial
    neutron (boolean)
    :param print_initials: print numbers of different initial particles, that cause neutron capture, but aren't neutrons
    :return:
    """
    # load the ROOT file:
    rfile = ROOT.TFile(input_file)
    # get the "geninfo"-Tree from the TFile:
    rtree_geninfo = rfile.Get("geninfo")
    # get the "secondaries"-TTree from the TFile:
    rtree_sec = rfile.Get("secondaries")
    # get the "nCapture"-TTree from the TFile:
    rtree_ncapture = rfile.Get("nCapture")

    # preallocate numbers of different initial particles, that cause neutron capture, but aren't neutrons:
    number_initial_proton = 0
    number_initial_piminus = 0
    number_initial_piplus = 0
    number_initial_pizero = 0
    # preallocate number of neutron, that cause neutron capture:
    number_initial_neutron = 0

    # loop over events, i.e. entries, in the TTree:
    for event in range(first_evt, last_evt+1, 1):

        if print_all:
            print("#####################    read event {0:d}:   ########################".format(event))

        """ read geninfo tree """
        # get the current event in the geninfo tree:
        rtree_geninfo.GetEntry(event)

        # get the value of the event ID:
        evt_id_geninfo = int(rtree_geninfo.GetBranch('evtID').GetLeaf('evtID').GetValue())
        if print_all:
            print("\n--------read geninfo Tree (evtID={0:d}):--------".format(evt_id_geninfo))

        # get the number of initial particles:
        n_init_particles = int(rtree_geninfo.GetBranch('nInitParticles').GetLeaf('nInitParticles').GetValue())

        # preallocate array for initial pdgid:
        init_pdgid_arr = np.array([])
        init_trkid_arr = np.array([])

        # preallocate number of initial neutrons:
        number_init_neutrons = 0

        flag_pion_plus = False

        # get initial PDG IDs:
        for index in range(n_init_particles):
            init_pdgid = int(rtree_geninfo.GetBranch('InitPDGID').GetLeaf('InitPDGID').GetValue(index))
            if print_all:
                print("Initial PDGID = {0:d}".format(init_pdgid))
            init_pdgid_arr = np.append(init_pdgid_arr, init_pdgid)

            if init_pdgid == 2112:
                number_init_neutrons += 1

            if init_pdgid == 211:
                flag_pion_plus = True

            init_trkid = int(rtree_geninfo.GetBranch('InitTRKID').GetLeaf('InitTRKID').GetValue(index))
            if print_all:
                print("Initial TRKID = {0:d}".format(init_trkid))
            init_trkid_arr = np.append(init_trkid_arr, init_trkid)

        # if not flag_pion_plus:
        #     continue

        """ read nCapture tree """
        # get the current event in the nCapture tree:
        rtree_ncapture.GetEntry(event)

        # get evtID of nCapture tree:
        evt_id_ncapture = int(rtree_ncapture.GetBranch('evtID').GetLeaf('evtID').GetValue())
        if print_all:
            print("\n--------read nCapture Tree (evtID={0:d}):--------".format(evt_id_ncapture))

        # get the number of nCaptures in the event:
        number_ncapture = int(rtree_ncapture.GetBranch('NeutronN').GetLeaf('NeutronN').GetValue())
        if print_all:
            print("Number of nCaptures from nCapture tree = {0:d}".format(number_ncapture))

        # preallocate trkid of neutron:
        neutron_trkid_arr = np.array([])
        # get TrkID of neutron:
        for index in range(number_ncapture):
            neutron_trkid = int(rtree_ncapture.GetBranch('NeutronTrkid').GetLeaf('NeutronTrkid').GetValue(index))
            if print_all:
                print("TrkID of neutron from nCapture Tree = {0:d}".format(neutron_trkid))
            neutron_trkid_arr = np.append(neutron_trkid_arr, neutron_trkid)

        """ read secondaries tree: """
        # set secondaries flag:
        flag_secondaries = False

        # preallocate arrays:
        pri_trkid_arr = np.array([])
        pri_pdgid_arr = np.array([])
        pri_kine_arr = np.array([])
        pri_time_arr = np.array([])
        sec_trkid_list = []
        sec_pdgid_list = []
        sec_process_name_list = []
        sec_kine_list = []
        sec_time_list = []

        # get the current event in the geninfo tree:
        rtree_sec.GetEntry(event)

        # get the value of the event ID:
        evt_id = int(rtree_sec.GetBranch('evtID').GetLeaf('evtID').GetValue())
        if print_all:
            print("\n--------read secondaries Tree (evtID={0:d}):--------".format(evt_id))

        # get the number of primaries:
        num_primaries = int(rtree_sec.GetBranch('numPrimaries').GetLeaf('numPrimaries').GetValue())

        # get the number of entries of the variable SecTrkID:
        num_secondaries = int(rtree_sec.GetBranch('numSec').GetLeaf('numSec').GetValue())

        if num_primaries == 0 and num_secondaries == 0:
            if print_all:
                print("NO particles in event\n")
            continue

        elif num_secondaries == 0:
            if print_all:
                print("NO secondaries in event\n")
            continue

        else:
            # set secondaries flag:
            flag_secondaries = True

        if flag_secondaries:
            # loop over primary particles:
            for num in range(num_primaries):
                # there are secondaries in the event:

                # get PrimaryTrkID:
                pri_trkid = int(rtree_sec.GetBranch("PrimaryTrkID").GetLeaf("PrimaryTrkID").GetValue(num))
                pri_trkid_arr = np.append(pri_trkid_arr, pri_trkid)

                # get PrimaryPDGID:
                pri_pdgid = int(rtree_sec.GetBranch("PrimaryPDGID").GetLeaf("PrimaryPDGID").GetValue(num))
                pri_pdgid_arr = np.append(pri_pdgid_arr, pri_pdgid)

                # get kinetic energy:
                pri_kine = float(rtree_sec.GetBranch("PrimaryKine").GetLeaf("PrimaryKine").GetValue(num))
                pri_kine_arr = np.append(pri_kine_arr, pri_kine)

                # get stop position:
                # pri_stopx = float(rtree_sec.GetBranch("PrimaryStopX").GetLeaf("PrimaryStopX").GetValue(num))
                # pri_stopy = float(rtree_sec.GetBranch("PrimaryStopY").GetLeaf("PrimaryStopY").GetValue(num))
                # pri_stopz = float(rtree_sec.GetBranch("PrimaryStopZ").GetLeaf("PrimaryStopZ").GetValue(num))

                # get time:
                pri_time = float(rtree_sec.GetBranch("PrimaryTime").GetLeaf("PrimaryTime").GetValue(num))
                pri_time_arr = np.append(pri_time_arr, pri_time)

                # get track length:
                # pri_trklength = float(rtree_sec.GetBranch("PrimaryTrkLen").GetLeaf("PrimaryTrkLen").GetValue(num))

                # preallocate arrays:
                sec_trkid_arr = np.array([])
                sec_pdgid_arr = np.array([])
                sec_process_name_arr = np.array([])
                sec_kine_arr = np.array([])
                sec_time_arr = np.array([])

                # loop over possible secondaries:
                for num1 in range(num_secondaries):

                    # get parent track ID ParentTrkID:
                    parent_trkid = int(rtree_sec.GetBranch("ParentTrkID").GetLeaf("ParentTrkID").GetValue(num1))

                    if pri_trkid != parent_trkid:
                        # secondary do NOT correspond to primary:
                        continue
                    else:
                        # secondary correspond to primary:
                        sec_trkid = int(rtree_sec.GetBranch("SecTrkID").GetLeaf("SecTrkID").GetValue(num1))
                        sec_trkid_arr = np.append(sec_trkid_arr, sec_trkid)

                        # get SecPDGID:
                        sec_pdgid = int(rtree_sec.GetBranch("SecPDGID").GetLeaf("SecPDGID").GetValue(num1))
                        sec_pdgid_arr = np.append(sec_pdgid_arr, sec_pdgid)

                        # get SecProcess (number of process):
                        sec_process = int(rtree_sec.GetBranch("SecProcess").GetLeaf("SecProcess").GetValue(num1))
                        # get name of process:
                        sec_process_name = get_process_name(sec_process)
                        sec_process_name_arr = np.append(sec_process_name_arr, sec_process_name)

                        # get SecKine:
                        sec_kine = float(rtree_sec.GetBranch("SecKine").GetLeaf("SecKine").GetValue(num1))
                        sec_kine_arr = np.append(sec_kine_arr, sec_kine)

                        # get position of secondary:
                        # sec_posx = float(rtree_sec.GetBranch("SecPosX").GetLeaf("SecPosX").GetValue(num1))
                        # sec_posy = float(rtree_sec.GetBranch("SecPosY").GetLeaf("SecPosY").GetValue(num1))
                        # sec_posz = float(rtree_sec.GetBranch("SecPosZ").GetLeaf("SecPosZ").GetValue(num1))

                        # get SecTime:
                        sec_time = float(rtree_sec.GetBranch("SecTime").GetLeaf("SecTime").GetValue(num1))
                        sec_time_arr = np.append(sec_time_arr, sec_time)

                # append arrays to lists:
                sec_trkid_list.append(sec_trkid_arr)
                sec_pdgid_list.append(sec_pdgid_arr)
                sec_process_name_list.append(sec_process_name_arr)
                sec_kine_list.append(sec_kine_arr)
                sec_time_list.append(sec_time_arr)

        else:
            continue

        # now the info about all primary particles is stored in arrays and the information about the secondaries is
        # stored in lists of arrays. For example: to get the pdgid of secondary 2 of primary 4: sec_pdgid_list[3][1]

        # first check, if there is nCapture in this event:
        num_n_capture = 0
        index_primary_n_capture = np.array([-1])

        # loop over all secondaries and search for nCapture:
        for index in range(num_primaries):
            for index1 in range(len(sec_trkid_list[index])):
                if sec_process_name_list[index][index1] != "nCapture":
                    continue
                else:
                    num_n_capture = num_n_capture + 1
                    # get index of corresponding primary particle:
                    # only append the index once, not twice for both secondaries (gamma and deuteron)
                    if index_primary_n_capture[-1] == index:
                        # if index is equal to last entry in the array, continue:
                        continue
                    else:
                        index_primary_n_capture = np.append(index_primary_n_capture, index)

        # divide num_n_capture by 2, because there are always TWO secondaries (gamma and deuteron) to ONE neutron
        # capture:
        num_n_capture = num_n_capture / 2
        if print_all:
            print("Number of nCaptures from secondaries tree = {0:d}\n".format(num_n_capture))

        # remove the first entry of index_primary_n_capture (the value -1), because it was just added as place holder:
        index_primary_n_capture = np.delete(index_primary_n_capture, 0)

        # check if number of nCaptures from nCapture tree agree with secondaries tree:
        if num_n_capture != number_ncapture:
            print("\nWARNING:             Different numbers of nCaptures: nCapture tree = {0:d}, secondaries tree = "
                  "{1:d}\n\n".format(number_ncapture, num_n_capture))

        # if num_n_capture == number_init_neutrons:
        #     print("**************   number of initial neutrons = number of neutron captures!    ***************\n")
        #     continue

        if num_n_capture == 0:
            if print_all:
                print("No neutron capture in the event!\n")
            continue

        # loop over number of nCaptures from secondaries tree:
        for index in range(num_n_capture):

            # set flag, that primary particle is equal to initial particle:
            flag_primary_equal_initial = False

            # set flag, that initial particle of neutron capture is no neutron:
            flag_initial_no_neutron = False

            # preallocate array, where indices of the primaries are saved:
            index_primary_array = np.array([])
            # preallocate list, where indices of secondaries corresponding to one primary are saved in arrays:
            index_secondary_list = []

            # check if pri_pdgid is neutron:
            if pri_pdgid_arr[index_primary_n_capture[index]] != 2112:
                print("WARNING: primary particles of nCapture is NO neutron")
            else:
                # get info of 'last' primary before nCapture (must be a neutron):
                if print_all:
                    print("nCapture of:")
                    print("PDG ID of ncapture pri. = {0:.0f}".format(pri_pdgid_arr[index_primary_n_capture[index]]))
                    print("TrkID of ncapture pri. = {0:.0f}".format(pri_trkid_arr[index_primary_n_capture[index]]))
                    print("KinE of ncapture pri. = {0:.3f} MeV".format(pri_kine_arr[index_primary_n_capture[index]]))
                    print("time of ncapture pri. = {0:.3f} ns".format(pri_time_arr[index_primary_n_capture[index]]))

                # set values of 'last' primary before nCapture (must be a neutron):
                trkid_to_be_check = pri_trkid_arr[index_primary_n_capture[index]]

                while not flag_primary_equal_initial:
                    # check if primary trkid of neutron before ncapture is equal to initial trkid:
                    for index1 in range(n_init_particles):
                        if pri_trkid_arr[index_primary_n_capture[index]] == init_trkid_arr[index1]:
                            # neutron is initial particle:
                            flag_primary_equal_initial = True
                            if print_all:
                                print("\nInitial Neutron is captured (TrkID = {0:.0f})\n\n"
                                      .format(init_trkid_arr[index1]))

                            if print_initials:
                                # increment number_initial_neutron:
                                number_initial_neutron += 1

                            # break from FOR loop
                            break

                    if flag_primary_equal_initial:
                        # break from WHILE loop
                        break

                    # loop over all particles to find the secondary with same trkid than primary:
                    for index2 in range(num_primaries):
                        for index3 in range(len(sec_trkid_list[index2])):
                            if sec_trkid_list[index2][index3] != trkid_to_be_check:
                                # check next secondary
                                continue
                            else:
                                # trk id's are equal:
                                # index of the primary to array:
                                index_primary_array = np.append(index_primary_array, index2)
                                # preallocate array, where indices of secondaries are saved:
                                index_secondary_array = np.array([])

                                # get the process name which created the secondary:
                                if print_all:
                                    print("\nProcess name = {0}".format(sec_process_name_list[index2][index3]))
                                    print("PrimaryPDGID of this process = {0:.0f}".format(pri_pdgid_arr[index2]))
                                    print("PrimaryTrkID of this process = {0:.0f}".format(pri_trkid_arr[index2]))
                                    print("PrimaryKine of this process = {0:.3f} MeV".format(pri_kine_arr[index2]))
                                    print("PrimaryTime of this process = {0:.3f} ns".format(pri_time_arr[index2]))
                                    print("Secondaries of this process:")
                                    print("SecPDGID = {0:.0f}".format(sec_pdgid_list[index2][index3]))
                                    print("SecTrkID = {0:.0f}".format(sec_trkid_list[index2][index3]))

                                # append index of secondary:
                                index_secondary_array = np.append(index_secondary_array, index3)

                                # get all other secondaries that are created from this primary:
                                if print_all:
                                    print("Other corresponding secondaries:")

                                for index5 in range(len(sec_trkid_list[index2])):
                                    if sec_trkid_list[index2][index5] == sec_trkid_list[index2][index3]:
                                        continue
                                    elif sec_pdgid_list[index2][index5] == 11:
                                        # do not print info about electrons:
                                        continue
                                    elif sec_pdgid_list[index2][index5] == 22:
                                        # do not print info about gammas:
                                        continue
                                    else:
                                        # append index of other secondary:
                                        index_secondary_array = np.append(index_secondary_array, index5)
                                        if print_all:
                                            print("SecPDGID = {0:.0f}".format(sec_pdgid_list[index2][index5]))
                                            print("SecTrkID = {0:.0f}".format(sec_trkid_list[index2][index5]))

                                # append array to list:
                                index_secondary_list.append(index_secondary_array)

                                # get the trkid of the primary corresponding to this secondary:
                                trkid_to_be_check = pri_trkid_arr[index2]

                                # check if this primary trkid is equal to initial trkid:
                                for index4 in range(n_init_particles):
                                    if trkid_to_be_check == init_trkid_arr[index4]:
                                        # primary is initial particle:
                                        flag_primary_equal_initial = True
                                        if print_all:
                                            print("\nPrimaryTrkID equal to InitTRKID = {0:.0f}"
                                                  .format(init_trkid_arr[index4]))
                                            print("InitPDGID = {0:.0f}\n\n".format(init_pdgid_arr[index4]))

                                        # check, if initial particle of nCapture is no neutron:
                                        if init_pdgid_arr[index4] != 2112:
                                            flag_initial_no_neutron = True

                                        else:
                                            if print_initials:
                                                # initial particle of nCapture is neutron:
                                                number_initial_neutron += 1

                                        # break from FOR loop
                                        break

                                if flag_primary_equal_initial:
                                    # break from For loop
                                    break

                        if flag_primary_equal_initial:
                            # break from FOR loop
                            break

                    if flag_primary_equal_initial:
                        # break from WHILE loop
                        break

            if flag_initial_no_neutron and print_process:
                print("************** event {0:d} ****************".format(event))

                # when initial is no neutron and print_process = True:
                for index6 in range(len(index_primary_array)):
                    print("\nneutron capture with no neutron as initial particle:")
                    print("TrkID of ncapture pri. = {0:.0f}".format(pri_trkid_arr[index_primary_n_capture[index]]))
                    print("PrimaryPDGID of this process = {0:.0f}"
                          .format(pri_pdgid_arr[int(index_primary_array[index6])]))
                    print("PrimaryTrkID of this process = {0:.0f}"
                          .format(pri_trkid_arr[int(index_primary_array[index6])]))
                    print("PrimaryKine of this process = {0:.3f} MeV"
                          .format(pri_kine_arr[int(index_primary_array[index6])]))
                    print("PrimaryTime of this process = {0:.3f} ns"
                          .format(pri_time_arr[int(index_primary_array[index6])]))

                    for index7 in range(len(index_secondary_list[index6])):

                        print("Process name = {0}"
                              .format(sec_process_name_list[int(index_primary_array[index6])][int(index_secondary_list[index6][index7])]))
                        print("Secondaries of this process:")
                        print("SecPDGID = {0:.0f}"
                              .format(sec_pdgid_list[int(index_primary_array[index6])][int(index_secondary_list[index6][index7])]))
                        print("SecTrkID = {0:.0f}"
                              .format(sec_trkid_list[int(index_primary_array[index6])][int(index_secondary_list[index6][index7])]))

            if flag_initial_no_neutron and print_initials:
                print("************** event {0:d} ****************".format(event))

                # preallocate variable, where PDG ID of initial particle is saved:
                initial_pdg_of_n_capture = 0

                # when initial is no neutron and print_process = True:
                for index6 in range(len(index_primary_array)):
                    # set variable, where PDG ID of initial particle is saved:
                    initial_pdg_of_n_capture = pri_pdgid_arr[int(index_primary_array[index6])]

                if initial_pdg_of_n_capture == 2212:
                    number_initial_proton += 1
                elif initial_pdg_of_n_capture == 211:
                    number_initial_piplus += 1
                elif initial_pdg_of_n_capture == -211:
                    number_initial_piminus += 1
                elif initial_pdg_of_n_capture == 111:
                    number_initial_pizero += 1
                else:
                    print("other initial cause nCapture:    PDG ID = {0:.0f}".format(initial_pdg_of_n_capture))

    return (number_initial_neutron, number_initial_proton, number_initial_piminus, number_initial_piplus,
            number_initial_pizero)


# get the date and time, when the script was run:
date = datetime.datetime.now()
now = date.strftime("%Y-%m-%d %H:%M")

# set the path of the input files (filename must be 'user_atmoNC_{}.root'):
# Input_path = "/home/astro/blum/juno/atmoNC/data_NC/secondary_anamgr/"
Input_path = "/local/scratch1/pipc51/astro/blum/detsim_output_data_noopt/"

# set path, where results should be saved:
# Output_path = "/home/astro/blum/juno/atmoNC/data_NC/secondary_anamgr/"
Output_path = "/local/scratch1/pipc51/astro/blum/detsim_output_data_noopt/"

# set the file number, the number of the first event and number of the last event that should be read:
file_number = 95
first_event = 326
last_event = 326

# set print flags:
Print_all = True
Print_process = True
Print_initials = True

# path to file:
# Input_file = Input_path + ".root"
Input_file = Input_path + "user_atmoNC_noopt_{0:d}.root".format(file_number)


print("########################## secondaries Tree #############################")

Number_initial_neutron, Number_initial_proton, Number_initial_piminus, Number_initial_piplus, Number_initial_pizero = \
    check_secondaries_tree(Input_file, first_event, last_event, Print_all, Print_process, Print_initials)

print("Initial particles that cause neutron capture:")
print("number of neutron = {0:.0f}".format(Number_initial_neutron))
print("number of proton = {0:.0f}".format(Number_initial_proton))
print("number of pi_minus = {0:.0f}".format(Number_initial_piminus))
print("number of pi_plus = {0:.0f}".format(Number_initial_piplus))
print("number of pi_zero = {0:.0f}".format(Number_initial_pizero))
