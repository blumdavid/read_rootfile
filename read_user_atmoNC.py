""" Script to read user_atmoNC_...root file from JUNO offline detector simulation (tut_detsim.py) and to
    calculate the visible spectrum of events that mimic IBD events:

    - file sample_detsim_user.root is generated with the JUNO detector simulation:
        python $TUTORIALROOT/share/tut_detsim.py
    - the hepevt file from the DSNB-NC.exe generator is used as input for the detector simulation

    This ROOT file is read and analyzed to get the spectrum of visible energy of the NC QEL events that mimic an IBD
    signal in the JUNO detector.

    Criteria for IBD signal (from JUNO PhysicsReport page 39):
        1.  fiducial volume cut: r < 17 m
        2.  prompt energy cut: 10 MeV < E_prompt < 105 MeV
        3.  delayed energy cut: 1.9 MeV < E_delayed < 2.5 MeV (neutron is captured on Carbon and releases gamma of
            around 2.2 MeV)
        4.  time cut: time interval between prompt and delayed signal: 600 ns < delta_T < 1.0 ms
            (Julia has taken 600 ns, in PhysicsReport there is no further specification)
        5.  distance cut between prompt and delayed signal: R_(prompt-delayed) < 1.5 m
        6.  neutron multiplicity cut: only 1 delayed signal (neutron capture with 1.9 MeV to 2.5 MeV) in the time window
            (600 ns to 1.0 ms)

        6. Muon veto criteria (83% efficiency by yellow book): -> Do I have to consider this??????????????

"""
import ROOT
from array import array
import sys
import numpy as np
from matplotlib import pyplot as plt
import xml.etree.ElementTree as ET
import datetime
import NC_background_functions
from scipy.optimize import curve_fit

# get the date and time, when the script was run:
date = datetime.datetime.now()
now = date.strftime("%Y-%m-%d %H:%M")

# set the path of the input files:
# input_path = "/local/scratch1/pipc51/astro/blum/detsim_output_data/"
input_path = "/home/astro/blum/juno/atmoNC/data_NC/output_neutron_cut/proton_B11_check/"

""" set the number of the first file and number of the last file that should be read: """
start_number = 0
stop_number = 0
# number of entries in the input files:
# number_entries_input = 100
number_entries_input = 10

# file number:
file_number = 11
# event number:
event_number = 2

# load the ROOT file:
# rfile = ROOT.TFile(input_path + "user_atmoNC_{0:d}.root".format(file_number))
rfile = ROOT.TFile(input_path + "user_p_atmoNC3_evt83.root")

# get the "evt"-TTree from the TFile:
rtree_evt = rfile.Get("evt")
# get the "geninfo"-TTree from the TFile:
rtree_geninfo = rfile.Get("geninfo")
# get the "prmtrkdep"-TTree from the TFile:
rtree_prmtrkdep = rfile.Get("prmtrkdep")

# get the number of events in the 'evt' Tree:
number_events_evt = rtree_evt.GetEntries()
# get the number of events in the geninfo Tree:
number_events_geninfo = rtree_geninfo.GetEntries()
# get the number of events in the prmtrkdep Tree:
number_events_prmtrkdep = rtree_prmtrkdep.GetEntries()
if number_events_geninfo == number_events_prmtrkdep and number_events_geninfo == number_events_evt:
    number_events = number_events_geninfo
else:
    sys.exit("ERROR: number of events in t Trees are NOT equal!!")

# check if number_events is equal to number_entries_input (if not, the detector simulation was incorrect!!):
if number_events != number_entries_input:
    sys.exit("ERROR: number of events are not equal to {0:d} -> Detector Simulation not correct!"
             .format(number_entries_input))

# loop over every event, i.e. every entry, in the TTree:
# for event in range(number_events):
for event in range(event_number, event_number + 1, 1):

    """ first read the "evt" Tree"""
    # get the current event in the TTree:
    rtree_evt.GetEntry(event)

    # get the value of the event ID:
    evt_id_evt = int(rtree_evt.GetBranch('evtID').GetLeaf('evtID').GetValue())

    # get number of photons of this event:
    n_photons = int(rtree_evt.GetBranch('nPhotons').GetLeaf('nPhotons').GetValue())
    print("number of photons in event = {0:d}".format(n_photons))

    # get edep in MeV of this event:
    edep = float(rtree_evt.GetBranch('edep').GetLeaf('edep').GetValue())

    # get edepX in mm of this event:
    edepX = float(rtree_evt.GetBranch('edepX').GetLeaf('edepX').GetValue())
    # get edepY in mm of this event:
    edepY = float(rtree_evt.GetBranch('edepY').GetLeaf('edepY').GetValue())
    # get edepZ in mm of this event:
    edepZ = float(rtree_evt.GetBranch('edepZ').GetLeaf('edepZ').GetValue())

    """ preallocate variables: """
    # number of pe in event:
    number_pe_event = 0
    # hittime of photons in event in ns:
    hittime_event = np.array([])

    # loop over every photon in the event:
    for index in range(n_photons):

        # get PMT ID, where photon is absorbed:
        pmt_id = int(rtree_evt.GetBranch('pmtID').GetLeaf('pmtID').GetValue(index))

        # only 20 inch PMTs (PMT ID of 20 inch PMTs are below 21000, PMT ID of 3 inch PMTs start at 290000):
        if pmt_id < 25000:
            # get nPE for this photon:
            n_pe = int(rtree_evt.GetBranch('nPE').GetLeaf('nPE').GetValue(index))
            # check, if photon produces only 1 PE:
            if n_pe != 1:
                print("{1:d} PE for 1 photon in event {0:d} in file {2}".format(evt_id_evt, n_pe, rootfile_input))

            # add n_pe to number_pe_event:
            number_pe_event = number_pe_event + n_pe

            # get hittime of this photon:
            hit_time = float(rtree_evt.GetBranch('hitTime').GetLeaf('hitTime').GetValue(index))
            # append hit_time to hittime_event:
            hittime_event = np.append(hittime_event, hit_time)

        else:
            continue

    print("evt {0:d} in file user_atmoNC_{1:d}.root:\n".format(evt_id_evt, file_number))
    print("number of PE = {0:0.1f}".format(number_pe_event))
    print("edep from evt-tree = {0:0.3f} MeV".format(edep))
    print("edepX from evt-tree = {0:0.3f} mm".format(edepX))
    print("edepY from evt-tree = {0:0.3f} mm".format(edepY))
    print("edepZ from evt-tree = {0:0.3f} mm\n".format(edepZ))

    """ read 'geninfo' tree: """
    rtree_geninfo.GetEntry(event)

    # get the value of the event ID:
    evt_id_geninfo = int(rtree_geninfo.GetBranch('evtID').GetLeaf('evtID').GetValue())

    # get number of particle in event:
    n_par = int(rtree_geninfo.GetBranch('nInitParticles').GetLeaf('nInitParticles').GetValue())

    for index in range(n_par):

        # get PDGID:
        PDG_ID = int(rtree_geninfo.GetBranch('InitPDGID').GetLeaf('InitPDGID').GetValue(index))
        print("PDG ID from geninfo = {0:d}".format(PDG_ID))

        # get initPx, initPy, initPz in MeV:
        InitPX = float(rtree_geninfo.GetBranch('InitPX').GetLeaf('InitPX').GetValue(index))
        InitPY = float(rtree_geninfo.GetBranch('InitPY').GetLeaf('InitPY').GetValue(index))
        InitPZ = float(rtree_geninfo.GetBranch('InitPZ').GetLeaf('InitPZ').GetValue(index))

        # calculate kinetic energy of initial particles in MeV:
        initE = np.sqrt(InitPX**2 + InitPY**2 + InitPZ**2)
        print("initial kinetic energy = {0:0.2f} MeV".format(initE))

    """ read 'prmtrkdep' tree: """
    rtree_prmtrkdep.GetEntry(event)

    # get number of particle in event:
    n_particle = int(rtree_prmtrkdep.GetBranch('nInitParticles').GetLeaf('nInitParticles').GetValue())

    for index in range(n_particle):

        # get edep in MeV of this event:
        edep_prm = float(rtree_prmtrkdep.GetBranch('edep').GetLeaf('edep').GetValue(index))

        # get edepX in mm of this event:
        edepX_prm = float(rtree_prmtrkdep.GetBranch('edepX').GetLeaf('edepX').GetValue(index))
        # get edepY in mm of this event:
        edepY_prm = float(rtree_prmtrkdep.GetBranch('edepY').GetLeaf('edepY').GetValue(index))
        # get edepZ in mm of this event:
        edepZ_prm = float(rtree_prmtrkdep.GetBranch('edepZ').GetLeaf('edepZ').GetValue(index))

        print("\nedep from prmtrkdep-tree = {0:0.3f} MeV".format(edep_prm))
        print("edepX from prmtrkdep-tree = {0:0.3f} mm".format(edepX_prm))
        print("edepY from prmtrkdep-tree = {0:0.3f} mm".format(edepY_prm))
        print("edepZ from prmtrkdep-tree = {0:0.3f} mm\n".format(edepZ_prm))

Bins = np.arange(0, 10000, 5)
Bins = np.arange(1000000, 100000000, 50)
plt.hist(hittime_event, Bins, align='mid', histtype='step')
plt.xlabel("hit-time in ns", fontsize=13)
# INFO-me: ylabel is only equal to number of PE, if nPE == 1 for all photons (1 PE each photon)
plt.ylabel("number of PE per bin", fontsize=13)
# plt.title("Example of PMT hit-time distribution of one 2.2 MeV $\\gamma$", fontsize=18)
plt.grid()
plt.show()
