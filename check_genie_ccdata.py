""" script to check the file genie_ccdata.root in folder /home/astro/blum/juno/JUNO-SOFT/data/Generator/DSNB/data/


"""
import ROOT
import numpy as np

# set the path of the inputs:
input_path = "/home/astro/blum/juno/JUNO-SOFT/data/Generator/DSNB/data/"

# file name of the input file:
input_name = input_path + "genie_ccdata.root"

# load the ROOT file:
rfile_input = ROOT.TFile(input_name)
# get the TTree from the TFile:
rtree_input = rfile_input.Get("particleT")

# get the number of entries in the ROOT-file:
number_entries = rtree_input.GetEntries()
# number_entries = 200

# set parameters:
pPdg_set = 12
tPdg_set = 1000060120
pdg_set = 11
E_min = 0.0
E_max = 100.0
number_interesting_events = 0
number_C11 = 0
number_B11 = 0
number_C10 = 0
number_Be10 = 0
number_B10 = 0

""" Read the data from the TTree: """
# loop over every entry, i.e. every event, in the TTree:
for event in range(0, number_entries):

    # get the current event in the TTree:
    rtree_input.GetEntry(event)

    # get PDG ID of neutrino:
    pPdg = int(rtree_input.GetBranch("pPdg").GetLeaf("pPdg").GetValue())
    if pPdg != pPdg_set:
        continue

    # get PDG ID of target:
    tPdg = int(rtree_input.GetBranch("tPdg").GetLeaf("tPdg").GetValue())
    if tPdg != tPdg_set:
        continue

    # get PDG ID of residual isotope:
    isoPdg = int(rtree_input.GetBranch("m_isoPdg").GetLeaf("m_isoPdg").GetValue())

    # get number of final particles:
    Npars = int(rtree_input.GetBranch("Npars").GetLeaf("Npars").GetValue())

    # preallocate arrays fof final particles:
    array_pdg = []
    array_E_kin = []
    flag_pions = False
    flag_kaons = False
    number_neutrons = 0

    # loop over Npars:
    for index in range(Npars):
        # get PDG ID of final particle:
        pdg = int(rtree_input.GetBranch("pdg").GetLeaf("pdg").GetValue(index))

        if pdg == -211 or pdg == 111 or pdg == 211:
            flag_pions = True

        if pdg == -321 or pdg == -311 or pdg == 311 or pdg == 321:
            flag_kaons = True

        if pdg == 2112:
            number_neutrons += 1

        # if final particle is electron, calculate kinetic energy:
        if pdg == pdg_set:
            # get momentum in x, y, z in MeV:
            px = float(rtree_input.GetBranch("px").GetLeaf("px").GetValue(index)) / 1000.0
            py = float(rtree_input.GetBranch("py").GetLeaf("py").GetValue(index)) / 1000.0
            pz = float(rtree_input.GetBranch("pz").GetLeaf("pz").GetValue(index)) / 1000.0

            # calculate kinetic energy in MeV:
            E_kin = np.sqrt(px**2 + py**2 + pz**2)
        else:
            # set E_kin to -1:
            E_kin = -1

        array_pdg.append(pdg)
        array_E_kin.append(E_kin)

    if flag_pions or flag_kaons:
        # only events without pions and without kaons
        continue

    if number_neutrons == 0 or number_neutrons > 2:
        # only events with 1 or 2 neutrons
        continue

    # only analyze event further that have 1 electron/positron between E_min and E_max:
    # flag, that electron energy is in energy window
    flag_electron = 0
    E_kin_total = 0.0

    for index in range(Npars):
        if array_pdg[index] == pdg_set:
            E_kin_total += array_E_kin[index]
        else:
            continue

    if E_kin_total < E_min or E_kin_total > E_max:
        # total electron energy not in energy window:
        continue
    else:
        number_interesting_events += 1

        if isoPdg == 1000060110:
            number_C11 += 1
        elif isoPdg == 1000050110:
            number_B11 += 1
        elif isoPdg == 1000060100:
            number_C10 += 1
        elif isoPdg == 1000050100:
            number_B10 += 1
        elif isoPdg == 1000040100:
            number_Be10 += 1
        else:
            print("m_isoPdg = {0:d}".format(isoPdg))
            print("pdg = {0}".format(array_pdg))

print("number C11 = {0:d}".format(number_C11))
print("number B11 = {0:d}".format(number_B11))
print("number C10 = {0:d}".format(number_C10))
print("number B10 = {0:d}".format(number_B10))
print("number Be10 = {0:d}".format(number_Be10))

print("number of events = {0:d}".format(number_entries))

print("number of events with nu_e + C12 -> electron + ... ({0:.0f} MeV <= E_kin_electron <= {1:.0f} MeV)"
      .format(E_min, E_max))
print(number_interesting_events)












