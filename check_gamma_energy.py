""" Script to check gammas (with energies of 1.9 MeV, 2.2 MeV and 2.5 MeV), which were simulated with tut_detsim.py of
    JUNO offline version J18v2r1-branch.

    For every simulated gamma, take totalPE, edep and Qedep and compare them

"""
import ROOT
import NC_background_functions
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit


def conversion(npe, radius):
    """

    :param npe: number pe
    :param radius: radial position in mm
    :return:
    """
    # if radius < 8000:
    #     evis = 0.0007524 * npe
    # elif 8000 <= radius < 12000:
    #     evis = 0.0007220 * npe
    # elif 12000 <= radius < 14000:
    #     evis = 0.0006927 * npe
    # elif 14000 <= radius < 16000:
    #     evis = 0.0006734 * npe
    # elif 16000 <= radius < 17000:
    #     evis = 0.0007227 * npe
    # else:
    #     evis = 0.0007877 * npe

    evis = 0.0007005 * npe

    return evis


# set the path of the input files:
input_path = "/local/scratch1/pipc51/astro/blum/gamma_2_2_MeV/"

# set path, where results should be saved:
output_path = "/home/astro/blum/juno/atmoNC/data_NC/output_gamma_2_2_MeV/"

# set the number of the first file and number of the last file (of 2.2 MeV gamma -> stop_number_2_2,
# of 1.9 MeV and 2.5 MeV gamma -> stop_number_rest) that should be read:
start_number = 0
stop_number_2_2 = 199
stop_number_rest = 99
# number of entries in the input files:
Number_entries_input = 10
# total number of events:
number_events_total_2_2 = (stop_number_2_2 - start_number + 1) * Number_entries_input
number_events_total_rest = (stop_number_rest - start_number + 1) * Number_entries_input

# set the radius for the volume cut in mm:
r_cut = 16000

energy_positron = []
for number in range(100):
    input_file = "/local/scratch1/pipc51/astro/blum/positron_output_CDcenter/user_positron_50MeV_{0:d}.root"\
        .format(number)
    # load ROOT file:
    rfile = ROOT.TFile(input_file)
    # get the "evt"-TTree from the TFile:
    rtree_evt = rfile.Get("evt")
    # get geninfo tree
    rtree_geninfo = rfile.Get("geninfo")

    # loop over every event, i.e. every entry, in the TTree:
    for event in range(50):
        # get the current event in the TTree:
        rtree_geninfo.GetEntry(event)

        # get position of 0th initial particle in x, y, z in mm (positions of the initial particles are equal):
        x_init = float(rtree_geninfo.GetBranch('InitX').GetLeaf('InitX').GetValue(0))
        y_init = float(rtree_geninfo.GetBranch('InitY').GetLeaf('InitY').GetValue(0))
        z_init = float(rtree_geninfo.GetBranch('InitZ').GetLeaf('InitZ').GetValue(0))

        # calculate the radius in mm:
        r_init = np.sqrt(x_init**2 + y_init**2 + z_init**2)

        # get current event in evt tree:
        rtree_evt.GetEntry(event)

        # get total number of pe of the event:
        totalPE = int(rtree_evt.GetBranch("totalPE").GetLeaf("totalPE").GetValue())
        # totalPE = NC_background_functions.conversion_npe_to_evis(totalPE)
        totalPE = conversion(totalPE, r_init)

        energy_positron.append(totalPE)

""" preallocate arrays for 2.2 MeV gamma: """
# number of PE of each event:
number_pe_2_2 = []
edep_2_2 = []
qedep_2_2 = []
print("\nstart reading 2.2 MeV gamma files...")

# loop over files of gamma = 2.2 MeV:
for number in range(start_number, stop_number_2_2+1, 1):
    # path to file:
    input_file = input_path + "user_gamma_2_2_MeV_{0:d}.root".format(number)

    # load ROOT file:
    rfile = ROOT.TFile(input_file)
    # get the "evt"-TTree from the TFile:
    rtree_evt = rfile.Get("evt")
    # get "geninfo" tree from root file:
    rtree_geninfo = rfile.Get("geninfo")
    # get "prmtrkdep" tree from root file:
    rtree_prmtrkdep = rfile.Get("prmtrkdep")

    # loop over every event, i.e. every entry, in the TTree:
    for event in range(Number_entries_input):

        """ check volume cut: """
        # get current event in geninfo tree:
        rtree_geninfo.GetEntry(event)

        # get initial x, y, z position:
        x_init = float(rtree_geninfo.GetBranch('InitX').GetLeaf('InitX').GetValue())
        y_init = float(rtree_geninfo.GetBranch('InitY').GetLeaf('InitY').GetValue())
        z_init = float(rtree_geninfo.GetBranch('InitZ').GetLeaf('InitZ').GetValue())

        # distance to detector center in mm:
        r_init = np.sqrt(x_init**2 + y_init**2 + z_init**2)

        # do volume cut:
        if r_init >= r_cut:
            continue

        # get current event in evt tree:
        rtree_evt.GetEntry(event)

        # get total number of pe of the event:
        totalPE = int(rtree_evt.GetBranch("totalPE").GetLeaf("totalPE").GetValue())
        # totalPE = NC_background_functions.conversion_npe_to_evis(totalPE)
        totalPE = conversion(totalPE, r_init)
        number_pe_2_2.append(totalPE)

        # get current event in prmtrkdep tree:
        rtree_prmtrkdep.GetEntry(event)

        # get quenched deposited energy in MeV:
        qedep = float(rtree_prmtrkdep.GetBranch("Qedep").GetLeaf("Qedep").GetValue())
        qedep_2_2.append(qedep)

        # get deposited energy in MeV
        edep = float(rtree_prmtrkdep.GetBranch("edep").GetLeaf("edep").GetValue())
        edep_2_2.append(edep)

""" preallocate arrays for 1.9 MeV gamma: """
# number of PE of each event:
number_pe_1_9 = []
edep_1_9 = []
qedep_1_9 = []
print("\nstart reading 1.9 MeV gamma files...")

# loop over files of gamma = 1.9 MeV:
for number in range(start_number, stop_number_rest+1, 1):
    # path to file:
    input_file = input_path + "user_gamma_1_9_MeV_{0:d}.root".format(number)

    # load ROOT file:
    rfile = ROOT.TFile(input_file)
    # get the "evt"-TTree from the TFile:
    rtree_evt = rfile.Get("evt")
    # get "geninfo" tree from root file:
    rtree_geninfo = rfile.Get("geninfo")
    # get "prmtrkdep" tree from root file:
    rtree_prmtrkdep = rfile.Get("prmtrkdep")

    # loop over every event, i.e. every entry, in the TTree:
    for event in range(Number_entries_input):

        """ check volume cut: """
        # get current event in geninfo tree:
        rtree_geninfo.GetEntry(event)

        # get initial x, y, z position:
        x_init = float(rtree_geninfo.GetBranch('InitX').GetLeaf('InitX').GetValue())
        y_init = float(rtree_geninfo.GetBranch('InitY').GetLeaf('InitY').GetValue())
        z_init = float(rtree_geninfo.GetBranch('InitZ').GetLeaf('InitZ').GetValue())

        # distance to detector center in mm:
        r_init = np.sqrt(x_init**2 + y_init**2 + z_init**2)

        # do volume cut:
        if r_init >= r_cut:
            continue

        # get current event in evt tree:
        rtree_evt.GetEntry(event)

        # get total number of pe of the event:
        totalPE = int(rtree_evt.GetBranch("totalPE").GetLeaf("totalPE").GetValue())
        # totalPE = NC_background_functions.conversion_npe_to_evis(totalPE)
        totalPE = conversion(totalPE, r_init)
        number_pe_1_9.append(totalPE)

        # get current event in prmtrkdep tree:
        rtree_prmtrkdep.GetEntry(event)

        # get quenched deposited energy in MeV:
        qedep = float(rtree_prmtrkdep.GetBranch("Qedep").GetLeaf("Qedep").GetValue())
        qedep_1_9.append(qedep)

        # get deposited energy in MeV
        edep = float(rtree_prmtrkdep.GetBranch("edep").GetLeaf("edep").GetValue())
        edep_1_9.append(edep)

""" preallocate arrays for 2.5 MeV gamma: """
# number of PE of each event:
number_pe_2_5 = []
edep_2_5 = []
qedep_2_5 = []
print("\nstart reading 2.5 MeV gamma files...")

# loop over files of gamma = 2.5 MeV:
for number in range(start_number, stop_number_rest+1, 1):
    # path to file:
    input_file = input_path + "user_gamma_2_5_MeV_{0:d}.root".format(number)

    # load ROOT file:
    rfile = ROOT.TFile(input_file)
    # get the "evt"-TTree from the TFile:
    rtree_evt = rfile.Get("evt")
    # get "geninfo" tree from root file:
    rtree_geninfo = rfile.Get("geninfo")
    # get "prmtrkdep" tree from root file:
    rtree_prmtrkdep = rfile.Get("prmtrkdep")

    # loop over every event, i.e. every entry, in the TTree:
    for event in range(Number_entries_input):

        """ check volume cut: """
        # get current event in geninfo tree:
        rtree_geninfo.GetEntry(event)

        # get initial x, y, z position:
        x_init = float(rtree_geninfo.GetBranch('InitX').GetLeaf('InitX').GetValue())
        y_init = float(rtree_geninfo.GetBranch('InitY').GetLeaf('InitY').GetValue())
        z_init = float(rtree_geninfo.GetBranch('InitZ').GetLeaf('InitZ').GetValue())

        # distance to detector center in mm:
        r_init = np.sqrt(x_init**2 + y_init**2 + z_init**2)

        # do volume cut:
        if r_init >= r_cut:
            continue

        # get current event in evt tree:
        rtree_evt.GetEntry(event)

        # get total number of pe of the event:
        totalPE = int(rtree_evt.GetBranch("totalPE").GetLeaf("totalPE").GetValue())
        # totalPE = NC_background_functions.conversion_npe_to_evis(totalPE)
        totalPE = conversion(totalPE, r_init)
        number_pe_2_5.append(totalPE)

        # get current event in prmtrkdep tree:
        rtree_prmtrkdep.GetEntry(event)

        # get quenched deposited energy in MeV:
        qedep = float(rtree_prmtrkdep.GetBranch("Qedep").GetLeaf("Qedep").GetValue())
        qedep_2_5.append(qedep)

        # get deposited energy in MeV
        edep = float(rtree_prmtrkdep.GetBranch("edep").GetLeaf("edep").GetValue())
        edep_2_5.append(edep)

# set energy range in MeV:
energy_MeV = np.arange(1.0, 3.0, 0.02)
# set energy range in PE:
energy_PE = np.arange(1300, 4000, 30)

h1 = plt.figure(1)
plt.hist(number_pe_1_9, energy_MeV, color='red', label='{0:.0f} gammas with 1.9 MeV\nmean = {1:.2f}, sigma = {2:.2f}'
         .format(len(number_pe_1_9), np.mean(number_pe_1_9), np.std(number_pe_1_9)))
plt.hist(number_pe_2_2, energy_MeV, color='blue', label='{0:.0f} gammas with 2.2 MeV\nmean = {1:.2f}, sigma = {2:.2f}'
         .format(len(number_pe_2_2), np.mean(number_pe_2_2), np.std(number_pe_2_2)))
plt.hist(number_pe_2_5, energy_MeV, color='green', label='{0:.0f} gammas with 2.5 MeV\nmean = {1:.2f}, sigma = {2:.2f}'
         .format(len(number_pe_2_5), np.mean(number_pe_2_5), np.std(number_pe_2_5)))
plt.xlabel("number of PE converted to MeV (per $\\gamma$)")
plt.ylabel("entries per bin")
plt.title("number of PE of 1.9 MeV, 2.2 MeV and 2.5 MeV gamma's")
plt.legend()
plt.grid()

h2 = plt.figure(2)
plt.hist(edep_1_9, energy_MeV, color='red', label='{0:.0f} gammas with 1.9 MeV\nmean = {1:.2f}, sigma = {2:.2f}'
         .format(len(edep_1_9), np.mean(edep_1_9), np.std(edep_1_9)))
plt.hist(edep_2_2, energy_MeV, color='blue', label='{0:.0f} gammas with 2.2 MeV\nmean = {1:.2f}, sigma = {2:.2f}'
         .format(len(edep_2_2), np.mean(edep_2_2), np.std(edep_2_2)))
plt.hist(edep_2_5, energy_MeV, color='green', label='{0:.0f} gammas with 2.5 MeV\nmean = {1:.2f}, sigma = {2:.2f}'
         .format(len(edep_2_5), np.mean(edep_2_5), np.std(edep_2_5)))
plt.xlabel("deposited energy in MeV (per $\\gamma$)")
plt.ylabel("entries per bin")
plt.title("Deposited energy of 1.9 MeV, 2.2 MeV and 2.5 MeV gamma's")
plt.legend()
plt.grid()

h3 = plt.figure(3)
plt.hist(qedep_1_9, energy_MeV, color='red', label='{0:.0f} gammas with 1.9 MeV\nmean = {1:.2f}, sigma = {2:.2f}'
         .format(len(qedep_1_9), np.mean(qedep_1_9), np.std(qedep_1_9)))
plt.hist(qedep_2_2, energy_MeV, color='blue', label='{0:.0f} gammas with 2.2 MeV\nmean = {1:.2f}, sigma = {2:.2f}'
         .format(len(qedep_2_2), np.mean(qedep_2_2), np.std(qedep_2_2)))
plt.hist(qedep_2_5, energy_MeV, color='green', label='{0:.0f} gammas with 2.5 MeV\nmean = {1:.2f}, sigma = {2:.2f}'
         .format(len(qedep_2_5), np.mean(qedep_2_5), np.std(qedep_2_5)))
plt.xlabel("Quenched deposited energy in MeV (per $\\gamma$)")
plt.ylabel("entries per bin")
plt.title("Quenched deposited energy of 1.9 MeV, 2.2 MeV and 2.5 MeV gamma's")
plt.legend()
plt.grid()

h4 = plt.figure(4)
plt.hist(energy_positron, np.arange(45, 55, 0.1), label='{0:.0f} positron with 50 MeV\nmean = {1:.2f}, sigma = {2:.2f}'
         .format(len(energy_positron), np.mean(energy_positron), np.std(energy_positron)))
plt.xlabel("number of PE converted to MeV ")
plt.ylabel("entries per bin")
plt.title("number of PE of positrons with kinetic energy of 50 MeV")
plt.legend()
plt.grid()

plt.show()



