""" Script to check conversion from nPE to MeV of neutrons/protons/positrons, respectively, which were simulated
    with tut_detsim.py of JUNO offline version J18v1r1-pre1.

    Results of this script are used to get the visible energy of a particle depending on the number of PE AND
    depending on the position of the energy deposition in the detector.

    DIFFERENCE to check_conversion_npe_mev.py:
    also the position of the particle inside the detector is used for the conversion from nPE to MeV
    (see script check_DMsignal_simulation.py).

    With this conversion the cut on the energy of a possible prompt signal can be made in the PE-regime and efficiency
    of this cut can be calculated.

"""
import datetime
import ROOT
import NC_background_functions
import numpy as np
from matplotlib import pyplot as plt
from decimal import Decimal
from matplotlib.colors import LogNorm

# get the date and time, when the script was run:
date = datetime.datetime.now()
now = date.strftime("%Y-%m-%d %H:%M")


def get_npe_redep_qedep_from_file(input_path, filename, first_file, last_file, e_min, e_max):
    """
    function to get total number of PE, distance to the detector center in mm and quenched deposited energy
    (visible energy) in MeV of each event of each file

    :param input_path: path, where the user-root files are saved
    :param filename: name of the files that are read (e.g. user_neutron_500_MeV_)
    :param first_file: number of the first file to be read
    :param last_file: number of the last file to be read
    :param e_min: minimum visible energy in MeV, that should be analyzed
    :param e_max: maximum visible energy in MeV, that should be analyzed
    :return:
    """
    # preallocate array, where total number of PE is stored:
    array_total_pe = []
    # preallocate array, where distance to detector center of position, where energy is deposited, is stored (in mm):
    array_r_edep = []
    # preallocate array, where quenched deposited energy is stored (in MeV):
    array_qedep = []
    # preallocate array, where deposited energy is stored (in MeV):
    array_edep = []

    # loop over user-root files:
    for filenumber in range(first_file, last_file+1, 1):

        # input user-root file:
        input_file = input_path + filename + "{0:d}.root".format(filenumber)
        # load the ROOT file:
        rfile = ROOT.TFile(input_file)
        # get the "evt"-TTree from the TFile:
        rtree_evt = rfile.Get("evt")
        # get 'prmtrkdep' tree from TFile:
        rtree_prmtrkdep = rfile.Get("prmtrkdep")

        # get number of event of the file:
        number_events_file = rtree_evt.GetEntries()

        # loop ever each event:
        for event in range(number_events_file):

            # get the current event in prmtrkdep tree:
            rtree_prmtrkdep.GetEntry(event)

            # get number of initial particles:
            n_init_particles = int(rtree_prmtrkdep.GetBranch("nInitParticles").GetLeaf("nInitParticles").GetValue())

            # preallocate total quenched deposited energy of event in MeV:
            qedep = 0.0

            # loop over number of initial particles:
            for index in range(n_init_particles):
                # get quenched deposited energy (visible energy) in MeV:
                qedep_value = float(rtree_prmtrkdep.GetBranch("Qedep").GetLeaf("Qedep").GetValue(index))

                # add it to qedep:
                qedep += qedep_value

            # check if qedep is in the correct energy window:
            if qedep > e_max or qedep < e_min:
                continue

            # append qedep to array:
            array_qedep.append(qedep)

            # get the current event in geninfo tree:
            rtree_evt.GetEntry(event)

            # get total number of PE:
            total_pe = int(rtree_evt.GetBranch("totalPE").GetLeaf("totalPE").GetValue())
            # append total_pe to array:
            array_total_pe.append(total_pe)

            # get deposited energy of the event in MeV:
            edep = float(rtree_evt.GetBranch("edep").GetLeaf("edep").GetValue())
            # append edep to array:
            array_edep.append(edep)

            # get x, y, z position, where energy is deposited, in mm:
            edepx = float(rtree_evt.GetBranch("edepX").GetLeaf("edepX").GetValue())
            edepy = float(rtree_evt.GetBranch("edepY").GetLeaf("edepY").GetValue())
            edepz = float(rtree_evt.GetBranch("edepZ").GetLeaf("edepZ").GetValue())

            # calculate distance to center in mm:
            r_edep = np.sqrt(edepx**2 + edepy**2 + edepz**2)

            # append r_edep to array:
            array_r_edep.append(r_edep)

    return array_total_pe, array_r_edep, array_qedep, array_edep


# set energy window, that should be analyzed in MeV:
E_min_window = 0.0
E_max_window = 120.0

# preallocate array, where total PE of all events are stored:
array_totalPE = []
# preallocate array, where distance to center of all events are stored in mm:
array_Redep = []
# preallocate array, where Qedep of all events are stored in MeV:
array_Qedep = []

""" read neutron files: """
print("read neutron events...")
path_neutron = "/local/scratch1/pipc51/astro/blum/conversion_nPE_MeV/neutron_output/"

first_file_neutron = 0
last_file_neutron = 99

# 10 MeV neutrons:
# filename_neutron_10MeV = "user_neutron_10_MeV_"
# nPE_neutron_10MeV, Redep_neutron_10MeV, Qedep_neutron_10MeV, Edep_neutron_10MeV = \
#     get_npe_redep_qedep_from_file(path_neutron, filename_neutron_10MeV, first_file_neutron, last_file_neutron,
#                                   E_min_window, E_max_window)

# 100 MeV neutrons:
# filename_neutron_100MeV = "user_neutron_100_MeV_"
# nPE_neutron_100MeV, Redep_neutron_100MeV, Qedep_neutron_100MeV, Edep_neutron_100MeV = \
#     get_npe_redep_qedep_from_file(path_neutron, filename_neutron_100MeV, first_file_neutron, last_file_neutron,
#                                   E_min_window, E_max_window)

# 300 MeV neutrons:
# filename_neutron_300MeV = "user_neutron_300_MeV_"
# nPE_neutron_300MeV, Redep_neutron_300MeV, Qedep_neutron_300MeV, Edep_neutron_300MeV = \
#     get_npe_redep_qedep_from_file(path_neutron, filename_neutron_300MeV, first_file_neutron, last_file_neutron,
#                                   E_min_window, E_max_window)

# 500 MeV neutron:
last_file_neutron_500MeV = 699
# filename_neutron_500MeV = "user_neutron_500_MeV_"
# nPE_neutron_500MeV, Redep_neutron_500MeV, Qedep_neutron_500MeV, Edep_neutron_500MeV = \
#     get_npe_redep_qedep_from_file(path_neutron, filename_neutron_500MeV, first_file_neutron, last_file_neutron_500MeV,
#                                   E_min_window, E_max_window)

# 1000 MeV neutron:
# filename_neutron_1000MeV = "user_neutron_1000_MeV_"
# nPE_neutron_1000MeV, Redep_neutron_1000MeV, Qedep_neutron_1000MeV, Edep_neutron_1000MeV = \
#     get_npe_redep_qedep_from_file(path_neutron, filename_neutron_1000MeV, first_file_neutron, last_file_neutron,
#                                   E_min_window, E_max_window)

# # total PE of all neutron events:
# totalPE_neutron = nPE_neutron_10MeV + nPE_neutron_100MeV + nPE_neutron_300MeV + nPE_neutron_500MeV + nPE_neutron_1000MeV
# # distance to center of all neutron events:
# Redep_neutron = (Redep_neutron_10MeV + Redep_neutron_100MeV + Redep_neutron_300MeV + Redep_neutron_500MeV +
#                  Redep_neutron_1000MeV)
# # Qedep of all neutron events:
# Qedep_neutron = (Qedep_neutron_10MeV + Qedep_neutron_100MeV + Qedep_neutron_300MeV + Qedep_neutron_500MeV +
#                  Qedep_neutron_1000MeV)
# Edep of all neutron events:
# Edep_neutron = (Edep_neutron_10MeV + Edep_neutron_100MeV + Edep_neutron_300MeV + Edep_neutron_500MeV +
#                 Edep_neutron_1000MeV)

""" read proton files: """
print("read proton events...")
path_proton = "/local/scratch1/pipc51/astro/blum/conversion_nPE_MeV/proton_output/"

first_file_proton = 0
last_file_proton = 99

# 10 MeV protons:
# filename_proton_10MeV = "user_proton_10_MeV_"
# nPE_proton_10MeV, Redep_proton_10MeV, Qedep_proton_10MeV, Edep_proton_10MeV = \
#     get_npe_redep_qedep_from_file(path_proton, filename_proton_10MeV, first_file_proton, last_file_proton,
#                                   E_min_window, E_max_window)

# 100 MeV proton:
# filename_proton_100MeV = "user_proton_100_MeV_"
# nPE_proton_100MeV, Redep_proton_100MeV, Qedep_proton_100MeV, Edep_proton_100MeV = \
#     get_npe_redep_qedep_from_file(path_proton, filename_proton_100MeV, first_file_proton, last_file_proton,
#                                   E_min_window, E_max_window)

# 1000 MeV proton:
# filename_proton_1000MeV = "user_proton_1000_MeV_"
# nPE_proton_1000MeV, Redep_proton_1000MeV, Qedep_proton_1000MeV, Edep_proton_1000MeV = \
#     get_npe_redep_qedep_from_file(path_proton, filename_proton_1000MeV, first_file_proton, last_file_proton,
#                                   E_min_window, E_max_window)

# # total PE of all proton events:
# totalPE_proton = nPE_proton_10MeV + nPE_proton_100MeV + nPE_proton_1000MeV
# # distance to center of all proton events:
# Redep_proton = Redep_proton_10MeV + Redep_proton_100MeV + Redep_proton_1000MeV
# # Qedep of all proton events:
# Qedep_proton = Qedep_proton_10MeV + Qedep_proton_100MeV + Qedep_proton_1000MeV
# # Edep of all proton events:
# Edep_proton = Edep_proton_10MeV + Edep_proton_100MeV + Edep_proton_1000MeV

""" read positron files: """
print("read positron events...")
path_positron = "/local/scratch1/pipc51/astro/blum/positron_output/"

first_file_positron = 0
last_file_positron = 99

# positrons from 10 MeV to 100 MeV:
# filename_positron = "user_positron_"
# totalPE_positron, Redep_positron, Qedep_positron, Edep_positron = \
#     get_npe_redep_qedep_from_file(path_positron, filename_positron, first_file_positron, last_file_positron,
#                                   E_min_window, E_max_window)

# positrons of 10 MeV:
filename_positron_10MeV = "user_positron_10_MeV_"
totalPE_positron_10MeV, Redep_positron_10MeV, Qedep_positron_10MeV, Edep_positron_10MeV = \
    get_npe_redep_qedep_from_file(path_positron, filename_positron_10MeV, first_file_positron, last_file_positron,
                                  E_min_window, E_max_window)

# positrons of 100 MeV:
filename_positron_100MeV = "user_positron_100_MeV_"
totalPE_positron_100MeV, Redep_positron_100MeV, Qedep_positron_100MeV, Edep_positron_100MeV = \
    get_npe_redep_qedep_from_file(path_positron, filename_positron_100MeV, first_file_positron, last_file_positron,
                                  E_min_window, E_max_window)

# positrons of 50 MeV in detector center:
path_positron_CDcenter = "/local/scratch1/pipc51/astro/blum/positron_output_CDcenter/"

filename_positron_50MeV = "user_positron_50MeV_"
totalPE_positron_50MeV, Redep_positron_50MeV, Qedep_positron_50MeV, Edep_positron_50MeV = \
    get_npe_redep_qedep_from_file(path_positron_CDcenter, filename_positron_50MeV, first_file_positron,
                                  last_file_positron, E_min_window, E_max_window)
# convert list Qedep_positron_50MeV to numpy array:
Qedep_positron_50MeV = np.asarray(Qedep_positron_50MeV)
# smear Qedep_positron_50MeV with the energy resolution of JUNO:
sigma = NC_background_functions.energy_resolution(Qedep_positron_50MeV)
Qedep_positron_50MeV_smeared = np.random.normal(Qedep_positron_50MeV, sigma)

""" read IBD files: """
print("read IBD events...")
path_IBD = "/local/scratch1/pipc51/astro/blum/IBD_hepevt/"

first_file_IBD = 0
last_file_IBD = 199

# IBD (positron and neutron) events from 10 MeV to 100 MeV:
# filename_IBD = "user_IBD_hepevt_"
# totalPE_IBD, Redep_IBD, Qedep_IBD, Edep_IBD = \
#     get_npe_redep_qedep_from_file(path_IBD, filename_IBD, first_file_IBD, last_file_IBD, E_min_window, E_max_window)

""" plot Qedep vs. nPE and Qedep vs. Redep: """
# h1 = plt.figure(1, figsize=(15, 8))
# plt.plot(totalPE_neutron, Qedep_neutron, "bx", label="neutron events")
# plt.plot(totalPE_proton, Qedep_proton, "rx", label="proton events")
# plt.plot(totalPE_positron, Qedep_positron, "gx", label="positron events")
# plt.plot(totalPE_IBD, Qedep_IBD, "mx", label="IBD events")
# plt.xlabel("number of p.e.")
# plt.ylabel("visible energy in MeV")
# plt.title("Visible energy vs. number of p.e.")
# plt.grid()
# plt.legend()
#
# h2 = plt.figure(2, figsize=(15, 8))
# plt.plot(Redep_neutron, Qedep_neutron, "bx", label="neutron events")
# plt.plot(Redep_proton, Qedep_proton, "rx", label="proton events")
# plt.plot(Redep_positron, Qedep_positron, "gx", label="positron events")
# plt.plot(Redep_IBD, Qedep_IBD, "mx", label="IBD events")
# plt.xlabel("distance to detector center in mm")
# plt.ylabel("visible energy in MeV")
# plt.title("Visible energy vs. distance to center")
# plt.grid()
# plt.legend()

h3 = plt.figure(3, figsize=(15, 8))
Bins = np.arange(E_min_window, E_max_window, 0.1)
plt.hist(Edep_positron_10MeV, bins=Bins, label="10 MeV positron")
plt.hist(Edep_positron_100MeV, bins=Bins, label="100 MeV positron")
plt.hist(Edep_positron_50MeV, bins=Bins, label="50 MeV positron")
plt.xlabel("Edep in MeV")
plt.ylabel("entries per bin")
plt.title("Deposited energy of 10 MeV, 50 MeV and 100 MeV positrons")
plt.grid()
plt.legend()

h4 = plt.figure(4, figsize=(15, 8))
plt.hist2d(Edep_positron_10MeV, Redep_positron_10MeV, bins=[Bins, np.arange(0.0, 17700+400, 1000)], cmap="Reds")
# plt.hist2d(Edep_positron_100MeV, Redep_positron_100MeV, bins=[Bins, np.arange(0.0, 17700+400, 1000)], cmap="Reds")
plt.xlabel("Edep in MeV")
plt.ylabel("distance to detector center in mm")
plt.title("10 MeV and 100 MeV positrons")
plt.colorbar()
plt.grid()

h5 = plt.figure(5, figsize=(15, 8))
Bins_PE = np.arange(min(totalPE_positron_10MeV), max(totalPE_positron_100MeV), 150)
plt.hist2d(totalPE_positron_10MeV, Redep_positron_10MeV, bins=[Bins_PE, np.arange(0.0, 17700+400, 1000)], cmap="Reds")
# plt.hist2d(totalPE_positron_100MeV, Redep_positron_100MeV, bins=[Bins_PE, np.arange(0.0, 17700+400, 1000)], cmap="Reds")
plt.xlabel("total number of p.e.")
plt.ylabel("distance to detector center in mm")
plt.title("10 MeV and 100 MeV positrons")
plt.colorbar()
plt.grid()

h6 = plt.figure(6, figsize=(15, 8))
plt.hist(totalPE_positron_10MeV, bins=Bins_PE, label="10 MeV positron")
plt.hist(totalPE_positron_100MeV, bins=Bins_PE, label="100 MeV positron")
plt.hist(totalPE_positron_50MeV, bins=Bins_PE, label="50 MeV positron in CD center")
plt.xlabel("total number of p.e.")
plt.ylabel("entries per bin")
plt.title("total number of p.e. of 10 MeV, 50 MeV and 100 MeV positrons")
plt.grid()
plt.legend()

h7 = plt.figure(7, figsize=(15, 8))
plt.hist(Qedep_positron_10MeV, bins=Bins, label="10 MeV positron")
plt.hist(Qedep_positron_100MeV, bins=Bins, label="100 MeV positron")
plt.hist(Qedep_positron_50MeV, bins=Bins, label="50 MeV positron in CD center")
plt.xlabel("Qedep in MeV")
plt.ylabel("entries per bin")
plt.title("Quenched deposited energy of 10 MeV, 50 MeV and 100 MeV positrons")
plt.grid()
plt.legend()

h8 = plt.figure(8, figsize=(15, 8))
plt.hist(Qedep_positron_50MeV_smeared, bins=Bins, label="50 MeV positron in CD center\n"
                                                        "(smeared with energy resolution)")
plt.xlabel("Qedep in MeV")
plt.ylabel("entries per bin")
plt.title("Quenched deposited energy of 50 MeV positrons in the detector center")
plt.grid()
plt.legend()

plt.show()




