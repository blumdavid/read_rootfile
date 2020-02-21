""" script to analyze the prompt signal of neutrons simulated with JUNO detsim:

    momentum of neutron = 15 MeV

    Do neutrons contribute to prompt signal of IBD event??

"""
import ROOT
import numpy as np
from matplotlib import pyplot as plt
from NC_background_functions import conversion_npe_to_evis

# path, where detsim output is saved:
input_path = "/home/astro/blum/PhD/work/MeVDM_JUNO/gen_spectrum_v2/neutron_signal_IBD_detsim/"

# path, where output should be saved:
output_path = input_path

# radius cut in mm:
radius_cut = 16000

# maximum hit-time in ns (for simplicity use hit-time and not the correct photon emission time. This is no problem,
# because neutrons are produced in detector center.
# Time-of-flight from detector center to PMT: 17700 mm / 194.67 mm / ns = 91 ns. Hit-time is corrected with 91 ns):
max_hittime_1 = 500
max_hittime_2 = 1000
max_hittime_3 = 1500

# array, where the number of pe of the prompt signal of neutrons is stored for the different maximum hit-times:
array_npe_1 = np.array([])
array_npe_2 = np.array([])
array_npe_3 = np.array([])

for filenumber in range(4):

    # input file:
    input_file = input_path + "user_neutron_15_MeV_{0:d}.root".format(filenumber)

    # load the ROOT file:
    rfile = ROOT.TFile(input_file)
    # get the "evt"-TTree from the TFile:
    rtree_evt = rfile.Get("evt")
    # get 'geninfo' tree from TFile:
    rtree_geninfo = rfile.Get("geninfo")
    # get 'prmtrkdep' tree from TFile:
    rtree_prmtrkdep = rfile.Get("prmtrkdep")

    # get the number of entries in the ROOT-file:
    number_entries = rtree_evt.GetEntries()

    for event in range(number_entries):

        # get the current event in geninfo tree:
        rtree_geninfo.GetEntry(event)

        # get InitX, InitY and InitZ from geninfo tree in mm:
        InitX = float(rtree_geninfo.GetBranch('InitX').GetLeaf('InitX').GetValue(0))
        InitY = float(rtree_geninfo.GetBranch('InitY').GetLeaf('InitY').GetValue(0))
        InitZ = float(rtree_geninfo.GetBranch('InitZ').GetLeaf('InitZ').GetValue(0))

        # calculate distance to detector center in mm
        r_init = np.sqrt(InitX**2 + InitY**2 + InitZ**2)

        # do fiducial volume cut:
        if r_init > radius_cut:
            continue

        # get event in evt tree:
        rtree_evt.GetEntry(event)

        # get total number of pe of the event:
        total_pe = int(rtree_evt.GetBranch("totalPE").GetLeaf("totalPE").GetValue())

        # preallocate number of pe of prompt signal:
        number_pe_prompt_1 = 0
        number_pe_prompt_2 = 0
        number_pe_prompt_3 = 0

        # loop over total number of pe:
        for index in range(total_pe):

            # get hittime of photon in ns:
            hittime = float(rtree_evt.GetBranch("hitTime").GetLeaf("hitTime").GetValue(index))

            # correct hittime to get photon emission time:
            photon_emission_time = hittime - 91.0

            # check if pe corresponds to prompt signal:
            if photon_emission_time < max_hittime_1:

                number_pe_prompt_1 += 1

            if photon_emission_time < max_hittime_2:

                number_pe_prompt_2 += 1

            if photon_emission_time < max_hittime_3:

                number_pe_prompt_3 += 1

        # check if there are pe in prompt signal and append value to array:
        if number_pe_prompt_1 != 0:
            array_npe_1 = np.append(array_npe_1, number_pe_prompt_1)

        if number_pe_prompt_2 != 0:
            array_npe_2 = np.append(array_npe_2, number_pe_prompt_2)

        if number_pe_prompt_3 != 0:
            array_npe_3 = np.append(array_npe_3, number_pe_prompt_3)

# convert number of pe into E_vis in MeV:
array_Evis_1 = conversion_npe_to_evis(array_npe_1)
array_Evis_2 = conversion_npe_to_evis(array_npe_2)
array_Evis_3 = conversion_npe_to_evis(array_npe_3)

# how many neutron deposit more than E_cut MeV:
E_cut = 0.04
number_1 = 0
number_2 = 0
number_3 = 0

for index in range(len(array_Evis_1)):
    if array_Evis_1[index] >= E_cut:
        number_1 += 1

for index in range(len(array_Evis_2)):
    if array_Evis_2[index] >= E_cut:
        number_2 += 1

for index in range(len(array_Evis_3)):
    if array_Evis_3[index] >= E_cut:
        number_3 += 1

# maximum number of pe from arrays:
print("max. E_vis for {0} ns time window = {1}, number of neutron with more than {3} MeV = {2}"
      .format(max_hittime_1, max(array_Evis_1), number_1, E_cut))
print("max. E_vis for {0} ns time window = {1}, number of neutron with more than {3} MeV = {2}"
      .format(max_hittime_2, max(array_Evis_2), number_2, E_cut))
print("max. E_vis for {0} ns time window = {1}, number of neutron with more than {3} MeV = {2}"
      .format(max_hittime_3, max(array_Evis_3), number_3, E_cut))

# display array in histogram:
bins = np.arange(0, max(array_Evis_3), 0.01)

h1 = plt.figure(1)
plt.hist(array_Evis_1, bins, color="red", histtype="step", label="entries = {1},\nmaximum hittime = {0} ns"
         .format(max_hittime_1, len(array_npe_1)))
plt.xlabel("visible energy in MeV of prompt signal")
plt.ylabel("entries")
plt.title("Number of p.e. of prompt signal (0 - {0} ns) of ".format(max_hittime_1) +"neutron with $E_{kin}$ = 15 MeV")
plt.legend()
plt.grid()

h2 = plt.figure(2)
plt.hist(array_Evis_2, bins, color="red", histtype="step", label="entries = {1},\nmaximum hittime = {0} ns"
         .format(max_hittime_2, len(array_npe_2)))
plt.xlabel("visible energy in MeV of prompt signal")
plt.ylabel("entries")
plt.title("Number of p.e. of prompt signal (0 - {0} ns) of ".format(max_hittime_2) +"neutron with $E_{kin}$ = 15 MeV")
plt.legend()
plt.grid()

h3 = plt.figure(3)
plt.hist(array_Evis_3, bins, color="red", histtype="step", label="entries = {1},\nmaximum hittime = {0} ns"
         .format(max_hittime_3, len(array_npe_3)))
plt.xlabel("visible energy in MeV of prompt signal")
plt.ylabel("entries")
plt.title("Number of p.e. of prompt signal (0 - {0} ns) of".format(max_hittime_2) +" neutron with $E_{kin}$ = 15 MeV")
plt.legend()
plt.grid()

plt.show()









