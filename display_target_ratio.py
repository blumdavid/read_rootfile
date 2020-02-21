""" script to display the energy spectrum of incoming neutrinos that interact via NC with the JUNO liquid scintillator
    for the different targets (e, p, C12, N14, O16, S32)

    This script is only used to display the target ratio.

    For all other plots (neutrino spectrum of the different neutrino flavour, ...) the output root file AFTER the
    DSNB-NC generator is used (see checkout_NCgen.py).

    Target ratio can not be displayed with checkout_NCgen.py, because there only NC interaction on C12 are considered
    (NO other targets)

"""
import ROOT
import datetime
import numpy as np
from matplotlib import pyplot as plt

# get the date and time, when the script was run:
date = datetime.datetime.now()
now = date.strftime("%Y-%m-%d %H:%M")

# set the path of the input root file (genie_data_NC.root consists NC interactions on all targets from
# gntp.101.gst.root):
input_file = "/home/astro/blum/juno/atmoNC/data_NC/genie_data_NC.root"

# path, where the output should be saved:
output_path = "/home/astro/blum/juno/atmoNC/data_NC/"

# bin-width of the array, which represents the incoming neutrino energy (in GeV) (float):
bin_width_incoming = 0.01

# load the ROOT file:
rfile = ROOT.TFile(input_file)
# get the TTree from the TFile:
rtree = rfile.Get("particleT")

# get the number of entries, i.e. events, in the ROOT-file:
number_entries = rtree.GetEntries()

""" preallocate the arrays: """
# energy of the incoming neutrino in GeV (1d array of float):
projectile_energy = np.array([])
# energy of the incoming neutrino interacting via NC with C12:
energy_nu_c12 = np.array([])
# energy of the incoming neutrino interacting via ES with free protons:
energy_nu_proton = np.array([])
# energy of the incoming neutrino interacting via NC with N14:
energy_nu_n14 = np.array([])
# energy of the incoming neutrino interacting via NC with O16:
energy_nu_o16 = np.array([])
# energy of the incoming neutrino interacting via ES with an electron:
energy_nu_electron = np.array([])
# energy of the incoming neutrino interacting via NC with S32:
energy_nu_s32 = np.array([])

# loop over every entry, i.e. every event, in the TTree:
for event in range(number_entries):

    # get the current event in the TTree:
    rtree.GetEntry(event)

    # get projectile energy and append it to the array:
    p_en = float(rtree.GetBranch('pEn').GetLeaf('pEn').GetValue())
    projectile_energy = np.append(projectile_energy, p_en)

    # get the target's pdgid:
    tgt_pdgid = int(rtree.GetBranch('tPdg').GetLeaf('tPdg').GetValue())

    if tgt_pdgid == 1000060120:
        # C12:
        energy_nu_c12 = np.append(energy_nu_c12, p_en)

    elif tgt_pdgid == 2212:
        # proton:
        energy_nu_proton = np.append(energy_nu_proton, p_en)

    elif tgt_pdgid == 1000070140:
        # N14:
        energy_nu_n14 = np.append(energy_nu_n14, p_en)

    elif tgt_pdgid == 1000080160:
        # O16:
        energy_nu_o16 = np.append(energy_nu_o16, p_en)

    elif tgt_pdgid == 11:
        # electron:
        energy_nu_electron = np.append(energy_nu_electron, p_en)

    elif tgt_pdgid == 1000160320:
        # S32:
        energy_nu_s32 = np.append(energy_nu_s32, p_en)

    else:
        print("other target than proton, electron, C12, N14, O16 or S32!!")

""" get number of events for the two different targets: """
# calculate the number of NC interaction events on C12 (integer):
n_c12 = len(energy_nu_c12)
# calculate the number of elastic scattering events on free protons (integer):
n_proton = len(energy_nu_proton)
# calculate the number of NC interaction events on N14 (integer):
n_n14 = len(energy_nu_n14)
# calculate the number of NC interaction events on O16 (integer):
n_o16 = len(energy_nu_o16)
# calculate the number of elastic scattering events on electrons (integer):
n_electron = len(energy_nu_electron)
# calculate the number of NC interaction events on S32 (integer):
n_s32 = len(energy_nu_s32)

""" get fraction of events for the two different targets (IN PERCENT): """
# calculate the fraction of NC interaction events on C12 in % (float):
fraction_c12 = float(n_c12)/float(number_entries)*100
# calculate the fraction of ES events on free protons in % (float):
fraction_proton = float(n_proton)/float(number_entries)*100
# calculate the fraction of NC interaction events on N14 in % (float):
fraction_n14 = float(n_n14)/float(number_entries)*100
# calculate the fraction of NC interaction events on O16 in % (float):
fraction_o16 = float(n_o16)/float(number_entries)*100
# calculate the fraction of ES events on electrons in % (float):
fraction_electron = float(n_electron)/float(number_entries)*100
# calculate the fraction of NC interaction events on S32 in % (float):
fraction_s32 = float(n_s32)/float(number_entries)*100

""" energy of the incoming neutrinos, defines the bin edges of the histogram """
energy_range = np.arange(0, np.max(projectile_energy)+2*bin_width_incoming, bin_width_incoming)

""" Display Output of 'get_target_ratio()' in plot: """
h1 = plt.figure(1, figsize=(15, 8))

# do not display the histogram, when there are no events:
plt.plot([], [], color='w', label="total number of entries = {0:d}".format(number_entries))
if n_c12 != 0:
    plt.hist(energy_nu_c12, energy_range, histtype='step', color='b',
             label="interactions on $^{12}$C: "+"fraction = {0:.3f}%".format(fraction_c12))
    plt.xscale("log")
    plt.yscale("log")

if n_proton != 0:
    plt.hist(energy_nu_proton, energy_range, histtype='step', color='r',
             label="interaction on protons: fraction = {0:.3f}%".format(fraction_proton))
    plt.xscale("log")
    plt.yscale("log")

plt.plot([], [], color='w', label="sum of fractions of $^{14}$N, $^{16}$O, electrons and $^{32}$S = "+"{0:.3f}%"
         .format(fraction_n14+fraction_o16+fraction_electron+fraction_s32))

plt.xlim(xmin=0.013, xmax=10)
plt.ylim(ymin=10)
plt.xlabel("Neutrino energy $E_{\\nu}$ in GeV", fontsize=15)
plt.ylabel("events", fontsize=15)
plt.title("Neutrino energy spectrum for NC interactions in JUNO liquid scintillator", fontsize=20)
plt.legend(fontsize=12)

plt.savefig(output_path + "target_ratio_NC_{0:.0f}evts.png".format(number_entries))

""" plot histogram in logarithmic scale: """
h2 = plt.figure(2, figsize=(15, 8))

# do not display the histogram, when there are no events:
plt.plot([], [], color='w', label="total number of entries = {0:d}".format(number_entries))
if n_c12 != 0:
    plt.hist(energy_nu_c12, energy_range, histtype='step', color='b', log=True,
             label="interactions on $^{12}$C: "+"fraction = {0:.3f}%".format(fraction_c12))

if n_proton != 0:
    plt.hist(energy_nu_proton, energy_range, histtype='step', color='r', log=True,
             label="interaction on protons: fraction = {0:.3f}%".format(fraction_proton))

if n_n14 != 0:
    plt.hist(energy_nu_n14, energy_range, histtype='step', color='g', log=True,
             label="interaction on $^{14}$N: "+"fraction = {0:.3f}%".format(fraction_n14))

if n_o16 != 0:
    plt.hist(energy_nu_o16, energy_range, histtype='step', color='k', log=True,
             label="interaction on $^{16}$O: "+"fraction = {0:.3f}%".format(fraction_o16))

if n_electron != 0:
    plt.hist(energy_nu_electron, energy_range, histtype='step', color='c', log=True,
             label="interaction on electrons: fraction = {0:.3f}%".format(fraction_electron))

if n_s32 != 0:
    plt.hist(energy_nu_s32, energy_range, histtype='step', color='m', log=True,
             label="interaction on $^{32}$S: "+"fraction = {0:.4f}%".format(fraction_s32))

plt.xlim(xmin=0, xmax=10)
plt.ylim(ymin=0.5)
plt.xlabel("Neutrino energy $E_{\\nu}$ in GeV", fontsize=15)
plt.ylabel("events", fontsize=15)
plt.title("Neutrino energy spectrum for NC interactions in JUNO liquid scintillator", fontsize=20)
plt.legend(fontsize=12)

plt.savefig(output_path + "target_ratio_NC_{0:.0f}evts_log.png".format(number_entries))

""" Save information about the target ratio into txt file: """
np.savetxt(output_path + "target_ratio_NC_{0:.0f}evts.txt".format(number_entries),
           np.array([number_entries, n_c12, n_proton, n_n14, n_o16, n_electron, n_s32, fraction_c12, fraction_proton,
                     fraction_n14, fraction_o16, fraction_electron, fraction_s32]), fmt='%4.5f',
           header="Information about the number and ratio of the target particles, where neutrinos interact via "
                  "NC in the JUNO liquid scintillator\n"
                  "(input file: {0}, script: display_target_ratio.py ({1})):\n"
                  "Number of events in the input file,\n"
                  "Number of interactions on C12,\n"
                  "Number of interactions on free protons,\n"
                  "Number of interactions on N14,\n"
                  "Number of interactions on O16,\n"
                  "Number of interactions on electrons,\n"
                  "Number of interactions on S32,\n"
                  "Fraction of interactions on C12,\n"
                  "Fraction of interactions on free protons,\n"
                  "Fraction of interactions on N14,\n"
                  "Fraction of interactions on O16,\n"
                  "Fraction of interactions electrons,\n"
                  "Fraction of interactions on S32:"
           .format(input_file, now))











