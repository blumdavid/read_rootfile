""" Script to check gammas (with energies of 1.9 MeV, 2.2 MeV and 2.5 MeV), which were simulated with tut_detsim.py of
    JUNO offline version J18v2r1-branch.

    Results of this script are used to convert the energy of a 2.2 MeV gamma, that is emitted by neutron capture on
    Hydrogen, to number of PE in the JUNO detector.

    With this conversion the cut on the energy of a possible delayed signal (from neutron capture in IBD) can be made
    in the PE-regime and efficiency of this cut can be calculated.

    More information: info_atmoNC.odt (/home/astro/blum/juno/atmoNC/)
"""
import datetime
import NC_background_functions
import numpy as np
from matplotlib import pyplot as plt

# get the date and time, when the script was run:
date = datetime.datetime.now()
now = date.strftime("%Y-%m-%d %H:%M")

# set the path of the input files:
input_path = "/local/scratch1/pipc51/astro/blum/gamma_2_2_MeV/"

# set path, where results should be saved:
output_path = "/home/astro/blum/juno/atmoNC/data_NC/output_gamma_2_2_MeV/"

# set the number of the first file and number of the last file that should be read:
start_number = 0
stop_number = 2
# number of entries in the input files:
Number_entries_input = 10
# total number of events:
number_events_total = (stop_number - start_number + 1) * Number_entries_input

""" preallocate arrays for 2.2 MeV gamma: """
# number of PE of each event:
number_pe_2_2 = np.array([])
# initial x position of each event in mm:
xposition_init_2_2 = np.array([])
# initial y position of each event in mm:
yposition_init_2_2 = np.array([])
# initial z position of each event in mm:
zposition_init_2_2 = np.array([])
# initial total momentum of each event in MeV:
momentum_init_2_2 = np.array([])

print("\nstart reading 2.2 MeV gamma files...")

# loop over files of gamma = 2.2 MeV:
for number in range(start_number, stop_number+1, 1):
    # path to file:
    input_file = input_path + "user_gamma_2_2_MeV_{0:d}.root".format(number)
    # print(input_file)

    # get number of PE per event (array of int), hit-times of the last event in ns (array of float),
    # initial x, y, z position per event in mm (arrays of float) and initial momentum per event in MeV (array of float):
    num_pe_2_2, hittimes_2_2, xpos_2_2, ypos_2_2, zpos_2_2, momentum_2_2 = \
        NC_background_functions.read_gamma_delayed_signal(input_file, Number_entries_input)

    # append arrays to array:
    number_pe_2_2 = np.append(number_pe_2_2, num_pe_2_2)
    xposition_init_2_2 = np.append(xposition_init_2_2, xpos_2_2)
    yposition_init_2_2 = np.append(yposition_init_2_2, ypos_2_2)
    zposition_init_2_2 = np.append(zposition_init_2_2, zpos_2_2)
    momentum_init_2_2 = np.append(momentum_init_2_2, momentum_2_2)


""" preallocate arrays for 1.9 MeV gamma: """
# number of PE of each event:
number_pe_1_9 = np.array([])
# initial x position of each event in mm:
xposition_init_1_9 = np.array([])
# initial y position of each event in mm:
yposition_init_1_9 = np.array([])
# initial z position of each event in mm:
zposition_init_1_9 = np.array([])
# initial total momentum of each event in MeV:
momentum_init_1_9 = np.array([])

print("\nstart reading 1.9 MeV gamma files...")

# loop over files of gamma = 1.9 MeV:
for number in range(start_number, stop_number+1, 1):
    # path to file:
    input_file = input_path + "user_gamma_1_9_MeV_{0:d}.root".format(number)
    # print(input_file)

    # get number of PE per event (array of int), hit-times of the last event in ns (array of float),
    # initial x, y, z position per event in mm (arrays of float) and initial momentum per event in MeV (array of float):
    num_pe_1_9, hittimes_1_9, xpos_1_9, ypos_1_9, zpos_1_9, momentum_1_9 = \
        NC_background_functions.read_gamma_delayed_signal(input_file, Number_entries_input)

    # append arrays to array:
    number_pe_1_9 = np.append(number_pe_1_9, num_pe_1_9)
    xposition_init_1_9 = np.append(xposition_init_1_9, xpos_1_9)
    yposition_init_1_9 = np.append(yposition_init_1_9, ypos_1_9)
    zposition_init_1_9 = np.append(zposition_init_1_9, zpos_1_9)
    momentum_init_1_9 = np.append(momentum_init_1_9, momentum_1_9)


""" preallocate arrays for 2.5 MeV gamma: """
# number of PE of each event:
number_pe_2_5 = np.array([])
# initial x position of each event in mm:
xposition_init_2_5 = np.array([])
# initial y position of each event in mm:
yposition_init_2_5 = np.array([])
# initial z position of each event in mm:
zposition_init_2_5 = np.array([])
# initial total momentum of each event in MeV:
momentum_init_2_5 = np.array([])

print("\nstart reading 2.5 MeV gamma files...")

# loop over files of gamma = 2.5 MeV:
for number in range(start_number, stop_number+1, 1):
    # path to file:
    input_file = input_path + "user_gamma_2_5_MeV_{0:d}.root".format(number)
    # print(input_file)

    # get number of PE per event (array of int), hit-times of the last event in ns (array of float),
    # initial x, y, z position per event in mm (arrays of float) and initial momentum per event in MeV (array of float):
    num_pe_2_5, hittimes_2_5, xpos_2_5, ypos_2_5, zpos_2_5, momentum_2_5 = \
        NC_background_functions.read_gamma_delayed_signal(input_file, Number_entries_input)

    # append arrays to array:
    number_pe_2_5 = np.append(number_pe_2_5, num_pe_2_5)
    xposition_init_2_5 = np.append(xposition_init_2_5, xpos_2_5)
    yposition_init_2_5 = np.append(yposition_init_2_5, ypos_2_5)
    zposition_init_2_5 = np.append(zposition_init_2_5, zpos_2_5)
    momentum_init_2_5 = np.append(momentum_init_2_5, momentum_2_5)

""" calculate R**2 with initial x, y, z-positions in mm**2: """
# r**2 of 2.2 MeV gamma:
r_squared_2_2 = xposition_init_2_2**2 + yposition_init_2_2**2 + zposition_init_2_2**2
# r**2 of 1_9 MeV gamma:
r_squared_1_9 = xposition_init_1_9**2 + yposition_init_1_9**2 + zposition_init_1_9**2
# r**2 of 2.5 MeV gamma:
r_squared_2_5 = xposition_init_2_5**2 + yposition_init_2_5**2 + zposition_init_2_5**2

h1 = plt.figure(1, figsize=(15, 8))
n_PE_2_2_MeV, bins_2_2_MeV, patches1 = plt.hist(number_pe_2_2, align='mid', bins=50,
                                                label="{0:d} $\\gamma'$s with E = 2.2 MeV".format(number_events_total))
plt.xlabel("number of PE", fontsize=13)
plt.ylabel("entries per bin", fontsize=13)
plt.title("Number of PE of 2.2 MeV $\\gamma$", fontsize=18)
plt.legend()
plt.grid()
plt.savefig(output_path + "hist_nPE_2_2_MeV.png")

h2 = plt.figure(2, figsize=(15, 8))
n_PE_1_9_MeV, bins_1_9_MeV, patches2 = plt.hist(number_pe_1_9, align='mid', bins=50,
                                                label="{0:d} $\\gamma'$s with E = 1.9 MeV".format(number_events_total))
plt.xlabel("number of PE", fontsize=13)
plt.ylabel("entries per bin", fontsize=13)
plt.title("Number of PE of 1.9 MeV $\\gamma$", fontsize=18)
plt.legend()
plt.grid()
plt.savefig(output_path + "hist_nPE_1_9_MeV.png")

h3 = plt.figure(3, figsize=(15, 8))
n_PE_2_5_MeV, bins_2_5_MeV, patches3 = plt.hist(number_pe_2_5, align='mid', bins=50,
                                                label="{0:d} $\\gamma'$s with E = 2.5 MeV".format(number_events_total))
plt.xlabel("number of PE", fontsize=13)
plt.ylabel("entries per bin", fontsize=13)
plt.title("Number of PE of 2.5 MeV $\\gamma$", fontsize=18)
plt.legend()
plt.grid()
plt.savefig(output_path + "hist_nPE_2_5_MeV.png")

h4 = plt.figure(4, figsize=(15, 8))
range_detector = (-17700, 17700)
plt.hist(xposition_init_2_2, bins=50, range=range_detector, align='mid', histtype='step', color='b', label="x position")
plt.hist(yposition_init_2_2, bins=50, range=range_detector, align='mid', histtype='step', color='r', label="y position")
plt.hist(zposition_init_2_2, bins=50, range=range_detector, align='mid', histtype='step', color='g', label="z position")
plt.xlim(xmin=-17700, xmax=17700)
plt.xlabel("initial position in mm", fontsize=13)
plt.ylabel("entries", fontsize=13)
plt.title("Initial positions of {0:d} $\\gamma'$s with E = 2.2 MeV".format(number_events_total), fontsize=18)
plt.legend()
plt.grid()
plt.savefig(output_path + "init_pos_2_2_MeV.png")

h5 = plt.figure(5, figsize=(15, 8))
plt.hist(xposition_init_1_9, bins=50, range=range_detector, align='mid', histtype='step', color='b', label="x position")
plt.hist(yposition_init_1_9, bins=50, range=range_detector, align='mid', histtype='step', color='r', label="y position")
plt.hist(zposition_init_1_9, bins=50, range=range_detector, align='mid', histtype='step', color='g', label="z position")
plt.xlim(xmin=-17700, xmax=17700)
plt.xlabel("initial position in mm", fontsize=13)
plt.ylabel("entries", fontsize=13)
plt.title("Initial positions of {0:d} $\\gamma'$s with E = 1.9 MeV".format(number_events_total), fontsize=18)
plt.legend()
plt.grid()
plt.savefig(output_path + "init_pos_1_9_MeV.png")

h6 = plt.figure(6, figsize=(15, 8))
plt.hist(xposition_init_2_5, bins=50, range=range_detector, align='mid', histtype='step', color='b', label="x position")
plt.hist(yposition_init_2_5, bins=50, range=range_detector, align='mid', histtype='step', color='r', label="y position")
plt.hist(zposition_init_2_5, bins=50, range=range_detector, align='mid', histtype='step', color='g', label="z position")
plt.xlim(xmin=-17700, xmax=17700)
plt.xlabel("initial position in mm", fontsize=13)
plt.ylabel("entries", fontsize=13)
plt.title("Initial positions of {0:d} $\\gamma'$s with E = 2.5 MeV".format(number_events_total), fontsize=18)
plt.legend()
plt.grid()
plt.savefig(output_path + "init_pos_2_5_MeV.png")

h9 = plt.figure(9, figsize=(15, 8))
range2 = (0, 313290000)
plt.hist(r_squared_2_2, bins=50, range=range2, align='mid', histtype='step', color='b',
         label="{0:d} gammas with E = 2.2 MeV".format(number_events_total))
plt.hist(r_squared_1_9, bins=50, range=range2, align='mid', histtype='step', color='r',
         label="{0:d} gammas with E = 2.2 MeV".format(number_events_total))
plt.hist(r_squared_2_5, bins=50, range=range2, align='mid', histtype='step', color='g',
         label="{0:d} gammas with E = 2.2 MeV".format(number_events_total))
plt.xlim(xmin=0, xmax=313290000)
plt.xlabel("initial R$^2$ in mm$^2$", fontsize=13)
plt.ylabel("entries", fontsize=13)
plt.title("Initial R$^2$ for 3 cases".format(number_events_total), fontsize=18)
plt.legend()
plt.grid()
plt.savefig(output_path + "init_R_square.png")

h7 = plt.figure(7, figsize=(15, 8))
bin_width = 0.005
Bins = np.arange(1.75, 2.65, bin_width)
plt.hist(momentum_init_2_2, bins=Bins, color='r', align='mid',
         label="{0:d} gammas with E = 2.2 MeV".format(number_events_total))
plt.hist(momentum_init_1_9, bins=Bins, color='b', align='mid',
         label="{0:d} gammas with E = 1.9 MeV".format(number_events_total))
plt.hist(momentum_init_2_5, bins=Bins, color='g', align='mid',
         label="{0:d} gammas with E = 2.5 MeV".format(number_events_total))
plt.xlabel("initial energy in MeV", fontsize=13)
plt.ylabel("entries per bin (bin-width = {0:0.3f} MeV)".format(bin_width), fontsize=13)
plt.title("Initial $\\gamma$-energy for 3 cases", fontsize=18)
plt.legend()
plt.grid()
plt.savefig(output_path + "init_energy.png")

h8 = plt.figure(8, figsize=(15, 8))
plt.hist(hittimes_2_2, bins=50, align='mid', histtype='step')
plt.xlabel("hit-time in ns", fontsize=13)
# INFO-me: ylabel is only equal to number of PE, if nPE == 1 for all photons (1 PE each photon)
plt.ylabel("number of PE per bin", fontsize=13)
plt.title("Example of PMT hit-time distribution of one 2.2 MeV $\\gamma$", fontsize=18)
plt.grid()
plt.savefig(output_path + "example_hittime_2_2_MeV.png")

# plt.show()
