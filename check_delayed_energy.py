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
from scipy.optimize import curve_fit

# get the date and time, when the script was run:
date = datetime.datetime.now()
now = date.strftime("%Y-%m-%d %H:%M")

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

# Set boolean variables:
PLOT_POSITION = True
PLOT_INITENERGY = True
PLOT_HITTIME = True

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
for number in range(start_number, stop_number_2_2+1, 1):
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
for number in range(start_number, stop_number_rest+1, 1):
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
for number in range(start_number, stop_number_rest+1, 1):
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

""" define gaussian function: """


def gaussian(x, a, b, c):
    return a * np.exp(- (x-b)**2 / (2*c**2))


h1 = plt.figure(1, figsize=(15, 8))
n_PE_2_2_MeV, bins_2_2_MeV, patches1 = plt.hist(number_pe_2_2, align='mid', bins=50,
                                                label="{0:d} $\\gamma'$s with E = 2.2 MeV"
                                                .format(number_events_total_2_2))
# Fit gaussian distribution to histogram:
# calculate mean of number of PE:
mean_2_2_MeV = np.mean(number_pe_2_2)
# calculate variance:
variance_2_2_MeV = np.var(number_pe_2_2)
# calculate sigma:
sigma_2_2_MeV = np.sqrt(variance_2_2_MeV)
# x-values for Gaussian distribution:
xvalues_2_2 = np.linspace(min(number_pe_2_2), max(number_pe_2_2), 250)
# fit gaussian to histogram (popt contains best-fit parameter [a, b, c]):
popt_2_2, pcov_2_2 = curve_fit(gaussian, bins_2_2_MeV[0:-1], n_PE_2_2_MeV, [100, mean_2_2_MeV, sigma_2_2_MeV])

plt.plot(xvalues_2_2, gaussian(xvalues_2_2, *popt_2_2), color='r',
         label='Gaussian fit: mean={0:0.2f}, $\\sigma$={1:0.2f}'.format(popt_2_2[1], popt_2_2[2]))
plt.xlabel("number of PE (per $\\gamma$)", fontsize=13)
plt.ylabel("entries per bin", fontsize=13)
plt.title("Number of PE of 2.2 MeV $\\gamma$", fontsize=18)
plt.legend()
plt.grid()
plt.savefig(output_path + "hist_nPE_2_2_MeV.png")


h2 = plt.figure(2, figsize=(15, 8))
n_PE_1_9_MeV, bins_1_9_MeV, patches2 = plt.hist(number_pe_1_9, align='mid', bins=50,
                                                label="{0:d} $\\gamma'$s with E = 1.9 MeV"
                                                .format(number_events_total_rest))
# Fit gaussian distribution to histogram:
# calculate mean of number of PE:
mean_1_9_MeV = np.mean(number_pe_1_9)
# calculate variance:
variance_1_9_MeV = np.var(number_pe_1_9)
# calculate sigma:
sigma_1_9_MeV = np.sqrt(variance_1_9_MeV)
# x-values for Gaussian distribution:
xvalues_1_9 = np.linspace(min(number_pe_1_9), max(number_pe_1_9), 250)
# fit gaussian to histogram:
popt_1_9, pcov_1_9 = curve_fit(gaussian, bins_1_9_MeV[0:-1], n_PE_1_9_MeV, [100, mean_1_9_MeV, sigma_1_9_MeV])

plt.plot(xvalues_1_9, gaussian(xvalues_1_9, *popt_1_9), color='r',
         label='Gaussian fit: mean={0:0.2f}, $\\sigma$={1:0.2f}'.format(popt_1_9[1], popt_1_9[2]))
plt.xlabel("number of PE (per $\\gamma$)", fontsize=13)
plt.ylabel("entries per bin", fontsize=13)
plt.title("Number of PE of 1.9 MeV $\\gamma$", fontsize=18)
plt.legend()
plt.grid()
plt.savefig(output_path + "hist_nPE_1_9_MeV.png")

h3 = plt.figure(3, figsize=(15, 8))
n_PE_2_5_MeV, bins_2_5_MeV, patches3 = plt.hist(number_pe_2_5, align='mid', bins=50,
                                                label="{0:d} $\\gamma'$s with E = 2.5 MeV"
                                                .format(number_events_total_rest))
# Fit gaussian distribution to histogram:
# calculate mean of number of PE:
mean_2_5_MeV = np.mean(number_pe_2_5)
# calculate variance:
variance_2_5_MeV = np.var(number_pe_2_5)
# calculate sigma:
sigma_2_5_MeV = np.sqrt(variance_2_5_MeV)
# x-values for Gaussian distribution:
xvalues_2_5 = np.linspace(min(number_pe_2_5), max(number_pe_2_5), 250)
# fit gaussian to histogram:
popt_2_5, pcov_2_5 = curve_fit(gaussian, bins_2_5_MeV[0:-1], n_PE_2_5_MeV, [100, mean_2_5_MeV, sigma_2_5_MeV])

plt.plot(xvalues_2_5, gaussian(xvalues_2_5, *popt_2_5), color='r',
         label='Gaussian fit: mean={0:0.2f}, $\\sigma$={1:0.2f}'.format(popt_2_5[1], popt_2_5[2]))
plt.xlabel("number of PE (per $\\gamma$)", fontsize=13)
plt.ylabel("entries per bin", fontsize=13)
plt.title("Number of PE of 2.5 MeV $\\gamma$", fontsize=18)
plt.legend()
plt.grid()
plt.savefig(output_path + "hist_nPE_2_5_MeV.png")

# plot initial position of gammas:
if PLOT_POSITION:
    h4 = plt.figure(4, figsize=(15, 8))
    range_detector = (-17700, 17700)
    plt.hist(xposition_init_2_2, bins=50, range=range_detector, align='mid', alpha=0.8, color='b', label="x position")
    plt.hist(yposition_init_2_2, bins=50, range=range_detector, align='mid', alpha=0.6, color='r', label="y position")
    plt.hist(zposition_init_2_2, bins=50, range=range_detector, align='mid', alpha=0.4, color='g', label="z position")
    plt.xlim(xmin=-17700, xmax=17700)
    plt.xlabel("initial position in mm", fontsize=13)
    plt.ylabel("entries", fontsize=13)
    plt.title("Initial positions of {0:d} $\\gamma'$s with E = 2.2 MeV".format(number_events_total_2_2), fontsize=18)
    plt.legend()
    plt.grid()
    plt.savefig(output_path + "init_pos_2_2_MeV.png")

    h5 = plt.figure(5, figsize=(15, 8))
    plt.hist(xposition_init_1_9, bins=50, range=range_detector, align='mid', alpha=0.8, color='b', label="x position")
    plt.hist(yposition_init_1_9, bins=50, range=range_detector, align='mid', alpha=0.6, color='r', label="y position")
    plt.hist(zposition_init_1_9, bins=50, range=range_detector, align='mid', alpha=0.4, color='g', label="z position")
    plt.xlim(xmin=-17700, xmax=17700)
    plt.xlabel("initial position in mm", fontsize=13)
    plt.ylabel("entries", fontsize=13)
    plt.title("Initial positions of {0:d} $\\gamma'$s with E = 1.9 MeV".format(number_events_total_rest), fontsize=18)
    plt.legend()
    plt.grid()
    plt.savefig(output_path + "init_pos_1_9_MeV.png")

    h6 = plt.figure(6, figsize=(15, 8))
    plt.hist(xposition_init_2_5, bins=50, range=range_detector, align='mid', alpha=0.8, color='b', label="x position")
    plt.hist(yposition_init_2_5, bins=50, range=range_detector, align='mid', alpha=0.6, color='r', label="y position")
    plt.hist(zposition_init_2_5, bins=50, range=range_detector, align='mid', alpha=0.4, color='g', label="z position")
    plt.xlim(xmin=-17700, xmax=17700)
    plt.xlabel("initial position in mm", fontsize=13)
    plt.ylabel("entries", fontsize=13)
    plt.title("Initial positions of {0:d} $\\gamma'$s with E = 2.5 MeV".format(number_events_total_rest), fontsize=18)
    plt.legend()
    plt.grid()
    plt.savefig(output_path + "init_pos_2_5_MeV.png")

    h9 = plt.figure(9, figsize=(15, 8))
    range2 = (0, 313290000)
    plt.hist(r_squared_2_2, bins=50, range=range2, align='mid', alpha=0.8, color='b',
             label="{0:d} gammas with E = 2.2 MeV".format(number_events_total_2_2))
    plt.hist(r_squared_1_9, bins=50, range=range2, align='mid', alpha=0.6, color='r',
             label="{0:d} gammas with E = 1.9 MeV".format(number_events_total_rest))
    plt.hist(r_squared_2_5, bins=50, range=range2, align='mid', alpha=0.4, color='g',
             label="{0:d} gammas with E = 2.5 MeV".format(number_events_total_rest))
    plt.xlim(xmin=0, xmax=313290000)
    plt.xlabel("initial R$^2$ in mm$^2$", fontsize=13)
    plt.ylabel("entries", fontsize=13)
    plt.title("Initial R$^2$ for 3 cases", fontsize=18)
    plt.legend()
    plt.grid()
    plt.savefig(output_path + "init_R_square.png")

# plot initial energy of gammas:
if PLOT_INITENERGY:
    h7 = plt.figure(7, figsize=(15, 8))
    bin_width = 0.005
    Bins = np.arange(1.75, 2.65, bin_width)
    plt.hist(momentum_init_2_2, bins=Bins, color='r', align='mid',
             label="{0:d} gammas with E = 2.2 MeV".format(number_events_total_2_2))
    plt.hist(momentum_init_1_9, bins=Bins, color='b', align='mid',
             label="{0:d} gammas with E = 1.9 MeV".format(number_events_total_rest))
    plt.hist(momentum_init_2_5, bins=Bins, color='g', align='mid',
             label="{0:d} gammas with E = 2.5 MeV".format(number_events_total_rest))
    plt.xlabel("initial energy in MeV", fontsize=13)
    plt.ylabel("entries per bin (bin-width = {0:0.3f} MeV)".format(bin_width), fontsize=13)
    plt.title("Initial $\\gamma$-energy for 3 cases", fontsize=18)
    plt.legend()
    plt.grid()
    plt.savefig(output_path + "init_energy.png")

# plot example of hittimes:
if PLOT_HITTIME:
    h8 = plt.figure(8, figsize=(15, 8))
    plt.hist(hittimes_2_2, bins=50, align='mid', histtype='step')
    plt.xlabel("hit-time in ns", fontsize=13)
    # INFO-me: ylabel is only equal to number of PE, if nPE == 1 for all photons (1 PE each photon)
    plt.ylabel("number of PE per bin", fontsize=13)
    plt.title("Example of PMT hit-time distribution of one 2.2 MeV $\\gamma$", fontsize=18)
    plt.grid()
    plt.savefig(output_path + "example_hittime_2_2_MeV.png")

""" Efficiency of energy cut of delayed signal (1.9 MeV < E_delayed < 2.5 MeV): """
# 1.9 MeV correspond to number of PE of popt_1_9[1] (mean of gaussian fit of number_pe_1_9):
PE_cut_min = popt_1_9[1]
# 2.5 MeV correspond to number of PE pf popt_2_5[1] (mean of gaussian fit of number_pe_2_5):
PE_cut_max = popt_2_5[1]
# preallocate number of entries in number_pe_2_2, that are outside of window (that are cut away):
entries_cutaway = 0
# preallocate number of entries in number_pe_2_2, that are inside of window (that are not cut away):
entries_notcutaway = 0
# how many entries of number_pe_2_2 are outside of the window specified by PE_cut_min and PE_cut_max:
for index in range(len(number_pe_2_2)):
    if number_pe_2_2[index] < PE_cut_min or number_pe_2_2[index] > PE_cut_max:
        entries_cutaway = entries_cutaway + 1
    else:
        entries_notcutaway = entries_notcutaway + 1

# efficiency of energy cut on delayed signal (in percent):
efficiency = float(entries_notcutaway) / float(number_events_total_2_2) * 100

h10 = plt.figure(10, figsize=(15, 8))
plt.hist(number_pe_2_2, align='mid', bins=50, label="{0:d} $\\gamma'$s with E = 2.2 MeV\ncut efficiency = {1:0.2f}%"
         .format(number_events_total_2_2, efficiency))
plt.axvline(PE_cut_min, ymin=0, ymax=max(n_PE_2_2_MeV), color='r', linestyle=":",
            label="lower cut parameter = {0:0.2f} PE (correspond to 1.9 MeV)".format(PE_cut_min))
plt.axvline(PE_cut_max, ymin=0, ymax=max(n_PE_2_2_MeV), color='r', linestyle="--",
            label="upper cut parameter = {0:0.2f} PE (correspond to 2.5 MeV)".format(PE_cut_max))
plt.xlabel("number of PE (per $\\gamma$)", fontsize=13)
plt.ylabel("entries per bin", fontsize=13)
plt.title("Number of PE of 2.2 MeV $\\gamma$", fontsize=18)
plt.legend()
plt.grid()
plt.savefig(output_path + "hist_nPE_2_2_MeV_efficiency.png")



