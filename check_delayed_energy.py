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

# set READ_DATA flag:
READ_DATA = True

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

""" define time window and bin width: """
# set time window of whole signal in ns:
min_time = -50
max_time = 10000
# Set bin-width of hittime histogram in ns:
binwidth = 5.0

""" thresholds for delayed signal: """
# Set threshold of number of PE per bin for possible delayed signal (bin-width = 5 ns):
threshold1 = 50
# set threshold2 of number of PEs per bin (signal peak is summed as long as nPE is above threshold2):
threshold2 = 0

# set the radius for the volume cut in mm:
r_cut = 16000

""" load position of the PMTs and corresponding PMT ID from file PMT_position.root: """
file_PMT_position = "/home/astro/blum/juno/atmoNC/PMT_information/PMT_position.root"
# array with PMT ID and corresponding x, y, z position in mm:
pmtID_pos_file, x_pos_pmt, y_pos_pmt, z_pos_pmt = NC_background_functions.get_pmt_position(file_PMT_position)

""" load 'time resolution' in ns of the 20 inch PMTs and corresponding PMT ID from file PmtData.root: """
file_PMT_time = "/home/astro/blum/juno/atmoNC/PMT_information/PmtData.root"
# array with PMT ID and corresponding sigma in ns:
pmtID_time_file, sigma_time_20inch = NC_background_functions.get_20inchpmt_tts(file_PMT_time)
# set TTS (FWHM) of the 3inch PMTs in ns:
tts_3inch = 5.0
# calculate time resolution (sigma) for the 3inch PMTs in ns:
sigma_time_3inch = tts_3inch / (2 * np.sqrt(2 * np.log(2)))
# set effective speed of light in the liquid scintillator in mm/ns (see page 7 of c_effective_JUNO-doc-3144-v2.pdf in
# folder /home/astro/blum/PhD/paper/Pulse_Shape_Discrimination/). Effective refraction index in LS n_eff = 1.54.
# c/n_eff = 299792458 m / 1.54 s ~ 194670427 m/s = 194670427 * 10**(-6) mm/ns ~ 194.67 mm/ns:
c_effective = 194.67


def save_array_to_file(arr, out_path, file_name, number_events):
    """
    function to save an array (either number_pe or qedep) to txt file to save time, because you must read the file only
    once
    :param arr: array that should be saved (array of float)
    :param out_path: path, where the txt file should be saved (string)
    :param file_name: file name of the txt file (string)
    :param number_events: number of events in the array/root-file
    :return:
    """
    np.savetxt(out_path + file_name + ".txt", arr, fmt='%1.5f',
               header="{0} of {1:d} events analyzed with script check_delayed_energy.py\n(number of photo-electron per "
                      "event OR quenched deposited energy/ visible energy per event in MeV).\n"
                      "({2})\n"
                      "(volume cut (R <= {3:d} mm) applied on initial position):"
               .format(file_name, number_events, now, r_cut))

    return


""" preallocate arrays for 2.2 MeV gamma: """
# number of PE of each event:
number_pe_2_2 = np.array([])

if READ_DATA:
    print("\nstart reading 2.2 MeV gamma files...")

    # loop over files of gamma = 2.2 MeV:
    for number in range(start_number, stop_number_2_2+1, 1):
        # path to file:
        input_file = input_path + "user_gamma_2_2_MeV_{0:d}.root".format(number)

        # get number of PE per event (array of int) from corrected hittime distribution:
        num_pe_2_2, num_analyzed_2_2 = \
            NC_background_functions.read_gamma_delayed_signal(input_file, Number_entries_input, r_cut, x_pos_pmt,
                                                              y_pos_pmt, z_pos_pmt, sigma_time_20inch, sigma_time_3inch,
                                                              c_effective, min_time, max_time, binwidth,
                                                              threshold1, threshold2)

        # append arrays to array:
        number_pe_2_2 = np.append(number_pe_2_2, num_pe_2_2)

    # save number of pe to txt file:
    save_array_to_file(number_pe_2_2, output_path, "number_pe_gamma_2_2", len(number_pe_2_2))
else:
    # load number of pe and qedep array from txt file:
    number_pe_2_2 = np.loadtxt(output_path + "number_pe_gamma_2_2.txt")

""" preallocate arrays for 1.9 MeV gamma: """
# number of PE of each event:
number_pe_1_9 = np.array([])

if READ_DATA:
    print("\nstart reading 1.9 MeV gamma files...")

    # loop over files of gamma = 1.9 MeV:
    for number in range(start_number, stop_number_2_2+1, 1):
        # path to file:
        input_file = input_path + "user_gamma_1_9_MeV_{0:d}.root".format(number)

        # get number of PE per event (array of int) from corrected hittime distribution:
        num_pe_1_9, num_analyzed_1_9 = \
            NC_background_functions.read_gamma_delayed_signal(input_file, Number_entries_input, r_cut, x_pos_pmt,
                                                              y_pos_pmt, z_pos_pmt, sigma_time_20inch, sigma_time_3inch,
                                                              c_effective, min_time, max_time, binwidth,
                                                              threshold1, threshold2)

        # append arrays to array:
        number_pe_1_9 = np.append(number_pe_1_9, num_pe_1_9)

    # save number of pe to txt file:
    save_array_to_file(number_pe_1_9, output_path, "number_pe_gamma_1_9", len(number_pe_1_9))
else:
    # load number of pe and qedep array from txt file:
    number_pe_1_9 = np.loadtxt(output_path + "number_pe_gamma_1_9.txt")

""" preallocate arrays for 2.5 MeV gamma: """
# number of PE of each event:
number_pe_2_5 = np.array([])

if READ_DATA:
    print("\nstart reading 2.5 MeV gamma files...")

    # loop over files of gamma = 2.5 MeV:
    for number in range(start_number, stop_number_rest+1, 1):
        # path to file:
        input_file = input_path + "user_gamma_2_5_MeV_{0:d}.root".format(number)

        # get number of PE per event (array of int) from corrected hittime distribution:
        num_pe_2_5, num_analyzed_2_5 = \
            NC_background_functions.read_gamma_delayed_signal(input_file, Number_entries_input, r_cut, x_pos_pmt,
                                                              y_pos_pmt, z_pos_pmt, sigma_time_20inch, sigma_time_3inch,
                                                              c_effective, min_time, max_time, binwidth,
                                                              threshold1, threshold2)

        # append arrays to array:
        number_pe_2_5 = np.append(number_pe_2_5, num_pe_2_5)

    # save number of pe to txt file:
    save_array_to_file(number_pe_2_5, output_path, "number_pe_gamma_2_5", len(number_pe_2_5))
else:
    # load number of pe and qedep array from txt file:
    number_pe_2_5 = np.loadtxt(output_path + "number_pe_gamma_2_5.txt")


""" define gaussian function: """


def gaussian(x, a, b, c):
    return a * np.exp(- (x-b)**2 / (2*c**2))


h1 = plt.figure(1, figsize=(15, 8))
n_PE_2_2_MeV, bins_2_2_MeV, patches1 = plt.hist(number_pe_2_2, align='mid', bins=100,
                                                label="{0:d} $\\gamma'$s with E = 2.2 MeV"
                                                .format(len(number_pe_2_2)))
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
n_PE_1_9_MeV, bins_1_9_MeV, patches2 = plt.hist(number_pe_1_9, align='mid', bins=100,
                                                label="{0:d} $\\gamma'$s with E = 1.9 MeV"
                                                .format(len(number_pe_1_9)))
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
n_PE_2_5_MeV, bins_2_5_MeV, patches3 = plt.hist(number_pe_2_5, align='mid', bins=100,
                                                label="{0:d} $\\gamma'$s with E = 2.5 MeV"
                                                .format(len(number_pe_2_5)))
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

""" Efficiency of energy cut of delayed signal (1.9 MeV < E_visible_delayed < 2.5 MeV): """
# 1.9 MeV correspond to number of PE of popt_n_10MeV[1] (mean of gaussian fit of number_pe_1_9):
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
efficiency = float(entries_notcutaway) / float(len(number_pe_2_2)) * 100

h11 = plt.figure(11, figsize=(15, 8))
plt.hist(number_pe_2_2, align='mid', bins=100, label="{0:d} $\\gamma'$s with E = 2.2 MeV\ncut efficiency = {1:0.2f}%"
         .format(len(number_pe_2_2), efficiency))
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



