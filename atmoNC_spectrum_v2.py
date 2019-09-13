""" Script to read the energy of the IBD-like events of atmospheric NC neutrino background for different cut parameters
    and calculate the spectrum of atmospheric NC neutrino background as function of the visible energy.

    1.  Read txt files where the filenumber and evtID of NC events, which pass the specific cut on the real
        (reconstructed) data, are saved. One file each cut: volume cut, prompt energy cut, delayed cut and PSD.

    2.  Read txt files where the filenumber and evtID of NC events, which pass the specific cut on the ideal
        (MC truth) data, are saved. One file each cut: volume cut, prompt energy cut, delayed cut and PSD.

    3.  Read also txt files where the filenumber and evtID of real IBD event, which pass the same cuts like the NC
        events (real/reconstructed data), are saved. One file each cut: volume cut, prompt energy cut, delayed cut
        and PSD.

    4.  Read also txt files where the filenumber and evtID of real IBD event, which pass the same cuts like the NC
        events (ideal/MC truth data), are saved. One file each cut: volume cut, prompt energy cut, delayed cut
        and PSD.

    3.  Compare the filenumber and evtID for NC events or IBD events, respectively, to get the events that pass all
        cuts and are therefore IBD-like events (do this for real and ideal data to be able to calculate the error of
        the efficiency of the total cuts)

    4.  For IBD events and NC events: calculate the efficiencies of the different cuts (eff_real, eff_ideal, ...) from
        files numbers_... with statistical errors and store them in txt files.

    5.  For NC events:

        5.1 take Evis of prompt signal of all IBD-like NC events and build energy spectrum (without PSD and with PSD)

        5.2 Consider the event rate of NC interactions on C12 inside the detector (calculated with cross-sections and
        neutrino fluxes) and calculate the 'real' spectrum of atmospheric NC background, JUNO will measure after 10
        years of data taking

        5.3 Consider the PSD efficiency and calculate the spectrum of atmospheric NC background, JUNO will measure
        after 10 years of data taking, after Pulse Shape Analysis

"""
import datetime
import sys
import NC_background_functions
import numpy as np
from matplotlib import pyplot as plt


def read_numbers_txt_files(filename_numbers):
    """
    function to read the file specified by filename_numbers and to calculate the efficiencies of this cut.

    :param filename_numbers: file name numbers_....txt of the cut to be analyzed
    :return:
    """
    # read txt file:
    array_numbers = np.loadtxt(filename_numbers)
    n_total_evts = float(array_numbers[0])
    n_analyzed = float(array_numbers[2])
    n_pass_real = float(array_numbers[3])
    n_rej_real = float(array_numbers[4])
    n_pass_ideal = float(array_numbers[5])
    n_rej_ideal = float(array_numbers[6])
    n_too_much = float(array_numbers[7])
    n_too_less = float(array_numbers[8])

    # check the numbers:
    if n_pass_real + n_rej_real != n_analyzed:
        sys.exit("n_pass_real + n_rej_real != n_analyzed")

    if n_pass_ideal + n_rej_ideal != n_analyzed:
        sys.exit("n_pass_ideal + n_rej_ideal != n_analyzed")

    # calculate efficiency for real case (with numbers from reconstruction):
    eff_real = n_pass_real / n_analyzed
    # calculate statistical error of eff_real:
    eff_real_stat = np.sqrt(n_analyzed) / n_analyzed

    # calculate efficiency for ideal case (with numbers from MC truth):
    eff_ideal = n_pass_ideal / n_analyzed
    # calculate statistical error of eff_ideal:
    eff_ideal_stat = np.sqrt(n_analyzed) / n_analyzed

    # calculate alpha = eff_ideal / eff_real:
    alpha = eff_ideal / eff_real

    # calculate statistical error of alpha ('efficiency of the efficiency') with Gaussian Error Propagation:
    alpha_stat = np.sqrt((eff_real_stat / eff_ideal) ** 2 + (eff_real * eff_ideal_stat / eff_ideal ** 2) ** 2)

    # as cross-check: calculate the 'efficiency of the efficiency' in another way considering n_too_much and n_too_less:
    eff_leak = 1 + (n_too_less - n_too_much) / n_pass_real

    # check alpha and eff_leak:
    if np.abs(alpha - eff_leak) > 0.000001:
        sys.exit("alpha != eff_leak")

    return eff_real, eff_real_stat, eff_ideal, eff_ideal_stat, alpha, alpha_stat, n_total_evts


# get the date and time, when the script was run:
date = datetime.datetime.now()
now = date.strftime("%Y-%m-%d %H:%M")

# path, where the information about the cuts on the NC events is stored:
input_path_NC = "/home/astro/blum/juno/atmoNC/data_NC/output_detsim_v2/"
# path, where the information about the cuts on the IBD events is stored:
input_path_IBD = "/home/astro/blum/juno/IBD_events/"
# path, where the information about the PSD cuts is stored:
input_path_PSD = "/home/astro/blum/juno/atmoNC/data_NC/output_PSD_v2/"

""" set values of the cut parameters: """
# fiducial volume cut on prompt signal in mm:
R_cut_prompt = 16000.0
# prompt energy cut in MeV:
E_min_prompt = 10.0
E_max_prompt = 100.0
# alpha_2 (error of the efficiency, if you convert the prompt energy from PE to MeV) (see
# info_conversion_proton_neutron.odt):
alpha_E_prompt_cut_2 = 1.0 + (62.0 - 122.0) / 5125.0
# statistical error of alpha_2 (statistical error of the error of the efficiency, if you convert the prompt
# energy from PE to MeV):
alpha_stat_E_prompt_cut_2 = np.sqrt(5125.0) / 5125.0
# time cut between prompt and delayed signal in ns:
time_cut_min = 500.0
time_cut_max = 1000000.0
# multiplicity cut:
multiplicity = 1
# delayed energy cut in PE:
E_min_delayed = 2500.0
E_max_delayed = 3400.0
# distance cut between prompt and delayed signal in mm:
distance_cut = 1500.0
# fiducial volume cut on delayed signal in mm:
R_cut_delayed = 16000.0
# TTR value from Pulse Shape Analysis:
# TODO-me: include PSD:
PSD_NC_suppression = 0.0
PSD_ttr_cut_value = 0.0

# path, where the output files will be saved:
output_path = (input_path_NC + "results_{0:.0f}mm_{1:.0f}MeVto{2:.0f}MeV_{3:.0f}nsto{4:.0f}ms_mult{5:d}_"
                               "{6:.0f}PEto{7:.0f}PE_dist{8:.0f}mm_R{9:.0f}mm_PSD{10:.0f}/"
               .format(R_cut_prompt, E_min_prompt, E_max_prompt, time_cut_min, time_cut_max/1000000.0, multiplicity,
                       E_min_delayed, E_max_delayed, distance_cut, R_cut_delayed, PSD_NC_suppression*100.0))
# path, where the numbers_.txt files of the delayed cut for NC events are stored:
input_path_NC_del = (input_path_NC + "delayed_cut_{0:.0f}nsto{1:.0f}ms_mult{2:d}_{3:.0f}PEto{4:.0f}PE_dist{5:.0f}mm_"
                                     "R{6:.0f}mm/".format(time_cut_min, time_cut_max/1000000.0, multiplicity,
                                                          E_min_delayed, E_max_delayed, distance_cut, R_cut_delayed))
# path, where the numbers_.txt file of the delayed cut for IBD events are stored:
input_path_IBD_del = (input_path_IBD + "delayed_cut_{0:.0f}nsto{1:.0f}ms_mult{2:d}_{3:.0f}PEto{4:.0f}PE_dist{5:.0f}mm_"
                                       "R{6:.0f}mm/".format(time_cut_min, time_cut_max/1000000.0, multiplicity,
                                                            E_min_delayed, E_max_delayed, distance_cut, R_cut_delayed))

""" prompt signal energy window: """
# bin-width of visible energy for histogram in MeV (must be the same like for DM signal, Reactor, CCatmo and DSNB;
# 100 keV = 0.1 MeV):
bin_width_energy = 0.5

# radius cut in m:
r_cut = R_cut_prompt / 1000.0
# time exposure in years:
time_in_years = 10
# time exposure in seconds:
time_seconds = time_in_years * 3.156 * 10 ** 7
# set booleans, that define, which plots are shown or saved (boolean):
PLOT_FLUX = False
SHOW_FLUXPLOT = True
SAVE_FLUXPLOT = False
PLOT_EVT_RATE = False
SHOW_EVT_RATE = True
SAVE_EVT_RATE = False

""" analyze the efficiencies of the cuts, define like above, for NC events: """
# volume cut on prompt signal:
string_numbers_volume_cut_NC = input_path_NC + "numbers_volume_cut_atmoNC_{0:.0f}mm.txt".format(R_cut_prompt)
(eff_real_volume_cut_NC, eff_real_stat_volume_cut_NC, eff_ideal_volume_cut_NC, eff_ideal_stat_volume_cut_NC,
 alpha_volume_cut_NC, alpha_stat_volume_cut_NC, number_events_total_NC) = \
    read_numbers_txt_files(string_numbers_volume_cut_NC)

# prompt energy cut:
string_numbers_E_prompt_NC = (input_path_NC_del + "numbers_prompt_energy_cut_atmoNC_{0:.0f}MeV_to_{1:.0f}MeV_0.txt"
                              .format(E_min_prompt, E_max_prompt))
(eff_real_E_prompt_cut_NC, eff_real_stat_E_prompt_cut_NC, eff_E_prompt_conversion_NC, eff_stat_E_prompt_conversion_NC,
 alpha_E_prompt_cut_1_NC, alpha_stat_E_prompt_cut_1_NC, number_events_total_NC) = \
    read_numbers_txt_files(string_numbers_E_prompt_NC)
# calculate eff_ideal (n_pass_ideal / n_analyzed) with alpha_2 * eff_conversion:
eff_ideal_E_prompt_cut_NC = alpha_E_prompt_cut_2 * eff_E_prompt_conversion_NC
# calculate statistical error of eff_ideal with Gaussian Error Propagation:
eff_ideal_stat_E_prompt_cut_NC = np.sqrt((eff_E_prompt_conversion_NC * alpha_stat_E_prompt_cut_2) ** 2 +
                                         (alpha_E_prompt_cut_2 * eff_stat_E_prompt_conversion_NC) ** 2)
# calculate alpha of prompt energy cut (alpha = alpha_1 * alpha_2, product of alpha of energy resolution and alpha of
# conversion):
alpha_E_prompt_cut_NC = alpha_E_prompt_cut_1_NC * alpha_E_prompt_cut_2
# calculate statistical error of alpha with Gaussian error propagation:
alpha_stat_E_prompt_cut_NC = np.sqrt((alpha_E_prompt_cut_2 * alpha_stat_E_prompt_cut_1_NC) ** 2 +
                                     (alpha_E_prompt_cut_1_NC * alpha_stat_E_prompt_cut_2) ** 2)

# time cut between prompt und delayed signal:
string_numbers_time_NC = (input_path_NC_del + "numbers_time_cut_atmoNC_{0:.0f}ns_to_{1:.0f}ms_0.txt"
                          .format(time_cut_min, time_cut_max/1000000.0))
(eff_real_time_cut_NC, eff_real_stat_time_cut_NC, eff_ideal_time_cut_NC, eff_ideal_stat_time_cut_NC, alpha_time_cut_NC,
 alpha_stat_time_cut_NC, number_events_total_NC) = read_numbers_txt_files(string_numbers_time_NC)

# multiplicity cut:
string_numbers_mult_NC = input_path_NC_del + "numbers_multiplicity_cut_atmoNC_mult{0:d}_0.txt".format(multiplicity)
(eff_real_mult_cut_NC, eff_real_stat_mult_cut_NC, eff_ideal_mult_cut_NC, eff_ideal_stat_mult_cut_NC, alpha_mult_cut_NC,
 alpha_stat_mult_cut_NC, number_events_total_NC) = read_numbers_txt_files(string_numbers_mult_NC)

# delayed energy cut:
string_numbers_E_delayed_NC = (input_path_NC_del + "numbers_delayed_energy_cut_atmoNC_{0:.0f}PE_to_{1:.0f}PE_0.txt"
                               .format(E_min_delayed, E_max_delayed))
(eff_real_E_delayed_cut_NC, eff_real_stat_E_delayed_cut_NC, eff_ideal_E_delayed_cut_NC, eff_ideal_stat_E_delayed_cut_NC,
 alpha_E_delayed_cut_NC, alpha_stat_E_delayed_cut_NC, number_events_total_NC) = \
    read_numbers_txt_files(string_numbers_E_delayed_NC)

# distance cut between prompt and delayed signal:
string_numbers_dist_NC = input_path_NC_del + "numbers_distance_cut_atmoNC_{0:.0f}mm_0.txt".format(distance_cut)
(eff_real_dist_cut_NC, eff_real_stat_dist_cut_NC, eff_ideal_dist_cut_NC, eff_ideal_stat_dist_cut_NC, alpha_dist_cut_NC,
 alpha_stat_dist_cut_NC, number_events_total_NC) = read_numbers_txt_files(string_numbers_dist_NC)

# volume cut on delayed signal:
string_numbers_vol_delayed_NC = (input_path_NC_del + "numbers_delayed_volume_cut_atmoNC_{0:.0f}mm_0.txt"
                                 .format(R_cut_delayed))
(eff_real_vol_delayed_NC, eff_real_stat_vol_delayed_NC, eff_ideal_vol_delayed_NC, eff_ideal_stat_vol_delayed_NC,
 alpha_vol_delayed_NC, alpha_stat_vol_delayed_NC, number_events_total_NC) = \
    read_numbers_txt_files(string_numbers_vol_delayed_NC)

# TODO-me: include PSD suppression


# efficiency of delayed cut (real events = reconstructed events):
eff_real_delayed_cut_NC = (eff_real_time_cut_NC * eff_real_mult_cut_NC * eff_real_E_delayed_cut_NC *
                           eff_real_dist_cut_NC * eff_real_vol_delayed_NC)
# statistical error of efficiency of delayed cut (real) with Gaussian error propagation:
eff_real_stat_delayed_cut_NC = (eff_real_delayed_cut_NC *
                                np.sqrt((eff_real_stat_time_cut_NC / eff_real_time_cut_NC)**2 +
                                        (eff_real_stat_mult_cut_NC / eff_real_mult_cut_NC)**2 +
                                        (eff_real_stat_E_delayed_cut_NC / eff_real_E_delayed_cut_NC)**2 +
                                        (eff_real_stat_dist_cut_NC / eff_real_dist_cut_NC)**2 +
                                        (eff_real_stat_vol_delayed_NC / eff_real_vol_delayed_NC)**2))

# efficiency of delayed cut (ideal events = MC truth):
eff_ideal_delayed_cut_NC = (eff_ideal_time_cut_NC * eff_ideal_mult_cut_NC * eff_ideal_E_delayed_cut_NC *
                            eff_ideal_dist_cut_NC * eff_ideal_vol_delayed_NC)
# statistical error of efficiency of delayed cut (ideal) with Gaussian error propagation:
eff_ideal_stat_delayed_cut_NC = (eff_ideal_delayed_cut_NC *
                                 np.sqrt((eff_ideal_stat_time_cut_NC / eff_ideal_time_cut_NC)**2 +
                                         (eff_ideal_stat_mult_cut_NC / eff_ideal_mult_cut_NC)**2 +
                                         (eff_ideal_stat_E_delayed_cut_NC / eff_ideal_E_delayed_cut_NC)**2 +
                                         (eff_ideal_stat_dist_cut_NC / eff_ideal_dist_cut_NC)**2 +
                                         (eff_ideal_stat_vol_delayed_NC / eff_ideal_vol_delayed_NC)**2))

# alpha of delayed cut:
alpha_delayed_cut_NC = (alpha_time_cut_NC * alpha_mult_cut_NC * alpha_E_delayed_cut_NC * alpha_dist_cut_NC *
                        alpha_vol_delayed_NC)
# statistical error of alpha of delayed cut:
alpha_stat_delayed_cut_NC = (alpha_delayed_cut_NC *
                             np.sqrt((alpha_stat_time_cut_NC / alpha_time_cut_NC)**2 +
                                     (alpha_stat_mult_cut_NC / alpha_mult_cut_NC)**2 +
                                     (alpha_stat_E_delayed_cut_NC / alpha_E_delayed_cut_NC)**2 +
                                     (alpha_stat_dist_cut_NC / alpha_dist_cut_NC)**2 +
                                     (alpha_stat_vol_delayed_NC / alpha_vol_delayed_NC)**2))

""" analyze the efficiencies of the cuts, define like above, for IBD events: """
# volume cut on prompt signal:
string_numbers_volume_cut_IBD = input_path_IBD + "numbers_volume_cut_IBD_{0:.0f}mm.txt".format(R_cut_prompt)
(eff_real_volume_cut_IBD, eff_real_stat_volume_cut_IBD, eff_ideal_volume_cut_IBD, eff_ideal_stat_volume_cut_IBD,
 alpha_volume_cut_IBD, alpha_stat_volume_cut_IBD, number_events_total_IBD) = \
    read_numbers_txt_files(string_numbers_volume_cut_IBD)

# prompt energy cut:
string_numbers_E_prompt_IBD = (input_path_IBD_del + "numbers_prompt_energy_cut_IBD_{0:.0f}MeV_to_{1:.0f}MeV_0.txt"
                               .format(E_min_prompt, E_max_prompt))
(eff_real_E_prompt_cut_IBD, eff_real_stat_E_prompt_cut_IBD, eff_E_prompt_conversion_IBD,
 eff_stat_E_prompt_conversion_IBD, alpha_E_prompt_cut_1_IBD, alpha_stat_E_prompt_cut_1_IBD, number_events_total_IBD) = \
    read_numbers_txt_files(string_numbers_E_prompt_IBD)
# calculate eff_ideal (n_pass_ideal / n_analyzed) with alpha_2 * eff_conversion:
eff_ideal_E_prompt_cut_IBD = alpha_E_prompt_cut_2 * eff_E_prompt_conversion_IBD
# calculate statistical error of eff_ideal with Gaussian Error Propagation:
eff_ideal_stat_E_prompt_cut_IBD = np.sqrt((eff_E_prompt_conversion_IBD * alpha_stat_E_prompt_cut_2) ** 2 +
                                          (alpha_E_prompt_cut_2 * eff_stat_E_prompt_conversion_IBD) ** 2)
# calculate alpha of prompt energy cut (alpha = alpha_1 * alpha_2, product of alpha of energy resolution and alpha of
# conversion):
alpha_E_prompt_cut_IBD = alpha_E_prompt_cut_1_IBD * alpha_E_prompt_cut_2
# calculate statistical error of alpha with Gaussian error propagation:
alpha_stat_E_prompt_cut_IBD = np.sqrt((alpha_E_prompt_cut_2 * alpha_stat_E_prompt_cut_1_IBD) ** 2 +
                                      (alpha_E_prompt_cut_1_IBD * alpha_stat_E_prompt_cut_2) ** 2)

# time cut between prompt und delayed signal:
string_numbers_time_IBD = (input_path_IBD_del + "numbers_time_cut_IBD_{0:.0f}ns_to_{1:.0f}ms_0.txt"
                          .format(time_cut_min, time_cut_max/1000000.0))
(eff_real_time_cut_IBD, eff_real_stat_time_cut_IBD, eff_ideal_time_cut_IBD, eff_ideal_stat_time_cut_IBD,
 alpha_time_cut_IBD, alpha_stat_time_cut_IBD, number_events_total_IBD) = \
    read_numbers_txt_files(string_numbers_time_IBD)

# multiplicity cut:
string_numbers_mult_IBD = input_path_IBD_del + "numbers_multiplicity_cut_IBD_mult{0:d}_0.txt".format(multiplicity)
(eff_real_mult_cut_IBD, eff_real_stat_mult_cut_IBD, eff_ideal_mult_cut_IBD, eff_ideal_stat_mult_cut_IBD,
 alpha_mult_cut_IBD, alpha_stat_mult_cut_IBD, number_events_total_IBD) = \
    read_numbers_txt_files(string_numbers_mult_IBD)

# delayed energy cut:
string_numbers_E_delayed_IBD = (input_path_IBD_del + "numbers_delayed_energy_cut_IBD_{0:.0f}PE_to_{1:.0f}PE_0.txt"
                                .format(E_min_delayed, E_max_delayed))
(eff_real_E_delayed_cut_IBD, eff_real_stat_E_delayed_cut_IBD, eff_ideal_E_delayed_cut_IBD,
 eff_ideal_stat_E_delayed_cut_IBD, alpha_E_delayed_cut_IBD, alpha_stat_E_delayed_cut_IBD, number_events_total_IBD) = \
    read_numbers_txt_files(string_numbers_E_delayed_IBD)

# distance cut between prompt and delayed signal:
string_numbers_dist_IBD = input_path_IBD_del + "numbers_distance_cut_IBD_{0:.0f}mm_0.txt".format(distance_cut)
(eff_real_dist_cut_IBD, eff_real_stat_dist_cut_IBD, eff_ideal_dist_cut_IBD, eff_ideal_stat_dist_cut_IBD,
 alpha_dist_cut_IBD, alpha_stat_dist_cut_IBD, number_events_total_IBD) = \
    read_numbers_txt_files(string_numbers_dist_IBD)

# volume cut on delayed signal:
string_numbers_vol_delayed_IBD = (input_path_IBD_del + "numbers_delayed_volume_cut_IBD_{0:.0f}mm_0.txt"
                                  .format(R_cut_delayed))
(eff_real_vol_delayed_IBD, eff_real_stat_vol_delayed_IBD, eff_ideal_vol_delayed_IBD, eff_ideal_stat_vol_delayed_IBD,
 alpha_vol_delayed_IBD, alpha_stat_vol_delayed_IBD, number_events_total_IBD) = \
    read_numbers_txt_files(string_numbers_vol_delayed_IBD)

# TODO-me: include PSD suppression

# efficiency of delayed cut (real events = reconstructed events):
eff_real_delayed_cut_IBD = (eff_real_time_cut_IBD * eff_real_mult_cut_IBD * eff_real_E_delayed_cut_IBD *
                            eff_real_dist_cut_IBD * eff_real_vol_delayed_IBD)
# statistical error of efficiency of delayed cut (real) with Gaussian error propagation:
eff_real_stat_delayed_cut_IBD = (eff_real_delayed_cut_IBD *
                                 np.sqrt((eff_real_stat_time_cut_IBD / eff_real_time_cut_IBD)**2 +
                                         (eff_real_stat_mult_cut_IBD / eff_real_mult_cut_IBD)**2 +
                                         (eff_real_stat_E_delayed_cut_IBD / eff_real_E_delayed_cut_IBD)**2 +
                                         (eff_real_stat_dist_cut_IBD / eff_real_dist_cut_IBD)**2 +
                                         (eff_real_stat_vol_delayed_IBD / eff_real_vol_delayed_IBD)**2))

# efficiency of delayed cut (ideal events = MC truth):
eff_ideal_delayed_cut_IBD = (eff_ideal_time_cut_IBD * eff_ideal_mult_cut_IBD * eff_ideal_E_delayed_cut_IBD *
                             eff_ideal_dist_cut_IBD * eff_ideal_vol_delayed_IBD)
# statistical error of efficiency of delayed cut (ideal) with Gaussian error propagation:
eff_ideal_stat_delayed_cut_IBD = (eff_ideal_delayed_cut_IBD *
                                  np.sqrt((eff_ideal_stat_time_cut_IBD / eff_ideal_time_cut_IBD)**2 +
                                          (eff_ideal_stat_mult_cut_IBD / eff_ideal_mult_cut_IBD)**2 +
                                          (eff_ideal_stat_E_delayed_cut_IBD / eff_ideal_E_delayed_cut_IBD)**2 +
                                          (eff_ideal_stat_dist_cut_IBD / eff_ideal_dist_cut_IBD)**2 +
                                          (eff_ideal_stat_vol_delayed_IBD / eff_ideal_vol_delayed_IBD)**2))

# alpha of delayed cut:
alpha_delayed_cut_IBD = (alpha_time_cut_IBD * alpha_mult_cut_IBD * alpha_E_delayed_cut_IBD * alpha_dist_cut_IBD *
                         alpha_vol_delayed_IBD)
# statistical error of alpha of delayed cut:
alpha_stat_delayed_cut_IBD = (alpha_delayed_cut_IBD *
                              np.sqrt((alpha_stat_time_cut_IBD / alpha_time_cut_IBD)**2 +
                                      (alpha_stat_mult_cut_IBD / alpha_mult_cut_IBD)**2 +
                                      (alpha_stat_E_delayed_cut_IBD / alpha_E_delayed_cut_IBD)**2 +
                                      (alpha_stat_dist_cut_IBD / alpha_dist_cut_IBD)**2 +
                                      (alpha_stat_vol_delayed_IBD / alpha_vol_delayed_IBD)**2))

# TODO-me: write information about efficiencies and alphas to file or print it!

""" Analyze IBD-like NC events (only information from filenumber_evtID_...txt files is needed) and create spectrum: """
# number of IBD-like events from real data before event rate (without PSD suppression and without alpha (MCtruth/real)):
number_IBDlike_real_simu_wo_PSD_wo_alpha = 0
# number of IBD-like events from real data before event rate (with PSD suppression, but without alpha (MCtruth/real)):
number_IBDlike_simu_w_PSD_wo_alpha = 0

# preallocate array, where Evis of IBD-like events (without PSD suppression and without alpha) is saved in MeV:
Evis_array_wo_PSD_wo_alpha = []
# preallocate array, where Evis of IBD-like events (with PSD suppression, but without alpha) is saved in MeV:
Evis_array_w_PSD_wo_alpha = []

# load files, where filenumber and evtID (and Evis) of events that pass the cut are stored:
array_pass_volume_cut = np.loadtxt(input_path_NC_del + "filenumber_evtID_volume_cut_atmoNC_{0:.0f}mm.txt"
                                   .format(R_cut_prompt))
array_pass_E_prompt_cut = np.loadtxt(input_path_NC_del + "filenumber_evtID_Evis_prompt_energy_cut_atmoNC_{0:.0f}MeV_"
                                                         "to_{1:.0f}MeV.txt".format(E_min_prompt, E_max_prompt))
array_pass_delayed_cut = np.loadtxt(input_path_NC_del + "filenumber_evtID_delayed_cut_atmoNC_{0:.0f}ns_to_{1:.0f}ms_"
                                                        "mult{2:d}_{3:.0f}PE_{4:.0f}PE_dist{5:.0f}mm_R{6:.0f}mm.txt"
                                    .format(time_cut_min, time_cut_max/1000000.0, multiplicity, E_min_delayed,
                                            E_max_delayed, distance_cut, R_cut_delayed))
# TODO-me: include PSD suppression:
array_pass_PSD = np.loadtxt()

# preallocate start indices of the arrays:
index_volume = 0
index_E_prompt = 0
index_PSD = 0

# loop over array_pass_delayed_cut:
for index in range(len(array_pass_delayed_cut)):

    # get filenumber and evtID of the event that passed delayed cut:
    filenumber_delayed_cut = array_pass_delayed_cut[index][0]
    evtID_delayed_cut = array_pass_delayed_cut[index][1]

    # loop over array_pass_volume_cut:
    for index1 in range(index_volume, len(array_pass_volume_cut), 1):

        # get filenumber and evtID of the event that passed volume cut:
        filenumber_volume_cut = array_pass_volume_cut[index1][0]
        evtID_volume_cut = array_pass_volume_cut[index1][1]

        # check if event also passed delayed cut:
        if filenumber_volume_cut == filenumber_delayed_cut and evtID_volume_cut == evtID_delayed_cut:

            # event passed volume cut AND delayed cut:

            # loop over array_pass_E_prompt_cut:
            for index2 in range(index_E_prompt, len(array_pass_E_prompt_cut), 1):

                # get filenumber, evtID and Evis of the event that passed prompt energy cut:
                filenumber_E_cut = array_pass_E_prompt_cut[index2][0]
                evtID_E_cut = array_pass_E_prompt_cut[index2][1]
                E_vis = array_pass_E_prompt_cut[index2][2]

                # check if event also passed volume and delayed cut:
                if filenumber_E_cut == filenumber_volume_cut and evtID_E_cut == evtID_volume_cut:

                    # event passed delayed, volume and prompt energy cut -> IBD-like event

                    # increment number_IBDlike_simu_wo_PSD_wo_alpha:
                    number_IBDlike_simu_wo_PSD_wo_alpha += 1
                    # append E_vis to Evis_array_wo_PSD:
                    Evis_array_wo_PSD.append(E_vis)

                    # loop over array_pass_PSD to check if event also passed PSD cut:
                    for index3 in range(index_PSD, len(array_pass_PSD), 1):

                        # get filenumber and evtID pf the event that passed PSD cut:
                        filenumber_PSD = array_pass_PSD[index3][0]
                        evtID_PSD = array_pass_PSD[index3][1]

                        # check if event also passed delayed, volume and prompt E cut:
                        if filenumber_PSD == filenumber_E_cut and evtID_PSD == filenumber_PSD:

                            # event passed all cuts:

                            # increment number_IBDlike_simu_w_PSD_wo_alpha:
                            number_IBDlike_simu_w_PSD_wo_alpha += 1
                            # append E_vis to Evis_array_w_PSD:
                            Evis_array_w_PSD.append(E_vis)

                            # set index_PSD = index3 -> start loop for next event at index_PSD:
                            index_PSD = index3
                            break

                        else:
                            # event passed delayed, volume and prompt energy cut, but is rejected by PSD cut.
                            # go to next entry of array_pass_PSD:
                            continue

                    # set index_E_prompt = index2 -> start loop for next event at index_E_prompt:
                    index_E_prompt = index2
                    break

                else:
                    # event passed delayed and volume cut, but is rejected by prompt energy cut.
                    # go to next entry of array_pass_E_prompt_cut:
                    continue

            # set index_volume = index1 -> start loop for next event at index_volume:
            index_volume = index1
            break

        else:
            # event passed delayed cut, but is rejected by volume cut.
            # go to next entry of array_pass_volume_cut:
            continue

""" Build histograms from E_vis_arrays: """
# set bin-edges of e_vis histogram in MeV:
bins_evis = np.arange(E_min_prompt, E_max_prompt + 2*bin_width_energy, bin_width_energy)

# build histogram from Evis_array_wo_PSD (without PSD suppression and alpha (IBD cuts and PSD)):
histo_Evis_wo_PSD_wo_alpha, bin_edges_evis = np.histogram(Evis_array_wo_PSD, bins_evis)

# build histogram from Evis_array_w_PSD (with PSD suppression, but without alpha (IBD cuts and PSD)):
histo_Evis_w_PSD_wo_alpha, bin_edges_evis = np.histogram(Evis_array_w_PSD, bins_evis)

print("number of NC events from simulation = {0:.0f}".format(number_events_total_NC))
print("number of IBD-like events from simulation (wo PSD, wo alpha) = {0:d}"
      .format(number_IBDlike_simu_wo_PSD_wo_alpha))
print("number of IBD-like events from simulation (w PSD, wo alpha) = {0:d}"
      .format(number_IBDlike_simu_w_PSD_wo_alpha))

""" consider alpha from cut efficiencies: """








# spectrum of all simulated IBD-like events (cut efficiencies are considered):
Evis_histo = Evis_histo_without_eff * cut_efficiency/100.0
# number of simulated IBD-like events (cut efficiencies are considered):
number_IBDlike_events_simu = number_IBDlike_events_simu_without_eff * cut_efficiency/100.0
print("number of IBD-like events from simulation (with cut efficiency) = {0:.2f}"
      .format(number_IBDlike_events_simu))

""" Event rate calculation: """
# calculate the theoretical event rate in NC events/sec in JUNO for neutrino energies from 0 MeV to 10 GeV (float)
# (event_rate = A * (flux_nue*xsec_nue + flux_nuebar*xsec_nuebar + flux_numu*xsec_numu + flux_numubar*xsec_numubar)):
event_rate = NC_background_functions.event_rate(bin_width_energy, r_cut, output_path, PLOT_FLUX, SHOW_FLUXPLOT,
                                                SAVE_FLUXPLOT, PLOT_EVT_RATE, SHOW_EVT_RATE, SAVE_EVT_RATE)

# number of NC events in JUNO after 10 years:
number_NC_events_JUNO = event_rate * time_seconds

# number of IBD-like events in JUNO after 10 years (cut efficiencies are considered):
number_IBDlike_events_JUNO = int(number_NC_events_JUNO * number_IBDlike_events_simu / number_NC_events_simu)

# normalize the spectrum of IBD-like events to the spectrum, JUNO will measure after 10 years (cut efficiencies are
# considered):
Evis_histo_JUNO = float(number_IBDlike_events_JUNO) / float(number_IBDlike_events_simu) * Evis_histo

""" display simulated spectrum: """
h1 = plt.figure(1, figsize=(15, 8))
plt.plot(bins_evis[:-1], Evis_histo, drawstyle="steps", linestyle="-", color="orange",
         label="atmospheric NC background\n(number of events = {0:.1f},\ncut efficiency = {1:.2f} % $\\pm$ {2:.2f} %)"
         .format(number_IBDlike_events_simu, cut_efficiency, error_cut_efficiency))
plt.xlim(xmin=min_energy, xmax=max_energy)
plt.ylim(ymin=0.0)
plt.xlabel("visible energy of prompt signal in MeV")
plt.ylabel("number of IBD-like events per bin (bin-width = {0:.2f} MeV)".format(bin_width_energy))
plt.title("Simulated spectrum of atmospheric NC neutrino events with IBD-like signature")
plt.legend()
plt.grid()
plt.savefig(output_path + "atmoNC_spectrum_simulated_bins{0:.0f}keV.png".format(bin_width_energy*1000))
# plt.show()
plt.close()

""" display visible spectrum in JUNO after 10 years: """
h2 = plt.figure(2, figsize=(15, 8))
plt.plot(bins_evis[:-1], Evis_histo_JUNO, drawstyle="steps", linestyle="-", color="orange",
         label="atmospheric NC background\n(number of events = {0:.1f},\ncut efficiency = {1:.2f} % $\\pm$ {2:.2f} %)"
         .format(number_IBDlike_events_JUNO, cut_efficiency, error_cut_efficiency))
plt.xlim(xmin=min_energy, xmax=max_energy)
plt.ylim(ymin=0.0)
plt.xlabel("visible energy of prompt signal in MeV")
plt.ylabel("number of IBD-like events per bin (bin-width = {0:.2f} MeV)".format(bin_width_energy))
plt.title("Expected spectrum of atmospheric NC neutrino events with IBD-like signature in JUNO after {0:.0f} years"
          .format(time_in_years))
plt.legend()
plt.grid()
plt.savefig(output_path + "atmoNC_spectrum_JUNO_bins{0:.0f}keV.png".format(bin_width_energy*1000))
# plt.show()
plt.close()

""" display visible spectrum in JUNO after 10 years after Pulse Shape Discrimination: """
# calculate spectrum after PSD:
Evis_histo_JUNO_PSD = Evis_histo_JUNO * (100.0-efficiency_PSD)/100.0
# number of events after PSD:
number_IBDlike_events_JUNO_PSD = number_IBDlike_events_JUNO * (100.0-efficiency_PSD)/100.0

h3 = plt.figure(3, figsize=(15, 8))
plt.plot(bins_evis[:-1], Evis_histo_JUNO_PSD, drawstyle="steps", linestyle="-", color="orange",
         label="atmospheric NC background after PSD\n(number of events = {0:0.1f})"
         .format(number_IBDlike_events_JUNO_PSD))
plt.xlim(xmin=min_energy, xmax=max_energy)
plt.ylim(ymin=0.0)
plt.xlabel("visible energy of prompt signal in MeV")
plt.ylabel("number of IBD-like events per bin (bin-width = {0:.2f} MeV)".format(bin_width_energy))
plt.title("Expected spectrum of atmospheric NC neutrino events with IBD-like signature in JUNO after {0:.0f} years\n"
          "after Pulse Shape Discrimination".format(time_in_years) + " (PSD efficiency $\\epsilon_{NC}$ = " +
          "{0:0.1f} %)"
          .format(efficiency_PSD))
plt.legend()
plt.grid()
plt.savefig(output_path + "atmoNC_spectrum_JUNO_afterPSD{1:.0f}_bins{0:.0f}keV.png".format(bin_width_energy*1000,
                                                                                           efficiency_PSD))
# plt.show()
plt.close()

""" save e_vis_array_JUNO to txt file (txt file must have same shape like files in folder 
    /home/astro/blum/PhD/work/MeVDM_JUNO/gen_spectrum_v2/).
    Save the array before Pulse Shape Discrimination. PSD efficiencies for real IBD signals and NC events is then 
    applied afterwards before calculating the Limits: """
# save Evis_histo_JUNO to txt-spectrum-file and information about simulation in txt-info-file:
print("... save data of spectrum to file...")
np.savetxt(output_path + 'NCatmo_onlyC12_bin{0:.0f}keV.txt'
           .format(bin_width_energy * 1000), Evis_histo_JUNO, fmt='%1.5e',
           header='Spectrum in IBD-like events/bin of atmospheric NC background events that mimic IBD signals '
                  '(calculated with atmoNC_spectrum_v1.py, {0}):'
                  '\n{3:d} NC events are simulated with JUNO detector software (tut_detsim.py).'
                  '\nNumber of IBD-like NC events in JUNO = {1:.2f}, binning of E_visible = {2:.3f} MeV,'
                  '\nNC interactions of nu_e, nu_e_bar, nu_mu and nu_mu_bar with C12 of liquid scintillator are '
                  'simulated with GENIE.'
                  '\nDeexcitation of residual isotopes are simulated with modified DSNB-NC.exe generator.'
                  '\nThen the final products are simulated with JUNO detector simulation and cuts are applied to get'
                  '\nthe number of NC events that mimic an IBD signal:'
                  '\ntotal cut efficiency = {4:.2f} %,'
                  '\nstatistical error of total cut efficiency = {11:.2f} %'
                  '\nefficiency volume cut = {5:.2f} %,'
                  '\nefficiency prompt energy cut = {6:.2f} %,'
                  '\nefficiency time cut = {7:.2f} %,'
                  '\nefficiency delayed energy cut = {8:.2f} %,'
                  '\nefficiency neutron multiplicity cut = {9:.2f} %,'
                  '\nefficiency distance cut = {10:.2f} %,'
           .format(now, number_IBDlike_events_JUNO, bin_width_energy, number_NC_events_simu, cut_efficiency,
                   efficiency_volume_cut, efficiency_prompt_energy_cut, efficiency_time_cut,
                   efficiency_delayed_energy_cut, efficiency_neutron_multiplicity_cut, efficiency_distance_cut,
                   error_cut_efficiency))
np.savetxt(output_path + 'NCatmo_info_onlyC12_bin{0:.0f}keV.txt'
           .format(bin_width_energy * 1000),
           np.array([min_energy, max_energy, bin_width_energy, time_in_years, r_cut, number_NC_events_simu,
                     number_IBDlike_events_JUNO, event_rate]),
           fmt='%1.9e',
           header='Information to simulation NCatmo_onlyC12_bin{0:.0f}keV.txt (analyzed files: user_atmoNC_{1:d}.root '
                  'to user_atmoNC_{2:d}.root):\n'
                  'values below: E_visible[0] in MeV, E_visible[-1] in MeV, interval_E_visible in MeV,'
                  '\nexposure time t_years in years, applied volume cut for radius in meter,'
                  '\nnumber of simulated NC events, number of IBD-like NC events in spectrum JUNO will measure, '
                  '\ntheoretical NC event rate in JUNO detector in NC events/sec,'
           .format(bin_width_energy * 1000, first_file, last_file))

