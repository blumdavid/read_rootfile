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

    5.  For IBD events and NC events that pass all cuts:
        take Evis of prompt signal and build energy spectrum (without PSD)

    6.  Consider the event rate of NC interactions on C12 inside the detector (calculated with cross-sections and
        neutrino fluxes) and calculate the 'real' spectrum of atmospheric NC background, JUNO will measure after 10
        years of data taking

    7.  Consider PSD efficiency:
        PSD efficiency depends on energy of prompt signal. To analyze PSD efficiency save filenumber, evtID and Evis
        of all IBD-like NC events in txt file and of all IBD events (that pass all IBD cuts) in txt file

    7.  Consider the PSD efficiency and calculate the spectrum of atmospheric NC background, JUNO will measure
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

    return eff_real, eff_real_stat, eff_ideal, eff_ideal_stat, alpha, alpha_stat, n_total_evts, n_analyzed


def write_info_to_file(out_path, event_type, cut_type, eff_real, stat_eff_real, eff_ideal, stat_eff_ideal, alpha,
                       stat_alpha, num_analyzed):
    """

    :param out_path: output path, where files are saved
    :param event_type: defines the events type (NC or IBD) (string)
    :param cut_type:    defines the cut (volume, prompt energy, time, multiplicity, delayed energy, distance, delayed
                        volume, delayed) with the cut parameters (string)
    :param eff_real: efficiency of real data
    :param stat_eff_real: statistical error of efficiency of real data
    :param eff_ideal: efficiency of ideal data
    :param stat_eff_ideal: statistical error of efficiency of ideal data
    :param alpha: 'error of the efficiency' (eff_ideal/eff_real)
    :param stat_alpha: statistical error of alpha
    :param num_analyzed: number of analyzed event (number events before the cut)
    :return:
    """
    np.savetxt(out_path + "efficiencies_{0}_{1}.txt".format(event_type, cut_type),
               np.array([eff_real, stat_eff_real, eff_ideal, stat_eff_ideal, alpha, stat_alpha]), fmt='%.5f',
               header="Efficiencies and alpha of the {0} cut of {1:.0f} {2} events (analyzed with "
                      "atmoNC_spectrum_v2.py):"
                      "\nefficiency of real data,"
                      "\nstatistical error of efficiency of real data,"
                      "\nefficiency of ideal data,"
                      "\nstatistical error of efficiency of ideal data,"
                      "\nalpha (eff_ideal/eff_real),"
                      "\nstatistical error of alpha, "
                      "\nnumber of analyzed events:"
               .format(cut_type, num_analyzed, event_type))

    return


# get the date and time, when the script was run:
date = datetime.datetime.now()
now = date.strftime("%Y-%m-%d %H:%M")

# path, where the information about the cuts on the NC events is stored:
input_path_NC = "/home/astro/blum/juno/atmoNC/data_NC/output_detsim_v2/"
# path, where the information about the cuts on the IBD events is stored:
input_path_IBD = "/home/astro/blum/juno/IBD_events/"
# path, where the information about the PSD cuts is stored:
input_path_PSD = "/home/astro/blum/juno/atmoNC/data_NC/output_PSD_v2/DCR/"
# input_path_PSD = "/home/astro/blum/juno/atmoNC/data_NC/output_PSD_v2/"

""" set values of the cut parameters: """
# fiducial volume cut on prompt signal in mm:
R_cut_prompt = 16000.0
# prompt energy cut in MeV:
E_min_prompt = 10.0
E_max_prompt = 100.0
# alpha_2 (error of the efficiency, if you convert the prompt energy from PE to MeV) (see
# info_conversion_proton_neutron.odt):
alpha_E_prompt_cut = 1.0 + (60.0 - 122.0) / 5124.0
# statistical error of alpha_2 (statistical error of the error of the efficiency, if you convert the prompt
# energy from PE to MeV):
alpha_stat_E_prompt_cut = np.sqrt(5124.0) / 5124.0
# time cut between prompt and delayed signal in ns:
time_cut_min = 1000.0
time_cut_max = 1000000.0
# multiplicity cut:
multiplicity = 1
# delayed energy cut in MeV:
E_min_delayed_MeV = 1.8
E_max_delayed_MeV = 2.55
# distance cut between prompt and delayed signal in mm:
distance_cut = 500.0
# fiducial volume cut on delayed signal in mm:
R_cut_delayed = 17700.0
# values from Pulse Shape Analysis:
# PSD NC suppression from /output_PSD_v2/:
# PSD_NC_suppression = 0.99
PSD_NC_suppression = 0.90

# start of tail, where least IBD events are cut away:
# tail_start = [275.0]
tail_start = [25.0, 25.0, 25.0, 75.0, 75.0, 75.0, 375.0, 375.0, 375.0, 400.0, 400.0, 400.0]

# end of tail, where least IBD events are cut away:
# tail_end = [600.0]
tail_end = [600.0, 800.0, 1000.0, 600.0, 800.0, 1000.0, 600.0, 800.0, 1000.0, 600.0, 800.0, 1000.0]

# tail to total ratio cut value:
# TTR_cut_value = [0.01662]
TTR_cut_value = [0.29019, 0.29438, 0.29654, 0.10727, 0.10901, 0.10977, 0.00431, 0.00544, 0.00588, 0.00347, 0.00460,
                 0.00503]

# set variable with the last file to be analyzed:
last_file_NC = 999
last_file_IBD = 199

# path, where the output files will be saved:
output_path = (input_path_NC + "DCR_results_{0:.0f}mm_{1:.0f}MeVto{2:.0f}MeV_{3:.0f}nsto{4:.0f}ms_mult{5:d}_"
                               "{6:.0f}keVto{7:.0f}keV_dist{8:.0f}mm_R{9:.0f}mm_PSD{10:.0f}/"
               .format(R_cut_prompt, E_min_prompt, E_max_prompt, time_cut_min, time_cut_max/1000000.0, multiplicity,
                       E_min_delayed_MeV*1000, E_max_delayed_MeV*1000, distance_cut, R_cut_delayed,
                       PSD_NC_suppression*100.0))
# path, where the numbers_.txt files of the delayed cut for NC events are stored:
input_path_NC_del = (input_path_NC + "delayed_cut_{0:.0f}nsto{1:.0f}ms_mult{2:d}_{3:.0f}keVto{4:.0f}keV_dist{5:.0f}mm_"
                                     "R{6:.0f}mm/".format(time_cut_min, time_cut_max/1000000.0, multiplicity,
                                                          E_min_delayed_MeV*1000, E_max_delayed_MeV*1000,
                                                          distance_cut, R_cut_delayed))
# path, where the numbers_.txt files of the delayed cut for IBD events are stored:
input_path_IBD_del = (input_path_IBD + "delayed_cut_{0:.0f}nsto{1:.0f}ms_mult{2:d}_{3:.0f}keVto{4:.0f}keV_"
                                       "dist{5:.0f}mm_R{6:.0f}mm/"
                      .format(time_cut_min, time_cut_max/1000000.0, multiplicity, E_min_delayed_MeV*1000,
                              E_max_delayed_MeV*1000, distance_cut, R_cut_delayed))

""" prompt signal energy window: """
# bin-width of visible energy for histogram in MeV (must be the same like for DM signal, Reactor, CCatmo and DSNB;
# 100 keV = 0.1 MeV):
bin_width_energy = 0.5

# number of events per root file:
n_events_per_rootfile = 100.0
# radius cut in m:
# INFO-me: take total volume here -> this is only used for event rate calculation
r_cut = 17.7
# time exposure in years:
time_in_years = 10
# time exposure in seconds:
time_seconds = time_in_years * 3.156 * 10 ** 7
# set booleans, that define, which plots are shown or saved (boolean):
PLOT_FLUX = False
SHOW_FLUXPLOT = False
SAVE_FLUXPLOT = False
PLOT_EVT_RATE = False
SHOW_EVT_RATE = False
SAVE_EVT_RATE = False

""" get information about PSD: """
# pulse shape discrimination cut on prompt signal:
string_numbers_PSD = input_path_PSD + "numbers_PSD.txt"
numbers_PSD = np.loadtxt(string_numbers_PSD)
number_events_analyzed_PSD_NC = numbers_PSD[1]
number_events_analyzed_PSD_IBD = numbers_PSD[3]

""" analyze the efficiencies of the cuts, defined like above, for NC events: """
# volume cut on prompt signal:
string_numbers_volume_cut_NC = input_path_NC + "numbers_volume_cut_atmoNC_{0:.0f}mm.txt".format(R_cut_prompt)
(eff_real_volume_cut_NC, eff_real_stat_volume_cut_NC, eff_ideal_volume_cut_NC, eff_ideal_stat_volume_cut_NC,
 alpha_volume_cut_NC, alpha_stat_volume_cut_NC, number_events_total_NC, number_analyzed_volume_cut_NC) = \
    read_numbers_txt_files(string_numbers_volume_cut_NC)

# prompt energy cut:
string_numbers_E_prompt_NC = (input_path_NC_del + "numbers_prompt_energy_cut_atmoNC_{0:.0f}MeV_to_{1:.0f}MeV_0.txt"
                              .format(E_min_prompt, E_max_prompt))
array_numbers_E_prompt_NC = np.loadtxt(string_numbers_E_prompt_NC)
number_events_total_NC = float(array_numbers_E_prompt_NC[0])
number_analyzed_E_prompt_NC = float(array_numbers_E_prompt_NC[2])
N_pass_real_E_prompt_NC = float(array_numbers_E_prompt_NC[3])
N_rej_real_E_prompt_NC = float(array_numbers_E_prompt_NC[4])
# check the numbers:
if N_pass_real_E_prompt_NC + N_rej_real_E_prompt_NC != number_analyzed_E_prompt_NC:
    sys.exit("N_pass_real_E_prompt_NC + N_rej_real_E_prompt_NC != number_analyzed_E_prompt_NC")
# calculate efficiency for real case (with numbers from conversion):
eff_real_E_prompt_cut_NC = N_pass_real_E_prompt_NC / number_analyzed_E_prompt_NC
# calculate statistical error of eff_real:
eff_real_stat_E_prompt_cut_NC = np.sqrt(number_analyzed_E_prompt_NC) / number_analyzed_E_prompt_NC
# get alpha of prompt energy cut:
alpha_E_prompt_cut_NC = alpha_E_prompt_cut
alpha_stat_E_prompt_cut_NC = alpha_stat_E_prompt_cut
# calculate efficiency for ideal case (without conversion):
eff_ideal_E_prompt_cut_NC = alpha_E_prompt_cut_NC * eff_real_E_prompt_cut_NC
# calculate statistical error of eff_ideal:
eff_ideal_stat_E_prompt_cut_NC = np.sqrt((alpha_E_prompt_cut_NC * eff_real_stat_E_prompt_cut_NC)**2 +
                                         (eff_real_E_prompt_cut_NC * alpha_stat_E_prompt_cut_NC)**2)

# time cut between prompt und delayed signal:
string_numbers_time_NC = (input_path_NC_del + "numbers_time_cut_atmoNC_{0:.0f}ns_to_{1:.0f}ms_0.txt"
                          .format(time_cut_min, time_cut_max/1000000.0))
(eff_real_time_cut_NC, eff_real_stat_time_cut_NC, eff_ideal_time_cut_NC, eff_ideal_stat_time_cut_NC, alpha_time_cut_NC,
 alpha_stat_time_cut_NC, number_events_total_NC, number_analyzed_time_cut_NC) \
    = read_numbers_txt_files(string_numbers_time_NC)

# multiplicity cut:
string_numbers_mult_NC = input_path_NC_del + "numbers_multiplicity_cut_atmoNC_mult{0:d}_0.txt".format(multiplicity)
(eff_real_mult_cut_NC, eff_real_stat_mult_cut_NC, eff_ideal_mult_cut_NC, eff_ideal_stat_mult_cut_NC, alpha_mult_cut_NC,
 alpha_stat_mult_cut_NC, number_events_total_NC, number_analyzed_mult_cut_NC) \
    = read_numbers_txt_files(string_numbers_mult_NC)

# delayed energy cut:
string_numbers_E_delayed_NC = (input_path_NC_del + "numbers_delayed_energy_cut_atmoNC_{0:.0f}keV_to_{1:.0f}keV_0.txt"
                               .format(E_min_delayed_MeV*1000, E_max_delayed_MeV*1000))
array_numbers_E_delayed_NC = np.loadtxt(string_numbers_E_delayed_NC)
number_events_total_NC = array_numbers_E_delayed_NC[0]
number_analyzed_E_delayed_NC = array_numbers_E_delayed_NC[2]
N_pass_real_E_delayed_NC = float(array_numbers_E_delayed_NC[3])
N_rej_real_E_delayed_NC = float(array_numbers_E_delayed_NC[4])
# check the numbers:
if N_pass_real_E_delayed_NC + N_rej_real_E_delayed_NC != number_analyzed_E_delayed_NC:
    sys.exit("N_pass_real_E_prompt_NC + N_rej_real_E_prompt_NC != number_analyzed_E_prompt_NC")
# calculate efficiency for real case (with numbers from conversion):
eff_real_E_delayed_cut_NC = N_pass_real_E_delayed_NC / number_analyzed_E_delayed_NC
# calculate statistical error of eff_real:
eff_real_stat_E_delayed_cut_NC = np.sqrt(number_analyzed_E_delayed_NC) / number_analyzed_E_delayed_NC
# get alpha of delayed energy cut:
alpha_E_delayed_cut_NC = alpha_E_prompt_cut
alpha_stat_E_delayed_cut_NC = alpha_stat_E_prompt_cut
# calculate efficiency for ideal case (without conversion):
eff_ideal_E_delayed_cut_NC = alpha_E_delayed_cut_NC * eff_real_E_delayed_cut_NC
# calculate statistical error of eff_ideal:
eff_ideal_stat_E_delayed_cut_NC = np.sqrt((alpha_E_delayed_cut_NC * eff_real_stat_E_delayed_cut_NC)**2 +
                                          (eff_real_E_delayed_cut_NC * alpha_stat_E_delayed_cut_NC)**2)

# distance cut between prompt and delayed signal:
string_numbers_dist_NC = input_path_NC_del + "numbers_distance_cut_atmoNC_{0:.0f}mm_0.txt".format(distance_cut)
(eff_real_dist_cut_NC, eff_real_stat_dist_cut_NC, eff_ideal_dist_cut_NC, eff_ideal_stat_dist_cut_NC, alpha_dist_cut_NC,
 alpha_stat_dist_cut_NC, number_events_total_NC, number_analyzed_dist_cut_NC) \
    = read_numbers_txt_files(string_numbers_dist_NC)

# volume cut on delayed signal:
string_numbers_vol_delayed_NC = (input_path_NC_del + "numbers_delayed_volume_cut_atmoNC_{0:.0f}mm_0.txt"
                                 .format(R_cut_delayed))
(eff_real_vol_delayed_NC, eff_real_stat_vol_delayed_NC, eff_ideal_vol_delayed_NC, eff_ideal_stat_vol_delayed_NC,
 alpha_vol_delayed_NC, alpha_stat_vol_delayed_NC, number_events_total_NC, number_analyzed_vol_delayed_NC) = \
    read_numbers_txt_files(string_numbers_vol_delayed_NC)

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
 alpha_volume_cut_IBD, alpha_stat_volume_cut_IBD, number_events_total_IBD, number_analyzed_volume_cut_IBD) = \
    read_numbers_txt_files(string_numbers_volume_cut_IBD)

# prompt energy cut:
string_numbers_E_prompt_IBD = (input_path_IBD_del + "numbers_prompt_energy_cut_IBD_{0:.0f}MeV_to_{1:.0f}MeV_0.txt"
                               .format(E_min_prompt, E_max_prompt))
array_numbers_E_prompt_IBD = np.loadtxt(string_numbers_E_prompt_IBD)
number_events_total_IBD = float(array_numbers_E_prompt_IBD[0])
number_analyzed_E_prompt_IBD = float(array_numbers_E_prompt_IBD[2])
N_pass_real_E_prompt_IBD = float(array_numbers_E_prompt_IBD[3])
N_rej_real_E_prompt_IBD = float(array_numbers_E_prompt_IBD[4])
# check the numbers:
if N_pass_real_E_prompt_IBD + N_rej_real_E_prompt_IBD != number_analyzed_E_prompt_IBD:
    sys.exit("N_pass_real_E_prompt_IBD + N_rej_real_E_prompt_IBD != number_analyzed_E_prompt_IBD")
# calculate efficiency for real case (with numbers from conversion):
eff_real_E_prompt_cut_IBD = N_pass_real_E_prompt_IBD / number_analyzed_E_prompt_IBD
# calculate statistical error of eff_real:
eff_real_stat_E_prompt_cut_IBD = np.sqrt(number_analyzed_E_prompt_IBD) / number_analyzed_E_prompt_IBD
# get alpha of prompt energy cut:
alpha_E_prompt_cut_IBD = alpha_E_prompt_cut
alpha_stat_E_prompt_cut_IBD = alpha_stat_E_prompt_cut
# calculate efficiency for ideal case (without conversion):
eff_ideal_E_prompt_cut_IBD = alpha_E_prompt_cut_IBD * eff_real_E_prompt_cut_IBD
# calculate statistical error of eff_ideal:
eff_ideal_stat_E_prompt_cut_IBD = np.sqrt((alpha_E_prompt_cut_IBD * eff_real_stat_E_prompt_cut_IBD)**2 +
                                          (eff_real_E_prompt_cut_IBD * alpha_stat_E_prompt_cut_IBD)**2)

# time cut between prompt und delayed signal:
string_numbers_time_IBD = (input_path_IBD_del + "numbers_time_cut_IBD_{0:.0f}ns_to_{1:.0f}ms_0.txt"
                           .format(time_cut_min, time_cut_max/1000000.0))
(eff_real_time_cut_IBD, eff_real_stat_time_cut_IBD, eff_ideal_time_cut_IBD, eff_ideal_stat_time_cut_IBD,
 alpha_time_cut_IBD, alpha_stat_time_cut_IBD, number_events_total_IBD, number_analyzed_time_cut_IBD) = \
    read_numbers_txt_files(string_numbers_time_IBD)

# multiplicity cut:
string_numbers_mult_IBD = input_path_IBD_del + "numbers_multiplicity_cut_IBD_mult{0:d}_0.txt".format(multiplicity)
(eff_real_mult_cut_IBD, eff_real_stat_mult_cut_IBD, eff_ideal_mult_cut_IBD, eff_ideal_stat_mult_cut_IBD,
 alpha_mult_cut_IBD, alpha_stat_mult_cut_IBD, number_events_total_IBD, number_analyzed_mult_cut_IBD) = \
    read_numbers_txt_files(string_numbers_mult_IBD)

# delayed energy cut:
string_numbers_E_delayed_IBD = (input_path_IBD_del + "numbers_delayed_energy_cut_IBD_{0:.0f}keV_to_{1:.0f}keV_0.txt"
                                .format(E_min_delayed_MeV*1000, E_max_delayed_MeV*1000))
array_numbers_E_delayed_IBD = np.loadtxt(string_numbers_E_delayed_IBD)
number_events_total_IBD = array_numbers_E_delayed_IBD[0]
number_analyzed_E_delayed_IBD = array_numbers_E_delayed_IBD[2]
N_pass_real_E_delayed_IBD = float(array_numbers_E_delayed_IBD[3])
N_rej_real_E_delayed_IBD = float(array_numbers_E_delayed_IBD[4])
# check the numbers:
if N_pass_real_E_delayed_IBD + N_rej_real_E_delayed_IBD != number_analyzed_E_delayed_IBD:
    sys.exit("N_pass_real_E_prompt_IBD + N_rej_real_E_prompt_IBD != number_analyzed_E_prompt_IBD")
# calculate efficiency for real case (with numbers from conversion):
eff_real_E_delayed_cut_IBD = N_pass_real_E_delayed_IBD / number_analyzed_E_delayed_IBD
# calculate statistical error of eff_real:
eff_real_stat_E_delayed_cut_IBD = np.sqrt(number_analyzed_E_delayed_IBD) / number_analyzed_E_delayed_IBD
# get alpha of delayed energy cut:
alpha_E_delayed_cut_IBD = alpha_E_prompt_cut
alpha_stat_E_delayed_cut_IBD = alpha_stat_E_prompt_cut
# calculate efficiency for ideal case (without conversion):
eff_ideal_E_delayed_cut_IBD = alpha_E_delayed_cut_IBD * eff_real_E_delayed_cut_IBD
# calculate statistical error of eff_ideal:
eff_ideal_stat_E_delayed_cut_IBD = np.sqrt((alpha_E_delayed_cut_IBD * eff_real_stat_E_delayed_cut_IBD)**2 +
                                           (eff_real_E_delayed_cut_IBD * alpha_stat_E_delayed_cut_IBD)**2)

# distance cut between prompt and delayed signal:
string_numbers_dist_IBD = input_path_IBD_del + "numbers_distance_cut_IBD_{0:.0f}mm_0.txt".format(distance_cut)
(eff_real_dist_cut_IBD, eff_real_stat_dist_cut_IBD, eff_ideal_dist_cut_IBD, eff_ideal_stat_dist_cut_IBD,
 alpha_dist_cut_IBD, alpha_stat_dist_cut_IBD, number_events_total_IBD, number_analyzed_dist_cut_IBD) = \
    read_numbers_txt_files(string_numbers_dist_IBD)

# volume cut on delayed signal:
string_numbers_vol_delayed_IBD = (input_path_IBD_del + "numbers_delayed_volume_cut_IBD_{0:.0f}mm_0.txt"
                                  .format(R_cut_delayed))
(eff_real_vol_delayed_IBD, eff_real_stat_vol_delayed_IBD, eff_ideal_vol_delayed_IBD, eff_ideal_stat_vol_delayed_IBD,
 alpha_vol_delayed_IBD, alpha_stat_vol_delayed_IBD, number_events_total_IBD, number_analyzed_vol_delayed_IBD) = \
    read_numbers_txt_files(string_numbers_vol_delayed_IBD)

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

""" write information about efficiencies and alphas for NC and IBD events of the single cuts to txt file: """
# set event type NC:
event_type_NC = "atmoNC"

# write info about volume cut of NC events to file:
vol_cut_NC = "volume_cut_R{0:.0f}mm".format(R_cut_prompt)
write_info_to_file(input_path_NC, event_type_NC, vol_cut_NC, eff_real_volume_cut_NC, eff_real_stat_volume_cut_NC,
                   eff_ideal_volume_cut_NC, eff_ideal_stat_volume_cut_NC, alpha_volume_cut_NC, alpha_stat_volume_cut_NC,
                   number_analyzed_volume_cut_NC)

# write info about prompt energy cut of NC events to file:
E_prompt_cut_NC = "prompt_energy_cut_{0:.0f}MeV_to_{1:.0f}MeV".format(E_min_prompt, E_max_prompt)
write_info_to_file(input_path_NC_del, event_type_NC, E_prompt_cut_NC, eff_real_E_prompt_cut_NC,
                   eff_real_stat_E_prompt_cut_NC, eff_ideal_E_prompt_cut_NC, eff_ideal_stat_E_prompt_cut_NC,
                   alpha_E_prompt_cut_NC, alpha_stat_E_prompt_cut_NC, number_analyzed_E_prompt_NC)

# write info about total delayed cut of NC events to file:
delayed_cut_NC = "delayed_cut_{0:.0f}ns_to_{1:.0f}ms_mult{2:d}_{3:.0f}keV_{4:.0f}keV_dist{5:.0f}mm_" \
                 "R{6:.0f}mm".format(time_cut_min, time_cut_max/1000000.0, multiplicity, E_min_delayed_MeV*1000,
                                     E_max_delayed_MeV*1000, distance_cut, R_cut_delayed)
write_info_to_file(input_path_NC_del, event_type_NC, delayed_cut_NC, eff_real_delayed_cut_NC,
                   eff_real_stat_delayed_cut_NC, eff_ideal_delayed_cut_NC, eff_ideal_stat_delayed_cut_NC,
                   alpha_delayed_cut_NC, alpha_stat_delayed_cut_NC, number_analyzed_time_cut_NC)

# write info about time cut of NC events to file:
time_cut_NC = "time_cut_{0:.0f}ns_to_{1:.0f}ms".format(time_cut_min, time_cut_max/1000000.0)
write_info_to_file(input_path_NC_del, event_type_NC, time_cut_NC, eff_real_time_cut_NC, eff_real_stat_time_cut_NC,
                   eff_ideal_time_cut_NC, eff_ideal_stat_time_cut_NC, alpha_time_cut_NC, alpha_stat_time_cut_NC,
                   number_analyzed_time_cut_NC)

# write info about multiplicity cut of NC events to file:
multiplicity_cut_NC = "mult{0:d}".format(multiplicity)
write_info_to_file(input_path_NC_del, event_type_NC, multiplicity_cut_NC, eff_real_mult_cut_NC,
                   eff_real_stat_mult_cut_NC, eff_ideal_mult_cut_NC, eff_ideal_stat_mult_cut_NC, alpha_mult_cut_NC,
                   alpha_stat_mult_cut_NC, number_analyzed_mult_cut_NC)

# write info about delayed energy cut of NC events to file:
delayed_energy_cut_NC = "delayed_energy_cut_{0:.0f}keV_{1:.0f}keV".format(E_min_delayed_MeV*1000,
                                                                          E_max_delayed_MeV*1000)
write_info_to_file(input_path_NC_del, event_type_NC, delayed_energy_cut_NC, eff_real_E_delayed_cut_NC,
                   eff_real_stat_E_delayed_cut_NC, eff_ideal_E_delayed_cut_NC, eff_ideal_stat_E_delayed_cut_NC,
                   alpha_E_delayed_cut_NC, alpha_stat_E_delayed_cut_NC, number_analyzed_E_delayed_NC)

# write info about distance cut of NC events to file:
distance_cut_NC = "distance_cut_dist{0:.0f}mm".format(distance_cut)
write_info_to_file(input_path_NC_del, event_type_NC, distance_cut_NC, eff_real_dist_cut_NC, eff_real_stat_dist_cut_NC,
                   eff_ideal_dist_cut_NC, eff_ideal_stat_dist_cut_NC, alpha_dist_cut_NC, alpha_stat_dist_cut_NC,
                   number_analyzed_dist_cut_NC)

# write info about delayed volume cut of NC events to file:
delayed_volume_cut_NC = "delayed_volume_cut_R{0:.0f}mm".format(R_cut_delayed)
write_info_to_file(input_path_NC_del, event_type_NC, delayed_volume_cut_NC, eff_real_vol_delayed_NC,
                   eff_real_stat_vol_delayed_NC, eff_ideal_vol_delayed_NC, eff_ideal_stat_vol_delayed_NC,
                   alpha_vol_delayed_NC, alpha_stat_vol_delayed_NC, number_analyzed_vol_delayed_NC)

# set event type IBD:
event_type_IBD = "IBD"

# write info about volume cut of IBD events to file:
vol_cut_IBD = "volume_cut_R{0:.0f}mm".format(R_cut_prompt)
write_info_to_file(input_path_IBD, event_type_IBD, vol_cut_IBD, eff_real_volume_cut_IBD, eff_real_stat_volume_cut_IBD,
                   eff_ideal_volume_cut_IBD, eff_ideal_stat_volume_cut_IBD, alpha_volume_cut_IBD,
                   alpha_stat_volume_cut_IBD, number_analyzed_volume_cut_IBD)

# write info about prompt energy cut of IBD events to file:
E_prompt_cut_IBD = "prompt_energy_cut_{0:.0f}MeV_to_{1:.0f}MeV".format(E_min_prompt, E_max_prompt)
write_info_to_file(input_path_IBD_del, event_type_IBD, E_prompt_cut_IBD, eff_real_E_prompt_cut_IBD,
                   eff_real_stat_E_prompt_cut_IBD, eff_ideal_E_prompt_cut_IBD, eff_ideal_stat_E_prompt_cut_IBD,
                   alpha_E_prompt_cut_IBD, alpha_stat_E_prompt_cut_IBD, number_analyzed_E_prompt_IBD)

# write info about total delayed cut of IBD events to file:
delayed_cut_IBD = "delayed_cut_{0:.0f}ns_to_{1:.0f}ms_mult{2:d}_{3:.0f}keV_{4:.0f}keV_dist{5:.0f}mm_" \
                  "R{6:.0f}mm".format(time_cut_min, time_cut_max/1000000.0, multiplicity, E_min_delayed_MeV*1000,
                                      E_max_delayed_MeV*1000, distance_cut, R_cut_delayed)
write_info_to_file(input_path_IBD_del, event_type_IBD, delayed_cut_IBD, eff_real_delayed_cut_IBD,
                   eff_real_stat_delayed_cut_IBD, eff_ideal_delayed_cut_IBD, eff_ideal_stat_delayed_cut_IBD,
                   alpha_delayed_cut_IBD, alpha_stat_delayed_cut_IBD, number_analyzed_time_cut_IBD)

# write info about time cut of IBD events to file:
time_cut_IBD = "time_cut_{0:.0f}ns_to_{1:.0f}ms".format(time_cut_min, time_cut_max/1000000.0)
write_info_to_file(input_path_IBD_del, event_type_IBD, time_cut_IBD, eff_real_time_cut_IBD, eff_real_stat_time_cut_IBD,
                   eff_ideal_time_cut_IBD, eff_ideal_stat_time_cut_IBD, alpha_time_cut_IBD, alpha_stat_time_cut_IBD,
                   number_analyzed_time_cut_IBD)

# write info about multiplicity cut of IBD events to file:
multiplicity_cut_IBD = "mult{0:d}".format(multiplicity)
write_info_to_file(input_path_IBD_del, event_type_IBD, multiplicity_cut_IBD, eff_real_mult_cut_IBD,
                   eff_real_stat_mult_cut_IBD, eff_ideal_mult_cut_IBD, eff_ideal_stat_mult_cut_IBD, alpha_mult_cut_IBD,
                   alpha_stat_mult_cut_IBD, number_analyzed_mult_cut_IBD)

# write info about delayed energy cut of IBD events to file:
delayed_energy_cut_IBD = "delayed_energy_cut_{0:.0f}keV_{1:.0f}keV".format(E_min_delayed_MeV*1000,
                                                                           E_max_delayed_MeV*1000)
write_info_to_file(input_path_IBD_del, event_type_IBD, delayed_energy_cut_IBD, eff_real_E_delayed_cut_IBD,
                   eff_real_stat_E_delayed_cut_IBD, eff_ideal_E_delayed_cut_IBD, eff_ideal_stat_E_delayed_cut_IBD,
                   alpha_E_delayed_cut_IBD, alpha_stat_E_delayed_cut_IBD, number_analyzed_E_delayed_IBD)

# write info about distance cut of IBD events to file:
distance_cut_IBD = "distance_cut_dist{0:.0f}mm".format(distance_cut)
write_info_to_file(input_path_IBD_del, event_type_IBD, distance_cut_IBD, eff_real_dist_cut_IBD,
                   eff_real_stat_dist_cut_IBD, eff_ideal_dist_cut_IBD, eff_ideal_stat_dist_cut_IBD, alpha_dist_cut_IBD,
                   alpha_stat_dist_cut_IBD, number_analyzed_dist_cut_IBD)

# write info about delayed volume cut of IBD events to file:
delayed_volume_cut_IBD = "delayed_volume_cut_R{0:.0f}mm".format(R_cut_delayed)
write_info_to_file(input_path_IBD_del, event_type_IBD, delayed_volume_cut_IBD, eff_real_vol_delayed_IBD,
                   eff_real_stat_vol_delayed_IBD, eff_ideal_vol_delayed_IBD, eff_ideal_stat_vol_delayed_IBD,
                   alpha_vol_delayed_IBD, alpha_stat_vol_delayed_IBD, number_analyzed_vol_delayed_IBD)

for index10 in range(len(tail_start)):
    """ Analyze NC events (only information from filenumber_evtID_...txt files is needed): """
    # number of IBD-like events from real data before event rate (without PSD suppression and
    # without alpha (MCtruth/real)):
    number_IBDlike_real_simu_wo_PSD_wo_alpha = 0
    # number of IBD-like events from real data before event rate (with PSD suppression, but without
    # alpha (MCtruth/real)):
    number_IBDlike_real_simu_w_PSD_wo_alpha = 0
    # number of IBD-like events from ideal data before event rate (without PSD suppression and
    # without alpha (MCtruth/real)):
    number_IBDlike_ideal_simu_wo_PSD_wo_alpha = 0
    # number of IBD-like events from ideal data before event rate (with PSD suppression, but without alpha
    # (MCtruth/real)):
    number_IBDlike_ideal_simu_w_PSD_wo_alpha = 0

    # preallocate array, where Evis of real (reconstructed) IBD-like events (without PSD suppression and without alpha)
    # is saved in MeV:
    Evis_array_real_wo_PSD_wo_alpha = []
    # preallocate array, where Evis of real (reconstructed) IBD-like events (with PSD suppression, but without alpha)
    # is saved in MeV:
    Evis_array_real_w_PSD_wo_alpha = []

    # The following arrays save filenumber and evtID of events that pass the cuts to check interaction channels,
    # deexcitation channels, ... of IBD-like events afterwards.
    # preallocate array, where filenumber of the events, that pass volume, prompt energy and delayed cut, are stored:
    filenumber_pass_vol_E_del = []
    # preallocate array, where evtID of the events, that pass volume, prompt energy and delayed cut, are stored:
    evtID_pass_vol_E_del = []
    # preallocate array, where filenumber of the events, that pass volume, prompt energy, delayed and PSD cut,
    # are stored:
    filenumber_pass_vol_E_del_PSD = []
    # preallocate array, where evtID of the events, that pass volume, prompt energy, delayed and PSD cut, are stored:
    evtID_pass_vol_E_del_PSD = []
    # preallocate array, where TTR values of NC events, that pass volume, prompt energy and delayed cut (before PSD),
    # are stored (this means the TTR value of IBD-like NC events):
    array_TTR_IBDlike = []

    # load files, where filenumber and evtID (and Evis) of events that pass the cut are stored (real/reconstructed
    # data):
    array_pass_volume_cut_real_NC = np.loadtxt(input_path_NC + "filenumber_evtID_volume_cut_atmoNC_{0:.0f}mm.txt"
                                               .format(R_cut_prompt))
    array_pass_E_prompt_cut_real_NC = np.loadtxt(input_path_NC_del + "filenumber_evtID_Evis_prompt_energy_cut_atmoNC_"
                                                                     "{0:.0f}MeV_to_{1:.0f}MeV_0.txt"
                                                 .format(E_min_prompt, E_max_prompt))
    array_pass_delayed_cut_real_NC = np.loadtxt(input_path_NC_del + "filenumber_evtID_delayed_cut_atmoNC_{0:.0f}ns_to_"
                                                                    "{1:.0f}ms_mult{2:d}_{3:.0f}keV_{4:.0f}keV_"
                                                                    "dist{5:.0f}mm_R{6:.0f}mm_0.txt"
                                                .format(time_cut_min, time_cut_max/1000000.0, multiplicity,
                                                        E_min_delayed_MeV*1000, E_max_delayed_MeV*1000, distance_cut,
                                                        R_cut_delayed))
    # load files, where filenumber and evtID (and Evis) of events that pass the cut are stored (ideal/MC truth data):
    array_pass_volume_cut_ideal_NC = np.loadtxt(input_path_NC + "filenumber_evtID_volume_cut_MCtruth_atmoNC_"
                                                                "{0:.0f}mm.txt".format(R_cut_prompt))
    # INFO-me: use filenumber and evtID of real data since no information about MCtruth data
    # INFO-me: alpha of prompt energy cut must be considered afterwards!
    array_pass_E_prompt_cut_ideal_NC = np.loadtxt(input_path_NC_del + "filenumber_evtID_Evis_prompt_energy_cut_"
                                                                      "atmoNC_{0:.0f}MeV_to_{1:.0f}MeV_0.txt"
                                                  .format(E_min_prompt, E_max_prompt))
    array_pass_delayed_cut_ideal_NC = np.loadtxt(input_path_NC_del + "filenumber_evtID_delayed_cut_MCtruth_atmoNC_"
                                                                     "{0:.0f}ns_to_{1:.0f}ms_mult{2:d}_{3:.0f}keV_"
                                                                     "{4:.0f}keV_dist{5:.0f}mm_R{6:.0f}mm_0.txt"
                                                 .format(time_cut_min, time_cut_max/1000000.0, multiplicity,
                                                         E_min_delayed_MeV*1000, E_max_delayed_MeV*1000, distance_cut,
                                                         R_cut_delayed))
    # load file, where filenumber, evtID and TTR values of all NC events are stored:
    array_TTR_NC = np.loadtxt(input_path_PSD + "TTR_atmoNC_{0:.0f}ns_{1:.0f}ns_0.txt".format(tail_start[index10],
                                                                                             tail_end[index10]))

    # preallocate start indices of the arrays:
    index_volume = 0
    index_E_prompt = 0
    index_TTR = 0

    # check real NC events:
    # loop over array_pass_delayed_cut_real_NC:
    for index in range(len(array_pass_delayed_cut_real_NC)):

        # get filenumber and evtID of the event that passed delayed cut:
        filenumber_delayed_cut = array_pass_delayed_cut_real_NC[index][0]
        evtID_delayed_cut = array_pass_delayed_cut_real_NC[index][1]

        if filenumber_delayed_cut > last_file_NC:
            continue

        # loop over array_pass_volume_cut_real_NC:
        for index1 in range(index_volume, len(array_pass_volume_cut_real_NC), 1):

            # get filenumber and evtID of the event that passed volume cut:
            filenumber_volume_cut = array_pass_volume_cut_real_NC[index1][0]
            evtID_volume_cut = array_pass_volume_cut_real_NC[index1][1]

            if filenumber_volume_cut > last_file_NC:
                continue

            # check if event also passed delayed cut:
            if filenumber_volume_cut == filenumber_delayed_cut and evtID_volume_cut == evtID_delayed_cut:

                # event passed volume cut AND delayed cut:

                # loop over array_pass_E_prompt_cut_real_NC:
                for index2 in range(index_E_prompt, len(array_pass_E_prompt_cut_real_NC), 1):

                    # get filenumber, evtID and Evis of the event that passed prompt energy cut:
                    filenumber_E_cut = array_pass_E_prompt_cut_real_NC[index2][0]
                    evtID_E_cut = array_pass_E_prompt_cut_real_NC[index2][1]
                    E_vis = array_pass_E_prompt_cut_real_NC[index2][2]

                    if filenumber_E_cut > last_file_NC:
                        continue

                    # check if event also passed volume and delayed cut:
                    if filenumber_E_cut == filenumber_volume_cut and evtID_E_cut == evtID_volume_cut:

                        # event passed delayed, volume and prompt energy cut -> IBD-like event!

                        # increment number_IBDlike_real_simu_wo_PSD_wo_alpha:
                        number_IBDlike_real_simu_wo_PSD_wo_alpha += 1
                        # append E_vis to Evis_array_real_wo_PSD_wo_alpha:
                        Evis_array_real_wo_PSD_wo_alpha.append(E_vis)

                        # append filenumber and evtID of event that pass the cuts to array:
                        filenumber_pass_vol_E_del.append(filenumber_E_cut)
                        evtID_pass_vol_E_del.append(evtID_E_cut)

                        # get TTR value of the IBD-like event corresponding to this filenumber and event:
                        # loop over array_TTR_NC:
                        for index4 in range(index_TTR, len(array_TTR_NC), 1):
                            filenumber_TTR = array_TTR_NC[index4][0]
                            evtID_TTR = array_TTR_NC[index4][1]
                            if filenumber_TTR == filenumber_E_cut and evtID_TTR == evtID_E_cut:
                                # get TTR value of this event:
                                TTR_value_NC = array_TTR_NC[index4][2]
                                # store TTR value to array_TTR_IBDlike:
                                array_TTR_IBDlike.append(TTR_value_NC)
                                # set index_TTR = index4 -> start loop for next event at index_TTR:
                                index_TTR = index4
                                break
                            else:
                                continue

                        # Check if event also passes PSD cut:
                        if TTR_value_NC <= TTR_cut_value[index10]:
                            # TTR value of the event is smaller than TTR cut value:
                            # -> event passes PSD cut!

                            # increment number_IBDlike_real_simu_w_PSD_wo_alpha:
                            number_IBDlike_real_simu_w_PSD_wo_alpha += 1
                            # append E_vis to Evis_array_real_w_PSD_wo_alpha:
                            Evis_array_real_w_PSD_wo_alpha.append(E_vis)

                            # append filenumber and evtID of event that pass the cuts to array:
                            filenumber_pass_vol_E_del_PSD.append(filenumber_E_cut)
                            evtID_pass_vol_E_del_PSD.append(evtID_E_cut)

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

    # preallocate start indices of the arrays:
    index_volume = 0
    index_E_prompt = 0
    index_TTR = 0

    # check ideal NC events:
    # loop over array_pass_delayed_cut_ideal_NC:
    for index in range(len(array_pass_delayed_cut_ideal_NC)):

        # get filenumber and evtID of the event that passed delayed cut:
        filenumber_delayed_cut = array_pass_delayed_cut_ideal_NC[index][0]
        evtID_delayed_cut = array_pass_delayed_cut_ideal_NC[index][1]

        if filenumber_delayed_cut > last_file_NC:
            continue

        # loop over array_pass_volume_cut_ideal_NC:
        for index1 in range(index_volume, len(array_pass_volume_cut_ideal_NC), 1):

            # get filenumber and evtID of the event that passed volume cut:
            filenumber_volume_cut = array_pass_volume_cut_ideal_NC[index1][0]
            evtID_volume_cut = array_pass_volume_cut_ideal_NC[index1][1]

            if filenumber_volume_cut > last_file_NC:
                continue

            # check if event also passed delayed cut:
            if filenumber_volume_cut == filenumber_delayed_cut and evtID_volume_cut == evtID_delayed_cut:

                # event passed volume cut AND delayed cut:

                # loop over array_pass_E_prompt_cut_ideal_NC:
                for index2 in range(index_E_prompt, len(array_pass_E_prompt_cut_ideal_NC), 1):

                    # get filenumber, evtID and Evis of the event that passed prompt energy cut:
                    filenumber_E_cut = array_pass_E_prompt_cut_ideal_NC[index2][0]
                    evtID_E_cut = array_pass_E_prompt_cut_ideal_NC[index2][1]
                    E_vis = array_pass_E_prompt_cut_ideal_NC[index2][2]

                    if filenumber_E_cut > last_file_NC:
                        continue

                    # check if event also passed volume and delayed cut:
                    if filenumber_E_cut == filenumber_volume_cut and evtID_E_cut == evtID_volume_cut:

                        # event passed delayed, volume and prompt energy cut -> IBD-like event

                        # increment number_IBDlike_ideal_simu_wo_PSD_wo_alpha:
                        number_IBDlike_ideal_simu_wo_PSD_wo_alpha += 1

                        # get TTR value of the IBD-like event corresponding to this filenumber and event:
                        # loop over array_TTR_NC:
                        for index4 in range(index_TTR, len(array_TTR_NC), 1):
                            filenumber_TTR = array_TTR_NC[index4][0]
                            evtID_TTR = array_TTR_NC[index4][1]
                            if filenumber_TTR == filenumber_E_cut and evtID_TTR == evtID_E_cut:
                                # get TTR value of this event:
                                TTR_value_NC = array_TTR_NC[index4][2]
                                # set index_TTR = index4 -> start loop for next event at index_TTR:
                                index_TTR = index4
                                break
                            else:
                                continue

                        # check if event also passes PSD cut:
                        if TTR_value_NC <= TTR_cut_value[index10]:
                            # TTR value of the event is smaller than TTR cut value:
                            # -> event passes PSD cut!

                            # increment number_IBDlike_real_simu_w_PSD_wo_alpha:
                            number_IBDlike_ideal_simu_w_PSD_wo_alpha += 1

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

    """ calculate the total NC efficiencies: """
    # total efficiency of real NC data (number IBD-like events / number total NC events):
    total_efficiency_NC_real_wo_PSD = float(number_IBDlike_real_simu_wo_PSD_wo_alpha) / float(number_events_total_NC)
    # statistical error of total_efficiency_NC_real_wo_PSD:
    stat_total_efficiency_NC_real_wo_PSD = np.sqrt(float(number_events_total_NC)) / float(number_events_total_NC)
    # total efficiency of ideal NC data without considering uncertainty of prompt energy cut:
    total_efficiency_NC_ideal_wo_PSD_1 = (float(number_IBDlike_ideal_simu_wo_PSD_wo_alpha) /
                                          float(number_events_total_NC))
    # statistical error of total_efficiency_NC_ideal_wo_PSD_1:
    stat_total_efficiency_NC_ideal_wo_PSD_1 = np.sqrt(float(number_events_total_NC)) / float(number_events_total_NC)
    # total efficiency of ideal NC data:
    # INFO-me: consider here alpha of prompt energy cut!
    total_efficiency_NC_ideal_wo_PSD = alpha_E_prompt_cut_NC * total_efficiency_NC_ideal_wo_PSD_1
    # statistical error of total_efficiency_NC_ideal_wo_PSD:
    stat_total_efficiency_NC_ideal_wo_PSD = np.sqrt((total_efficiency_NC_ideal_wo_PSD_1 * alpha_stat_E_prompt_cut_NC)**2
                                                    + (alpha_E_prompt_cut_NC *
                                                       stat_total_efficiency_NC_ideal_wo_PSD_1)**2)
    # total alpha for NC events ('error of efficiency'):
    total_alpha_NC_wo_PSD = total_efficiency_NC_ideal_wo_PSD / total_efficiency_NC_real_wo_PSD
    # statistical error of total_alpha_wo_PSD:
    stat_total_alpha_NC_wo_PSD = np.sqrt((stat_total_efficiency_NC_real_wo_PSD / total_efficiency_NC_ideal_wo_PSD) ** 2
                                         +
                                         (total_efficiency_NC_real_wo_PSD * stat_total_efficiency_NC_ideal_wo_PSD /
                                          total_efficiency_NC_ideal_wo_PSD) ** 2)

    # total efficiency of real NC data with PSD:
    total_efficiency_NC_real_w_PSD = float(number_IBDlike_real_simu_w_PSD_wo_alpha) / float(number_events_total_NC)
    # statistical error of total_efficiency_NC_real_w_PSD:
    stat_total_efficiency_NC_real_w_PSD = np.sqrt(float(number_events_total_NC)) / float(number_events_total_NC)
    # total efficiency of ideal NC data without considering uncertainty of prompt energy cut:
    total_efficiency_NC_ideal_w_PSD_1 = float(number_IBDlike_ideal_simu_w_PSD_wo_alpha) / float(number_events_total_NC)
    # statistical error of total_efficiency_NC_ideal_w_PSD_1:
    stat_total_efficiency_NC_ideal_w_PSD_1 = np.sqrt(float(number_events_total_NC)) / float(number_events_total_NC)
    # total efficiency of ideal NC data with PSD:
    # INFO-me: consider here alpha of prompt energy cut!
    total_efficiency_NC_ideal_w_PSD = alpha_E_prompt_cut_NC * total_efficiency_NC_ideal_w_PSD_1
    # statistical error of total_efficiency_NC_ideal_w_PSD:
    stat_total_efficiency_NC_ideal_w_PSD = np.sqrt((total_efficiency_NC_ideal_w_PSD_1 * alpha_stat_E_prompt_cut_NC)**2
                                                   + (alpha_E_prompt_cut_NC *
                                                      stat_total_efficiency_NC_ideal_w_PSD_1)**2)
    # total alpha for NC events with PSD:
    total_alpha_NC_w_PSD = total_efficiency_NC_ideal_w_PSD / total_efficiency_NC_real_w_PSD
    # statistical error of total_alpha_w_PSD:
    stat_total_alpha_NC_w_PSD = np.sqrt((stat_total_efficiency_NC_real_w_PSD / total_efficiency_NC_ideal_w_PSD) ** 2 +
                                        (total_efficiency_NC_real_w_PSD * stat_total_efficiency_NC_ideal_w_PSD /
                                         total_efficiency_NC_ideal_w_PSD) ** 2)

    """ Build histograms from E_vis_arrays: """
    # set bin-edges of e_vis histogram in MeV:
    bins_evis = np.arange(E_min_prompt, E_max_prompt + 2*bin_width_energy, bin_width_energy)

    # build histogram from Evis_array_real_wo_PSD_wo_alpha (without PSD suppression and alpha):
    histo_Evis_wo_PSD_wo_alpha, bin_edges_evis = np.histogram(Evis_array_real_wo_PSD_wo_alpha, bins_evis)

    # build histogram from Evis_array_real_w_PSD_wo_alpha (with PSD suppression, but without alpha (IBD cuts and PSD)):
    histo_Evis_w_PSD_wo_alpha, bin_edges_evis = np.histogram(Evis_array_real_w_PSD_wo_alpha, bins_evis)

    """ Event rate calculation: """
    # calculate the theoretical event rate in NC events/sec in JUNO for neutrino energies from 0 MeV to 10 GeV (float)
    # (event_rate = A * (flux_nue*xsec_nue + flux_nuebar*xsec_nuebar + flux_numu*xsec_numu +
    # flux_numubar*xsec_numubar)):
    event_rate = NC_background_functions.event_rate(bin_width_energy, r_cut, output_path, PLOT_FLUX, SHOW_FLUXPLOT,
                                                    SAVE_FLUXPLOT, PLOT_EVT_RATE, SHOW_EVT_RATE, SAVE_EVT_RATE)

    # number of NC events in JUNO after 10 years:
    number_NC_events_JUNO = event_rate * time_seconds

    # number of IBD-like events in JUNO after 10 years without PSD and without alpha:
    number_IBDlike_JUNO_wo_PSD = int(number_NC_events_JUNO * number_IBDlike_real_simu_wo_PSD_wo_alpha /
                                     number_events_total_NC)
    # number of IBD-like events in JUNO after 10 years with PSD , but without alpha:
    number_IBDlike_JUNO_w_PSD = int(number_NC_events_JUNO * number_IBDlike_real_simu_w_PSD_wo_alpha /
                                    number_events_total_NC)

    # normalize the spectrum of IBD-like events to the spectrum, JUNO will measure after 10 years (without PSD and
    # without alpha):
    Evis_histo_JUNO_wo_PSD = (float(number_IBDlike_JUNO_wo_PSD) / float(number_IBDlike_real_simu_wo_PSD_wo_alpha) *
                              histo_Evis_wo_PSD_wo_alpha)
    # normalize the spectrum of IBD-like events to the spectrum, JUNO will measure after 10 years (with PSD, but without
    # alpha):
    Evis_histo_JUNO_w_PSD = (float(number_IBDlike_JUNO_w_PSD) / float(number_IBDlike_real_simu_w_PSD_wo_alpha) *
                             histo_Evis_w_PSD_wo_alpha)

    """ Analyze IBD events (only information from filenumber_evtID_...txt files is needed): """
    # number of IBD events from real data before event rate (without PSD suppression and without alpha (MCtruth/real)):
    number_IBD_real_simu_wo_PSD_wo_alpha = 0
    # number of IBD events from real data before event rate (with PSD suppression, but without alpha (MCtruth/real)):
    number_IBD_real_simu_w_PSD_wo_alpha = 0
    # number of IBD-like events from ideal data before event rate (without PSD suppression and
    # without alpha (MCtruth/real)):
    number_IBD_ideal_simu_wo_PSD_wo_alpha = 0
    # number of IBD-like events from ideal data before event rate (with PSD suppression, but without alpha
    # (MCtruth/real)):
    number_IBD_ideal_simu_w_PSD_wo_alpha = 0

    # preallocate array, where Evis of real (reconstructed) IBD events that pass all cuts
    # (without PSD suppression and without alpha) is saved in MeV:
    Evis_array_real_wo_PSD_wo_alpha_IBD = []
    # preallocate array, where Evis of real (reconstructed) IBD events that pass all cuts
    # (with PSD suppression, but without alpha) is saved in MeV:
    Evis_array_real_w_PSD_wo_alpha_IBD = []

    # The following arrays save filenumber and evtID of events that pass the cuts to check interaction channels,
    # deexcitation channels, ... of IBD-like events afterwards.
    # preallocate array, where filenumber of the events, that pass volume, prompt energy and delayed cut, are stored:
    filenumber_pass_vol_E_del_IBD = []
    # preallocate array, where evtID of the events, that pass volume, prompt energy and delayed cut, are stored:
    evtID_pass_vol_E_del_IBD = []
    # preallocate array, where filenumber of the events, that pass volume, prompt energy, delayed and PSD cut,
    # are stored:
    filenumber_pass_vol_E_del_PSD_IBD = []
    # preallocate array, where evtID of the events, that pass volume, prompt energy, delayed and PSD cut, are stored:
    evtID_pass_vol_E_del_PSD_IBD = []
    # preallocate array, where TTR values of IBD events, that pass volume, prompt energy and delayed cut
    # (before PSD), are
    # stored:
    array_TTR_IBD_beforePSD = []

    # load files, where filenumber and evtID of events that pass the cut are stored (real/reconstructed data):
    array_pass_volume_cut_real_IBD = np.loadtxt(input_path_IBD + "filenumber_evtID_volume_cut_IBD_{0:.0f}mm.txt"
                                                .format(R_cut_prompt))
    array_pass_E_prompt_cut_real_IBD = np.loadtxt(input_path_IBD_del + "filenumber_evtID_Evis_prompt_energy_cut_IBD_"
                                                                       "{0:.0f}MeV_to_{1:.0f}MeV_0.txt"
                                                  .format(E_min_prompt, E_max_prompt))
    array_pass_delayed_cut_real_IBD = np.loadtxt(input_path_IBD_del + "filenumber_evtID_delayed_cut_IBD_{0:.0f}ns_to_"
                                                                      "{1:.0f}ms_mult{2:d}_{3:.0f}keV_{4:.0f}keV_"
                                                                      "dist{5:.0f}mm_R{6:.0f}mm_0.txt"
                                                 .format(time_cut_min, time_cut_max/1000000.0, multiplicity,
                                                         E_min_delayed_MeV*1000, E_max_delayed_MeV*1000, distance_cut,
                                                         R_cut_delayed))
    # load files, where filenumber and evtID (and Evis) of events that pass the cut are stored (ideal/MC truth data):
    array_pass_volume_cut_ideal_IBD = np.loadtxt(input_path_IBD + "filenumber_evtID_volume_cut_MCtruth_IBD_"
                                                                  "{0:.0f}mm.txt".format(R_cut_prompt))
    # INFO-me: use filenumber and evtID of real data since no information about MCtruth data
    # INFO-me: alpha of prompt energy cut must be considered afterwards!
    array_pass_E_prompt_cut_ideal_IBD = np.loadtxt(input_path_IBD_del + "filenumber_evtID_Evis_prompt_energy_cut_"
                                                                        "IBD_{0:.0f}MeV_to_{1:.0f}MeV_0.txt"
                                                   .format(E_min_prompt, E_max_prompt))
    array_pass_delayed_cut_ideal_IBD = np.loadtxt(input_path_IBD_del + "filenumber_evtID_delayed_cut_MCtruth_IBD_"
                                                                       "{0:.0f}ns_"
                                                                       "to_{1:.0f}ms_mult{2:d}_{3:.0f}keV_{4:.0f}keV_"
                                                                       "dist{5:.0f}mm_R{6:.0f}mm_0.txt"
                                                  .format(time_cut_min, time_cut_max/1000000.0, multiplicity,
                                                          E_min_delayed_MeV*1000, E_max_delayed_MeV*1000, distance_cut,
                                                          R_cut_delayed))
    # load file, where filenumber, evtID and TTR values of all IBD events are stored:
    array_TTR_IBD = np.loadtxt(input_path_PSD + "TTR_IBD_{0:.0f}ns_{1:.0f}ns_0.txt".format(tail_start[index10],
                                                                                           tail_end[index10]))

    # preallocate start indices of the arrays:
    index_volume = 0
    index_E_prompt = 0
    index_TTR = 0

    # check real IBD events:
    # loop over array_pass_delayed_cut_real_IBD:
    for index in range(len(array_pass_delayed_cut_real_IBD)):

        # get filenumber and evtID of the event that passed delayed cut:
        filenumber_delayed_cut = array_pass_delayed_cut_real_IBD[index][0]
        evtID_delayed_cut = array_pass_delayed_cut_real_IBD[index][1]

        if filenumber_delayed_cut > last_file_IBD:
            continue

        # loop over array_pass_volume_cut_real_IBD:
        for index1 in range(index_volume, len(array_pass_volume_cut_real_IBD), 1):

            # get filenumber and evtID of the event that passed volume cut:
            filenumber_volume_cut = array_pass_volume_cut_real_IBD[index1][0]
            evtID_volume_cut = array_pass_volume_cut_real_IBD[index1][1]

            if filenumber_volume_cut > last_file_IBD:
                continue

            # check if event also passed delayed cut:
            if filenumber_volume_cut == filenumber_delayed_cut and evtID_volume_cut == evtID_delayed_cut:

                # event passed volume cut AND delayed cut:

                # loop over array_pass_E_prompt_cut_real_IBD:
                for index2 in range(index_E_prompt, len(array_pass_E_prompt_cut_real_IBD), 1):

                    # get filenumber, evtID and Evis of the event that passed prompt energy cut:
                    filenumber_E_cut = array_pass_E_prompt_cut_real_IBD[index2][0]
                    evtID_E_cut = array_pass_E_prompt_cut_real_IBD[index2][1]
                    E_vis = array_pass_E_prompt_cut_real_IBD[index2][2]

                    if filenumber_E_cut > last_file_IBD:
                        continue

                    # check if event also passed volume and delayed cut:
                    if filenumber_E_cut == filenumber_volume_cut and evtID_E_cut == evtID_volume_cut:

                        # event passed delayed, volume and prompt energy cut!

                        # increment number_IBD_real_simu_wo_PSD_wo_alpha:
                        number_IBD_real_simu_wo_PSD_wo_alpha += 1
                        # append E_vis to Evis_array_real_wo_PSD_wo_alpha:
                        Evis_array_real_wo_PSD_wo_alpha_IBD.append(E_vis)

                        # append filenumber and evtID of event that pass the cuts to array:
                        filenumber_pass_vol_E_del_IBD.append(filenumber_E_cut)
                        evtID_pass_vol_E_del_IBD.append(evtID_E_cut)

                        # get TTR value of the IBD event corresponding to this filenumber and event:
                        # loop over array_TTR_IBD:
                        for index4 in range(index_TTR, len(array_TTR_IBD), 1):
                            filenumber_TTR = array_TTR_IBD[index4][0]
                            evtID_TTR = array_TTR_IBD[index4][1]
                            if filenumber_TTR == filenumber_E_cut and evtID_TTR == evtID_E_cut:
                                # get TTR value of this event:
                                TTR_value_IBD = array_TTR_IBD[index4][2]
                                # store TTR value to array_TTR_IBD_beforePSD:
                                array_TTR_IBD_beforePSD.append(TTR_value_IBD)
                                # set index_TTR = index4 -> start loop for next event at index_TTR:
                                index_TTR = index4
                                break
                            else:
                                continue

                        # check if event also pass PSD cut:
                        if TTR_value_IBD <= TTR_cut_value[index10]:
                            # TTR value of the event is smaller than TTR cut value:
                            # -> event passes PSD cut!

                            # increment number_IBD_real_simu_w_PSD_wo_alpha:
                            number_IBD_real_simu_w_PSD_wo_alpha += 1
                            # append E_vis to Evis_array_real_w_PSD_wo_alpha_IBD:
                            Evis_array_real_w_PSD_wo_alpha_IBD.append(E_vis)

                            # append filenumber and evtID of event that pass the cuts to array:
                            filenumber_pass_vol_E_del_PSD_IBD.append(filenumber_E_cut)
                            evtID_pass_vol_E_del_PSD_IBD.append(evtID_E_cut)

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

    # preallocate start indices of the arrays:
    index_volume = 0
    index_E_prompt = 0
    index_TTR = 0

    # check ideal IBD events:
    # loop over array_pass_delayed_cut_ideal_IBD:
    for index in range(len(array_pass_delayed_cut_ideal_IBD)):

        # get filenumber and evtID of the event that passed delayed cut:
        filenumber_delayed_cut = array_pass_delayed_cut_ideal_IBD[index][0]
        evtID_delayed_cut = array_pass_delayed_cut_ideal_IBD[index][1]

        if filenumber_delayed_cut > last_file_IBD:
            continue

        # loop over array_pass_volume_cut_ideal_IBD:
        for index1 in range(index_volume, len(array_pass_volume_cut_ideal_IBD), 1):

            # get filenumber and evtID of the event that passed volume cut:
            filenumber_volume_cut = array_pass_volume_cut_ideal_IBD[index1][0]
            evtID_volume_cut = array_pass_volume_cut_ideal_IBD[index1][1]

            if filenumber_volume_cut > last_file_IBD:
                continue

            # check if event also passed delayed cut:
            if filenumber_volume_cut == filenumber_delayed_cut and evtID_volume_cut == evtID_delayed_cut:

                # event passed volume cut AND delayed cut:

                # loop over array_pass_E_prompt_cut_ideal_IBD:
                for index2 in range(index_E_prompt, len(array_pass_E_prompt_cut_ideal_IBD), 1):

                    # get filenumber, evtID and Evis of the event that passed prompt energy cut:
                    filenumber_E_cut = array_pass_E_prompt_cut_ideal_IBD[index2][0]
                    evtID_E_cut = array_pass_E_prompt_cut_ideal_IBD[index2][1]
                    E_vis = array_pass_E_prompt_cut_ideal_IBD[index2][2]

                    if filenumber_E_cut > last_file_IBD:
                        continue

                    # check if event also passed volume and delayed cut:
                    if filenumber_E_cut == filenumber_volume_cut and evtID_E_cut == evtID_volume_cut:

                        # event passed delayed, volume and prompt energy cut:

                        # increment number_IBD_ideal_simu_wo_PSD_wo_alpha:
                        number_IBD_ideal_simu_wo_PSD_wo_alpha += 1

                        # get TTR value of the IBD event corresponding to this filenumber and event:
                        # loop over array_TTR_IBD:
                        for index4 in range(index_TTR, len(array_TTR_IBD), 1):
                            filenumber_TTR = array_TTR_IBD[index4][0]
                            evtID_TTR = array_TTR_IBD[index4][1]
                            if filenumber_TTR == filenumber_E_cut and evtID_TTR == evtID_E_cut:
                                # get TTR value of this event:
                                TTR_value_IBD = array_TTR_IBD[index4][2]
                                # set index_TTR = index4 -> start loop for next event at index_TTR:
                                index_TTR = index4
                                break
                            else:
                                continue

                        # check if event also pass PSD cut:
                        if TTR_value_IBD <= TTR_cut_value[index10]:
                            # TTR value of the event is smaller than TTR cut value:
                            # -> event passes PSD cut!

                            # increment number_IBD_real_simu_w_PSD_wo_alpha:
                            number_IBD_ideal_simu_w_PSD_wo_alpha += 1

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

    """ calculate the total IBD efficiencies: """
    # total efficiency of real IBD data (number IBD events pass cut / number total IBD events):
    total_efficiency_IBD_real_wo_PSD = float(number_IBD_real_simu_wo_PSD_wo_alpha) / float(number_events_total_IBD)
    # statistical error of total_efficiency_IBD_real_wo_PSD:
    stat_total_efficiency_IBD_real_wo_PSD = np.sqrt(float(number_events_total_IBD)) / float(number_events_total_IBD)
    # total efficiency of ideal IBD data without considering uncertainty of prompt energy cut:
    total_efficiency_IBD_ideal_wo_PSD_1 = float(number_IBD_ideal_simu_wo_PSD_wo_alpha) / float(number_events_total_IBD)
    # statistical error of total_efficiency_IBD_ideal_wo_PSD without considering uncertainty of prompt energy cut:
    stat_total_efficiency_IBD_ideal_wo_PSD_1 = np.sqrt(float(number_events_total_IBD)) / float(number_events_total_IBD)
    # total efficiency of ideal IBD data:
    # INFO-me: consider here alpha of prompt energy cut!
    total_efficiency_IBD_ideal_wo_PSD = alpha_E_prompt_cut_IBD * total_efficiency_IBD_ideal_wo_PSD_1
    # statistical error of total_efficiency_IBD_ideal_wo_PSD:
    stat_total_efficiency_IBD_ideal_wo_PSD = np.sqrt((total_efficiency_IBD_ideal_wo_PSD_1 *
                                                      alpha_stat_E_prompt_cut_IBD) ** 2 +
                                                     (alpha_E_prompt_cut_IBD *
                                                      stat_total_efficiency_IBD_ideal_wo_PSD_1) ** 2)
    # total alpha for IBD events ('error of efficiency'):
    total_alpha_IBD_wo_PSD = total_efficiency_IBD_ideal_wo_PSD / total_efficiency_IBD_real_wo_PSD
    # statistical error of total_alpha_IBD_wo_PSD:
    stat_total_alpha_IBD_wo_PSD = np.sqrt((stat_total_efficiency_IBD_real_wo_PSD /
                                           total_efficiency_IBD_ideal_wo_PSD) ** 2 +
                                          (total_efficiency_IBD_real_wo_PSD * stat_total_efficiency_IBD_ideal_wo_PSD /
                                           total_efficiency_IBD_ideal_wo_PSD) ** 2)

    # total efficiency of real IBD data with PSD:
    total_efficiency_IBD_real_w_PSD = float(number_IBD_real_simu_w_PSD_wo_alpha) / float(number_events_total_IBD)
    # statistical error of total_efficiency_IBD_real_w_PSD:
    stat_total_efficiency_IBD_real_w_PSD = np.sqrt(float(number_events_total_IBD)) / float(number_events_total_IBD)
    # total efficiency of ideal IBD data with PSD without considering uncertainty of prompt energy cut:
    total_efficiency_IBD_ideal_w_PSD_1 = float(number_IBD_ideal_simu_w_PSD_wo_alpha) / float(number_events_total_IBD)
    # statistical error of total_efficiency_IBD_ideal_w_PSD_1 without considering uncertainty of prompt energy cut:
    stat_total_efficiency_IBD_ideal_w_PSD_1 = np.sqrt(float(number_events_total_IBD)) / float(number_events_total_IBD)
    # total efficiency of ideal IBD data with PSD:
    # INFO-me: consider here alpha of prompt energy cut!
    total_efficiency_IBD_ideal_w_PSD = alpha_E_prompt_cut_IBD * total_efficiency_IBD_ideal_w_PSD_1
    # statistical error of total_efficiency_IBD_ideal_w_PSD:
    stat_total_efficiency_IBD_ideal_w_PSD = np.sqrt((total_efficiency_IBD_ideal_w_PSD_1 *
                                                     alpha_stat_E_prompt_cut_IBD) ** 2 +
                                                    (alpha_E_prompt_cut_IBD *
                                                     stat_total_efficiency_IBD_ideal_w_PSD_1) ** 2)
    # total alpha for IBD events with PSD:
    total_alpha_IBD_w_PSD = total_efficiency_IBD_ideal_w_PSD / total_efficiency_IBD_real_w_PSD
    # statistical error of total_alpha_IBD_w_PSD:
    stat_total_alpha_IBD_w_PSD = np.sqrt((stat_total_efficiency_IBD_real_w_PSD /
                                          total_efficiency_IBD_ideal_w_PSD) ** 2 +
                                         (total_efficiency_IBD_real_w_PSD * stat_total_efficiency_IBD_ideal_w_PSD /
                                          total_efficiency_IBD_ideal_w_PSD) ** 2)

    """ Build histograms from E_vis_arrays: """
    # build histogram from Evis_array_real_wo_PSD_wo_alpha_IBD (without PSD suppression and alpha):
    histo_Evis_wo_PSD_wo_alpha_IBD, bin_edges_evis = np.histogram(Evis_array_real_wo_PSD_wo_alpha_IBD, bins_evis)

    # build histogram from Evis_array_real_w_PSD_wo_alpha_IBD (with PSD suppression, but without alpha
    # (IBD cuts and PSD)):
    histo_Evis_w_PSD_wo_alpha_IBD, bin_edges_evis = np.histogram(Evis_array_real_w_PSD_wo_alpha_IBD, bins_evis)

    """ save information about total efficiencies and total alphas to txt file: """
    np.savetxt(output_path + "total_efficiencies_atmoNC_wo_PSD.txt",
               np.array([total_efficiency_NC_real_wo_PSD, stat_total_efficiency_NC_real_wo_PSD,
                         total_efficiency_NC_ideal_wo_PSD, stat_total_efficiency_NC_ideal_wo_PSD,
                         total_alpha_NC_wo_PSD, stat_total_alpha_NC_wo_PSD, number_IBDlike_real_simu_wo_PSD_wo_alpha,
                         number_IBDlike_JUNO_wo_PSD, number_events_total_NC]), fmt='%.5f',
               header="Total efficiency and total alpha of {0:.0f} atmoNC events without PSD analyzed with "
                      "atmoNC_spectrum_v2.py:"
                      "\ntotal efficiency of real data without PSD,"
                      "\nstatistical error of total efficiency of real data,"
                      "\ntotal efficiency of ideal data without PSD,"
                      "\nstatistical error of total efficiency of ideal data,"
                      "\ntotal alpha without PSD,"
                      "\nstatistical error of total alpha,"
                      "\nnumber of IBD-like events after all cuts (except of PSD, alpha not considered),"
                      "\nnumber of IBD-like events in JUNO after {1:.0f} years after all cuts (except of PSD, alpha "
                      "not considered),"
                      "\nnumber of analyzed NC events:".format(number_events_total_NC, time_in_years))

    np.savetxt(output_path + "total_efficiencies_atmoNC_w_PSD{0:.0f}.txt".format(PSD_NC_suppression*100.0),
               np.array([total_efficiency_NC_real_w_PSD, stat_total_efficiency_NC_real_w_PSD,
                         total_efficiency_NC_ideal_w_PSD, stat_total_efficiency_NC_ideal_w_PSD,
                         total_alpha_NC_w_PSD, stat_total_alpha_NC_w_PSD, number_IBDlike_real_simu_w_PSD_wo_alpha,
                         number_IBDlike_JUNO_w_PSD, number_events_total_NC]), fmt='%.5f',
               header="Total efficiency and total alpha of {0:.0f} atmoNC events with PSD analyzed with "
                      "atmoNC_spectrum_v2.py:"
                      "\ntotal efficiency of real data with PSD,"
                      "\nstatistical error of total efficiency of real data,"
                      "\ntotal efficiency of ideal data with PSD,"
                      "\nstatistical error of total efficiency of ideal data,"
                      "\ntotal alpha with PSD,"
                      "\nstatistical error of total alpha,"
                      "\nnumber of IBD-like events after all cuts (including PSD, alpha not considered),"
                      "\nnumber of IBD-like events in JUNO after {1:.0f} years after all cuts (including PSD, alpha "
                      "not considered),"
                      "\nnumber of analyzed NC events:".format(number_events_total_NC, time_in_years))

    np.savetxt(output_path + "total_efficiencies_IBD_wo_PSD.txt",
               np.array([total_efficiency_IBD_real_wo_PSD, stat_total_efficiency_IBD_real_wo_PSD,
                         total_efficiency_IBD_ideal_wo_PSD, stat_total_efficiency_IBD_ideal_wo_PSD,
                         total_alpha_IBD_wo_PSD, stat_total_alpha_IBD_wo_PSD, number_IBD_real_simu_wo_PSD_wo_alpha,
                         number_events_total_IBD]), fmt='%.5f',
               header="Total efficiency and total alpha of {0:.0f} IBD events without PSD analyzed with "
                      "atmoNC_spectrum_v2.py:"
                      "\ntotal efficiency of real data without PSD,"
                      "\nstatistical error of total efficiency of real data,"
                      "\ntotal efficiency of ideal data without PSD,"
                      "\nstatistical error of total efficiency of ideal data,"
                      "\ntotal alpha without PSD,"
                      "\nstatistical error of total alpha,"
                      "\nnumber of IBD events after all cuts (except of PSD, alpha not considered),"
                      "\nnumber of analyzed IBD events:".format(number_events_total_IBD))

    np.savetxt(output_path + "total_efficiencies_IBD_w_PSD{0:.0f}.txt".format(PSD_NC_suppression*100.0),
               np.array([total_efficiency_IBD_real_w_PSD, stat_total_efficiency_IBD_real_w_PSD,
                         total_efficiency_IBD_ideal_w_PSD, stat_total_efficiency_IBD_ideal_w_PSD,
                         total_alpha_IBD_w_PSD, stat_total_alpha_IBD_w_PSD, number_IBD_real_simu_w_PSD_wo_alpha,
                         number_events_total_IBD]), fmt='%.5f',
               header="Total efficiency and total alpha of {0:.0f} IBD events with PSD analyzed with "
                      "atmoNC_spectrum_v2.py:"
                      "\ntotal efficiency of real data with PSD,"
                      "\nstatistical error of total efficiency of real data,"
                      "\ntotal efficiency of ideal data with PSD,"
                      "\nstatistical error of total efficiency of ideal data,"
                      "\ntotal alpha with PSD,"
                      "\nstatistical error of total alpha,"
                      "\nnumber of IBD events after all cuts (including PSD, alpha not considered),"
                      "\nnumber of analyzed IBD events:".format(number_events_total_IBD))

    print("PSD NC suppression = {0:.2f}".format(PSD_NC_suppression))
    print("tail {0:.0f} ns to {1:.0f} ns, TTR = {2:.5f}".format(tail_start[index10], tail_end[index10],
                                                                TTR_cut_value[index10]))
    print("NC events:")
    print(number_IBDlike_real_simu_wo_PSD_wo_alpha)
    print(number_IBDlike_real_simu_w_PSD_wo_alpha)
    # calculate PSD NC suppression in %:
    PSD_NC_suppression_real = (100.0 - float(number_IBDlike_real_simu_w_PSD_wo_alpha) /
                               float(number_IBDlike_real_simu_wo_PSD_wo_alpha) * 100.0)
    print(PSD_NC_suppression_real)
    print("IBD events:")
    print(number_IBD_real_simu_wo_PSD_wo_alpha)
    print(number_IBD_real_simu_w_PSD_wo_alpha)
    # calculate PSD IBD suppression in %:
    PSD_IBD_suppression_real = (100.0 - float(number_IBD_real_simu_w_PSD_wo_alpha) /
                                float(number_IBD_real_simu_wo_PSD_wo_alpha) * 100.0)
    print(PSD_IBD_suppression_real)

    """ display simulated spectrum: """
    h1 = plt.figure(1, figsize=(11, 6))
    plt.plot(bins_evis[:-1], histo_Evis_wo_PSD_wo_alpha, drawstyle="steps", linestyle="-", color="orange",
             label="atmospheric NC background:\nnumber of events = {0:.0f},\ncut efficiency = {1:.1f} % $\\pm$ {2:.1f} "
                   "%"
             .format(number_IBDlike_real_simu_wo_PSD_wo_alpha, total_alpha_NC_wo_PSD * 100.0,
                     stat_total_alpha_NC_wo_PSD * 100.0))
    plt.xlim(xmin=E_min_prompt, xmax=E_max_prompt)
    plt.ylim(ymin=0.0)
    plt.xlabel("visible energy of prompt signal in MeV")
    plt.ylabel("number of IBD-like events per bin (bin-width = {0:.2f} MeV)".format(bin_width_energy))
    plt.title("Simulated spectrum of atmospheric NC neutrino events with IBD-like signature")
    plt.legend()
    plt.grid()
    plt.savefig(output_path + "atmoNC_spectrum_simu_woPSD_bins{0:.0f}keV.png".format(bin_width_energy*1000))
    # plt.show()
    plt.close()

    """ display simulated spectrum with PSD: """
    h2 = plt.figure(2, figsize=(11, 6))
    plt.plot(bins_evis[:-1], histo_Evis_w_PSD_wo_alpha, drawstyle="steps", linestyle="-", color="orange",
             label="atmospheric NC background:\nnumber of events = {0:.0f},\nPSD suppression of atmo. NC = {1:.2f} %,\n"
                   "PSD suppression of real IBD = {2:.2f} %"
             .format(number_IBDlike_real_simu_w_PSD_wo_alpha, PSD_NC_suppression_real, PSD_IBD_suppression_real))
    plt.xlim(xmin=E_min_prompt, xmax=E_max_prompt)
    plt.ylim(ymin=0.0)
    plt.xlabel("visible energy of prompt signal in MeV")
    plt.ylabel("number of IBD-like events per bin (bin-width = {0:.2f} MeV)".format(bin_width_energy))
    plt.title("Simulated spectrum of atmospheric NC neutrino events with IBD-like signature\n"
              "after Pulse Shape Discrimination")
    plt.legend()
    plt.grid()
    plt.savefig(output_path + "atmoNC_spectrum_simu_wPSD{1:.0f}_bins{0:.0f}keV.png".format(bin_width_energy*1000,
                                                                                           PSD_NC_suppression*100))
    # plt.show()
    plt.close()

    """ display simulated spectrum without and with PSD in JUNO after 10 years: """
    h9 = plt.figure(9, figsize=(11, 6))
    plt.semilogy(bins_evis[:-1], Evis_histo_JUNO_wo_PSD, drawstyle="steps", linestyle="--", color="orange",
                 label="without PSD: number of events = {0:.1f}"
                 .format(number_IBDlike_JUNO_wo_PSD))
    plt.semilogy(bins_evis[:-1], Evis_histo_JUNO_w_PSD, drawstyle="steps", linestyle="-", color="orange",
                 label="with PSD: number of events = {0:.0f}, PSD suppression of atmo. NC = {1:.2f} %"
                 .format(number_IBDlike_JUNO_w_PSD, PSD_NC_suppression_real))
    plt.xlim(xmin=E_min_prompt, xmax=E_max_prompt)
    plt.ylim(ymin=0.05)
    plt.xlabel("visible energy of prompt signal in MeV")
    plt.ylabel("number of IBD-like events per bin (bin-width = {0:.2f} MeV)".format(bin_width_energy))
    plt.title("Expected spectrum of atmospheric NC neutrino events with IBD-like signature in JUNO after {0:.0f} years"
              "\n(before and after Pulse Shape Discrimination)".format(time_in_years))
    plt.legend()
    plt.grid()
    plt.savefig(output_path + "atmoNC_spectrum_simu_wo_wPSD{1:.0f}_bins{0:.0f}keV.png".format(bin_width_energy*1000,
                                                                                              PSD_NC_suppression*100))
    # plt.show()
    plt.close()

    """ display visible spectrum in JUNO after 10 years: """
    h3 = plt.figure(3, figsize=(11, 6))
    plt.plot(bins_evis[:-1], Evis_histo_JUNO_wo_PSD, drawstyle="steps", linestyle="-", color="orange",
             label="atmospheric NC background:\nnumber of events = {0:.0f},\ncut efficiency = {1:.1f} % $\\pm$ {2:.1f} "
                   "%"
             .format(number_IBDlike_JUNO_wo_PSD, total_alpha_NC_wo_PSD*100.0, stat_total_alpha_NC_wo_PSD*100.0))
    plt.xlim(xmin=E_min_prompt, xmax=E_max_prompt)
    plt.ylim(ymin=0.0)
    plt.xlabel("visible energy of prompt signal in MeV")
    plt.ylabel("number of IBD-like events per bin (bin-width = {0:.2f} MeV)".format(bin_width_energy))
    plt.title("Expected spectrum of atmospheric NC neutrino events with IBD-like signature in JUNO after {0:.0f} years"
              .format(time_in_years))
    plt.legend()
    plt.grid()
    plt.savefig(output_path + "atmoNC_spectrum_JUNO_woPSD_bins{0:.0f}keV.png".format(bin_width_energy*1000))
    # plt.show()
    plt.close()

    """ display visible spectrum in JUNO after 10 years with PSD: """
    h4 = plt.figure(4, figsize=(11, 6))
    plt.plot(bins_evis[:-1], Evis_histo_JUNO_w_PSD, drawstyle="steps", linestyle="-", color="orange",
             label="atmospheric NC background:\nnumber of events = {0:.0f},\nPSD suppression of atmo. NC = {1:.2f} %,\n"
                   "PSD suppression of real IBD = {2:.2f} %"
             .format(number_IBDlike_JUNO_w_PSD, PSD_NC_suppression_real, PSD_IBD_suppression_real))
    plt.xlim(xmin=E_min_prompt, xmax=E_max_prompt)
    plt.ylim(ymin=0.0)
    plt.xlabel("visible energy of prompt signal in MeV")
    plt.ylabel("number of IBD-like events per bin (bin-width = {0:.2f} MeV)".format(bin_width_energy))
    plt.title("Expected spectrum of atmospheric NC neutrino events with IBD-like signature in JUNO after {0:.0f} "
              "years\nafter Pulse Shape Discrimination"
              .format(time_in_years))
    plt.legend()
    plt.grid()
    plt.savefig(output_path + "atmoNC_spectrum_JUNO_wPSD{1:.0f}_bins{0:.0f}keV.png".format(bin_width_energy*1000,
                                                                                           PSD_NC_suppression*100))
    # plt.show()
    plt.close()

    """ display simulated IBD and atmoNC spectrum: """
    h5 = plt.figure(5, figsize=(11, 6))
    plt.plot(bins_evis[:-1], histo_Evis_wo_PSD_wo_alpha, drawstyle="steps", linestyle="-", color="orange",
             label="atmospheric NC background: number of events = {0:.0f},\ncut efficiency = {1:.1f} % $\\pm$ {2:.1f} "
                   "%)"
             .format(number_IBDlike_real_simu_wo_PSD_wo_alpha, total_alpha_NC_wo_PSD * 100.0,
                     stat_total_alpha_NC_wo_PSD * 100.0))
    plt.plot(bins_evis[:-1], histo_Evis_wo_PSD_wo_alpha_IBD, drawstyle="steps", linestyle="-", color="blue",
             label="IBD signal: number of events = {0:.0f},\ncut efficiency = {1:.1f} % $\\pm$ {2:.1f} "
                   "%)"
             .format(number_IBD_real_simu_wo_PSD_wo_alpha, total_alpha_IBD_wo_PSD * 100.0,
                     stat_total_alpha_IBD_wo_PSD * 100.0))
    plt.xlim(xmin=E_min_prompt, xmax=E_max_prompt)
    plt.ylim(ymin=0.0)
    plt.xlabel("visible energy of prompt signal in MeV")
    plt.ylabel("number of events per bin (bin-width = {0:.2f} MeV)".format(bin_width_energy))
    plt.title("Simulated spectra of IBD and atmospheric NC events that pass the IBD selection criteria")
    plt.legend()
    plt.grid()
    plt.savefig(output_path + "IBD_atmoNC_spectra_simu_woPSD_bins{0:.0f}keV.png".format(bin_width_energy * 1000))
    # plt.show()
    plt.close()

    """ display simulated IBD and atmoNC spectrum with PSD: """
    h6 = plt.figure(6, figsize=(11, 6))
    plt.plot(bins_evis[:-1], histo_Evis_w_PSD_wo_alpha, drawstyle="steps", linestyle="-", color="orange",
             label="atmospheric NC background: number of events = {0:.0f},\nPSD suppression of atmo. NC = {1:.2f} %)"
             .format(number_IBDlike_real_simu_w_PSD_wo_alpha, PSD_NC_suppression_real))
    plt.plot(bins_evis[:-1], histo_Evis_w_PSD_wo_alpha_IBD, drawstyle="steps", linestyle="-", color="blue",
             label="IBD signal: number of events = {0:.0f},\nPSD suppression of real IBD = {1:.2f} %)"
             .format(number_IBD_real_simu_w_PSD_wo_alpha, PSD_IBD_suppression_real))
    plt.xlim(xmin=E_min_prompt, xmax=E_max_prompt)
    plt.ylim(ymin=0.0)
    plt.xlabel("visible energy of prompt signal in MeV")
    plt.ylabel("number of events per bin (bin-width = {0:.2f} MeV)".format(bin_width_energy))
    plt.title("Simulated spectra of IBD and atmospheric NC events that pass the IBD selection criteria\n"
              "(after Pulse Shape Discrimination)")
    plt.legend()
    plt.grid()
    plt.savefig(output_path + "IBD_atmoNC_spectrum_simu_wPSD{1:.0f}_bins{0:.0f}keV.png"
                .format(bin_width_energy * 1000, PSD_NC_suppression * 100))
    # plt.show()
    plt.close()

    """ display simulated IBD spectrum without and with PSD: """
    h10 = plt.figure(10, figsize=(11, 6))
    plt.plot(bins_evis[:-1], histo_Evis_wo_PSD_wo_alpha_IBD, drawstyle="steps", linestyle="--", color="blue",
             label="without PSD: number of events = {0:.0f}"
             .format(number_IBD_real_simu_wo_PSD_wo_alpha))
    plt.plot(bins_evis[:-1], histo_Evis_w_PSD_wo_alpha_IBD, drawstyle="steps", linestyle="-", color="blue",
             label="with PSD: number of events = {0:.0f},\nPSD suppression of real IBD = {1:.2f} %"
             .format(number_IBD_real_simu_w_PSD_wo_alpha, PSD_IBD_suppression_real))
    plt.xlim(xmin=E_min_prompt, xmax=E_max_prompt)
    plt.ylim(ymin=0.0)
    plt.xlabel("visible energy of prompt signal in MeV")
    plt.ylabel("number of events per bin (bin-width = {0:.2f} MeV)".format(bin_width_energy))
    plt.title("Simulated spectra of IBD events that pass the IBD selection criteria\n"
              "(before and after Pulse Shape Discrimination)")
    plt.legend()
    plt.grid()
    plt.savefig(output_path + "IBD_spectrum_simu_wo_wPSD{1:.0f}_bins{0:.0f}keV.png"
                .format(bin_width_energy * 1000, PSD_NC_suppression * 100))
    # plt.show()
    plt.close()

    """ save Evis_histo_JUNO_wo_PSD to txt file (txt file must have same shape like files in folder 
        /home/astro/blum/PhD/work/MeVDM_JUNO/gen_spectrum_v2/).
        Save the array before Pulse Shape Discrimination: """
    # save Evis_histo_JUNO_wo_PSD to txt-spectrum-file and information about simulation in txt-info-file:
    print("... save data of spectrum to file...")
    np.savetxt(output_path + 'NCatmo_onlyC12_woPSD_bin{0:.0f}keV.txt'
               .format(bin_width_energy * 1000), Evis_histo_JUNO_wo_PSD, fmt='%1.5e',
               header='Spectrum in IBD-like events/bin of atmospheric NC background events that mimic IBD signals '
                      'without PSD (calculated with atmoNC_spectrum_v2.py, {0}):'
                      '\n{3:.0f} NC events are simulated with JUNO detector software (tut_detsim.py).'
                      '\nNumber of IBD-like NC events in JUNO = {1:.2f}, binning of E_visible = {2:.3f} MeV,'
                      '\nNC interactions of nu_e, nu_e_bar, nu_mu and nu_mu_bar with C12 of liquid scintillator are '
                      'simulated with GENIE.'
                      '\nDeexcitation of residual isotopes are simulated with modified DSNB-NC.exe generator.'
                      '\nThen the final products are simulated with JUNO detector simulation and cuts are applied to '
                      'get'
                      '\nthe number of NC events that mimic an IBD signal:'
                      '\ntotal cut efficiency = {4:.2f} %,'
                      '\nstatistical error of total cut efficiency = {5:.2f} %:'
               .format(now, number_IBDlike_JUNO_wo_PSD, bin_width_energy, number_events_total_NC,
                       total_alpha_NC_wo_PSD*100.0, stat_total_alpha_NC_wo_PSD*100.0))
    np.savetxt(output_path + 'NCatmo_info_onlyC12_woPSD_bin{0:.0f}keV.txt'
               .format(bin_width_energy * 1000),
               np.array([E_min_prompt, E_max_prompt, bin_width_energy, time_in_years, R_cut_prompt,
                         number_events_total_NC,
                         number_IBDlike_JUNO_wo_PSD, event_rate, total_efficiency_NC_real_wo_PSD,
                         stat_total_efficiency_NC_real_wo_PSD, total_efficiency_NC_ideal_wo_PSD,
                         stat_total_efficiency_NC_ideal_wo_PSD, total_alpha_NC_wo_PSD, stat_total_alpha_NC_wo_PSD]),
               fmt='%1.9e',
               header='Information to simulation NCatmo_onlyC12_bin{0:.0f}keV.txt (analyzed files: user_atmoNC_0.root '
                      'to user_atmoNC_999.root):\n'
                      'values below: E_visible[0] in MeV, E_visible[-1] in MeV, interval_E_visible in MeV,'
                      '\nexposure time t_years in years, applied volume cut for radius in mm,'
                      '\nnumber of simulated NC events, '
                      '\nnumber of IBD-like NC events in spectrum JUNO will measure, '
                      '\ntheoretical NC event rate in JUNO detector in NC events/sec,'
                      '\ntotal efficiency of real NC data wo PSD,'
                      '\nstatistical error of total efficiency of real NC data wo PSD,'
                      '\ntotal efficiency of ideal NC data wo PSD,'
                      '\nstatistical error of total efficiency of ideal NC data wo PSD,'
                      '\ntotal alpha for NC events wo PSD (error of efficiency),'
                      '\nstatistical error of total alpha for NC events wo PSD:'
               .format(bin_width_energy * 1000))

    """ save Evis_histo_JUNO_w_PSD to txt file (txt file must have same shape like files in folder 
        /home/astro/blum/PhD/work/MeVDM_JUNO/gen_spectrum_v2/).
        Save the array before Pulse Shape Discrimination: """
    # save Evis_histo_JUNO_w_PSD to txt-spectrum-file and information about simulation in txt-info-file:
    print("... save data of spectrum to file...")
    np.savetxt(output_path + 'NCatmo_onlyC12_wPSD{1:.0f}_bin{0:.0f}keV.txt'
               .format(bin_width_energy * 1000, PSD_NC_suppression*100.0), Evis_histo_JUNO_w_PSD, fmt='%1.5e',
               header='Spectrum in IBD-like events/bin of atmospheric NC background events that mimic IBD signals '
                      'with PSD (calculated with atmoNC_spectrum_v2.py, {0}):'
                      '\n{3:.0f} NC events are simulated with JUNO detector software (tut_detsim.py).'
                      '\nNumber of IBD-like NC events in JUNO wo PSD = {1:.2f}, number of IBD-like NC events in JUNO w '
                      'PSD = {6:.2f}'
                      '\nbinning of E_visible = {2:.3f} MeV,'
                      '\nNC interactions of nu_e, nu_e_bar, nu_mu and nu_mu_bar with C12 of liquid scintillator are '
                      'simulated with GENIE.'
                      '\nDeexcitation of residual isotopes are simulated with modified DSNB-NC.exe generator.'
                      '\nThen the final products are simulated with JUNO detector simulation and cuts are applied to '
                      'get'
                      '\nthe number of NC events that mimic an IBD signal:'
                      '\nPSD NC suppression = {7:.2f} %'
                      '\nPSD IBD suppression = {8:.2f} %'
                      '\ntotal cut efficiency = {4:.2f} %,'
                      '\nstatistical error of total cut efficiency = {5:.2f} %:'
               .format(now, number_IBDlike_JUNO_wo_PSD, bin_width_energy, number_events_total_NC,
                       total_alpha_NC_w_PSD*100.0, stat_total_alpha_NC_w_PSD*100.0, number_IBDlike_JUNO_w_PSD,
                       PSD_NC_suppression_real, PSD_IBD_suppression_real))
    np.savetxt(output_path + 'NCatmo_info_onlyC12_wPSD{1:.0f}_bin{0:.0f}keV.txt'
               .format(bin_width_energy * 1000, PSD_NC_suppression*100.0),
               np.array([E_min_prompt, E_max_prompt, bin_width_energy, time_in_years, R_cut_prompt,
                         number_events_total_NC,
                         number_IBDlike_JUNO_wo_PSD, number_IBDlike_JUNO_w_PSD, event_rate,
                         total_efficiency_NC_real_w_PSD,
                         stat_total_efficiency_NC_real_w_PSD, total_efficiency_NC_ideal_w_PSD,
                         stat_total_efficiency_NC_ideal_w_PSD, total_alpha_NC_w_PSD, stat_total_alpha_NC_w_PSD]),
               fmt='%1.9e',
               header='Information to simulation NCatmo_onlyC12_bin{0:.0f}keV.txt (analyzed files: user_atmoNC_0.root '
                      'to user_atmoNC_999.root):\n'
                      'values below: E_visible[0] in MeV, E_visible[-1] in MeV, interval_E_visible in MeV,'
                      '\nexposure time t_years in years, applied volume cut for radius in mm,'
                      '\nnumber of simulated NC events, '
                      '\nnumber of IBD-like NC events in spectrum JUNO will measure wo PSD, '
                      '\nnumber of IBD-like events after PSD,'
                      '\ntheoretical NC event rate in JUNO detector in NC events/sec,'
                      '\ntotal efficiency of real NC data w PSD,'
                      '\nstatistical error of total efficiency of real NC data w PSD,'
                      '\ntotal efficiency of ideal NC data w PSD,'
                      '\nstatistical error of total efficiency of ideal NC data w PSD,'
                      '\ntotal alpha for NC events w PSD (error of efficiency),'
                      '\nstatistical error of total alpha for NC events w PSD:'
               .format(bin_width_energy * 1000))

    """ save filenumber_pass_vol_E_del, evtID_pass_vol_E_del and Evis_array_real_wo_PSD_wo_alpha to txt file: """
    np.savetxt(output_path + 'atmoNC_filenumber_evtID_Evis_pass_all_cuts_wo_PSD.txt',
               np.c_[filenumber_pass_vol_E_del, evtID_pass_vol_E_del, Evis_array_real_wo_PSD_wo_alpha], fmt='%i',
               header='filenumber | evtID | E_vis in MeV of events that pass all cuts (without PSD)')

    """ save filenumber_pass_vol_E_del_PSD, evtID_pass_vol_E_del_PSD and Evis_array_real_w_PSD_wo_alpha to txt file: """
    np.savetxt(output_path + 'atmoNC_filenumber_evtID_Evis_pass_all_cuts_w_PSD.txt',
               np.c_[filenumber_pass_vol_E_del_PSD, evtID_pass_vol_E_del_PSD, Evis_array_real_w_PSD_wo_alpha], fmt='%i',
               header='filenumber | evtID | E_vis in MeV of events that pass all cuts (with PSD)')

    """ save filenumber_pass_vol_E_del_IBD, evtID_pass_vol_E_del_IBD and Evis_array_real_wo_PSD_wo_alpha_IBD to txt 
    file: """
    np.savetxt(output_path + 'IBD_filenumber_evtID_Evis_pass_all_cuts_wo_PSD.txt',
               np.c_[filenumber_pass_vol_E_del_IBD, evtID_pass_vol_E_del_IBD, Evis_array_real_wo_PSD_wo_alpha_IBD],
               fmt='%i',
               header='filenumber | evtID | E_vis in MeV of events that pass all cuts (without PSD)')

    """ save filenumber_pass_vol_E_del_PSD_IBD, evtID_pass_vol_E_del_PSD_IBD and Evis_array_real_w_PSD_wo_alpha_IBD to 
    txt file: """
    np.savetxt(output_path + 'IBD_filenumber_evtID_Evis_pass_all_cuts_w_PSD.txt',
               np.c_[filenumber_pass_vol_E_del_PSD_IBD, evtID_pass_vol_E_del_PSD_IBD,
                     Evis_array_real_w_PSD_wo_alpha_IBD], fmt='%i',
               header='filenumber | evtID | E_vis in MeV of events that pass all cuts (with PSD)')

    """ display TTR value of events, that pass volume, prompt energy and delayed cut (before PSD) in histogram: """
    h7 = plt.figure(7, figsize=(11, 6))
    First_bin = 0.0

    if max(array_TTR_IBD_beforePSD) >= max(array_TTR_IBDlike):
        maximum_tot_value = max(array_TTR_IBD_beforePSD)
    else:
        maximum_tot_value = max(array_TTR_IBDlike)

    maximum_tot_value = 0.07

    # Last_bin = maximum_tot_value
    Last_bin = maximum_tot_value
    Bin_width = (Last_bin-First_bin) / 200
    Bins = np.arange(First_bin, Last_bin+Bin_width, Bin_width)
    plt.hist(array_TTR_IBD_beforePSD, bins=Bins, histtype="step", align='mid', color="r", linewidth=1.5,
             label="prompt signal of IBD events (entries = {0:d})".format(len(array_TTR_IBD_beforePSD)))
    plt.hist(array_TTR_IBDlike, bins=Bins, histtype="step", align='mid', color="b", linewidth=1.5,
             label="prompt signal of IBD-like NC events (entries = {0:d})".format(len(array_TTR_IBDlike)))
    plt.xlim(xmin=0.0, xmax=maximum_tot_value)
    plt.xlabel("tail-to-total ratio")
    plt.ylabel("events")
    plt.title("Tail-to-total ratio of prompt signals of IBD and NC events" +
              "\n(tail window {0:0.1f} ns to {1:0.1f} ns)".format(tail_start[index10], tail_end[index10]))
    plt.legend()
    plt.grid()
    plt.savefig(output_path + "TTR_beforePSD_{0:0.0f}ns_to_{1:0.0f}ns.png".format(tail_start[index10],
                                                                                  tail_end[index10]))
    plt.close()

    # with efficiency:
    h8 = plt.figure(8, figsize=(11, 6))
    n_1, bins_1, patches_1 = plt.hist(array_TTR_IBD_beforePSD, bins=Bins, histtype="step", align='mid', color="r",
                                      linewidth=1.5, label="prompt signal of IBD events (entries = {0:d})"
                                      .format(len(array_TTR_IBD_beforePSD)))
    n_2, bins_2, patches_2 = plt.hist(array_TTR_IBDlike, bins=Bins, histtype="step", align='mid', color="b",
                                      linewidth=1.5,
                                      label="prompt signal of IBD-like NC events (entries = {0:d})"
                                      .format(len(array_TTR_IBDlike)))
    plt.vlines(TTR_cut_value[index10], 0.0, max(n_1)+max(n_1)/10, colors="k", linestyles="-",
               label="$\\epsilon_{NC}$ = "+"{0:0.2f} %\n".format(PSD_NC_suppression_real) +
                     "$\\epsilon_{IBD}$"+" = {0:0.2f} %\n".format(PSD_IBD_suppression_real) +
                     "ttr cut value = {0:.5f}".format(TTR_cut_value[index10]))
    plt.xlim(xmin=0.0, xmax=maximum_tot_value)
    plt.xlabel("tail-to-total ratio")
    plt.ylabel("events")
    plt.title("Tail-to-total ratio of prompt signals of IBD and NC events" +
              "\n(tail window {0:0.1f} ns to {1:0.1f} ns)".format(tail_start[index10], tail_end[index10]))
    plt.legend()
    plt.grid()
    plt.savefig(output_path + "TTR_beforePSD_{0:0.0f}ns_to_{1:0.0f}ns_with_Eff.png".format(tail_start[index10],
                                                                                           tail_end[index10]))
    plt.close()

    """ save TTR values of NC events and IBD events that pass all cuts (before PSD) to txt file: """
    np.savetxt(output_path + "TTR_IBDlike_NCevents_{0:.0f}ns_to_{1:.0f}ns.txt".format(tail_start[index10],
                                                                                      tail_end[index10]),
               array_TTR_IBDlike, fmt='%.5f',
               header="TTR values of {0:d} IBDlike NC events (NC events that pass all cuts):"
               .format(len(array_TTR_IBDlike)))

    np.savetxt(output_path + "TTR_beforePSD_IBDevents_{0:.0f}ns_to_{1:.0f}ns.txt".format(tail_start[index10],
                                                                                         tail_end[index10]),
               array_TTR_IBD_beforePSD, fmt='%.5f',
               header="TTR values of {0:d} IBD events, that pass all cuts (before PSD):"
               .format(len(array_TTR_IBD_beforePSD)))

    """ 2D histogram with TTR values and Evis for IBD events before PSD: """
    h11 = plt.figure(11)
    plt.hist2d(Evis_array_real_wo_PSD_wo_alpha_IBD, array_TTR_IBD_beforePSD, [bins_evis, Bins])
    plt.hlines(TTR_cut_value[index10], xmin=min(Evis_array_real_wo_PSD_wo_alpha_IBD),
               xmax=max(Evis_array_real_wo_PSD_wo_alpha_IBD), colors="k", linestyles="-",
               label="TTR cut value = {0:.5f}".format(TTR_cut_value[index10]))
    plt.ylim(ymin=0.0, ymax=maximum_tot_value)
    plt.xlabel("visible energy in MeV")
    plt.ylabel("TTR values for tail window from {0:.0f} ns to {1:.0f}".format(tail_start[index10], tail_end[index10]))
    plt.title("TTR vs visible energy of IBD events that pass IBD selection criteria")
    plt.legend()
    plt.grid()
    plt.savefig(output_path + "2D_IBD_TTR_vs_Evis_{0:0.0f}ns_to_{1:0.0f}ns.png".format(tail_start[index10],
                                                                                       tail_end[index10]))
    plt.close()

    """ 2D histogram with TTR values and Evis for IBD events before PSD: """
    h12 = plt.figure(12)
    plt.hist2d(Evis_array_real_wo_PSD_wo_alpha, array_TTR_IBDlike, [bins_evis, Bins])
    plt.hlines(TTR_cut_value[index10], xmin=min(Evis_array_real_wo_PSD_wo_alpha),
               xmax=max(Evis_array_real_wo_PSD_wo_alpha), colors="k", linestyles="-",
               label="TTR cut value = {0:.5f}".format(TTR_cut_value[index10]))
    plt.ylim(ymin=0.0, ymax=maximum_tot_value)
    plt.xlabel("visible energy in MeV")
    plt.ylabel("TTR values for tail window from {0:.0f} ns to {1:.0f}".format(tail_start[index10], tail_end[index10]))
    plt.title("TTR vs visible energy of IBD-like atmospheric NC events")
    plt.legend()
    plt.grid()
    plt.savefig(output_path + "2D_atmoNC_TTR_vs_Evis_{0:0.0f}ns_to_{1:0.0f}ns.png".format(tail_start[index10],
                                                                                          tail_end[index10]))
    plt.close()


