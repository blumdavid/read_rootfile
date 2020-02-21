""" Script to get the NC interaction channel and particles produced via NC, the deexcitation channel and particles
    produced via deexcitation, and the total (NC + deex.) channels for the NC events, that mimic an IBD-like signal.

    This script is very similar to checkout_NCgen.py, BUT not all events of gen_NC_onlyC12_250000evts_seed1.root are
    read, but only atmospheric NC events, that can mimic an IBD-like signal in JUNO detector.

    For all IBD-like events (events, that pass all cuts like volume, prompt energy, delayed energy, time, neutron
    multiplicity and distance cut), the evtID and the file number of user_atmoNC_{}.root are stored in folder
    /home/astro/blum/juno/atmoNC/data_NC/output_detsim/ in the name of the txt- and png-files of the prompt signals.

"""
import datetime
import NC_background_functions
import os
import numpy as np
from matplotlib import pyplot as plt
import re
import sys


def get_numbers_from_filename(filename):
    """
    function to get number as integer out of a string. For example: filename='file235' -> num = 235 of type integer

    :param filename: string of one part of the filename 'file{}_evt{}_prompt_signal.txt'
    :return:
    """
    # get the number out of filename and convert it into integer:
    num = int(re.search(r'\d+', filename).group(0))

    return num


# get the date and time, when the script was run:
date = datetime.datetime.now()
now = date.strftime("%Y-%m-%d %H:%M")

# set SAVE_FIG, defines if figures are saved:
SAVE_FIG = True

# set SAVE_TXT, defines if txt files are saved:
SAVE_TXT = True

# set SHOW_PLOT, defines if the figures are shown:
SHOW_PLOT = True

# set the path of the input root file:
input_path = "/home/astro/blum/juno/atmoNC/data_NC/output_generator/"
# file name of the input rootfile:
input_name = input_path + "gen_NC_onlyC12_250000evts_seed1.root"

# set the path of the folder, where prompt signal of IBD-like signals are stored:
path_ibdlike_event = "/home/astro/blum/juno/atmoNC/data_NC/output_detsim/"

# set the path, where the outputs are saved:
output_path = "/home/astro/blum/juno/atmoNC/data_NC/output_detsim/interaction_channels_IBDlike_signal/"

# bin-width of the array, which represents the incoming neutrino energy (in GeV) (float):
bin_width_incoming = 0.1

# number of events per user_atmoNC_{}.root file:
number_per_user_atmoNC_file = 100

# set the NC interaction event rate. Set evt_rate_Genie = 1, because the event rate is not yet included into function
# get_neutrino_energy():
evt_rate_Genie = 1

# total exposure time. Set time = 1, because the exposure time is not yet included into function get_neutrino_energy():
time = 1

""" get the evtID and file number of IBD-like events to calculate the event number that corresponds to t_evtID in 
gen_NC_onlyC12_250000evts_seed1.root: """
# preallocate array, where the event number is stored:
event_number = []
# preallocate number of file, that are read:
number_of_files_read = 0
# loop over all files in folder path_ibdlike_event, that start with 'file' and end with '_prompt_signal.txt' (the name
# of these files contains the file number and evtID):
for file_ibdlike in os.listdir(path_ibdlike_event):

    if file_ibdlike.startswith("file") and file_ibdlike.endswith("_prompt_signal.txt"):

        # increment number_of_file_read:
        number_of_files_read += 1

        # split string file_ibdlike into two parts: 1. part: 'file{}', 2. part: 'vt{}_prompt_signal.txt'
        x = file_ibdlike.split("_e")

        # x[0] is string 'file{}':
        file_string = x[0]
        # x[1] is string 'vt{}_prompt_signal.txt':
        event_string = x[1]

        # get file_number of file_string:
        file_number = get_numbers_from_filename(file_string)
        # get evtID of event_string:
        evtID = get_numbers_from_filename(event_string)

        # calculate the event_num with file_number, evtID and number_per_user_atmoNC_file:
        event_num = file_number * number_per_user_atmoNC_file + evtID

        # append event_num to event_number array:
        event_number.append(event_num)

# event numbers in event_number are not sorted. Therefore sort event_number array in ascending:
event_number.sort()

# read NC generator data only of the IBD-like NC events to arrays:
(event_ID, projectile_PDG, projectile_E, target_PDG, NC_inter_ch_ID, deexcitation_ID, isotope_PDG, Nparticles,
 final_PDG, final_Px, final_Py, final_Pz) = NC_background_functions.read_nc_data_ibdlike_signal(input_name,
                                                                                                event_number)

# check, if number_of_files_read is equal to len(event_ID):
if number_of_files_read != len(event_ID):
    sys.exit("ERROR: number of files, that are read ({0:d}), is NOT equal to len(event_ID) ({1:d})"
             .format(number_of_files_read, len(event_ID)))
else:
    number_ibdlike_events = number_of_files_read

""" get the number of events as function of the energy of the incoming neutrinos for each neutrino type: """
(Energy_nu_incoming,
 Event_nu_e_IN, Event_nu_e_bar_IN, Event_nu_mu_IN, Event_nu_mu_bar_IN, Event_nu_tau_IN, Event_nu_tau_bar_IN,
 Number_nu_e_IN, Number_nu_e_bar_IN, Number_nu_mu_IN, Number_nu_mu_bar_IN, Number_nu_tau_IN, Number_nu_tau_bar_IN,
 Frac_nu_e_IN, Frac_nu_e_bar_IN, Frac_nu_mu_IN, Frac_nu_mu_bar_IN, Frac_nu_tau_IN, Frac_nu_tau_bar_IN) \
    = NC_background_functions.get_neutrino_energy(projectile_PDG, projectile_E, bin_width_incoming, evt_rate_Genie,
                                                  time)

""" get the fractions of events of the different types on NC interaction channels: """
(Frac_C12_B11_p, Frac_C12_B11_n_piplus, Frac_C12_B11_n_piminus_2piplus, Frac_C12_B11_p_piminus_piplus,
 Frac_C12_B11_p_2piminus_2piplus, Frac_C12_B11_piplus, Frac_C12_B11_other,
 Frac_C12_C11_n, Frac_C12_C11_p_piminus, Frac_C12_C11_n_piminus_piplus, Frac_C12_C11_p_2piminus_piplus,
 Frac_C12_C11_p_3piminus_2piplus, Frac_C12_C11_n_2piminus_2piplus, Frac_C12_C11_other,
 Frac_C12_B10_p_n, Frac_C12_B10_2p_piminus, Frac_C12_B10_p_n_piminus_piplus,
 Frac_C12_B10_2n_piplus, Frac_C12_B10_2n_piminus_2piplus, Frac_C12_B10_2p_2piminus_piplus,
 Frac_C12_B10_2p_3piminus_2piplus, Frac_C12_B10_p_n_2piminus_2piplus, Frac_C12_B10_other,
 Frac_C12_C10_2n, Frac_C12_C10_p_n_piminus, Frac_C12_C10_p_n_2piminus_piplus, Frac_C12_C10_2n_piminus_piplus,
 Frac_C12_C10_2p_2piminus, Frac_C12_C10_other,
 Frac_C12_Be10_2p, Frac_C12_Be10_p_n_piplus, Frac_C12_Be10_p_n_piminus_2piplus, Frac_C12_Be10_2p_piminus_piplus,
 Frac_C12_Be10_2n_2piplus, Frac_C12_Be10_p_n_2piminus_3piplus, Frac_C12_Be10_2p_2piminus_2piplus,
 Frac_C12_Be10_2p_3piminus_3piplus, Frac_C12_Be10_other,
 Frac_C12_B9_p_2n, Frac_C12_B9_p_2n_piminus_piplus, Frac_C12_B9_2p_n_3piminus_2piplus, Frac_C12_B9_2p_n_piminus,
 Frac_C12_B9_3n_piplus, Frac_C12_B9_p_2n_2piminus_2piplus, Frac_C12_B9_2p_n_2piminus_piplus, Frac_C12_B9_other,
 Frac_C12_Be9_2p_n, Frac_C12_Be9_p_2n_piplus, Frac_C12_Be9_3p_piminus, Frac_C12_Be9_p_2n_piminus_2piplus,
 Frac_C12_Be9_2p_n_piminus_piplus, Frac_C12_Be9_2p_n_3piminus_3piplus, Frac_C12_Be9_2p_n_2piminus_2piplus,
 Frac_C12_Be9_3n_2piplus, Frac_C12_Be9_3p_2piminus_piplus, Frac_C12_Be9_other,
 Frac_C12_Be8_2p_2n, Frac_C12_Be8_3p_n_piminus, Frac_C12_Be8_p_3n_piplus, Frac_C12_Be8_2p_2n_2piminus_2piplus,
 Frac_C12_Be8_4n_2piplus, Frac_C12_Be8_2p_2n_piminus_piplus, Frac_C12_Be8_3p_n_2piminus_piplus,
 Frac_C12_Be8_4p_2piminus, Frac_C12_Be8_other,
 Frac_C12_C9_p_2n_piminus, Frac_C12_C9_3n, Frac_C12_C9_2p_n_2piminus, Frac_C12_C9_3n_2piminus_2piplus,
 Frac_C12_C9_other,
 Frac_C12_Be7_2p_3n, Frac_C12_Be7_p_4n_piplus, Frac_C12_Be7_2p_3n_2piminus_2piplus, Frac_C12_Be7_3p_2n_piminus,
 Frac_C12_Be7_4p_n_2piminus, Frac_C12_Be7_3p_2n_2piminus_piplus, Frac_C12_Be7_other,
 Frac_C12_Li6_3p_3n, Frac_C12_Li6_2p_4n_piplus, Frac_C12_Li6_5p_n_2piminus, Frac_C12_Li6_2p_4n_piminus_2piplus,
 Frac_C12_Li6_4p_2n_piminus, Frac_C12_Li6_3p_3n_piminus_piplus, Frac_C12_Li6_other,
 Frac_C12_Li8_3p_n, Frac_C12_Li8_4p_piminus, Frac_C12_Li8_4p_2piminus_piplus, Frac_C12_Li8_2p_2n_piplus,
 Frac_C12_Li8_3p_n_piminus_piplus, Frac_C12_Li8_other,
 Frac_C12_Li7_2p_3n_piplus, Frac_C12_Li7_4p_n_piminus, Frac_C12_Li7_3p_2n, Frac_C12_Li7_3p_2n_piminus_piplus,
 Frac_C12_Li7_4p_n_2piminus_piplus, Frac_C12_Li7_2p_3n_piminus_2piplus, Frac_C12_Li7_other,
 Frac_C12_B8_p_3n, Frac_C12_B8_p_3n_piminus_piplus, Frac_C12_B8_2p_2n_2piminus_piplus, Frac_C12_B8_2p_2n_piminus,
 Frac_C12_B8_4n_piplus, Frac_C12_B8_other,
 Frac_C12_Li9_2p_n_piplus, Frac_C12_Li9_3p, Frac_C12_Li9_3p_piminus_piplus, Frac_C12_Li9_2p_n_piminus_2piplus,
 Frac_C12_Li9_p_2n_piminus_3piplus, Frac_C12_Li9_other,
 Frac_C12_C8_4n, Frac_C12_C8_4n_other, Frac_C12_He8_4p, Frac_C12_He8_4p_other, Frac_C12_B7_p_4n, Frac_C12_B7_p_4n_other,
 Frac_C12_He7_4p_n, Frac_C12_He7_4p_n_other, Frac_C12_H7_5p, Frac_C12_H7_5p_other,
 Frac_C12_Be6_2p_4n, Frac_C12_Be6_2p_4n_other, Frac_C12_He6_4p_2n, Frac_C12_He6_4p_2n_other,
 Frac_C12_H6_5p_n, Frac_C12_H6_5p_n_other,
 Frac_C12_mass11u, Frac_C12_mass10u, Frac_C12_mass9u, Frac_C12_mass8u, Frac_C12_mass7u, Frac_C12_mass6u,
 Frac_C12_mass5orless,
 Frac_C12_C12, Frac_C12_NoIso, Frac_C12_NoIso_5p_6n,
 Frac_no_C12, Frac_ES_proton_chID, Frac_ES_electron_chID, Frac_ES_O16_chID, Frac_ES_N14_chID, Frac_ES_S32_chID,
 Frac_C12_faulty) \
    = NC_background_functions.get_interaction_channel(NC_inter_ch_ID, isotope_PDG, target_PDG)

""" get the numbers of events of the different deexcitation channels of the different produced isotopes: """
(number_Entries, number_target_C12, number_no_C12, number_light_Iso,
 number_c11_notex, number_c11_deex, number_c11_li6_p_alpha, number_c11_be7_alpha, number_c11_b10_p, number_c11_b9_n_p,
 number_c11_be8_p_d, number_c11_be9_2p, number_c11_b9_d, number_c11_be8_he3, number_c11_c10_n, number_c11_li5_d_alpha,
 number_c11_li5_n_p_alpha, number_c11_missing,
 number_b11_notex, number_b11_deex, number_b11_li6_n_alpha, number_b11_b9_2n, number_b11_be8_n_d, number_b11_be9_d,
 number_b11_be10_p, number_b11_b10_n, number_b11_be9_n_p, number_b11_li7_alpha, number_b11_be8_t,
 number_b11_he5_d_alpha, number_b11_he6_p_alpha, number_b11_be8_2n_p, number_b11_missing,
 number_c10_notex, number_c10_deex, number_c10_b9_p, number_c10_be7_p_d, number_c10_li6_p_he3, number_c10_he4_p_d_he3,
 number_c10_li6_2p_d, number_c10_be8_2p, number_c10_be7_n_2p, number_c10_li6_n_3p, number_c10_be6_n_p_d,
 number_c10_li5_n_2p_d, number_c10_he3_p_d_alpha, number_c10_li5_d_he3, number_c10_li5_p_2d, number_c10_he3_n_2p_alpha,
 number_c10_li4_n_p_alpha, number_c10_b8_n_p, number_c10_b8_d, number_c10_be6_p_t, number_c10_he4_n_2p_he3,
 number_c10_li5_n_p_he3, number_c10_missing,
 number_b10_notex, number_b10_deex, number_b10_be9_p, number_b10_be8_d, number_b10_b9_n, number_b10_be7_t,
 number_b10_be8_n_p, number_b10_li7_he3, number_b10_he5_p_alpha, number_b10_li6_alpha, number_b10_li5_n_alpha,
 number_b10_li7_p_d, number_b10_missing,
 number_be10_notex, number_be10_deex, number_be10_li6_2n_d, number_be10_li7_2n_p, number_be10_he7_n_2p,
 number_be10_li6_n_t, number_be10_li6_3n_p, number_be10_he5_n_2d, number_be10_be8_2n, number_be10_he5_n_alpha,
 number_be10_t_n_d_alpha, number_be10_li8_n_p, number_be10_he4_n_d_t, number_be10_he5_n_p_t, number_be10_li7_n_d,
 number_be10_h4_n_p_alpha, number_be10_he5_2n_p_d, number_be10_be9_n, number_be10_he6_n_p_d, number_be10_he5_d_t,
 number_be10_h4_d_alpha, number_be10_he4_2n_p_t, number_be10_missing,
 number_c9_notex, number_c9_deex, number_c9_be7_2p, number_c9_missing,
 number_b9_notex, number_b9_deex, number_b9_be8_p, number_b9_missing,
 number_be9_notex, number_be9_deex, number_be9_li8_p, number_be9_missing,
 number_li9_notex, number_li9_deex, number_li9_h4_n_alpha, number_li9_he7_d, number_li9_li8_n, number_li9_missing,
 number_b8_notex, number_b8_deex, number_b8_li6_2p, number_b8_missing,
 number_li8_notex, number_li8_deex, number_li8_li7_n, number_li8_li6_2n, number_li8_missing,
 number_be7_notex, number_be7_deex, number_be7_li5_d, number_be7_li6_p, number_be7_missing,
 number_li7_notex, number_li7_deex, number_li7_li6_n, number_li7_missing) \
    = NC_background_functions.get_deex_channel(deexcitation_ID, isotope_PDG, target_PDG)

""" get the number of events as function of the incoming neutrino energy for the different residual isotopes BEFORE
    deexcitation (C12, C11, B11, C10, B10, Be10, C9, B9, Be9, Li9, C8, B8, Be8, Li8, He8, B7, Be7, Li7, He7, H7, Be6, 
    Li6, He6, H6): """
(Energy_nu_incoming_2, Events_c12, Events_c11, Events_b11, Events_c10, Events_b10, Events_be10, Events_c9, Events_b9,
 Events_be9, Events_li9, Events_c8, Events_b8, Events_be8, Events_li8, Events_he8, Events_b7, Events_be7, Events_li7,
 Events_he7, Events_h7, Events_be6, Events_li6, Events_he6, Events_h6,
 Number_c12, Number_c11, Number_b11, Number_c10, Number_b10, Number_be10, Number_c9, Number_b9, Number_be9, Number_li9,
 Number_c8, Number_b8, Number_be8, Number_li8, Number_he8, Number_b7, Number_be7, Number_li7, Number_he7, Number_h7,
 Number_be6, Number_li6, Number_he6, Number_h6, Number_rest,
 Fraction_c12, Fraction_c11, Fraction_b11, Fraction_c10, Fraction_b10, Fraction_be10, Fraction_c9, Fraction_b9,
 Fraction_be9, Fraction_li9, Fraction_c8, Fraction_b8, Fraction_be8, Fraction_li8, Fraction_he8, Fraction_b7,
 Fraction_be7, Fraction_li7, Fraction_he7, Fraction_h7, Fraction_be6, Fraction_li6, Fraction_he6, Fraction_h6,
 Fraction_rest) = NC_background_functions.get_residual_isotopes_before_deex(projectile_E, isotope_PDG, target_PDG,
                                                                            bin_width_incoming)

""" ge the number of events of the different channels, that mimic IBD-like signal. (combine NC interaction channel and 
    deexcitation channel to get to know, which final particles produce IBD-like signal inside JUNO detector): """
(num_c12, num_no_c12,
 num_b11, num_b11_n_piplus, num_b11_p,
 num_c11, num_c11_n, num_c11_p_piminus,
 num_c10, num_c10_2n, num_c10_n_p_piminus,
 num_b10, num_b10_n_p,
 num_be10, num_be10_2p, num_be10_n_p_piplus,
 num_c9, num_c9_3n, num_c9_2n_p_piminus,
 num_b9, num_b9_n_d, num_b9_2n_p, num_b9_n_2p_piminus, num_b9_3n_piplus,
 num_be9, num_be9_n_2p, num_be9_p_d, num_be9_2n_p_piplus,
 num_li9, num_li9_3p, num_li9_n_2p_piplus,
 num_b8, num_b8_3n_p, num_b8_2n_d,
 num_be8, num_be8_n_p_d, num_be8_2n_2p, num_be8_n_3p_piminus, num_be8_n_he3, num_be8_p_t,
 num_li8, num_li8_n_3p,
 num_he8, num_he8_4p,
 num_b7,
 num_be7, num_be7_n_p_t, num_be7_3n_2p, num_be7_n_alpha, num_be7_2n_p_d,
 num_li7, num_li7_n_p_he3, num_li7_2n_3p, num_li7_n_2p_d, num_li7_p_alpha, num_li7_n_alpha_piplus,
 num_he7, num_he7_n_4p,
 num_be6, num_be6_3n_p_d, num_be6_2n_2d, num_be6_4n_2p, num_be6_2n_p_t,
 num_li6, num_li6_n_p_alpha, num_li6_2n_2p_d, num_li6_n_2p_t, num_li6_3n_3p, num_li6_2n_p_he3, num_li6_2n_alpha_piplus,
 num_li6_2n_4p_piminus, num_li6_n_2p_he3_piminus,
 num_he6, num_he6_2n_4p, num_he6_n_3p_d, num_he6_n_2p_he3,
 num_h6, num_h6_n_5p,
 num_li5, num_li5_2n_p_2d, num_li5_2n_p_alpha, num_li5_n_alpha_d, num_li5_3n_p_he3, num_li5_3n_2p_d, num_li5_2n_he3_d,
 num_he5, num_he5_n_3p_t, num_he5_n_2p_alpha, num_he5_n_2p_2d, num_he5_3n_4p,
 num_h5, num_h5_2n_5p,
 num_li4, num_li4_3n_p_alpha,
 num_he4, num_he4_n_2p_t_d, num_he4_3n_2p_he3, num_he4_n_alpha_he3, num_he4_2n_p_he3_d, num_he4_4n_4p, num_he4_2n_2he3,
 num_h4, num_h4_n_3p_alpha, num_h4_3n_5p,
 number_false_channel,
 interesting_index) = NC_background_functions.get_combined_channel(Nparticles, final_PDG, target_PDG)

print("event number {0:.0f} of interesting index {1:d}".format(event_ID[interesting_index], interesting_index))

""" Display Output of 'get_neutrino_energy()' in plot: """
h1 = plt.figure(1, figsize=(15, 8))

# do not display the histogram, when there are no events:
if Number_nu_e_IN != 0:
    plt.plot(Energy_nu_incoming[:-1], Event_nu_e_IN, drawstyle='steps', color='b',
             label="interactions with $\\nu_e$: $N_{events} = $"+"{0:.2f}; fraction = {1:.2f}%"
             .format(Number_nu_e_IN, Frac_nu_e_IN))
if Number_nu_e_bar_IN != 0:
    plt.plot(Energy_nu_incoming[:-1], Event_nu_e_bar_IN, drawstyle='steps', color='b', linestyle='--',
             label="interactions with $\\bar{\\nu}_e$: $N_{events} = $"+"{0:.2f}; fraction = {1:.2f}%"
             .format(Number_nu_e_bar_IN, Frac_nu_e_bar_IN))
if Number_nu_mu_IN != 0:
    plt.plot(Energy_nu_incoming[:-1], Event_nu_mu_IN, drawstyle='steps', color='r',
             label="interactions with $\\nu_\\mu$: $N_{events} = $"+"{0:.2f}; fraction = {1:.2f}%"
             .format(Number_nu_mu_IN, Frac_nu_mu_IN))
if Number_nu_mu_bar_IN != 0:
    plt.plot(Energy_nu_incoming[:-1], Event_nu_mu_bar_IN, drawstyle='steps', color='r', linestyle='--',
             label="interactions with $\\bar{\\nu}_\\mu$: $N_{events} = $"+"{0:.2f}; fraction = {1:.2f}%"
             .format(Number_nu_mu_bar_IN, Frac_nu_mu_bar_IN))
if Number_nu_tau_IN != 0:
    plt.plot(Energy_nu_incoming[:-1], Event_nu_tau_IN, drawstyle='steps', color='g',
             label="interactions with $\\nu_\\tau$: $N_{events} = $"+"{0:.2f}; fraction = {1:.2f}%"
             .format(Number_nu_tau_IN, Frac_nu_tau_IN))
if Number_nu_tau_bar_IN != 0:
    plt.plot(Energy_nu_incoming[:-1], Event_nu_tau_bar_IN, drawstyle='steps', color='g', linestyle='--',
             label="interactions with $\\bar{\\nu}_\\tau$: $N_{events} = $"+"{0:.2f}; fraction = {1:.2f}%"
             .format(Number_nu_tau_bar_IN, Frac_nu_tau_bar_IN))

plt.xlim(xmin=0, xmax=10.0)
plt.ylim(ymin=0)
plt.xlabel("Neutrino energy $E_{\\nu}$ in GeV", fontsize=15)
plt.ylabel("events per bin (bin-width = {0:.2f} GeV)".format(bin_width_incoming), fontsize=15)
plt.title("Energy spectrum of atmospheric neutrinos \n"
          "(interacting via NC on $^{12}C$ and mimicking IBD-like signal in the JUNO detector)", fontsize=20)
plt.legend(fontsize=12)

if SAVE_FIG:
    plt.savefig(output_path + "IBDlike_signals_incoming_neutrino_spectrum_{0:.0f}evts.png"
                .format(number_ibdlike_events))

""" Save information about the ratio of incoming neutrinos into txt file: """
if SAVE_TXT:
    np.savetxt(output_path + "IBDlike_signals_incoming_neutrino_spectrum_{0:.0f}evts.txt".format(number_ibdlike_events),
               np.array([number_ibdlike_events, Number_nu_e_IN, Number_nu_e_bar_IN, Number_nu_mu_IN,
                         Number_nu_mu_bar_IN,
                         Number_nu_tau_IN, Number_nu_tau_bar_IN, Frac_nu_e_IN, Frac_nu_e_bar_IN, Frac_nu_mu_IN,
                         Frac_nu_mu_bar_IN, Frac_nu_tau_IN, Frac_nu_tau_bar_IN]), fmt='%4.5f',
               header="Information about the ratio of the incoming atmospheric neutrinos that interact with the\n"
                      "JUNO liquid scintillator via NC interactions on C12 and mimic an IBD-like signal in the detector"
                      "\n"
                      "(input file: {0}, script: OLD_IBDlike_events_channels_v1.py ({1})):\n"
                      "Number of events in the input file,\n"
                      "Number of interactions with nu_e,\n"
                      "Number of interactions with nu_e_bar,\n"
                      "Number of interactions with nu_mu,\n"
                      "Number of interactions with nu_mu_bar,\n"
                      "Number of interactions with nu_tau,\n"
                      "Number of interactions with nu_tau_bar,\n"
                      "Fraction of interactions with nu_e,\n"
                      "Fraction of interactions with nu_e_bar,\n"
                      "Fraction of interactions with nu_mu,\n"
                      "Fraction of interactions with nu_mu_bar,\n"
                      "Fraction of interactions with nu_tau,\n"
                      "Fraction of interactions with nu_tau_bar:"
               .format(input_name, now))


""" Save information from get_interaction_channel() into txt file: """
if SAVE_TXT:
    np.savetxt(output_path + "IBDlike_signals_NC_onlyC12_interaction_channels_{0:.0f}evts.txt"
               .format(number_ibdlike_events),
               np.array([Frac_C12_B11_p, Frac_C12_B11_n_piplus, Frac_C12_B11_n_piminus_2piplus,
                         Frac_C12_B11_p_piminus_piplus, Frac_C12_B11_p_2piminus_2piplus, Frac_C12_B11_piplus,
                         Frac_C12_B11_other,
                         Frac_C12_C11_n, Frac_C12_C11_p_piminus, Frac_C12_C11_n_piminus_piplus,
                         Frac_C12_C11_p_2piminus_piplus, Frac_C12_C11_p_3piminus_2piplus,
                         Frac_C12_C11_n_2piminus_2piplus, Frac_C12_C11_other,
                         Frac_C12_B10_p_n, Frac_C12_B10_2p_piminus, Frac_C12_B10_p_n_piminus_piplus,
                         Frac_C12_B10_2n_piplus, Frac_C12_B10_2n_piminus_2piplus, Frac_C12_B10_2p_2piminus_piplus,
                         Frac_C12_B10_2p_3piminus_2piplus, Frac_C12_B10_p_n_2piminus_2piplus, Frac_C12_B10_other,
                         Frac_C12_C10_2n, Frac_C12_C10_p_n_piminus, Frac_C12_C10_p_n_2piminus_piplus,
                         Frac_C12_C10_2n_piminus_piplus, Frac_C12_C10_2p_2piminus, Frac_C12_C10_other,
                         Frac_C12_Be10_2p, Frac_C12_Be10_p_n_piplus, Frac_C12_Be10_p_n_piminus_2piplus,
                         Frac_C12_Be10_2p_piminus_piplus, Frac_C12_Be10_2n_2piplus, Frac_C12_Be10_p_n_2piminus_3piplus,
                         Frac_C12_Be10_2p_2piminus_2piplus, Frac_C12_Be10_2p_3piminus_3piplus, Frac_C12_Be10_other,
                         Frac_C12_B9_p_2n, Frac_C12_B9_p_2n_piminus_piplus, Frac_C12_B9_2p_n_3piminus_2piplus,
                         Frac_C12_B9_2p_n_piminus, Frac_C12_B9_3n_piplus, Frac_C12_B9_p_2n_2piminus_2piplus,
                         Frac_C12_B9_2p_n_2piminus_piplus, Frac_C12_B9_other,
                         Frac_C12_Be9_2p_n, Frac_C12_Be9_p_2n_piplus, Frac_C12_Be9_3p_piminus,
                         Frac_C12_Be9_p_2n_piminus_2piplus, Frac_C12_Be9_2p_n_piminus_piplus,
                         Frac_C12_Be9_2p_n_3piminus_3piplus, Frac_C12_Be9_2p_n_2piminus_2piplus,
                         Frac_C12_Be9_3n_2piplus, Frac_C12_Be9_3p_2piminus_piplus, Frac_C12_Be9_other,
                         Frac_C12_Be8_2p_2n, Frac_C12_Be8_3p_n_piminus, Frac_C12_Be8_p_3n_piplus,
                         Frac_C12_Be8_2p_2n_2piminus_2piplus, Frac_C12_Be8_4n_2piplus,
                         Frac_C12_Be8_2p_2n_piminus_piplus, Frac_C12_Be8_3p_n_2piminus_piplus, Frac_C12_Be8_4p_2piminus,
                         Frac_C12_Be8_other,
                         Frac_C12_C9_p_2n_piminus, Frac_C12_C9_3n, Frac_C12_C9_2p_n_2piminus,
                         Frac_C12_C9_3n_2piminus_2piplus, Frac_C12_C9_other,
                         Frac_C12_Be7_2p_3n, Frac_C12_Be7_p_4n_piplus, Frac_C12_Be7_2p_3n_2piminus_2piplus,
                         Frac_C12_Be7_3p_2n_piminus, Frac_C12_Be7_4p_n_2piminus, Frac_C12_Be7_3p_2n_2piminus_piplus,
                         Frac_C12_Be7_other,
                         Frac_C12_Li6_3p_3n, Frac_C12_Li6_2p_4n_piplus, Frac_C12_Li6_5p_n_2piminus,
                         Frac_C12_Li6_2p_4n_piminus_2piplus, Frac_C12_Li6_4p_2n_piminus,
                         Frac_C12_Li6_3p_3n_piminus_piplus, Frac_C12_Li6_other,
                         Frac_C12_Li8_3p_n, Frac_C12_Li8_4p_piminus, Frac_C12_Li8_4p_2piminus_piplus,
                         Frac_C12_Li8_2p_2n_piplus, Frac_C12_Li8_3p_n_piminus_piplus, Frac_C12_Li8_other,
                         Frac_C12_Li7_2p_3n_piplus, Frac_C12_Li7_4p_n_piminus, Frac_C12_Li7_3p_2n,
                         Frac_C12_Li7_3p_2n_piminus_piplus, Frac_C12_Li7_4p_n_2piminus_piplus,
                         Frac_C12_Li7_2p_3n_piminus_2piplus, Frac_C12_Li7_other,
                         Frac_C12_B8_p_3n, Frac_C12_B8_p_3n_piminus_piplus, Frac_C12_B8_2p_2n_2piminus_piplus,
                         Frac_C12_B8_2p_2n_piminus, Frac_C12_B8_4n_piplus, Frac_C12_B8_other,
                         Frac_C12_Li9_2p_n_piplus, Frac_C12_Li9_3p, Frac_C12_Li9_3p_piminus_piplus,
                         Frac_C12_Li9_2p_n_piminus_2piplus, Frac_C12_Li9_p_2n_piminus_3piplus, Frac_C12_Li9_other,
                         Frac_C12_C8_4n, Frac_C12_C8_4n_other, Frac_C12_He8_4p, Frac_C12_He8_4p_other,
                         Frac_C12_B7_p_4n, Frac_C12_B7_p_4n_other, Frac_C12_He7_4p_n, Frac_C12_He7_4p_n_other,
                         Frac_C12_H7_5p, Frac_C12_H7_5p_other, Frac_C12_Be6_2p_4n, Frac_C12_Be6_2p_4n_other,
                         Frac_C12_He6_4p_2n, Frac_C12_He6_4p_2n_other, Frac_C12_H6_5p_n, Frac_C12_H6_5p_n_other,
                         Frac_C12_mass11u, Frac_C12_mass10u, Frac_C12_mass9u, Frac_C12_mass8u, Frac_C12_mass7u,
                         Frac_C12_mass6u, Frac_C12_mass5orless,
                         Frac_C12_C12, Frac_C12_NoIso, Frac_C12_NoIso_5p_6n,
                         Frac_no_C12, Frac_ES_proton_chID, Frac_ES_electron_chID, Frac_ES_O16_chID, Frac_ES_N14_chID,
                         Frac_ES_S32_chID, Frac_C12_faulty]), fmt='%4.5f',
               header="Information about the fractions of the different NC interaction channels on C12, that produce an "
                      "IBD-like signal in the detector\n"
                      "(input file: {0}, script: OLD_IBDlike_events_channels_v1.py ({1})):\n"
                      "Fraction of NC channel nu + C12 -> B11 + p,\n"
                      "Fraction of NC channel nu + C12 -> B11 + n + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> B11 + n + pi_minus + 2*pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> B11 + p + pi_minus + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> B11 + p + 2*pi_minus + 2*pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> B11 + pi_plus,\n"
                      "Fraction of other NC channels of nu + C12 -> B11 + ...,\n"
                      "Fraction of NC channel nu + C12 -> C11 + n,\n"
                      "Fraction of NC channel nu + C12 -> C11 + p + pi_minus,\n"
                      "Fraction of NC channel nu + C12 -> C11 + n + pi_minus + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> C11 + p + 2*pi_minus + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> C11 + p + 3*pi_minus + 2*pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> C11 + n + 2*pi_minus + 2*pi_plus,\n"
                      "Fraction of other NC channels of nu + C12 -> C11 + ...,\n"
                      "Fraction of NC channel nu + C12 -> B10 + p + n,\n"
                      "Fraction of NC channel nu + C12 -> B10 + 2p + pi_minus,\n"
                      "Fraction of NC channel nu + C12 -> B10 + p + n + pi_minus + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> B10 + 2n + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> B10 + 2n + pi_minus + 2*pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> B10 + 2p + 2*pi_minus + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> B10 + 2p + 3*pi_minus + 2*pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> B10 + p + n + 2*pi_minus + 2*pi_plus,\n"
                      "Fraction of other NC channels of nu + C12 -> B10 + ...,\n"
                      "Fraction of NC channel nu + C12 -> C10 + 2n,\n"
                      "Fraction of NC channel nu + C12 -> C10 + p + n + pi_minus,\n"
                      "Fraction of NC channel nu + C12 -> C10 + p + n + 2*pi_minus + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> C10 + 2n + pi_minus + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> C10 + 2p + 2*pi_minus,\n"
                      "Fraction of other NC channels of nu + C12 -> C10 + ...,\n"
                      "Fraction of NC channel nu + C12 -> Be10 + 2p,\n"
                      "Fraction of NC channel nu + C12 -> Be10 + p + n + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Be10 + p + n + pi_minus + 2*pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Be10 + 2p + pi_minus + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Be10 + 2n + 2*pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Be10 + p + n + 2*pi_minus + 3*pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Be10 + 2p + 2*pi_minus + 2*pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Be10 + 2p + 3*pi_minus + 3*pi_plus,\n"
                      "Fraction of other NC channels of nu + C12 -> Be10 + ...,\n"
                      "Fraction of NC channel nu + C12 -> B9 + p + 2n,\n"
                      "Fraction of NC channel nu + C12 -> B9 + p + 2n + pi_minus + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> B9 + 2p + n + 3*pi_minus + 2*pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> B9 + 2p + n + pi_minus,\n"
                      "Fraction of NC channel nu + C12 -> B9 + 3n + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> B9 + p + 2n + 2*pi_minus + 2*pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> B9 + 2p + n + 2*pi_minus + pi_plus,\n"
                      "Fraction of other NC channels of nu + C12 -> B9 + ...,\n"               
                      "Fraction of NC channel nu + C12 -> Be9 + 2p + n,\n"     
                      "Fraction of NC channel nu + C12 -> Be9 + p + 2n + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Be9 + 3p + pi_minus,\n"
                      "Fraction of NC channel nu + C12 -> Be9 + p + 2n + pi_minus + 2*pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Be9 + 2p + n + pi_minus + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Be9 + 2p + n + 3*pi_minus + 3*pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Be9 + 2p + n + 2*pi_minus + 2*pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Be9 + 3n + 2*pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Be9 + 3p + 2*pi_minus + pi_plus,\n"
                      "Fraction of other NC channels of nu + C12 -> Be9 + ..,\n"
                      "Fraction of NC channel nu + C12 -> Be8 + 2p + 2n,\n"
                      "Fraction of NC channel nu + C12 -> Be8 + 3p + n + pi_minus,\n"
                      "Fraction of NC channel nu + C12 -> Be8 + p + 3n + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Be8 + 2p + 2n + 2*pi_minus + 2*pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Be8 + 4n + 2*pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Be8 + 2p + 2n + pi_minus + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Be8 + 3p + n + 2*pi_minus + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Be8 + 4p + 2*pi_minus,\n"
                      "Fraction of other NC channels of nu + C12 -> Be8 + ...,\n"
                      "Fraction of NC channel nu + C12 -> C9 + p + 2n + pi_minus,\n"
                      "Fraction of NC channel nu + C12 -> C9 + 3n,\n"
                      "Fraction of NC channel nu + C12 -> C9 + 2p + n + 2*pi_minus,\n"
                      "Fraction of NC channel nu + C12 -> C9 + 3n + 2*pi_minus + 2*pi_plus,\n"
                      "Fraction of other NC channels of nu + C12 -> C9 + ...,\n"
                      "Fraction of NC channel nu + C12 -> Be7 + 2p + 3n,\n"
                      "Fraction of NC channel nu + C12 -> Be7 + p + 4n + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Be7 + 2p + 3n + 2*pi_minus + 2*pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Be7 + 3p + 2n + pi_minus,\n"
                      "Fraction of NC channel nu + C12 -> Be7 + 4p + n + 2*pi_minus,\n"
                      "Fraction of NC channel nu + C12 -> Be7 + 3p + 2n + 2*pi_minus + pi_plus,\n"
                      "Fraction of other NC channels of nu + C12 -> Be7 + ..,\n"
                      "Fraction of NC channel nu + C12 -> Li6 + 3p + 3n,\n"
                      "Fraction of NC channel nu + C12 -> Li6 + 2p + 4n + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Li6 + 5p + n + 2*pi_minus,\n"
                      "Fraction of NC channel nu + C12 -> Li6 + 2p + 4n + pi_minus + 2*pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Li6 + 4p + 2n + pi_minus,\n"
                      "Fraction of NC channel nu + C12 -> Li6 + 3p + 3n + pi_minus + pi_plus,\n"
                      "Fraction of other NC channels of nu + C12 -> Li6 + ...,\n"
                      "Fraction of NC channel nu + C12 -> Li8 + 3p + n,\n"
                      "Fraction of NC channel nu + C12 -> Li8 + 4p + pi_minus,\n"
                      "Fraction of NC channel nu + C12 -> Li8 + 4p + 2*pi_minus + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Li8 + 2p + 2n + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Li8 + 3p + n + pi_minus + pi_plus,\n"
                      "Fraction of other NC channels of nu + C12 -> Li8 + ...,\n"
                      "Fraction of NC channel nu + C12 -> Li7 + 2p + 3n + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Li7 + 4p + n + pi_minus,\n"
                      "Fraction of NC channel nu + C12 -> Li7 + 3p + 2n,\n"
                      "Fraction of NC channel nu + C12 -> Li7 + 3p + 2n + pi_minus + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Li7 + 4p + n + 2*pi_minus + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Li7 + 2p + 3n + pi_minus + 2*pi_plus,\n"
                      "Fraction of other NC channels of nu + C12 -> Li7 + ...,\n"
                      "Fraction of NC channel nu + C12 -> B8 + p + 3n,\n"
                      "Fraction of NC channel nu + C12 -> B8 + p + 3n + pi_minus + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> B8 + 2p + 2n + 2*pi_minus + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> B8 + 2p + 2n + pi_minus,\n"
                      "Fraction of NC channel nu + C12 -> B8 + 4n + pi_plus,\n"
                      "Fraction of other NC channels of nu + C12 -> B8 + ...,\n"
                      "Fraction of NC channel nu + C12 -> Li9 + 2p + n + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Li9 + 3p,\n"
                      "Fraction of NC channel nu + C12 -> Li9 + 3p + pi_minus + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Li9 + 2p + n + pi_minus + 2*pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Li9 + p + 2n + pi_minus + 3*pi_plus,\n"
                      "Fraction of other NC channels of nu + C12 -> Li9 + ...,\n"
                      "Fraction of NC channel nu + C12 -> C8 + 4n,\n"
                      "Fraction of NC channel nu + C12 -> C8 + 4n + N*pi_minus + N*pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> He8 + 4p,\n"
                      "Fraction of NC channel nu + C12 -> He8 + 4p + N*pi_minus + N*pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> B7 + p + 4n,\n"
                      "Fraction of NC channel nu + C12 -> B7 + p + 4n + N*pi_minus + N*pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> He7 + 4p + n,\n"
                      "Fraction of NC channel nu + C12 -> He7 + 4p + n + N*pi_minus + N*pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> H7 + 5p,\n"
                      "Fraction of NC channel nu + C12 -> H7 + 5p + N*pi_minus + N*pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Be6 + 2p + 4n,\n"
                      "Fraction of NC channel nu + C12 -> Be6 + 2p + 4n + N*pi_minus + N*pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> He6 + 4p + 2n,\n"
                      "Fraction of NC channel nu + C12 -> He6 + 4p + 2n + N*pi_minus + N*pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> H6 + 5p + n,\n"
                      "Fraction of NC channel nu + C12 -> H6 + 5p + n + N*pi_minus + N*pi_plus,\n"
                      "Fraction of other NC channels nu + C12 -> C11/B11 + ...\n"
                      "Fraction of other NC channels nu + C12 -> C10/B10/Be10 + ...\n"
                      "Fraction of other NC channels nu + C12 -> C9/B9/Be9/Li9 + ...\n"
                      "Fraction of other NC channels nu + C12 -> C8/B8/Be8/Li8/He8 + ...\n"
                      "Fraction of other NC channels nu + C12 -> B7/Be7/Li7/He7/H7 + ...\n"
                      "Fraction of other NC channels nu + C12 -> Be6/Li6/He6/H6 + ...\n"
                      "Fraction of NC channels with isotopes with mass <=5: nu + C12 -> X + ...\n"
                      "Fraction of NC channels nu + C12 -> nu + C12 + ...,\n"
                      "Fraction of NC channels without isotope nu + C12 -> nu + X*p + Y*n + Z*pion,\n"
                      "Fraction of NC channels without isotope nu + C12 -> nu + 5p + 6n + 1pi_plus,\n"
                      "Fraction of NC channels WITHOUT C12 as target,\n"
                      "Fraction of ES channel nu + p -> nu + p + ...,\n"
                      "Fraction of ES channel nu + electron -> nu + electron + ...,\n"
                      "Fraction of ES channel nu + O16 -> nu + O16 + ...,\n"
                      "Fraction of ES channel nu + N14 -> nu + N14 + ...,\n"
                      "Fraction of ES channel nu + S32 -> nu + S32 + ...,\n"
                      "Fraction of faulty interactions (isotopes or particles are missing):"
               .format(input_name, now))


""" Save information from get_deex_channel() into txt file: """
if SAVE_TXT:
    np.savetxt(output_path + "IBDlike_signals_deex_channels_NC_onlyC12_from_generator_{0:.0f}evts.txt"
               .format(number_ibdlike_events),
               np.array([number_Entries, number_target_C12, number_no_C12, number_light_Iso,
                         number_c11_notex, number_c11_deex, number_c11_li6_p_alpha, number_c11_be7_alpha,
                         number_c11_b10_p, number_c11_b9_n_p, number_c11_be8_p_d, number_c11_be9_2p, number_c11_b9_d,
                         number_c11_be8_he3, number_c11_c10_n, number_c11_li5_d_alpha, number_c11_li5_n_p_alpha,
                         number_c11_missing,
                         number_b11_notex, number_b11_deex, number_b11_li6_n_alpha, number_b11_b9_2n,
                         number_b11_be8_n_d, number_b11_be9_d, number_b11_be10_p, number_b11_b10_n, number_b11_be9_n_p,
                         number_b11_li7_alpha, number_b11_be8_t, number_b11_he5_d_alpha, number_b11_he6_p_alpha,
                         number_b11_be8_2n_p, number_b11_missing,
                         number_c10_notex, number_c10_deex, number_c10_b9_p, number_c10_be7_p_d, number_c10_li6_p_he3,
                         number_c10_he4_p_d_he3, number_c10_li6_2p_d, number_c10_be8_2p, number_c10_be7_n_2p,
                         number_c10_li6_n_3p, number_c10_be6_n_p_d, number_c10_li5_n_2p_d, number_c10_he3_p_d_alpha,
                         number_c10_li5_d_he3, number_c10_li5_p_2d, number_c10_he3_n_2p_alpha, number_c10_li4_n_p_alpha,
                         number_c10_b8_n_p, number_c10_b8_d, number_c10_be6_p_t, number_c10_he4_n_2p_he3,
                         number_c10_li5_n_p_he3, number_c10_missing,
                         number_b10_notex, number_b10_deex, number_b10_be9_p, number_b10_be8_d, number_b10_b9_n,
                         number_b10_be7_t, number_b10_be8_n_p, number_b10_li7_he3, number_b10_he5_p_alpha,
                         number_b10_li6_alpha, number_b10_li5_n_alpha, number_b10_li7_p_d, number_b10_missing,
                         number_be10_notex, number_be10_deex, number_be10_li6_2n_d, number_be10_li7_2n_p,
                         number_be10_he7_n_2p, number_be10_li6_n_t, number_be10_li6_3n_p, number_be10_he5_n_2d,
                         number_be10_be8_2n, number_be10_he5_n_alpha, number_be10_t_n_d_alpha, number_be10_li8_n_p,
                         number_be10_he4_n_d_t, number_be10_he5_n_p_t, number_be10_li7_n_d, number_be10_h4_n_p_alpha,
                         number_be10_he5_2n_p_d, number_be10_be9_n, number_be10_he6_n_p_d, number_be10_he5_d_t,
                         number_be10_h4_d_alpha, number_be10_he4_2n_p_t, number_be10_missing,
                         number_c9_notex, number_c9_deex, number_c9_be7_2p, number_c9_missing,
                         number_b9_notex, number_b9_deex, number_b9_be8_p, number_b9_missing,
                         number_be9_notex, number_be9_deex, number_be9_li8_p, number_be9_missing,
                         number_li9_notex, number_li9_deex, number_li9_h4_n_alpha, number_li9_he7_d, number_li9_li8_n,
                         number_li9_missing,
                         number_b8_notex, number_b8_deex, number_b8_li6_2p, number_b8_missing,
                         number_li8_notex, number_li8_deex, number_li8_li7_n, number_li8_li6_2n, number_li8_missing,
                         number_be7_notex, number_be7_deex, number_be7_li5_d, number_be7_li6_p, number_be7_missing,
                         number_li7_notex, number_li7_deex, number_li7_li6_n, number_li7_missing]), fmt='%4.5f',
               header="Information about the number of the different deexcitation channels, which occur in the DSNB-NC "
                      "generator (for NC interactions on C12), for events that mimic an IBD-like signal in JUNO "
                      "detector\n"
                      "(input file: {0}, script: OLD_IBDlike_events_channels_v1.py ({1})):\n"
                      "number of entries in the input file,\n"
                      "number of events with C12 as target,\n"
                      "number of events without C12 as target,\n"
                      "number of events with 'light' isotopes (isotopes, where no TALYS deex. root file exists),\n"
                      "C11 is not excited,\n"
                      "C11 is excited and de-excites,\n"
                      "C11* -> p + alpha + Li6,\n"
                      "C11* -> alpha + Be7,\n"
                      "C11* -> p + B10,\n"
                      "C11* -> n + p + B9,\n"
                      "C11* -> p + d + Be8,\n"
                      "C11* -> 2p + Be9,\n"
                      "C11* -> d + B9,\n"
                      "C11* -> He3 + Be8,\n"
                      "C11* -> n + C10,\n"
                      "C11* -> d + alpha + Li5,\n"
                      "C11* -> n + p + alpha + Li5,\n"
                      "deexcitations of C11 not yet included,\n"
                      "B11 is not excited,\n"
                      "B11 is excited and de-excites,\n"
                      "B11* -> n + alpha + Li6,\n"
                      "B11* -> 2n + B9,\n"
                      "B11* -> n + d + Be8,\n"
                      "B11* -> d + Be9,\n"
                      "B11* -> p + Be10,\n"
                      "B11* -> n + B10,\n"
                      "B11* -> n + p + Be9,\n"
                      "B11* -> alpha + Li7,\n"
                      "B11* -> t + Be8,\n"
                      "B11* -> d + alpha + He5,\n"
                      "B11* -> p + alpha + He6,\n"
                      "B11* -> 2n + p + Be8,\n"
                      "deexcitations of B11 not yet included,\n"
                      "C10 is not excited,\n"
                      "C10 is excited and de-excites,\n"
                      "C10* -> p + B9,\n"
                      "C10* -> p + d + Be7,\n"
                      "C10* -> p + He3 + Li6,\n"
                      "C10* -> p + d + He3 + He4,\n"
                      "C10* -> 2p + d + Li6,\n"
                      "C10* -> 2p + Be8,\n"
                      "C10* -> n + 2p + Be7,\n"
                      "C10* -> n + 3p + Li6,\n"
                      "C10* -> n + p + d + Be6,\n"
                      "C10* -> n + 2p + d + Li5,\n"
                      "C10* -> p + d + alpha + He3,\n"
                      "C10* -> d + He3 + Li5,\n"
                      "C10* -> p + 2d + Li5,\n"
                      "C10* -> n + 2p + alpha + He3,\n"
                      "C10* -> n + p + alpha + Li4,\n"
                      "C10* -> n + p + B8,\n"
                      "C10* -> d + B8,\n"
                      "C10* -> p + t + Be6,\n"
                      "C10* -> n + 2p + He3 + He4,\n"
                      "C10* -> n + p + He3 + Li5,\n"
                      "deexcitations of C10 not yet included,\n"
                      "B10 is not excited,\n"
                      "B10 is excited and de-excites,\n"
                      "B10* -> p + Be9,\n"
                      "B10* -> d + Be8,\n"
                      "B10* -> n + B9,\n"
                      "B10* -> t + Be7,\n"
                      "B10* -> n + p + Be8,\n"
                      "B10* -> He3 + Li7,\n"
                      "B10* -> p + alpha + He5,\n"
                      "B10* -> alpha + Li6,\n"
                      "B10* -> n + alpha + Li5,\n"
                      "B10* -> p + d + Li7,\n"
                      "deexcitations of B10 not yet included,\n"
                      "Be10 is not excited,\n"
                      "Be10 is excited and de-excites,\n"
                      "Be10* -> 2n + d + Li6,\n"
                      "Be10* -> 2n + p + Li7,\n"
                      "Be10* -> n + 2p + He7,\n"
                      "Be10* -> n + t + Li6,\n"
                      "Be10* -> 3n + p + Li6,\n"
                      "Be10* -> n + 2d + He5,\n"
                      "Be10* -> 2n + Be8,\n"
                      "Be10* -> n + alpha + He5,\n"
                      "Be10* -> n + d + alpha + tritium,\n"
                      "Be10* -> n + p + Li8,\n"
                      "Be10* -> n + d + t + He4,\n"
                      "Be10* -> n + p + t + He5,\n"
                      "Be10* -> n + d + Li7,\n"
                      "Be10* -> n + p + alpha + H4,\n"
                      "Be10* -> 2n + p + d + He5,\n"
                      "Be10* -> n + Be9,\n"
                      "Be10* -> n + p + d + He6,\n"
                      "Be10* -> d + t + He5,\n"
                      "Be10* -> d + alpha + H4,\n"
                      "Be10* -> 2n + p + t + He4,\n"
                      "deexcitations of Be10 not yet included,\n"
                      "C9 is not excited,\n"
                      "C9 is excited an de-excites,\n"
                      "C9* -> 2p + Be7,\n"
                      "deexcitations of C9 not yet included,\n"
                      "B9 is not excited,\n"
                      "B9 is excited and de-excites,\n"
                      "B9* -> p + Be8,\n"
                      "deexcitations of B9 not yet included,\n"
                      "Be9 is not excited,\n"
                      "Be9 is excite and de-excites,\n"
                      "Be9* -> p + Li8,\n"
                      "deexcitations of Be9 not yet included,\n"
                      "Li9 is not excited,\n"
                      "Li9 is excited and de-excites,\n"
                      "Li9* -> n + alpha + H4,\n"
                      "Li9* -> d + He7,\n"
                      "Li9* -> n + Li8,\n"
                      "deexcitations of Li9 not yet included,\n"
                      "B8 is not excited,\n"
                      "B8 is excited and de-excites,\n"
                      "B8* -> 2p + Li6,\n"
                      "deexcitations of B8 not yet included,\n"
                      "Li8 is not excited,\n"
                      "Li8 is excited and de-excites,\n"
                      "Li8* -> n + Li7,\n"
                      "Li8* -> 2n + Li6,\n"
                      "deexcitations of Li8 not yet included,\n"
                      "Be7 is not excited,\n"
                      "Be7 is excited and de-excites,\n"
                      "Be7* -> d + Li5,\n"
                      "Be7* -> p + Li6,\n"
                      "deexcitations of Be7 not yet included,\n"
                      "Li7 is not excited,\n"
                      "Li7 is excited and de-excites,\n"
                      "Li7* -> n + Li6,\n"
                      "deexcitations of Li7 not yet included:"
               .format(input_name, now))

""" Display Output of 'get_residual_isotopes_before_deex()' in plot: """
h3 = plt.figure(3, figsize=(15, 8))

plt.plot([], [], color='w', label="total number of entries = {0:d}".format(number_ibdlike_events))
# do not display the histogram, when fraction is to small:
if Fraction_c12 > 2.0:
    plt.plot(Energy_nu_incoming_2[:-1], Events_c12, drawstyle='steps', color='m',
             label="$^{12}$C:"+" fraction = {0:.2f}%".format(Fraction_c12))
if Fraction_c11 > 2.0:
    plt.plot(Energy_nu_incoming_2[:-1], Events_c11, drawstyle='steps', color='b', linestyle='solid',
             label="$^{11}$C:"+" fraction = {0:.2f}%".format(Fraction_c11))
if Fraction_b11 > 2.0:
    plt.plot(Energy_nu_incoming_2[:-1], Events_b11, drawstyle='steps', color='r', linestyle='solid',
             label="$^{11}$B:"+" fraction = {0:.2f}%".format(Fraction_b11))
if Fraction_c10 > 2.0:
    plt.plot(Energy_nu_incoming_2[:-1], Events_c10, drawstyle='steps', color='b', linestyle='dashed',
             label="$^{10}$C:"+" fraction = {0:.2f}%".format(Fraction_c10))
if Fraction_b10 > 2.0:
    plt.plot(Energy_nu_incoming_2[:-1], Events_b10, drawstyle='steps', color='r', linestyle='dashed',
             label="$^{10}$B:"+" fraction = {0:.2f}%".format(Fraction_b10))
if Fraction_be10 > 2.0:
    plt.plot(Energy_nu_incoming_2[:-1], Events_be10, drawstyle='steps', color='g', linestyle='dashed',
             label="$^{10}$Be:"+" fraction = {0:.2f}%".format(Fraction_be10))
if Fraction_c9 > 2.0:
    plt.plot(Energy_nu_incoming_2[:-1], Events_c9, drawstyle='steps', color='b', linestyle='dotted',
             label="$^{9}$C:"+" fraction = {0:.2f}%".format(Fraction_c9))
if Fraction_b9 > 2.0:
    plt.plot(Energy_nu_incoming_2[:-1], Events_b9, drawstyle='steps', color='r', linestyle='dotted',
             label="$^{9}$B:"+" fraction = {0:.2f}%".format(Fraction_b9))
if Fraction_be9 > 2.0:
    plt.plot(Energy_nu_incoming_2[:-1], Events_be9, drawstyle='steps', color='g', linestyle='dotted',
             label="$^{9}$Be:"+" fraction = {0:.2f}%".format(Fraction_be9))
if Fraction_li9 > 2.0:
    plt.plot(Energy_nu_incoming_2[:-1], Events_li9, drawstyle='steps', color='k', linestyle='dotted',
             label="$^{9}$Li:"+" fraction = {0:.2f}%".format(Fraction_li9))
if Fraction_c8 > 2.0:
    plt.plot(Energy_nu_incoming_2[:-1], Events_c8, drawstyle='steps', color='m',
             label="$^{8}$C:"+" fraction = {0:.2f}%".format(Fraction_c8))
if Fraction_b8 > 2.0:
    plt.plot(Energy_nu_incoming_2[:-1], Events_b8, drawstyle='steps', color='m',
             label="$^{8}$B:"+" fraction = {0:.2f}%".format(Fraction_b8))
if Fraction_be8 > 2.0:
    plt.plot(Energy_nu_incoming_2[:-1], Events_be8, drawstyle='steps', color='m',
             label="$^{8}$Be:"+" fraction = {0:.2f}%".format(Fraction_be8))
if Fraction_li8 > 2.0:
    plt.plot(Energy_nu_incoming_2[:-1], Events_li8, drawstyle='steps', color='m',
             label="$^{8}$Li:"+" fraction = {0:.2f}%".format(Fraction_li8))
if Fraction_he8 > 2.0:
    plt.plot(Energy_nu_incoming_2[:-1], Events_he8, drawstyle='steps', color='m',
             label="$^{8}$He:"+" fraction = {0:.2f}%".format(Fraction_he8))
if Fraction_b7 > 2.0:
    plt.plot(Energy_nu_incoming_2[:-1], Events_b7, drawstyle='steps', color='r', linestyle='dashdot',
             label="$^{7}$B:"+" fraction = {0:.2f}%".format(Fraction_b7))
if Fraction_be7 > 2.0:
    plt.plot(Energy_nu_incoming_2[:-1], Events_be7, drawstyle='steps', color='g', linestyle='dashdot',
             label="$^{7}$Be:"+" fraction = {0:.2f}%".format(Fraction_be7))
if Fraction_li7 > 2.0:
    plt.plot(Energy_nu_incoming_2[:-1], Events_li7, drawstyle='steps', color='k', linestyle='dashdot',
             label="$^{7}$Li:"+" fraction = {0:.2f}%".format(Fraction_li7))
if Fraction_he7 > 2.0:
    plt.plot(Energy_nu_incoming_2[:-1], Events_he7, drawstyle='steps', color='m', linestyle='dashdot',
             label="$^{7}$He:"+" fraction = {0:.2f}%".format(Fraction_he7))
if Fraction_h7 > 2.0:
    plt.plot(Energy_nu_incoming_2[:-1], Events_h7, drawstyle='steps', color='c', linestyle='dashdot',
             label="$^{7}$H:"+" fraction = {0:.2f}%".format(Fraction_h7))
if Fraction_be6 > 2.0:
    plt.plot(Energy_nu_incoming_2[:-1], Events_be6, drawstyle='steps', color='m',
             label="$^{6}$Be:"+" fraction = {0:.2f}%".format(Fraction_be6))
if Fraction_li6 > 2.0:
    plt.plot(Energy_nu_incoming_2[:-1], Events_li6, drawstyle='steps', color='m',
             label="$^{6}$Li:"+" fraction = {0:.2f}%".format(Fraction_li6))
if Fraction_he6 > 2.0:
    plt.plot(Energy_nu_incoming_2[:-1], Events_he6, drawstyle='steps', color='m',
             label="$^{6}$He:"+" fraction = {0:.2f}%".format(Fraction_he6))
if Fraction_h6 > 2.0:
    plt.plot(Energy_nu_incoming_2[:-1], Events_h6, drawstyle='steps', color='m',
             label="$^{6}$H:"+" fraction = {0:.2f}%".format(Fraction_h6))
plt.xlim(xmin=0, xmax=10.0)
plt.ylim(ymin=0)
plt.xlabel("Neutrino energy $E_{\\nu}$ in GeV", fontsize=15)
plt.ylabel("events per bin (bin-width = {0:.2f} GeV)".format(bin_width_incoming), fontsize=15)
plt.title("Neutrino energy spectrum for NC interactions on $^{12}C$, that mimic IBD-like signal, "
          "\nfor different residual isotopes "
          "($\\nu_x$+$^{12}C$ $\\rightarrow$ $\\nu_x$ + ...)", fontsize=20)
plt.legend(fontsize=12)

if SAVE_FIG:
    plt.savefig(output_path + "IBDlike_signals_residual_isotopes_before_deex_{0:.0f}evts.png"
                .format(number_ibdlike_events))

""" Save information about residual isotopes before deexcitation into txt file: """
if SAVE_TXT:
    np.savetxt(output_path + "IBDlike_signals_residual_isotopes_before_deex_{0:.0f}evts.txt"
               .format(number_ibdlike_events),
               np.array([number_ibdlike_events, Fraction_c12, Fraction_c11, Fraction_b11, Fraction_c10, Fraction_b10,
                         Fraction_be10,
                         Fraction_c9, Fraction_b9, Fraction_be9, Fraction_li9, Fraction_c8, Fraction_b8, Fraction_be8,
                         Fraction_li8, Fraction_he8, Fraction_b7, Fraction_be7, Fraction_li7, Fraction_he7,
                         Fraction_h7, Fraction_be6, Fraction_li6, Fraction_he6, Fraction_h6, Fraction_rest]),
               fmt='%4.5f',
               header="Information about the fraction (in %) of residual isotopes before deexcitation, which are"
                      " produced, when neutrinos interact via NC on C12 in the JUNO liquid scintillator\n,"
                      "but only events, that mimic an IBD-like signal."
                      "(input file: {0}, script: OLD_IBDlike_events_channels_v1.py ({1})):\n"
                      "Number of events in the input file,\n"
                      "Fraction of C12 in %,\n"
                      "Fraction of C11 in %,\n"
                      "Fraction of B11 in %,\n"
                      "Fraction of C10 in %,\n"
                      "Fraction of B10 in %,\n"
                      "Fraction of Be10 in %,\n"
                      "Fraction of C9 in %,\n"
                      "Fraction of B9 in %,\n"
                      "Fraction of Be9 in %,\n"
                      "Fraction of Li9 in %,\n"
                      "Fraction of C8 in %,\n"
                      "Fraction of B8 in %,\n"
                      "Fraction of Be8 in %,\n"
                      "Fraction of Li8 in %,\n"
                      "Fraction of He8 in %,\n"
                      "Fraction of B7 in %,\n"
                      "Fraction of Be7 in %,\n"
                      "Fraction of Li7 in %,\n"
                      "Fraction of He7 in %,\n"
                      "Fraction of H7 in %,\n"
                      "Fraction of Be6 in %,\n"
                      "Fraction of Li6 in %,\n"
                      "Fraction of He6 in %,\n"
                      "Fraction of H6 in %,\n"
                      "Fraction of isotopes not yet included:"
               .format(input_name, now))

""" Save information from get_combined_channel() into txt file: """
if SAVE_TXT:
    np.savetxt(output_path + "IBDlike_signals_combined_channels_NC_onlyC12_from_generator_{0:.0f}evts.txt"
               .format(number_ibdlike_events),
               np.array([number_ibdlike_events, num_c12, num_no_c12,
                         num_b11, num_b11_n_piplus, num_b11_p,
                         num_c11, num_c11_n, num_c11_p_piminus,
                         num_c10, num_c10_2n, num_c10_n_p_piminus,
                         num_b10, num_b10_n_p,
                         num_be10, num_be10_2p, num_be10_n_p_piplus,
                         num_c9, num_c9_3n, num_c9_2n_p_piminus,
                         num_b9, num_b9_n_d, num_b9_2n_p, num_b9_n_2p_piminus, num_b9_3n_piplus,
                         num_be9, num_be9_n_2p, num_be9_p_d, num_be9_2n_p_piplus,
                         num_li9, num_li9_3p, num_li9_n_2p_piplus,
                         num_b8, num_b8_3n_p, num_b8_2n_d,
                         num_be8, num_be8_n_p_d, num_be8_2n_2p, num_be8_n_3p_piminus, num_be8_n_he3, num_be8_p_t,
                         num_li8, num_li8_n_3p,
                         num_he8, num_he8_4p,
                         num_b7,
                         num_be7, num_be7_n_p_t, num_be7_3n_2p, num_be7_n_alpha, num_be7_2n_p_d,
                         num_li7, num_li7_n_p_he3, num_li7_2n_3p, num_li7_n_2p_d, num_li7_p_alpha,
                         num_li7_n_alpha_piplus,
                         num_he7, num_he7_n_4p,
                         num_be6, num_be6_3n_p_d, num_be6_2n_2d, num_be6_4n_2p, num_be6_2n_p_t,
                         num_li6, num_li6_n_p_alpha, num_li6_2n_2p_d, num_li6_n_2p_t, num_li6_3n_3p, num_li6_2n_p_he3,
                         num_li6_2n_alpha_piplus, num_li6_2n_4p_piminus, num_li6_n_2p_he3_piminus,
                         num_he6, num_he6_2n_4p, num_he6_n_3p_d, num_he6_n_2p_he3,
                         num_h6, num_h6_n_5p,
                         num_li5, num_li5_2n_p_2d, num_li5_2n_p_alpha, num_li5_n_alpha_d, num_li5_3n_p_he3,
                         num_li5_3n_2p_d, num_li5_2n_he3_d,
                         num_he5, num_he5_n_3p_t, num_he5_n_2p_alpha, num_he5_n_2p_2d, num_he5_3n_4p,
                         num_h5, num_h5_2n_5p,
                         num_li4, num_li4_3n_p_alpha,
                         num_he4, num_he4_n_2p_t_d, num_he4_3n_2p_he3, num_he4_n_alpha_he3, num_he4_2n_p_he3_d,
                         num_he4_4n_4p, num_he4_2n_2he3,
                         num_h4, num_h4_n_3p_alpha, num_h4_3n_5p,
                         number_false_channel]), fmt='%4.5f',
               header="Information about the number of the different combined channels (NC interaction + deexcitation) "
                      "(final particles after NC interaction and after deexcitation), for events that mimic an "
                      "IBD-like signal in JUNO detector\n"
                      "(input file: {0}, script: OLD_IBDlike_events_channels_v1.py ({1})):\n"
                      "number of entries in the input file,\n"
                      "number of events with C12 as target,\n"
                      "number of events without C12 as target,\n"
                      "num_b11,\n"
                      "num_b11_n_piplus,\n"
                      "num_b11_p,\n"
                      "num_c11,\n"
                      "num_c11_n,\n"
                      "num_c11_p_piminus,\n"
                      "num_c10,\n"
                      "num_c10_2n,\n"
                      "num_c10_n_p_piminus,\n"
                      "num_b10,\n"
                      "num_b10_n_p,\n"
                      "num_be10,\n"
                      "num_be10_2p,\n"
                      "num_be10_n_p_piplus,\n"
                      "num_c9,\n"
                      "num_c9_3n,\n"
                      "num_c9_2n_p_piminus,\n"
                      "num_b9,\n"
                      "num_b9_n_d,\n"
                      "num_b9_2n_p,\n"
                      "num_b9_n_2p_piminus,\n"
                      "num_b9_3n_piplus,\n"
                      "num_be9,\n"
                      "num_be9_n_2p,\n"
                      "num_be9_p_d,\n"
                      "num_be9_2n_p_piplus,\n"
                      "num_li9,\n"
                      "num_li9_3p,\n"
                      "num_li9_n_2p_piplus,\n"
                      "num_b8,\n"
                      "num_b8_3n_p,\n"
                      "num_b8_2n_d,\n"
                      "num_be8,\n"
                      "num_be8_n_p_d,\n"
                      "num_be8_2n_2p,\n"
                      "num_be8_n_3p_piminus,\n"
                      "num_be8_n_he3,\n"
                      "num_be8_p_t,\n"
                      "num_li8,\n"
                      "num_li8_n_3p,\n"
                      "num_he8,\n"
                      "num_he8_4p,\n"
                      "num_b7,\n"
                      "num_be7,\n"
                      "num_be7_n_p_t,\n"
                      "num_be7_3n_2p,\n"
                      "num_be7_n_alpha,\n"
                      "num_be7_2n_p_d,\n"
                      "num_li7,\n"
                      "num_li7_n_p_he3,\n"
                      "num_li7_2n_3p,\n"
                      "num_li7_n_2p_d,\n"
                      "num_li7_p_alpha,\n"
                      "num_li7_n_alpha_piplus,\n"
                      "num_he7,\n"
                      "num_he7_n_4p,\n"
                      "num_be6,\n"
                      "num_be6_3n_p_d,\n"
                      "num_be6_2n_2d,\n"
                      "num_be6_4n_2p,\n"
                      "num_be6_2n_p_t,\n"
                      "num_li6,\n"
                      "num_li6_n_p_alpha,\n"
                      "num_li6_2n_2p_d,\n"
                      "num_li6_n_2p_t,\n"
                      "num_li6_3n_3p,\n"
                      "num_li6_2n_p_he3,\n"
                      "num_li6_2n_alpha_piplus,\n"
                      "num_li6_2n_4p_piminus,\n"
                      "num_li6_n_2p_he3_piminus,\n"
                      "num_he6,\n"
                      "num_he6_2n_4p,\n"
                      "num_he6_n_3p_d,\n"
                      "num_he6_n_2p_he3,\n"
                      "num_h6,\n"
                      "num_h6_n_5p,\n"
                      "num_li5,\n"
                      "num_li5_2n_p_2d,\n"
                      "num_li5_2n_p_alpha,\n"
                      "num_li5_n_alpha_d,\n"
                      "num_li5_3n_p_he3,\n"
                      "num_li5_3n_2p_d,\n"
                      "num_li5_2n_he3_d,\n"
                      "num_he5,\n"
                      "num_he5_n_3p_t,\n"
                      "num_he5_n_2p_alpha,\n"
                      "num_he5_n_2p_2d,\n"
                      "num_he5_3n_4p,\n"
                      "num_h5,\n"
                      "num_h5_2n_5p,\n"
                      "num_li4,\n"
                      "num_li4_3n_p_alpha,\n"
                      "num_he4,\n"
                      "num_he4_n_2p_t_d,\n"
                      "num_he4_3n_2p_he3,\n"
                      "num_he4_n_alpha_he3,\n"
                      "num_he4_2n_p_he3_d,\n"
                      "num_he4_4n_4p,\n"
                      "num_he4_2n_2he3,\n"
                      "num_h4,\n"
                      "num_h4_n_3p_alpha,\n"
                      "num_h4_3n_5p,\n"
                      "number event without isotope (only 3n, 4n or 5n as 'isotope'):"
               .format(input_name, now))


