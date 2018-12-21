""" Script to check out the results from the DSNB-NC generator of Jie Cheng with the original GENIE file of Julia.

    The ROOT-file generated with DSNB-NC.exe is read and analyzed with this script.

    DSNB-NC.exe is the generator of the Neutral Current background from atmospheric neutrinos and antineutrinos
    (build by Jie Cheng).

    The interactions of the neutrino and antineutrinos with the liquid scintillator are simulated with GENIE by Julia.
    For the flux of the atmospheric neutrinos (nu_e, nu_e_bar, nu_mu, nu_mu_bar) the flux of Honda for JUNO site and
    solar maximum is used (same flux like in my simulation of atmospheric CC background but for solar maximum and not
    for mean of solar maximum and minimum; energy range from 0 GeV to 10 GeV).

    Is it really from 0 to 10 GeV OR from 100 MeV to 10 GeV????

    For the interactions with the LS, the interactions on C12, protons, electron, N14, O16 and S32 are considered
    (composition of LS: C12=0.87924, proton=0.1201, O16=0.00034, N14=0.00027, S32=5e-5).

    The output of the GENIE simulation is saved in "genie_data.root". Only QEL AND NC events are taken into account
    (but all targets)

    Then the deexcitation of the products from the NC interactions is simulated:
    The status of the residual nucleus is unknown. Therefore the deexcitation of the residual nucleus is assumed from
    a advanced shell model (modified DSNB-NC generator -> see info_atmoNC.odt)

    The deexitation of the residual nucleus is simulated with the TALYS software and saved in "*deexcitation.root" for
    all nuclei (Li7, Li8, Li9, C9, C10, C11, Be7, Be9, Be10, B8, B9, B10, B11)

    With these root-files, the deexcitation is calculated and all final particles with their PDG ID, momentum and
    mass printed in the terminal.
    These information together with the target ID, the channel ID, the deexcitation ID, the isotope ID and the energy
    of the incoming neutrino are saved in the root output file.

"""

# import ROOT
import datetime
# import glob
import NC_background_functions
import numpy as np
from matplotlib import pyplot as plt

# get the date and time, when the script was run:
date = datetime.datetime.now()
now = date.strftime("%Y-%m-%d %H:%M")


# set SAVE_FIG, defines if figures are saved:
SAVE_FIG = True

# set SAVE_TXT, defines if txt files are saved:
SAVE_TXT = True

# set SHOW_PLOT, defines if the figures are shown:
SHOW_PLOT = True


# set the path of the inputs:
# input_path = "/home/astro/blum/juno/test_output_DSNB_gen/input/"
input_path = "/home/astro/blum/juno/atmoNC/data_NC/output_generator/"

# file name of the input file:
input_name = input_path + "gen_NC_onlyC12_250000evts_seed1.root"
# input_name = input_path + "test_11330events_newVersion_seed11330.root"

# set the path, where the outputs are saved:
output_path = "/home/astro/blum/juno/atmoNC/data_NC/output_checkoutNC/"
# output_path = "/home/astro/blum/juno/test_output_DSNB_gen/output/"


# bin-width of the array, which represents the incoming neutrino energy (in GeV) (float):
bin_width_incoming = 0.1


# set the NC interaction event rate for GENIE in units of events/(s*20ktons) (see NC_generator.pdf from Jie Cheng)
# (float):
# TODO: include the event rate to get "real" spectra
# TODO-me: Get the Event Rate for GENIE simulations of Julia!!!
# evt_rate_Genie = 3.59E-5
evt_rate_Genie = 1

# total exposure time in years (float):
t_years = 10
# total time-exposure in seconds (1yr = 365.2425 * 24 * 60 * 60 sec = 3.1556952 * 10^7 sec), 10 years (float):
# TODO: include the total exposure time to get "real" spectra
# time = t_years * 3.156 * 10 ** 7
time = 1


# read NC generator data to arrays:
(event_ID, projectile_PDG, projectile_E, target_PDG, NC_inter_ch_ID, deexcitation_ID, isotope_PDG, Nparticles,
 final_PDG, final_Px, final_Py, final_Pz) = NC_background_functions.read_nc_data(input_name)


# get the number of events in the root file:
# from the maximum of event_ID +1 (float):
evtnumber_fromMax = int(np.max(event_ID) + 1)
# from the length of the array:
evtnumber_fromLength = len(event_ID)
if evtnumber_fromMax != evtnumber_fromLength:
    print("WARNING: different values for max(event_ID) and len(event_ID)!!!")
else:
    NumEvent = evtnumber_fromLength


# get the number of event as function of the energy of the incoming neutrinos for each neutrino type:
(Energy_nu_incoming,
 Event_nu_e_IN, Event_nu_e_bar_IN, Event_nu_mu_IN, Event_nu_mu_bar_IN, Event_nu_tau_IN, Event_nu_tau_bar_IN,
 Number_nu_e_IN, Number_nu_e_bar_IN, Number_nu_mu_IN, Number_nu_mu_bar_IN, Number_nu_tau_IN, Number_nu_tau_bar_IN,
 Frac_nu_e_IN, Frac_nu_e_bar_IN, Frac_nu_mu_IN, Frac_nu_mu_bar_IN, Frac_nu_tau_IN, Frac_nu_tau_bar_IN) \
    = NC_background_functions.get_neutrino_energy(projectile_PDG, projectile_E, bin_width_incoming, evt_rate_Genie,
                                                  time)


# get the number of events as function of the energy of the incoming neutrino for the two types of target
# (C12 and proton):
(Energy_nu_incoming_1, Event_NC_C12, Event_ES_proton, Event_NC_N14, Event_NC_O16, Event_ES_electron, Event_NC_S32,
 Number_NC_C12, Number_ES_proton, Number_NC_N14, Number_NC_O16, Number_ES_electron, Number_NC_S32,
 Frac_NC_C12, Frac_ES_proton_target, Frac_NC_N14_target, Frac_NC_O16_target, Frac_ES_electron_target,
 Frac_NC_S32_target) \
    = NC_background_functions.get_target_ratio(projectile_E, target_PDG, bin_width_incoming)


# get the fractions of events of the different types on NC interaction channels:
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


# get the numbers of events of the different deexcitation channels of the different produced isotopes:
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
             label="interactions with $\\nu_\mu$: $N_{events} = $"+"{0:.2f}; fraction = {1:.2f}%"
             .format(Number_nu_mu_IN, Frac_nu_mu_IN))
if Number_nu_mu_bar_IN != 0:
    plt.plot(Energy_nu_incoming[:-1], Event_nu_mu_bar_IN, drawstyle='steps', color='r', linestyle='--',
             label="interactions with $\\bar{\\nu}_\mu$: $N_{events} = $"+"{0:.2f}; fraction = {1:.2f}%"
             .format(Number_nu_mu_bar_IN, Frac_nu_mu_bar_IN))
if Number_nu_tau_IN != 0:
    plt.plot(Energy_nu_incoming[:-1], Event_nu_tau_IN, drawstyle='steps', color='g',
             label="interactions with $\\nu_\\tau$: $N_{events} = $"+"{0:.2f}; fraction = {1:.2f}%"
             .format(Number_nu_tau_IN, Frac_nu_tau_IN))
if Number_nu_tau_bar_IN != 0:
    plt.plot(Energy_nu_incoming[:-1], Event_nu_tau_bar_IN, drawstyle='steps', color='g', linestyle='--',
             label="interactions with $\\bar{\\nu}_\\tau$: $N_{events} = $"+"{0:.2f}; fraction = {1:.2f}%"
             .format(Number_nu_tau_bar_IN, Frac_nu_tau_bar_IN))

plt.xlim(xmin=0)
plt.ylim(ymin=0)
plt.xlabel("Neutrino energy $E_{\\nu}$ in GeV", fontsize=15)
plt.ylabel("events", fontsize=15)
plt.title("Energy spectrum of atmospheric neutrinos (interacting via NC on $^{12}C$ in the JUNO detector)", fontsize=20)
plt.legend(fontsize=12)

if SAVE_FIG:
    plt.savefig(output_path + "incoming_neutrino_spectrum_{0:.0f}evts.png".format(NumEvent))

""" Save information about the ratio of incoming neutrinos into txt file: """
if SAVE_TXT:
    np.savetxt(output_path + "incoming_neutrino_spectrum_{0:.0f}evts.txt".format(NumEvent),
               np.array([NumEvent, Number_nu_e_IN, Number_nu_e_bar_IN, Number_nu_mu_IN, Number_nu_mu_bar_IN,
                         Number_nu_tau_IN, Number_nu_tau_bar_IN, Frac_nu_e_IN, Frac_nu_e_bar_IN, Frac_nu_mu_IN,
                         Frac_nu_mu_bar_IN, Frac_nu_tau_IN, Frac_nu_tau_bar_IN]), fmt='%4.5f',
               header="Information about the ratio of the incoming atmospheric neutrinos that interact with the\n"
                      "JUNO liquid scintillator via NC interactions on C12\n"
                      "(input file: {0}, script: checkout_NCgen.py ({1})):\n"
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


""" Display Output of 'get_target_ratio()' in plot: """
h2 = plt.figure(2, figsize=(15, 8))

# do not display the histogram, when there are no events:
if Number_NC_C12 != 0:
    plt.plot(Energy_nu_incoming_1[:-1], Event_NC_C12, drawstyle='steps', color='b',
             label="interactions on $^{12}$C: $N_{events} = $"+"{0:.2f}; fraction = {1:.2f}%"
             .format(Number_NC_C12, Frac_NC_C12))
if Number_ES_proton != 0:
    plt.plot(Energy_nu_incoming_1[:-1], Event_ES_proton, drawstyle='steps', color='r',
             label="interaction on protons: $N_{events} = $"+"{0:.2f}; fraction = {1:.2f}%"
             .format(Number_ES_proton, Frac_ES_proton_target))

plt.xlim(xmin=0)
plt.ylim(ymin=0)
plt.xlabel("Neutrino energy $E_{\\nu}$ in GeV", fontsize=15)
plt.ylabel("events", fontsize=15)
plt.title("Neutrino energy spectrum for NC interactions on $^{12}C$", fontsize=20)
plt.legend(fontsize=12)

if SAVE_FIG:
    plt.savefig(output_path + "target_ratio_{0:.0f}evts.png".format(NumEvent))

""" Save information about the target ratio into txt file: """
if SAVE_TXT:
    np.savetxt(output_path + "target_ratio_{0:.0f}evts.txt".format(NumEvent),
               np.array([NumEvent, Number_NC_C12, Number_ES_proton, Number_NC_N14, Number_NC_O16, Number_ES_electron,
                         Number_NC_S32, Frac_NC_C12, Frac_ES_proton_target, Frac_NC_N14_target, Frac_NC_O16_target,
                         Frac_ES_electron_target, Frac_NC_S32_target]), fmt='%4.5f',
               header="Information about the number and ratio of the target particles, where neutrinos interact via "
                      "NC on C12, in the JUNO liquid scintillator\n"
                      "(input file: {0}, script: checkout_NCgen.py ({1})):\n"
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
               .format(input_name, now))


""" Save information from get_interaction_channel() into txt file: """
if SAVE_TXT:
    np.savetxt(output_path + "NC_onlyC12_interaction_channels_{0:.0f}evts.txt".format(NumEvent),
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
               header="Information about the fractions of the different NC interaction channels on C12\n"
                      "(input file: {0}, script: checkout_NCgen.py ({1})):\n"
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
    np.savetxt(output_path + "deex_channels_NC_onlyC12_from_generator_{0:.0f}evts.txt".format(NumEvent),
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
                      "generator (for NC interactions on C12)\n"
                      "(input file: {0}, script: checkout_NCgen.py ({1})):\n"
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



if SHOW_PLOT:
    plt.show()
