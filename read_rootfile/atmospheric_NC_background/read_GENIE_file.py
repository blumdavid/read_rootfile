""" Script to read the 'original' ROOT file from Julia's GENIE simulation to investigate the different NC interaction
    channels.

    The output txt file with the fraction of the different interaction channels is saved in
    "/home/astro/blum/juno/atmoNC/interaction_channels_qel_NC_{0:d}evts.txt".

    In "interaction_channels_edit.ods" the channels are sorted descending.

    In python script "display_interaction_channels.py", the main interaction channels are displayed with their
    fractions.

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

# set SAVE_TXT, defines if txt files are saved:
SAVE_TXT = True

# set the path of the inputs:
input_path = "/home/astro/blum/juno/atmoNC/data_Julia/"

# file name of the input file:
input_name = input_path + "gntp.101.gst.root"

# set the path, where the outputs are saved:
output_path = "/home/astro/blum/juno/atmoNC/"


# get the fraction of the different types on NC interaction channels:
(Number_Events,
 Frac_C12_B11_p, Frac_C12_B11_n_piplus, Frac_C12_B11_n_piminus_2piplus, Frac_C12_B11_p_piminus_piplus,
 Frac_C12_B11_p_2piminus_2piplus, Frac_C12_B11_piplus,
 Frac_C12_C11_n, Frac_C12_C11_p_piminus, Frac_C12_C11_n_piminus_piplus, Frac_C12_C11_p_2piminus_piplus,
 Frac_C12_C11_p_3piminus_2piplus, Frac_C12_C11_n_2piminus_2piplus,
 Frac_C12_B10_p_n, Frac_C12_B10_2p_piminus, Frac_C12_B10_p_n_piminus_piplus,
 Frac_C12_B10_2n_piplus, Frac_C12_B10_2n_piminus_2piplus, Frac_C12_B10_2p_2piminus_piplus,
 Frac_C12_B10_2p_3piminus_2piplus, Frac_C12_B10_p_n_2piminus_2piplus,
 Frac_C12_C10_2n, Frac_C12_C10_p_n_piminus, Frac_C12_C10_p_n_2piminus_piplus, Frac_C12_C10_2n_piminus_piplus,
 Frac_C12_C10_2p_2piminus,
 Frac_C12_Be10_2p, Frac_C12_Be10_p_n_piplus, Frac_C12_Be10_p_n_piminus_2piplus, Frac_C12_Be10_2p_piminus_piplus,
 Frac_C12_Be10_2n_2piplus, Frac_C12_Be10_p_n_2piminus_3piplus, Frac_C12_Be10_2p_2piminus_2piplus,
 Frac_C12_Be10_2p_3piminus_3piplus,
 Frac_C12_B9_p_2n, Frac_C12_B9_p_2n_piminus_piplus, Frac_C12_B9_2p_n_3piminus_2piplus, Frac_C12_B9_2p_n_piminus,
 Frac_C12_B9_3n_piplus, Frac_C12_B9_p_2n_2piminus_2piplus, Frac_C12_B9_2p_n_2piminus_piplus,
 Frac_C12_Be9_2p_n, Frac_C12_Be9_p_2n_piplus, Frac_C12_Be9_3p_piminus, Frac_C12_Be9_p_2n_piminus_2piplus,
 Frac_C12_Be9_2p_n_piminus_piplus, Frac_C12_Be9_2p_n_3piminus_3piplus, Frac_C12_Be9_2p_n_2piminus_2piplus,
 Frac_C12_Be9_3n_2piplus, Frac_C12_Be9_3p_2piminus_piplus,
 Frac_C12_Be8_2p_2n, Frac_C12_Be8_3p_n_piminus, Frac_C12_Be8_p_3n_piplus, Frac_C12_Be8_2p_2n_2piminus_2piplus,
 Frac_C12_Be8_4n_2piplus, Frac_C12_Be8_2p_2n_piminus_piplus, Frac_C12_Be8_3p_n_2piminus_piplus,
 Frac_C12_Be8_4p_2piminus,
 Frac_C12_C9_p_2n_piminus, Frac_C12_C9_3n, Frac_C12_C9_2p_n_2piminus, Frac_C12_C9_3n_2piminus_2piplus,
 Frac_C12_Be7_2p_3n, Frac_C12_Be7_p_4n_piplus, Frac_C12_Be7_2p_3n_2piminus_2piplus, Frac_C12_Be7_3p_2n_piminus,
 Frac_C12_Be7_4p_n_2piminus, Frac_C12_Be7_3p_2n_2piminus_piplus,
 Frac_C12_Li6_3p_3n, Frac_C12_Li6_2p_4n_piplus, Frac_C12_Li6_5p_n_2piminus, Frac_C12_Li6_2p_4n_piminus_2piplus,
 Frac_C12_Li6_4p_2n_piminus, Frac_C12_Li6_3p_3n_piminus_piplus,
 Frac_C12_Li8_3p_n, Frac_C12_Li8_4p_piminus, Frac_C12_Li8_4p_2piminus_piplus, Frac_C12_Li8_2p_2n_piplus,
 Frac_C12_Li8_3p_n_piminus_piplus,
 Frac_C12_Li7_2p_3n_piplus, Frac_C12_Li7_4p_n_piminus, Frac_C12_Li7_3p_2n, Frac_C12_Li7_3p_2n_piminus_piplus,
 Frac_C12_Li7_4p_n_2piminus_piplus, Frac_C12_Li7_2p_3n_piminus_2piplus,
 Frac_C12_B8_p_3n, Frac_C12_B8_p_3n_piminus_piplus, Frac_C12_B8_2p_2n_2piminus_piplus, Frac_C12_B8_2p_2n_piminus,
 Frac_C12_B8_4n_piplus,
 Frac_C12_Li9_2p_n_piplus, Frac_C12_Li9_3p, Frac_C12_Li9_3p_piminus_piplus, Frac_C12_Li9_2p_n_piminus_2piplus,
 Frac_C12_Li9_p_2n_piminus_3piplus,
 Frac_C12_C8_4n, Frac_C12_He8_4p, Frac_C12_B7_p_4n, Frac_C12_He7_4p_n, Frac_C12_H7_5p, Frac_C12_Be6_2p_4n,
 Frac_C12_Li5_3p_4n, Frac_C12_Li4_3p_5n, Frac_C12_He6_4p_2n, Frac_C12_He5_4p_3n, Frac_C12_He4_4p_4n, Frac_C12_He3_4p_5n,
 Frac_C12_H6_5p_n, Frac_C12_H5_5p_2n, Frac_C12_H4_5p_3n, Frac_C12_H3_5p_4n, Frac_C12_H2_5p_5n,
 Frac_C12_C12, Frac_C12_NoIso,
 Frac_no_C12, Frac_ES_proton_chID, Frac_ES_electron_chID, Frac_ES_O16_chID, Frac_ES_N14_chID, Frac_ES_S32_chID,
 Frac_C12_missing) \
    = NC_background_functions.get_channels_from_original_genie_file(input_name)


""" Save information from get_interaction_channel() into txt file: """
if SAVE_TXT:
    np.savetxt(output_path + "interaction_channels_qel_NC_{0:d}evts.txt".format(Number_Events),
               np.array([Frac_C12_B11_p, Frac_C12_B11_n_piplus, Frac_C12_B11_n_piminus_2piplus,
                         Frac_C12_B11_p_piminus_piplus, Frac_C12_B11_p_2piminus_2piplus, Frac_C12_B11_piplus,
                         Frac_C12_C11_n, Frac_C12_C11_p_piminus, Frac_C12_C11_n_piminus_piplus,
                         Frac_C12_C11_p_2piminus_piplus, Frac_C12_C11_p_3piminus_2piplus,
                         Frac_C12_C11_n_2piminus_2piplus,
                         Frac_C12_B10_p_n, Frac_C12_B10_2p_piminus, Frac_C12_B10_p_n_piminus_piplus,
                         Frac_C12_B10_2n_piplus, Frac_C12_B10_2n_piminus_2piplus, Frac_C12_B10_2p_2piminus_piplus,
                         Frac_C12_B10_2p_3piminus_2piplus, Frac_C12_B10_p_n_2piminus_2piplus,
                         Frac_C12_C10_2n, Frac_C12_C10_p_n_piminus, Frac_C12_C10_p_n_2piminus_piplus,
                         Frac_C12_C10_2n_piminus_piplus, Frac_C12_C10_2p_2piminus,
                         Frac_C12_Be10_2p, Frac_C12_Be10_p_n_piplus, Frac_C12_Be10_p_n_piminus_2piplus,
                         Frac_C12_Be10_2p_piminus_piplus, Frac_C12_Be10_2n_2piplus, Frac_C12_Be10_p_n_2piminus_3piplus,
                         Frac_C12_Be10_2p_2piminus_2piplus, Frac_C12_Be10_2p_3piminus_3piplus,
                         Frac_C12_B9_p_2n, Frac_C12_B9_p_2n_piminus_piplus, Frac_C12_B9_2p_n_3piminus_2piplus,
                         Frac_C12_B9_2p_n_piminus, Frac_C12_B9_3n_piplus, Frac_C12_B9_p_2n_2piminus_2piplus,
                         Frac_C12_B9_2p_n_2piminus_piplus,
                         Frac_C12_Be9_2p_n, Frac_C12_Be9_p_2n_piplus, Frac_C12_Be9_3p_piminus,
                         Frac_C12_Be9_p_2n_piminus_2piplus, Frac_C12_Be9_2p_n_piminus_piplus,
                         Frac_C12_Be9_2p_n_3piminus_3piplus, Frac_C12_Be9_2p_n_2piminus_2piplus,
                         Frac_C12_Be9_3n_2piplus, Frac_C12_Be9_3p_2piminus_piplus,
                         Frac_C12_Be8_2p_2n, Frac_C12_Be8_3p_n_piminus, Frac_C12_Be8_p_3n_piplus,
                         Frac_C12_Be8_2p_2n_2piminus_2piplus, Frac_C12_Be8_4n_2piplus,
                         Frac_C12_Be8_2p_2n_piminus_piplus, Frac_C12_Be8_3p_n_2piminus_piplus, Frac_C12_Be8_4p_2piminus,
                         Frac_C12_C9_p_2n_piminus, Frac_C12_C9_3n, Frac_C12_C9_2p_n_2piminus,
                         Frac_C12_C9_3n_2piminus_2piplus,
                         Frac_C12_Be7_2p_3n, Frac_C12_Be7_p_4n_piplus, Frac_C12_Be7_2p_3n_2piminus_2piplus,
                         Frac_C12_Be7_3p_2n_piminus, Frac_C12_Be7_4p_n_2piminus, Frac_C12_Be7_3p_2n_2piminus_piplus,
                         Frac_C12_Li6_3p_3n, Frac_C12_Li6_2p_4n_piplus, Frac_C12_Li6_5p_n_2piminus,
                         Frac_C12_Li6_2p_4n_piminus_2piplus, Frac_C12_Li6_4p_2n_piminus,
                         Frac_C12_Li6_3p_3n_piminus_piplus,
                         Frac_C12_Li8_3p_n, Frac_C12_Li8_4p_piminus, Frac_C12_Li8_4p_2piminus_piplus,
                         Frac_C12_Li8_2p_2n_piplus, Frac_C12_Li8_3p_n_piminus_piplus,
                         Frac_C12_Li7_2p_3n_piplus, Frac_C12_Li7_4p_n_piminus, Frac_C12_Li7_3p_2n,
                         Frac_C12_Li7_3p_2n_piminus_piplus, Frac_C12_Li7_4p_n_2piminus_piplus,
                         Frac_C12_Li7_2p_3n_piminus_2piplus,
                         Frac_C12_B8_p_3n, Frac_C12_B8_p_3n_piminus_piplus, Frac_C12_B8_2p_2n_2piminus_piplus,
                         Frac_C12_B8_2p_2n_piminus, Frac_C12_B8_4n_piplus,
                         Frac_C12_Li9_2p_n_piplus, Frac_C12_Li9_3p, Frac_C12_Li9_3p_piminus_piplus,
                         Frac_C12_Li9_2p_n_piminus_2piplus, Frac_C12_Li9_p_2n_piminus_3piplus,
                         Frac_C12_C8_4n, Frac_C12_He8_4p, Frac_C12_B7_p_4n, Frac_C12_He7_4p_n, Frac_C12_H7_5p,
                         Frac_C12_Be6_2p_4n, Frac_C12_Li5_3p_4n, Frac_C12_Li4_3p_5n, Frac_C12_He6_4p_2n,
                         Frac_C12_He5_4p_3n, Frac_C12_He4_4p_4n, Frac_C12_He3_4p_5n,
                         Frac_C12_H6_5p_n, Frac_C12_H5_5p_2n, Frac_C12_H4_5p_3n, Frac_C12_H3_5p_4n, Frac_C12_H2_5p_5n,
                         Frac_C12_C12, Frac_C12_NoIso,
                         Frac_no_C12, Frac_ES_proton_chID, Frac_ES_electron_chID, Frac_ES_O16_chID, Frac_ES_N14_chID,
                         Frac_ES_S32_chID, Frac_C12_missing]), fmt='%4.5f',
               header="Information about the number of the different NC interaction channels of Julia's GENIE file\n"
                      "(input file: {0}, script: read_GENIE_file.py ({1})):\n"
                      "Fraction of NC channel nu + C12 -> B11 + p,\n"
                      "Fraction of NC channel nu + C12 -> B11 + n + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> B11 + n + pi_minus + 2*pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> B11 + p + pi_minus + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> B11 + p + 2*pi_minus + 2*pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> B11 + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> C11 + n,\n"
                      "Fraction of NC channel nu + C12 -> C11 + p + pi_minus,\n"
                      "Fraction of NC channel nu + C12 -> C11 + n + pi_minus + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> C11 + p + 2*pi_minus + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> C11 + p + 3*pi_minus + 2*pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> C11 + n + 2*pi_minus + 2*pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> B10 + p + n,\n"
                      "Fraction of NC channel nu + C12 -> B10 + 2p + pi_minus,\n"
                      "Fraction of NC channel nu + C12 -> B10 + p + n + pi_minus + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> B10 + 2n + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> B10 + 2n + pi_minus + 2*pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> B10 + 2p + 2*pi_minus + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> B10 + 2p + 3*pi_minus + 2*pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> B10 + p + n + 2*pi_minus + 2*pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> C10 + 2n,\n"
                      "Fraction of NC channel nu + C12 -> C10 + p + n + pi_minus,\n"
                      "Fraction of NC channel nu + C12 -> C10 + p + n + 2*pi_minus + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> C10 + 2n + pi_minus + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> C10 + 2p + 2*pi_minus,\n"
                      "Fraction of NC channel nu + C12 -> Be10 + 2p,\n"
                      "Fraction of NC channel nu + C12 -> Be10 + p + n + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Be10 + p + n + pi_minus + 2*pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Be10 + 2p + pi_minus + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Be10 + 2n + 2*pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Be10 + p + n + 2*pi_minus + 3*pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Be10 + 2p + 2*pi_minus + 2*pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Be10 + 2p + 3*pi_minus + 3*pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> B9 + p + 2n,\n"
                      "Fraction of NC channel nu + C12 -> B9 + p + 2n + pi_minus + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> B9 + 2p + n + 3*pi_minus + 2*pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> B9 + 2p + n + pi_minus,\n"
                      "Fraction of NC channel nu + C12 -> B9 + 3n + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> B9 + p + 2n + 2*pi_minus + 2*pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> B9 + 2p + n + 2*pi_minus + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Be9 + 2p + n,\n"     
                      "Fraction of NC channel nu + C12 -> Be9 + p + 2n + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Be9 + 3p + pi_minus,\n"
                      "Fraction of NC channel nu + C12 -> Be9 + p + 2n + pi_minus + 2*pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Be9 + 2p + n + pi_minus + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Be9 + 2p + n + 3*pi_minus + 3*pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Be9 + 2p + n + 2*pi_minus + 2*pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Be9 + 3n + 2*pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Be9 + 3p + 2*pi_minus + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Be8 + 2p + 2n,\n"
                      "Fraction of NC channel nu + C12 -> Be8 + 3p + n + pi_minus,\n"
                      "Fraction of NC channel nu + C12 -> Be8 + p + 3n + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Be8 + 2p + 2n + 2*pi_minus + 2*pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Be8 + 4n + 2*pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Be8 + 2p + 2n + pi_minus + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Be8 + 3p + n + 2*pi_minus + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Be8 + 4p + 2*pi_minus,\n"
                      "Fraction of NC channel nu + C12 -> C9 + p + 2n + pi_minus,\n"
                      "Fraction of NC channel nu + C12 -> C9 + 3n,\n"
                      "Fraction of NC channel nu + C12 -> C9 + 2p + n + 2*pi_minus,\n"
                      "Fraction of NC channel nu + C12 -> C9 + 3n + 2*pi_minus + 2*pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Be7 + 2p + 3n,\n"
                      "Fraction of NC channel nu + C12 -> Be7 + p + 4n + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Be7 + 2p + 3n + 2*pi_minus + 2*pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Be7 + 3p + 2n + pi_minus,\n"
                      "Fraction of NC channel nu + C12 -> Be7 + 4p + n + 2*pi_minus,\n"
                      "Fraction of NC channel nu + C12 -> Be7 + 3p + 2n + 2*pi_minus + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Li6 + 3p + 3n,\n"
                      "Fraction of NC channel nu + C12 -> Li6 + 2p + 4n + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Li6 + 5p + n + 2*pi_minus,\n"
                      "Fraction of NC channel nu + C12 -> Li6 + 2p + 4n + pi_minus + 2*pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Li6 + 4p + 2n + pi_minus,\n"
                      "Fraction of NC channel nu + C12 -> Li6 + 3p + 3n + pi_minus + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Li8 + 3p + n,\n"
                      "Fraction of NC channel nu + C12 -> Li8 + 4p + pi_minus,\n"
                      "Fraction of NC channel nu + C12 -> Li8 + 4p + 2*pi_minus + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Li8 + 2p + 2n + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Li8 + 3p + n + pi_minus + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Li7 + 2p + 3n + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Li7 + 4p + n + pi_minus,\n"
                      "Fraction of NC channel nu + C12 -> Li7 + 3p + 2n,\n"
                      "Fraction of NC channel nu + C12 -> Li7 + 3p + 2n + pi_minus + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Li7 + 4p + n + 2*pi_minus + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Li7 + 2p + 3n + pi_minus + 2*pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> B8 + p + 3n,\n"
                      "Fraction of NC channel nu + C12 -> B8 + p + 3n + pi_minus + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> B8 + 2p + 2n + 2*pi_minus + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> B8 + 2p + 2n + pi_minus,\n"
                      "Fraction of NC channel nu + C12 -> B8 + 4n + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Li9 + 2p + n + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Li9 + 3p,\n"
                      "Fraction of NC channel nu + C12 -> Li9 + 3p + pi_minus + pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Li9 + 2p + n + pi_minus + 2*pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> Li9 + p + 2n + pi_minus + 3*pi_plus,\n"
                      "Fraction of NC channel nu + C12 -> C8 + 4n,\n"
                      "Fraction of NC channel nu + C12 -> He8 + 4p,\n"
                      "Fraction of NC channel nu + C12 -> B7 + p + 4n,\n"
                      "Fraction of NC channel nu + C12 -> He7 + 4p + n,\n"
                      "Fraction of NC channel nu + C12 -> H7 + 5p,\n"
                      "Fraction of NC channel nu + C12 -> Be6 + 2p + 4n,\n"
                      "Fraction of NC channel nu + C12 -> Li5 + 3p + 4n,\n"
                      "Fraction of NC channel nu + C12 -> Li4 + 3p + 5n,\n"
                      "Fraction of NC channel nu + C12 -> He6 + 4p + 2n,\n"
                      "Fraction of NC channel nu + C12 -> He5 + 4p + 3n,\n"
                      "Fraction of NC channel nu + C12 -> He4 + 4p + 4n,\n"
                      "Fraction of NC channel nu + C12 -> He3 + 4p + 5n,\n"
                      "Fraction of NC channel nu + C12 -> H6 + 5p + n,\n"
                      "Fraction of NC channel nu + C12 -> H5 + 5p + 2n,\n"
                      "Fraction of NC channel nu + C12 -> H4 + 5p + 3n (+ maybe pions and Kaons),\n"
                      "Fraction of NC channel nu + C12 -> H3 + 5p + 4n (+ maybe pions and Kaons),\n"
                      "Fraction of NC channel nu + C12 -> H2 + 5p + 5n (+ maybe pions and Kaons),\n"
                      "Fraction of NC channels nu + C12 -> nu + C12 + ...,\n"
                      "Fraction of NC channels without isotope nu + C12 -> nu + X*p + Y*n + Z*pion,\n"
                      "Fraction of NC channels WITHOUT C12 as target,\n"
                      "Fraction of ES channel nu + p -> nu + p + ...,\n"
                      "Fraction of ES channel nu + electron -> nu + electron + ...,\n"
                      "Fraction of ES channel nu + O16 -> nu + O16 + ...,\n"
                      "Fraction of ES channel nu + N14 -> nu + N14 + ...,\n"
                      "Fraction of ES channel nu + S32 -> nu + S32 + ...,\n"
                      "Fraction of interactions not yet included in function get_channels_from_original_genie_file():"
               .format(input_name, now))


