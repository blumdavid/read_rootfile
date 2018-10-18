""" Script to check out the results from the DSNB-NC generator of Jie Cheng.

    The ROOT-file generated with DSNB-NC.exe is read and analyzed with this script.

    DSNB-NC.exe is the generator of the Neutral Current background from atmospheric neutrinos and antineutrinos
    (build by Jie Cheng).

    The interactions of the neutrino and antineutrinos with the liquid scintillator are simulated with GENIE (or NuWro).
    For the flux of the atmospheric neutrinos (nu_e, nu_e_bar, nu_mu, nu_mu_bar) the flux of Honda for JUNO site is used
    (same flux like in my simulation of atmospheric CC background, energy range from 100 MeV to 10 GeV).

    -> is the solar average used? or solar minimum? or solar maximum?

    For the interactions with the LS, the interactions on C12 (-> neutral current interactions) and on free protons
    (-> elastic scattering) are considered.

    The output of the GENIE simulation is saved in "genie_data.root".

    Then the deexcitation of the products from the NC interactions is simulated:
    The status of the residual nucleus is unknown. Therefor the deexcitation of the residual nucleus is assumed from
    a simple shell model
    -> 1/3 of the residual nuclei are excited
    -> 2/3 of the residual nuclei are NOT excited, but in the ground state
    (Jie is working on the implementation of a more complicated nuclear model, but it is not implemented in the DSNB-NC
    generator yet (10.10.2018))

    The deexitation of the residual nucleus is simulated with the TALYS software and saved in "*deexcitation.root" for
    all nuclei (Li7, Li8, Li9, C9, C10, C11, Be7, Be9, Be10, B8, B9, B10, B11)

    With these root-files, the deexcitation is calculated and all final particles with their PDG ID, momentum and
    mass printed in the terminal.
    These information together with the target ID, the channel ID, the deexcitation ID, the isotope ID and the energy
    of the incoming neutrino are saved in the root output file.

"""

import ROOT
import glob
import NC_background_functions
import numpy as np
from matplotlib import pyplot as plt


""" set the path of the inputs: """
input_path = "/home/astro/blum/juno/test_output_DSNB_gen/"

""" file name of the input file: """
input_name = input_path + "test_100events.root"

# set the NC interaction event rate for GENIE in units of events/(s*20ktons) (see NC_generator.pdf from Jie Cheng)
# (float):
evt_rate_Genie = 3.59E-5

# total exposure time in years (float):
t_years = 10
# total time-exposure in seconds (1yr = 365.2425 * 24 * 60 * 60 sec = 3.1556952 * 10^7 sec), 10 years (float):
time = t_years * 3.156 * 10 ** 7

# read NC generator data to arrays:
(event_ID, projectile_PDG, projectile_E, target_PDG, NC_inter_ch_ID, deexcitation_ID, isotope_PDG, Nparticles,
 final_PDG, final_Px, final_Py, final_Pz) = NC_background_functions.read_nc_data(input_name)



