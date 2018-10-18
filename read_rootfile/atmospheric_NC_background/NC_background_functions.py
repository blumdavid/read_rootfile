""" Script contains different functions which can be used to check the root-output file from DSNB-NC.exe:

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
import numpy as np
from matplotlib import pyplot as plt

"""
class NCData:
    def __init__(self):  # this method creates the class object.
        # event ID (starts from 0):
        self.eventID = 0
        # PDG ID of the projectile (i.e. which neutrino is interacting):
        self.projectilePDG = 0
        # energy of the incoming neutrino in GeV:
        self.projectileEnergy = 0
        # PDG ID of the target particle (either C12 or proton):
        self.targetPDG = 0
        # Channel ID of the NC interaction, represents which particles are produced via the NC interaction:
        self.NCinteration_ch_ID = 0
        # Channel ID of the deexcitation, represents which particles are produced via the deexication of the produced
        # excited isotope:
        self.deexcitation_ID = 0
        # PDG ID of the isotope after the NC interaction, BUT before the deexcitation:
        self.isotopePDG = 0
        # number of final particles after NC interactions and deexcitation:
        self.Nparticles = 0
        # PDG ID of the final particles. It is an array of length "Nparticles":
        self.finalPDG = 0
        # momentum in x-direction of the final particles. It is an array of length "Nparticles":
        self.finalPx = 0
        # momentum in y-direction of the final particles. It is an array of length "Nparticles":
        self.finalPy = 0
        # momentum in z-direction of the final particles. It is an array of length "Nparticles":
        self.finalPz = 0
"""


def get_mass_from_pdg(pdg):
    """
    function to get the mass in GeV of a particle from its PDG ID (Monte Carlo Particle Number Scheme)
    (values taken from NCGenerator.cc, )

    :param pdg: PDG ID of the particle (integer)
    :return: mass[pdg]: mass of the particle in GeV (float)
    """
    mass = dict()
    # mass of gamma:
    mass[22] = 0
    # mass of pion_0:
    mass[111] = 0.13957
    # mass of pion_plus:
    mass[211] = 0.13957
    # mass of pion_minus:
    mass[-211] = 0.13957
    # mass of neutron:
    mass[2112] = 0.93957
    # proton:
    mass[2212] = 0.93827
    # deuterium H2 (stable):
    mass[1000010020] = 1.8756
    # tritium H3:
    mass[1000010030] = 2.8089
    # He3 (stable):
    mass[1000020030] = 2.8084
    # He4 or alpha (stable):
    mass[1000020040] = 3.7274
    # Li6 (stable):
    mass[1000030060] = 5.6015
    # Li7 (stable):
    mass[1000030070] = 6.5335
    # Li8:
    mass[1000030080] = 7.4708
    # Li9:
    mass[1000030090] = 8.4061
    # Be7:
    mass[1000040070] = 6.5344
    # Be8:
    mass[1000040080] = 7.4548
    # Be9 (stable):
    mass[1000040090] = 8.3925
    # Be10:
    mass[1000040100] = 9.3249
    # B8:
    mass[1000050080] = 7.4728
    # B9:
    mass[1000050090] = 8.3935
    # B10 (stable):
    mass[1000050100] = 9.3244
    # B11 (stable):
    mass[1000050110] = 10.2522
    # C9:
    mass[1000060090] = 8.4100
    # C10:
    mass[1000060100] = 9.3280
    # C11:
    mass[1000060110] = 10.2542
    # C12 (stable):
    mass[1000060120] = 11.1748

    return mass[pdg]


def read_nc_data(rootfile):
    """
    function reads a ROOT-file and saves the values from the root-tree to numpy arrays.

    :param rootfile: path to the input root file (string)

    :return:
    """
    # load the ROOT file:
    rfile = ROOT.TFile(rootfile)
    # get the TTree from the TFile:
    rtree = rfile.Get("genEvt")

    # get the number of entries, i.e. events, in the ROOT-file:
    number_entries = rtree.GetEntries()

    "preallocate all arrays: "
    # event ID (starts from 0) (1d array of integer):
    event_id = np.array([])
    # PDG ID of the projectile (i.e. which neutrino is interacting) (1d array integer):
    projectile_pdg = np.array([])
    # energy of the incoming neutrino in GeV (1d array of float):
    projectile_energy = np.array([])
    # PDG ID of the target particle (either C12 or proton) (1d array integer):
    target_pdg = np.array([])
    # Channel ID of the NC interaction, represents which particles are produced via the NC interaction
    # (1d array integer):
    nc_interaction_ch_id = np.array([])
    # Channel ID of the deexcitation, represents which particles are produced via the deexication of the produced
    # excited isotope (1d array integer):
    deexcitation_id = np.array([])
    # PDG ID of the isotope after the NC interaction, BUT before the deexcitation (1d array integer):
    isotope_pdg = np.array([])
    # number of final particles after NC interactions and deexcitation (1d array integer):
    n_particles = np.array([])

    # PDG ID of the final particles (list of np.arrays of integers):
    final_pdg = []
    # momentum in x-direction of the final particles (list of np.arrays of floats):
    final_px = []
    # momentum in y-direction of the final particles (list of np.arrays of floats):
    final_py = []
    # momentum in z-direction of the final particles (list of np.arrays of floats):
    final_pz = []

    # loop over every entry, i.e. every event, in the TTree:
    for event in range(number_entries):

        # get the current event in the TTree:
        rtree.GetEntry(event)

        # get the value of the event ID and append it to the array:
        evt_id = rtree.GetBranch('t_evtID').GetLeaf('t_evtID').GetValue()
        event_id = np.append(event_id, evt_id)

        # get the value of the projectile PDG and append it to the array:
        pjt_pdg = rtree.GetBranch('t_pPdg').GetLeaf('t_pPdg').GetValue()
        projectile_pdg = np.append(projectile_pdg, pjt_pdg)

        # get the value of the neutrino energy and append it to the array:
        pjt_e = rtree.GetBranch('t_pEn').GetLeaf('t_pEn').GetValue()
        projectile_energy = np.append(projectile_energy, pjt_e)

        # get the value of the target PDG and append it to the array:
        trt_pdg = rtree.GetBranch('t_tPdg').GetLeaf('t_tPdg').GetValue()
        target_pdg = np.append(target_pdg, trt_pdg)

        # get the value of the channel ID and append it to the array:
        ch_id = rtree.GetBranch('t_channelID').GetLeaf('t_channelID').GetValue()
        nc_interaction_ch_id = np.append(nc_interaction_ch_id, ch_id)

        # get the value of the deexcitation ID  and append it to the array:
        deex_id = rtree.GetBranch('t_deexID').GetLeaf('t_deexID').GetValue()
        deexcitation_id = np.append(deexcitation_id, deex_id)

        # get the value of the PDG of the produced isotope and append it to the array:
        iso_pdg = rtree.GetBranch('t_isoPdg').GetLeaf('t_isoPdg').GetValue()
        isotope_pdg = np.append(isotope_pdg, iso_pdg)

        # get the value of the number of particles and append it to the array:
        n_par = rtree.GetBranch('t_Npars').GetLeaf('t_Npars').GetValue()
        n_particles = np.append(n_particles, n_par)

        # get final PDGs of all final particles
        # preallocate an array, where all "n_par" values are stored:
        f_pdg_array = np.array([])

        # loop over the number of particles, get the final PDG and append it to the array:
        for index in range(int(n_par)):
            # get the value of the final PDG and append it to the array:
            f_pdg = rtree.GetBranch('t_pdg').GetLeaf('t_pdg').GetValue(index)
            f_pdg_array = np.append(f_pdg_array, f_pdg)

        # append the np.array to the list:
        final_pdg.append(f_pdg_array)

        # get final momentum Px of all final particles
        # preallocate an array, where all "n_par" values are stored:
        f_px_array = np.array([])

        # loop over the number of particles, get the final momentum and append it to the array:
        for index in range(int(n_par)):
            # get the value of the final momentum Px and append it to the array:
            f_px = rtree.GetBranch('t_px').GetLeaf('t_px').GetValue(index)
            f_px_array = np.append(f_px_array, f_px)

        # append the np.array to the list:
        final_px.append(f_px_array)

        # get final momentum Py of all final particles
        # preallocate an array, where all "n_par" values are stored:
        f_py_array = np.array([])

        # loop over the number of particles, get the final momentum and append it to the array:
        for index in range(int(n_par)):
            # get the value of the final momentum Py and append it to the array:
            f_py = rtree.GetBranch('t_py').GetLeaf('t_py').GetValue(index)
            f_py_array = np.append(f_py_array, f_py)

        # append the np.array to the list:
        final_py.append(f_py_array)

        # get final momentum Pz of all final particles
        # preallocate an array, where all "n_par" values are stored:
        f_pz_array = np.array([])

        # loop over the number of particles, get the final momentum and append it to the array:
        for index in range(int(n_par)):
            # get the value of the final momentum Pz and append it to the array:
            f_pz = rtree.GetBranch('t_pz').GetLeaf('t_pz').GetValue(index)
            f_pz_array = np.append(f_pz_array, f_pz)

        # append the np.array to the list:
        final_pz.append(f_py_array)

    return event_id, projectile_pdg, projectile_energy, target_pdg, nc_interaction_ch_id, deexcitation_id, \
           isotope_pdg, n_particles, final_pdg, final_px, final_py, final_pz


def get_neutrino_energy(projectile_pdg, projectile_energy):
    """
    function to get the energies of the different types of neutrinos (electron-neutrino, electron-antineutrino,
    muon-neutrino, muon-antineutrino, tau-neutrino, tau-antineutrino)

    :param projectile_pdg: PDG ID of the projectile, i.e. of the incoming neutrino (array of integers)
    :param projectile_energy: energy of the projectile, i.e. of the incoming neutrinos, in GeV (array of integers)

    :return:
    """

    """ Preallocate all arrays: """
    # energy of incoming electron neutrinos in GeV (np.array of float):
    energy_nu_e_incoming = np.array([])
    # energy of incoming electron antineutrinos in GeV (np.array of float):
    energy_nu_e_bar_incoming = np.array([])
    # energy of incoming muon neutrinos in GeV (np.array of float):
    energy_nu_mu_incoming = np.array([])
    # energy of incoming muon antineutrinos in GeV (np.array of float):
    energy_nu_mu_bar_incoming = np.array([])
    # energy of incoming tau neutrinos in GeV (np.array of float):
    energy_nu_tau_incoming = np.array([])
    # energy of incoming tau antineutrinos in GeV (np.array of float):
    energy_nu_tau_bar_incoming = np.array([])



