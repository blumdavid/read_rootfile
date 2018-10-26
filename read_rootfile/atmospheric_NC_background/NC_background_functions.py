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


def get_number_of_p_and_n_of_isotope(pdg_id):
    """
    function to get the number of protons and neutrons of a nulcei from its PDG ID

    :param pdg_id: PDG ID of the nuclei (float)
    :return:
    """
    # get the number of protons (integer)
    number_p = int((pdg_id - 1000000000) / 10000)

    # get the sum of protons and neutrons (integer):
    number_p_and_n = int((pdg_id - 1000000000 - number_p*10000) / 10)

    # calculate the number of neutrons (integer):
    number_n = number_p_and_n - number_p

    return number_p, number_n


def get_number_of_particles_of_channelid(channel_id):
    """
    Function to calculate the number of neutrinos, protons, neutron, pion_minus and pion_plus from the channel ID

    :param channel_id: Channel ID of the NC interaction, represents which particles are produced via the NC
    interaction (float)
    :return:
    """
    # preallocate variables:
    number_p = 0
    number_n = 0
    number_pion_minus = 0
    number_pion_plus = 0

    if channel_id > 10000:
        # channel_id > 10000 -> at least one pion is created:
        # number of protons from channel ID (integer):
        number_p = int((channel_id - 10000) / 1000)
        # number of neutrons from channel ID (integer):
        number_n = int((channel_id - 10000 - number_p*1000) / 100)
        # number of pion_minus from channel ID (integer):
        number_pion_minus = int((channel_id - 10000 - number_p*1000 - number_n*100) / 10)
        # number of pion_plus from channel ID (integer):
        number_pion_plus = int(channel_id - 10000 - number_p*1000 - number_n*100 - number_pion_minus*10)

    elif channel_id > 3:
        # channel_id > 3 -> only protons and neutrons are created:
        # number of protons from channel ID (integer):
        number_p = int((channel_id - 100) / 10)
        # number of neutrons from channel ID (integer):
        number_n = int(channel_id - 100 - number_p*10)

    elif channel_id == 2:
        print("channel ID = 2 in get_number_of_particles_of_channelid()")

    elif channel_id == 3:
        print("channel ID = 3 in get_number_of_particles_of_channelid()")

    return number_p, number_n, number_pion_minus, number_pion_plus


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


def get_neutrino_energy(projectile_pdg, projectile_energy, bin_width, event_rate, time):
    """
    function to get the number of events (number of NC interactions of atmospheric neutrinos) as function of the
    neutrino energy of the different types of neutrinos (electron-neutrino, electron-antineutrino, muon-neutrino,
    muon-antineutrino, tau-neutrino, tau-antineutrino).
    Also the total number of events for each neutrino type and the fraction of each neutrino type to the total number
    of events is calcualted.

    :param projectile_pdg: PDG ID of the projectile, i.e. of the incoming neutrino (array float)
    :param projectile_energy: energy of the projectile, i.e. of the incoming neutrinos, in GeV (array of float)
    :param bin_width: bin width of the array, which represents the energy of incoming neutrinos in GeV (float)
    :param event_rate: NC interaction event rate in units of 1/(sec * 20 kton) (float)
    :param time: total exposure time in seconds (float)

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

    # check, if input arrays have same length:
    if len(projectile_pdg) != len(projectile_energy):
        print("WARNING (in function get_neutrino_energy()): input arrays don't have same length!!")

    # get number of entries in the arrays:
    number_entries = len(projectile_pdg)

    # loop over all entries of the arrays:
    for index in range(number_entries):

        # check the PDG ID:
        if projectile_pdg[index] == 12:
            # for electron-neutrino PDG = 12:
            energy = projectile_energy[index]
            energy_nu_e_incoming = np.append(energy_nu_e_incoming, energy)

        elif projectile_pdg[index] == -12:
            # for electron-antineutrino PDG = -12:
            energy = projectile_energy[index]
            energy_nu_e_bar_incoming = np.append(energy_nu_e_bar_incoming, energy)

        elif projectile_pdg[index] == 14:
            # for muon-neutrino PDG = 14:
            energy = projectile_energy[index]
            energy_nu_mu_incoming = np.append(energy_nu_mu_incoming, energy)

        elif projectile_pdg[index] == -14:
            # for muon-antineutrino PDG == -14:
            energy = projectile_energy[index]
            energy_nu_mu_bar_incoming = np.append(energy_nu_mu_bar_incoming, energy)

        elif projectile_pdg[index] == 16:
            # for tau-neutrino PDG = 16:
            energy = projectile_energy[index]
            energy_nu_tau_incoming = np.append(energy_nu_tau_incoming, energy)

        elif projectile_pdg[index] == -16:
            # for tau-antineutrino PDG = -16:
            energy = projectile_energy[index]
            energy_nu_tau_bar_incoming = np.append(energy_nu_tau_bar_incoming, energy)

        else:
            print("WARNING (in function get_neutrino_energy(): NOT only neutrinos as projectile!")

    """ get number of events from NC interactions of different neutrino types: """
    # calculate the number of events from NC interaction of electron-neutrinos (integer):
    n_nu_e = len(energy_nu_e_incoming)
    # calculate the number of events from NC interaction of electron-antineutrinos (integer):
    n_nu_e_bar = len(energy_nu_e_bar_incoming)
    # calculate the number of events from NC interaction of muon-neutrinos (integer):
    n_nu_mu = len(energy_nu_mu_incoming)
    # calculate the number of events from NC interaction of muon-antineutrinos (integer):
    n_nu_mu_bar = len(energy_nu_mu_bar_incoming)
    # calculate the number of events from NC interaction of tau-neutrinos (integer):
    n_nu_tau = len(energy_nu_tau_incoming)
    # calculate the number of events from NC interaction of tau-antineutrinos (integer):
    n_nu_tau_bar = len(energy_nu_tau_bar_incoming)

    """ get fraction of events from NC interactions of different neutrino types (IN PERCENT): """
    # calculate the fraction of events from electron-neutrinos of the all events in % (float):
    fraction_nu_e = float(n_nu_e)/float(number_entries)*100
    # calculate the fraction of events from electron-antineutrinos of the all events in % (float):
    fraction_nu_e_bar = float(n_nu_e_bar)/float(number_entries)*100
    # calculate the fraction of events from muon-neutrinos of the all events in % (float):
    fraction_nu_mu = float(n_nu_mu)/float(number_entries)*100
    # calculate the fraction of events from muon-antineutrinos of the all events in % (float):
    fraction_nu_mu_bar = float(n_nu_mu_bar)/float(number_entries)*100
    # calculate the fraction of events from tau-neutrinos of the all events in % (float):
    fraction_nu_tau = float(n_nu_tau)/float(number_entries)*100
    # calculate the fraction of events from tau-antineutrinos of the all events in % (float):
    fraction_nu_tau_bar = float(n_nu_tau_bar)/float(number_entries)*100

    """ create histograms with the energy arrays from above to get the number of events per bin: """
    energy_range = np.arange(0, np.max(projectile_energy)+2, bin_width)
    events_nu_e_in, bins1 = np.histogram(energy_nu_e_incoming, energy_range)
    events_nu_e_bar_in, bins1 = np.histogram(energy_nu_e_bar_incoming, energy_range)
    events_nu_mu_in, bins1 = np.histogram(energy_nu_mu_incoming, energy_range)
    events_nu_mu_bar_in, bins1 = np.histogram(energy_nu_mu_bar_incoming, energy_range)
    events_nu_tau_in, bins1 = np.histogram(energy_nu_tau_incoming, energy_range)
    events_nu_tau_bar_in, bins1 = np.histogram(energy_nu_tau_bar_incoming, energy_range)

    """ calculate the number of neutrino NC interactions as function of energy per "time" per 20 ktons: """
    # TODO-me: event rate and exposure time are NOT included correctly!!!
    # INFO-me: Calculation is only correct, when event_rate=1 and time=1!!!
    event_nu_e_incoming = events_nu_e_in * event_rate * time
    event_nu_e_bar_incoming = events_nu_e_bar_in * event_rate * time
    event_nu_mu_incoming = events_nu_mu_in * event_rate * time
    event_nu_mu_bar_incoming = events_nu_mu_bar_in * event_rate * time
    event_nu_tau_incoming = events_nu_tau_in * event_rate * time
    event_nu_tau_bar_incoming = events_nu_tau_bar_in * event_rate * time

    """ calculate the total number of neutrino NC interactions per "time" and 20 ktons: """
    number_nu_e_incoming = np.sum(event_nu_e_incoming)
    number_nu_e_bar_incoming = np.sum(event_nu_e_bar_incoming)
    number_nu_mu_incoming = np.sum(event_nu_mu_incoming)
    number_nu_mu_bar_incoming = np.sum(event_nu_mu_bar_incoming)
    number_nu_tau_incoming = np.sum(event_nu_tau_incoming)
    number_nu_tau_bar_incoming = np.sum(event_nu_tau_bar_incoming)

    return energy_range, event_nu_e_incoming, event_nu_e_bar_incoming, event_nu_mu_incoming, event_nu_mu_bar_incoming, \
           event_nu_tau_incoming, event_nu_tau_bar_incoming, number_nu_e_incoming, number_nu_e_bar_incoming, \
           number_nu_mu_incoming, number_nu_mu_bar_incoming, number_nu_tau_incoming, number_nu_tau_bar_incoming, \
           fraction_nu_e, fraction_nu_e_bar, fraction_nu_mu, fraction_nu_mu_bar, fraction_nu_tau, \
           fraction_nu_tau_bar


def get_target_ratio(projectile_energy, target_pdg, bin_width):
    """
    function to get the fraction of NC interaction events on target particles (C12, N14, O16, S32) and of elastic
    scattering events on free protons and electrons, and to get the number of events of the different target particles
    as function of the energy of the incoming neutrinos.

    :param projectile_energy: energy of the projectile, i.e. of the incoming neutrinos, in GeV (array of float)
    :param target_pdg: PDG ID of the target particle (array of float)
    :param bin_width: bin width of the array, which represents the energy of incoming neutrinos in GeV (float)

    :return:
    """

    """ preallocate the arrays: """
    # energy of the incoming neutrino interacting via NC with C12:
    energy_nu_c12 = np.array([])
    # energy of the incoming neutrino interacting via ES with free protons:
    energy_nu_proton = np.array([])
    # energy of the incoming neutrino interacting via NC with N14:
    energy_nu_n14 = np.array([])
    # energy of the incoming neutrino interacting via NC with O16:
    energy_nu_o16 = np.array([])
    # energy of the incoming neutrino interacting via ES with an electron:
    energy_nu_electron = np.array([])
    # energy of the incoming neutrino interacting via NC with S32:
    energy_nu_s32 = np.array([])

    # check, if input arrays have same length:
    if len(projectile_energy) != len(target_pdg):
        print("WARNING (in function get_target_ratio()): input arrays don't have same length!!")

    # get number of entries in the arrays:
    number_entries = len(projectile_energy)

    # loop over all entries in the arrays:
    for index in range(number_entries):

        # check the PDG ID:
        if target_pdg[index] == 1000060120:
            # PDG ID of C12 = 1000060120:
            energy = projectile_energy[index]
            energy_nu_c12 = np.append(energy_nu_c12, energy)

        elif target_pdg[index] == 2212:
            # PDG ID of proton = 2212:
            energy = projectile_energy[index]
            energy_nu_proton = np.append(energy_nu_proton, energy)

        elif target_pdg[index] == 1000070140:
            # PDG ID of N14 = 1000070140:
            energy = projectile_energy[index]
            energy_nu_n14 = np.append(energy_nu_n14, energy)

        elif target_pdg[index] == 1000080160:
            # PDG ID of O16 = 1000080160:
            energy = projectile_energy[index]
            energy_nu_o16 = np.append(energy_nu_o16, energy)

        elif target_pdg[index] == 11:
            # PDG ID of electron = 11:
            energy = projectile_energy[index]
            energy_nu_electron = np.append(energy_nu_electron, energy)

        elif target_pdg[index] == 1000160320:
            # PDG ID of S32 = 1000160320:
            energy = projectile_energy[index]
            energy_nu_s32 = np.append(energy_nu_s32, energy)

        else:
            print("WARNING (in function get_target_ratio(): NOT only C12 or protons as targets!")
            print(target_pdg[index])


    """ get number of events for the two different targets: """
    # calculate the number of NC interaction events on C12 (integer):
    n_c12 = len(energy_nu_c12)
    # calculate the number of elastic scattering events on free protons (integer):
    n_proton = len(energy_nu_proton)
    # calculate the number of NC interaction events on N14 (integer):
    n_n14 = len(energy_nu_n14)
    # calculate the number of NC interaction events on O16 (integer):
    n_o16 = len(energy_nu_o16)
    # calculate the number of elastic scattering events on electrons (integer):
    n_electron = len(energy_nu_electron)
    # calculate the number of NC interaction events on S32 (integer):
    n_s32 = len(energy_nu_s32)

    """ get fraction of events for the two different targets (IN PERCENT): """
    # calculate the fraction of NC interaction events on C12 in % (float):
    fraction_c12 = float(n_c12)/float(number_entries)*100
    # calculate the fraction of ES events on free protons in % (float):
    fraction_proton = float(n_proton)/float(number_entries)*100
    # calculate the fraction of NC interaction events on N14 in % (float):
    fraction_n14 = float(n_n14)/float(number_entries)*100
    # calculate the fraction of NC interaction events on O16 in % (float):
    fraction_o16 = float(n_o16)/float(number_entries)*100
    # calculate the fraction of ES events on electrons in % (float):
    fraction_electron = float(n_electron)/float(number_entries)*100
    # calculate the fraction of NC interaction events on S32 in % (float):
    fraction_s32 = float(n_s32)/float(number_entries)*100

    """ create histograms with the energy arrays from above to get the number of events per bin: """
    energy_range = np.arange(0, np.max(projectile_energy)+2, bin_width)
    events_c12, bins1 = np.histogram(energy_nu_c12, energy_range)
    events_proton, bins1 = np.histogram(energy_nu_proton, energy_range)
    events_n14, bins1 = np.histogram(energy_nu_n14, energy_range)
    events_o16, bins1 = np.histogram(energy_nu_o16, energy_range)
    events_electron, bins1 = np.histogram(energy_nu_electron, energy_range)
    events_s32, bins1 = np.histogram(energy_nu_s32, energy_range)

    # TODO: the event rate and exposure time is NOT included yet!!!

    return energy_range, events_c12, events_proton, events_n14, events_o16, events_electron, events_s32, \
           n_c12, n_proton, n_n14, n_o16, n_electron, n_s32, \
           fraction_c12, fraction_proton, fraction_n14, fraction_o16, fraction_electron, fraction_s32


def get_interaction_channel(channel_id, isotope_pdg, target_pdg):
    """
    function to get the different NC interaction channels from the interaction of neutrinos with the target particles

    :param channel_id: ID of the NC interaction channel, defines the product particles of the NC interactions
    (array of float)
    :param isotope_pdg: PDG ID of the isotope, that is created via the NC interaction (before deexcitation)
    (array of float)
    :param target_pdg: PDG ID of the target particle (array of float)

    :return:
    """

    # get the number of entries of the array (integer):
    number_entries = len(channel_id)

    """ preallocate arrays and variables: """
    number_no_c12 = 0

    """ B11 """
    # number of interaction channel: nu + C12 -> B11 + p (integer):
    number_c12_b11_p = 0
    # number of interaction channel: nu + C12 -> B11 + n + pi_plus (integer):
    number_c12_b11_n_piplus = 0
    # number of interaction channel: nu + C12 -> B11 + n + pi_minus + 2*pi_plus (integer):
    number_c12_b11_n_piminus_2piplus = 0
    # number of interaction channel: nu + C12 -> B11 + p + pi_minus + pi_plus (integer):
    number_c12_b11_p_piminus_piplus = 0
    # number of interaction channel: nu + C12 -> B11 + p + 2*pi_minus + 2*pi_plus (integer):
    number_c12_b11_p_2piminus_2piplus = 0
    # number of interaction channel: nu + C12 -> B11 + pi_plus (integer):
    number_c12_b11_piplus = 0
    # number of OTHER interaction channels: nu + C12 -> B11 + ...:
    number_c12_b11_other = 0

    """ C11 """
    # number of interaction channel: nu + C12 -> C11 + n (integer):
    number_c12_c11_n = 0
    # number of interaction channel: nu + C12 -> C11 + p + pi_minus (integer):
    number_c12_c11_p_piminus = 0
    # number of interaction channel: nu + C12 -> C11 + n + pi_minus + pi_plus (integer):
    number_c12_c11_n_piminus_piplus = 0
    # number of interaction channel: nu + C12 -> C11 + p + 2*pi_minus + pi_plus (integer):
    number_c12_c11_p_2piminus_piplus = 0
    # number of interaction channel: nu + C12 -> C11 + p + 3*pi_minus + 2*pi_plus (integer):
    number_c12_c11_p_3piminus_2piplus = 0
    # number of interaction channel: nu + C12 -> C11 + n + 2*pi_minus + 2*pi_plus (integer):
    number_c12_c11_n_2piminus_2piplus = 0
    # number of OTHER interaction channels: nu + C12 -> C11 + ...:
    number_c12_c11_other = 0

    """ B10 """
    # number of interaction channel: nu + C12 -> B10 + p + n (integer):
    number_c12_b10_p_n = 0
    # number of interaction channel: nu + C12 -> B10 + 2p + pi_minus (integer):
    number_c12_b10_2p_piminus = 0
    # number of interaction channel: nu + C12 -> B10 + p + n + pi_minus + pi_plus (integer):
    number_c12_b10_p_n_piminus_piplus = 0
    # number of interaction channel: nu + C12 -> B10 + 2n + pi_plus (integer):
    number_c12_b10_2n_piplus = 0
    # number of interaction channel: nu + C12 -> B10 + 2n + pi_minus + 2*pi_plus (integer):
    number_c12_b10_2n_piminus_2piplus = 0
    # number of interaction channel: nu + C12 -> B10 + 2p + 2*pi_minus + pi_plus (integer):
    number_c12_b10_2p_2piminus_piplus = 0
    # number of interaction channel: nu + C12 -> B10 + 2p + 3*pi_minus + 2*pi_plus (integer):
    number_c12_b10_2p_3piminus_2piplus = 0
    # number of interaction channel: nu + C12 -> B10 + p + n + 2*pi_minus + 2*pi_plus (integer):
    number_c12_b10_p_n_2piminus_2piplus = 0
    # number of OTHER interaction channels: nu + C12 -> B10 + ...:
    number_c12_b10_other = 0

    """ C10 """
    # number of interaction channel: nu + C12 -> C10 + 2n (integer):
    number_c12_c10_2n = 0
    # number of interaction channel: nu + C12 -> C10 + p + n + pi_minus (integer):
    number_c12_c10_p_n_piminus = 0
    # number of interaction channel: nu + C12 -> C10 + p + n + 2*pi_minus + pi_plus (integer):
    number_c12_c10_p_n_2piminus_piplus = 0
    # number of interaction channel: nu + C12 -> C10 + 2n + pi_minus + pi_plus (integer):
    number_c12_c10_2n_piminus_piplus = 0
    # number of interaction channel: nu + C12 -> C10 + 2p + 2*pi_minus (integer):
    number_c12_c10_2p_2piminus = 0
    # number of OTHER interaction channels: nu + C12 -> C10 + ...:
    number_c12_c10_other = 0

    """ Be10 """
    # number of interaction channel: nu + C12 -> Be10 + 2*p (integer):
    number_c12_be10_2p = 0
    # number of interaction channel: nu + C12 -> Be10 + p + n + pi_plus (integer):
    number_c12_be10_p_n_piplus = 0
    # number of interaction channel: nu + C12 -> Be10 + p + n + pi_minus + 2*pi_plus (integer):
    number_c12_be10_p_n_piminus_2piplus = 0
    # number of interaction channel: nu + C12 -> Be10 + 2*p + pi_minus + pi_plus (integer):
    number_c12_be10_2p_piminus_piplus = 0
    # number of interaction channel: nu + C12 -> Be10 + 2*n + 2*pi_plus (integer):
    number_c12_be10_2n_2piplus = 0
    # number of interaction channel: nu + C12 -> Be10 + p + n + 2*pi_minus + 3*pi_plus (integer):
    number_c12_be10_p_n_2piminus_3piplus = 0
    # number of interaction channel: nu + C12 -> Be10 + 2*p + 2*pi_minus + 2*pi_plus (integer):
    number_c12_be10_2p_2piminus_2piplus = 0
    # number of interaction channel: nu + C12 -> Be10 + 2*p + 3*pi_minus + 3*pi_plus (integer):
    number_c12_be10_2p_3piminus_3piplus = 0
    # number of OTHER interaction channels: nu + C12 -> Be10 + ...:
    number_c12_be10_other = 0

    """ B9 """
    # number of interaction channel: nu + C12 -> B9 + p + 2*n (integer):
    number_c12_b9_p_2n = 0
    # number of interaction channel: nu + C12 -> B9 + p + 2n + pi_minus + pi_plus (integer):
    number_c12_b9_p_2n_piminus_piplus = 0
    # number of interaction channel: nu + C12 -> B9 + 2p + n + 3*pi_minus + 2*pi_plus (integer):
    number_c12_b9_2p_n_3piminus_2piplus = 0
    # number of interaction channel: nu + C12 -> B9 + 2p + n + pi_minus (integer):
    number_c12_b9_2p_n_piminus = 0
    # number of interaction channel: nu + C12 -> B9 + 3n + pi_plus (integer):
    number_c12_b9_3n_piplus = 0
    # number of interaction channel: nu + C12 -> B9 + p + 2n + 2*pi_minus + 2*pi_plus (integer):
    number_c12_b9_p_2n_2piminus_2piplus = 0
    # number of interaction channel: nu + C12 -> B9 + 2p + n + 2*pi_minus+ pi_plus (integer):
    number_c12_b9_2p_n_2piminus_piplus = 0
    # number of OTHER interaction channels: nu + C12 -> B9 + ...
    number_c12_b9_other = 0

    """ Be9 """
    # number of interaction channel: nu + C12 -> Be9 + 2*p + n (integer):
    number_c12_be9_2p_n = 0
    # number of interaction channel: nu + C12 -> Be9 + p + 2n + pi_plus (integer):
    number_c12_be9_p_2n_piplus = 0
    # number of interaction channel: nu + C12 -> Be9 + 3p + pi_minus (integer):
    number_c12_be9_3p_piminus = 0
    # number of interaction channel: nu + C12 -> Be9 + p + 2n + pi_minus + 2*pi_plus (integer):
    number_c12_be9_p_2n_piminus_2piplus = 0
    # number of interaction channel: nu + C12 -> Be9 + 2p + n + pi_minus + pi_plus (integer):
    number_c12_be9_2p_n_piminus_piplus = 0
    # number of interaction channel: nu + C12 -> Be9 + 2p + n + 3*pi_minus + 3*pi_plus (integer):
    number_c12_be9_2p_n_3piminus_3piplus = 0
    # number of interaction channel: nu + C12 -> Be9 + 2p + n + 2*pi_minus + 2*pi_plus (integer):
    number_c12_be9_2p_n_2piminus_2piplus = 0
    # number of interaction channel: nu + C12 -> Be9 + 3n + 2*pi_plus (integer):
    number_c12_be9_3n_2piplus = 0
    # number of interaction channel: nu + C12 -> Be9 + 3p + 2*pi_minus + pi_plus (integer):
    number_c12_be9_3p_2piminus_piplus = 0
    # number of OTHER interaction channels: nu + C12 -> Be9 + ...:
    number_c12_be9_other = 0

    """ Be8 """
    # number of interaction channel: nu + C12 -> Be8 + 2p + 2n (integer):
    number_c12_be8_2p_2n = 0
    # number of interaction channel: nu + C12 -> Be8 + 3p + n + pi_minus (integer):
    number_c12_be8_3p_n_piminus = 0
    # number of interaction channel: nu + C12 -> Be8 + p + 3n + pi_plus (integer):
    number_c12_be8_p_3n_piplus = 0
    # number of interaction channel: nu + C12 -> Be8 + 2p + 2n + 2*pi_minus + 2*pi_plus (integer):
    number_c12_be8_2p_2n_2piminus_2piplus = 0
    # number of interaction channel: nu + C12 -> Be8 + 4n + 2*pi_plus (integer):
    number_c12_be8_4n_2piplus = 0
    # number of interaction channel: nu + C12 -> Be8 + 2p + 2n + pi_minus * pi_plus (integer):
    number_c12_be8_2p_2n_piminus_piplus = 0
    # number of interaction channel: nu + C12 -> Be8 + 3p + n + 2*pi_minus + pi_plus (integer):
    number_c12_be8_3p_n_2piminus_piplus = 0
    # number of interaction channel: nu + C12 -> Be8 + 4p + 2*pi_minus (integer):
    number_c12_be8_4p_2piminus = 0
    # number of OTHER interaction channels: nu + C12 -> Be8 + ...:
    number_c12_be8_other = 0

    """ C9 """
    # number of interaction channel: nu + C12 -> C9 + p + 2n + pi_minus (integer):
    number_c12_c9_p_2n_piminus = 0
    # number of interaction channel: nu + C12 -> C9 + 3n (integer):
    number_c12_c9_3n = 0
    # number of interaction channel: nu + C12 -> C9 + 2p + n + 2*pi_minus (integer):
    number_c12_c9_2p_n_2piminus = 0
    # number of interaction channel: nu + C12 -> C9 + 3n + 2*pi_minus + 2*pi_plus (integer):
    number_c12_c9_3n_2piminus_2piplus = 0
    # number of OTHER interaction channels: nu + C12 -> C9 + ...:
    number_c12_c9_other = 0

    """ Be7 """
    # number of interaction channel: nu + C12 -> Be7 + 2p + 3n (integer):
    number_c12_be7_2p_3n = 0
    # number of interaction channel: nu + C12 -> Be7 + p + 4n + pi_plus (integer):
    number_c12_be7_p_4n_piplus = 0
    # number of interaction channel: nu + C12 -> Be7 + 2p + 3n + 2*pi_minus + 2*pi_plus (integer):
    number_c12_be7_2p_3n_2piminus_2piplus = 0
    # number of interaction channel: nu + C12 -> Be7 + 3p + 2n + pi_minus (integer):
    number_c12_be7_3p_2n_piminus = 0
    # number of interaction channel: nu + C12 -> Be7 + 4p + n + 2*pi_minus (integer):
    number_c12_be7_4p_n_2piminus = 0
    # number of interaction channel: nu + C12 -> Be7 + 3p + 2n + 2*pi_minus + pi_plus (integer):
    number_c12_be7_3p_2n_2piminus_piplus = 0
    # number of OTHER interaction channels: nu + C12 -> Be7 + ...:
    number_c12_be7_other = 0

    """ Li6 """
    # number of interaction channel: nu + C12 -> Li6 + 3p + 3n (integer):
    number_c12_li6_3p_3n = 0
    # number of interaction channel: nu + C12 -> Li6 + 2p + 4n + pi_plus (integer):
    number_c12_li6_2p_4n_piplus = 0
    # number of interaction channel: nu + C12 -> Li6 + 5p + n + 2*pi_minus (integer):
    number_c12_li6_5p_n_2piminus = 0
    # number of interaction channel: nu + C12 -> Li6 + 2p + 4n + pi_minus + 2*pi_plus (integer):
    number_c12_li6_2p_4n_piminus_2piplus = 0
    # number of interaction channel: nu + C12 -> Li6 + 4p + 2n + pi_minus (integer):
    number_c12_li6_4p_2n_piminus = 0
    # number of interaction channel: nu + C12 -> Li6 + 3p + 3n + pi_minus + pi_plus (integer):
    number_c12_li6_3p_3n_piminus_piplus = 0
    # number of OTHER interaction channels: nu + C12 -> Li6 + ...:
    number_c12_li6_other = 0

    """ Li8 """
    # number of interaction channel: nu + C12 -> Li8 + 3p + n (integer):
    number_c12_li8_3p_n = 0
    # number of interaction channel: nu + C12 -> Li8 + 4p + pi_minus (integer):
    number_c12_li8_4p_piminus = 0
    # number of interaction channel: nu + C12 -> Li8 + 4p + 2*pi_minus + pi_plus (integer):
    number_c12_li8_4p_2piminus_piplus = 0
    # number of interaction channel: nu + C12 -> Li8 + 2p + 2n + pi_plus (integer):
    number_c12_li8_2p_2n_piplus = 0
    # number of interaction channel: nu + C12 -> Li8 + 3p + n + pi_minus + pi_plus (integer):
    number_c12_li8_3p_n_piminus_piplus = 0
    # number of OTHER interaction channels: nu + C12 -> Li8 + ...:
    number_c12_li8_other = 0

    """ Li7 """
    # number of interaction channel: nu + C12 -> Li7 + 2p + 3n + pi_plus (integer):
    number_c12_li7_2p_3n_piplus = 0
    # number of interaction channel: nu + C12 -> Li7 + 4p + n + pi_minus (integer):
    number_c12_li7_4p_n_piminus = 0
    # number of interaction channel: nu + C12 -> Li7 + 3p + 2n (integer):
    number_c12_li7_3p_2n = 0
    # number of interaction channel: nu + C12 -> Li7 + 3p + 2n + pi_minus + pi_plus (integer):
    number_c12_li7_3p_2n_piminus_piplus = 0
    # number of interaction channel: nu + C12 -> Li7 + 4p + n + 2*pi_minus + pi_plus (integer):
    number_c12_li7_4p_n_2piminus_piplus = 0
    # number of interaction channel: nu + C12 -> Li7 + 2p + 3n + pi_minus + 2*pi_plus (integer):
    number_c12_li7_2p_3n_piminus_2piplus = 0
    # number of OTHER interaction channels: nu + C12 -> Li7 + ...:
    number_c12_li7_other = 0

    """ B8 """
    # number of interaction channel: nu + C12 -> B8 + p + 3n (integer):
    number_c12_b8_p_3n = 0
    # number of interaction channel: nu + C12 -> B8 + p + 3n + pi_minus + pi_plus (integer):
    number_c12_b8_p_3n_piminus_piplus = 0
    # number of interaction channel: nu + C12 -> B8 + 2p + 2n + 2*pi_minus + pi_plus (integer):
    number_c12_b8_2p_2n_2piminus_piplus = 0
    # number of interaction channel: nu + C12 -> B8 + 2p + 2n + pi_minus (integer):
    number_c12_b8_2p_2n_piminus = 0
    # number of interaction channel: nu + C12 -> B8 + 4n + pi_plus (integer):
    number_c12_b8_4n_piplus = 0
    # number of OTHER interaction channels: nu + C12 -> B8 + ...:
    number_c12_b8_other = 0

    """ Li9 """
    # number of interaction channel: nu + C12 -> Li9 + 2p + n + pi_plus (integer):
    number_c12_li9_2p_n_piplus = 0
    # number of interaction channel: nu + C12 -> Li9 + 3p (integer).
    number_c12_li9_3p = 0
    # number of interaction channel: nu + C12 -> Li9 + 3p + pi_minus + pi_plus (integer):
    number_c12_li9_3p_piminus_piplus = 0
    # number of interaction channel: nu + C12 -> Li9 + 2p + n + pi_minus + 2*pi_plus (integer):
    number_c12_li9_2p_n_piminus_2piplus = 0
    # number of interaction channel: nu + C12 -> Li9 + p + 2n + pi_minus + 3*pi_plus (integer):
    number_c12_li9_p_2n_piminus_3piplus = 0
    # number of interaction channel: nu + C12 -> Li9 + ...:
    number_c12_li9_other = 0

    # loop over all entries of the array:
    for index in range(number_entries):

        # check, if target is C12 (PDG ID = 1000060120):
        if target_pdg[index] == 1000060120:

            # check the PDG ID of the created isotopes:
            if isotope_pdg[index] == 1000050110:
                # for B11:
                num_p, num_n, num_pi_minus, num_pi_plus = get_number_of_particles_of_channelid(channel_id[index])

                if num_p == 1 and num_n == 0 and num_pi_minus == 0 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> B11 + proton:
                    number_c12_b11_p = number_c12_b11_p + 1

                elif num_p == 0 and num_n == 1 and num_pi_minus == 0 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> B11 + n + pi_plus:
                    number_c12_b11_n_piplus = number_c12_b11_n_piplus + 1

                elif num_p == 0 and num_n == 1 and num_pi_minus == 1 and num_pi_plus == 2:
                    # interaction channel: nu + C12 -> B11 + n + pi_minus * 2*pi_plus:
                    number_c12_b11_n_piminus_2piplus = number_c12_b11_n_piminus_2piplus + 1

                elif num_p == 1 and num_n == 0 and num_pi_minus == 1 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> B11 + p + pi_minus + pi_plus:
                    number_c12_b11_p_piminus_piplus = number_c12_b11_p_piminus_piplus + 1

                elif num_p == 1 and num_n == 0 and num_pi_minus == 2 and num_pi_plus == 2:
                    # interaction channel: nu + C12 -> B11 + p + 2*pi_minus + 2*pi_plus:
                    number_c12_b11_p_2piminus_2piplus = number_c12_b11_p_2piminus_2piplus + 1

                elif num_p == 0 and num_n == 0 and num_pi_minus == 0 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> B11 + pi_plus:
                    number_c12_b11_piplus = number_c12_b11_piplus + 1

                else:
                    number_c12_b11_other = number_c12_b11_other + 1
                    # print("new interaction channel with nu + C12 -> B11: channel ID = {0:.0f}"
                    #       .format(channel_id[index]))


            elif isotope_pdg[index] == 1000060110:
                # for C11:
                num_p, num_n, num_pi_minus, num_pi_plus = get_number_of_particles_of_channelid(channel_id[index])

                if num_p == 0 and num_n == 1 and num_pi_minus == 0 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> C11 + n:
                    number_c12_c11_n = number_c12_c11_n + 1

                elif num_p == 1 and num_n == 0 and num_pi_minus == 1 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> C11 + p + pi_minus:
                    number_c12_c11_p_piminus = number_c12_c11_p_piminus + 1

                elif num_p == 0 and num_n == 1 and num_pi_minus == 1 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> C11 + n + pi_minus + pi_plus:
                    number_c12_c11_n_piminus_piplus = number_c12_c11_n_piminus_piplus + 1

                elif num_p == 1 and num_n == 0 and num_pi_minus == 2 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> C11 + p + 2*pi_minus + pi_plus:
                    number_c12_c11_p_2piminus_piplus = number_c12_c11_p_2piminus_piplus + 1

                elif num_p == 1 and num_n == 0 and num_pi_minus == 3 and num_pi_plus == 2:
                    # interaction channel: nu + C12 -> C11 + p + 3*pi_minus + 2*pi_plus:
                    number_c12_c11_p_3piminus_2piplus = number_c12_c11_p_3piminus_2piplus + 1

                elif num_p == 0 and num_n == 1 and num_pi_minus == 2 and num_pi_plus == 2:
                    # interaction channel: nu + C12 -> C11 + n + 2*pi_minus + 2*pi_plus:
                    number_c12_c11_n_2piminus_2piplus = number_c12_c11_n_2piminus_2piplus + 1

                else:
                    number_c12_c11_other = number_c12_c11_other + 1
                    # print("new interaction channel with nu + C12 -> C11: channel ID = {0:.0f}"
                    #       .format(channel_id[index]))


            elif isotope_pdg[index] == 1000050100:
                # for B10:
                num_p, num_n, num_pi_minus, num_pi_plus = get_number_of_particles_of_channelid(channel_id[index])

                if num_p == 1 and num_n == 1 and num_pi_minus == 0 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> B10 + p + n:
                    number_c12_b10_p_n = number_c12_b10_p_n + 1

                elif num_p == 2 and num_n == 0 and num_pi_minus == 1 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> B10 + 2*p + pi_minus:
                    number_c12_b10_2p_piminus = number_c12_b10_2p_piminus + 1

                elif num_p == 1 and num_n == 1 and num_pi_minus == 1 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> B10 + p + n + pi_minus + pi_plus:
                    number_c12_b10_p_n_piminus_piplus = number_c12_b10_p_n_piminus_piplus + 1

                elif num_p == 0 and num_n == 2 and num_pi_minus == 0 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> B10 + 2*n + pi_plus:
                    number_c12_b10_2n_piplus = number_c12_b10_2n_piplus + 1

                elif num_p == 0 and num_n == 2 and num_pi_minus == 1 and num_pi_plus == 2:
                    # interaction channel: nu + C12 -> B10 + 2*n + pi_minus + 2*pi_plus:
                    number_c12_b10_2n_piminus_2piplus = number_c12_b10_2n_piminus_2piplus + 1

                elif num_p == 2 and num_n == 0 and num_pi_minus == 2 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> B10 + 2*p + 2*pi_minus + pi_plus:
                    number_c12_b10_2p_2piminus_piplus = number_c12_b10_2p_2piminus_piplus + 1

                elif num_p == 2 and num_n == 0 and num_pi_minus == 3 and num_pi_plus == 2:
                    # interaction channel: nu + C12 -> B10 + 2*p + 3*pi_minus + 2*pi_plus:
                    number_c12_b10_2p_3piminus_2piplus = number_c12_b10_2p_3piminus_2piplus + 1

                elif num_p == 1 and num_n == 1 and num_pi_minus == 2 and num_pi_plus == 2:
                    # interaction channel: nu + C12 -> B10 + p + n + 2*pi_minus + 2*pi_plus:
                    number_c12_b10_p_n_2piminus_2piplus = number_c12_b10_p_n_2piminus_2piplus + 1

                else:
                    number_c12_b10_other = number_c12_b10_other + 1
                    # print("new interaction channel with nu + C12 -> B10: channel ID = {0:.0f}"
                    #       .format(channel_id[index]))


            elif isotope_pdg[index] == 1000060100:
                # for C10:
                num_p, num_n, num_pi_minus, num_pi_plus = get_number_of_particles_of_channelid(channel_id[index])

                if num_p == 0 and num_n == 2 and num_pi_minus == 0 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> C10 + 2n:
                    number_c12_c10_2n = number_c12_c10_2n + 1

                elif num_p == 1 and num_n == 1 and num_pi_minus == 1 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> C10 + p + n + pi_minus:
                    number_c12_c10_p_n_piminus = number_c12_c10_p_n_piminus + 1

                elif num_p == 1 and num_n == 1 and num_pi_minus == 2 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> C10 + p + n + 2*pi_minus + pi_plus:
                    number_c12_c10_p_n_2piminus_piplus = number_c12_c10_p_n_2piminus_piplus + 1

                elif num_p == 0 and num_n == 2 and num_pi_minus == 1 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> C10 + 2*n + pi_minus + pi_plus:
                    number_c12_c10_2n_piminus_piplus = number_c12_c10_2n_piminus_piplus + 1

                elif num_p == 2 and num_n == 0 and num_pi_minus == 2 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> C10 + 2*p + 2*pi_minus:
                    number_c12_c10_2p_2piminus = number_c12_c10_2p_2piminus + 1

                else:
                    number_c12_c10_other = number_c12_c10_other + 1
                    # print("new interaction channel with nu + C12 -> C10: channel ID = {0:.0f}"
                    #       .format(channel_id[index]))


            elif isotope_pdg[index] == 1000040100:
                # for Be10:
                num_p, num_n, num_pi_minus, num_pi_plus = get_number_of_particles_of_channelid(channel_id[index])

                if num_p == 2 and num_n == 0 and num_pi_minus == 0 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> Be10 + 2*p:
                    number_c12_be10_2p = number_c12_be10_2p + 1

                elif num_p == 1 and num_n == 1 and num_pi_minus == 0 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> Be10 + p + n + pi_plus:
                    number_c12_be10_p_n_piplus = number_c12_be10_p_n_piplus + 1

                elif num_p == 1 and num_n == 1 and num_pi_minus == 1 and num_pi_plus == 2:
                    # interaction channel: nu + C12 -> Be10 + p + n + pi_minus + 2*pi_plus:
                    number_c12_be10_p_n_piminus_2piplus = number_c12_be10_p_n_piminus_2piplus + 1

                elif num_p == 2 and num_n == 0 and num_pi_minus == 1 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> Be10 + 2*p + pi_minus + pi_plus:
                    number_c12_be10_2p_piminus_piplus = number_c12_be10_2p_piminus_piplus + 1

                elif num_p == 0 and num_n == 2 and num_pi_minus == 0 and num_pi_plus == 2:
                    # interaction channel: nu + C12 -> Be10 + 2n + 2*pi_plus:
                    number_c12_be10_2n_2piplus = number_c12_be10_2n_2piplus + 1

                elif num_p == 1 and num_n == 1 and num_pi_minus == 2 and num_pi_plus == 3:
                    # interaction channel: nu + C12 -> Be10 + p + n + 2*pi_minus + 3*pi_plus:
                    number_c12_be10_p_n_2piminus_3piplus = number_c12_be10_p_n_2piminus_3piplus + 1

                elif num_p == 2 and num_n == 0 and num_pi_minus == 2 and num_pi_plus == 2:
                    # interaction channel: nu + C12 -> Be10 + 2p + 2*pi_minus + 2*pi_plus:
                    number_c12_be10_2p_2piminus_2piplus = number_c12_be10_2p_2piminus_2piplus + 1

                elif num_p == 2 and num_n == 0 and num_pi_minus == 3 and num_pi_plus == 3:
                    # interaction channel: nu + C12 -> Be10 + 2p + 3*pi_minus + 3*pi_plus:
                    number_c12_be10_2p_3piminus_3piplus = number_c12_be10_2p_3piminus_3piplus + 1

                else:
                    number_c12_be10_other = number_c12_be10_other + 1
                    # print("new interaction channel with nu + C12 -> Be10: channel ID = {0:.0f}"
                    #       .format(channel_id[index]))


            elif isotope_pdg[index] == 1000050090:
                # for B9:
                num_p, num_n, num_pi_minus, num_pi_plus = get_number_of_particles_of_channelid(channel_id[index])

                if num_p == 1 and num_n == 2 and num_pi_minus == 0 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> B9 + p + 2*n:
                    number_c12_b9_p_2n = number_c12_b9_p_2n + 1

                elif num_p == 1 and num_n == 2 and num_pi_minus == 1 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> B9 + p + 2n + pi_minus + pi_plus:
                    number_c12_b9_p_2n_piminus_piplus = number_c12_b9_p_2n_piminus_piplus + 0

                elif num_p == 2 and num_n == 1 and num_pi_minus == 3 and num_pi_plus == 2:
                    # interaction channel: nu + C12 -> B9 + 2p + n + 3*pi_minus + 2*pi_plus:
                    number_c12_b9_2p_n_3piminus_2piplus = number_c12_b9_2p_n_3piminus_2piplus + 0

                elif num_p == 2 and num_n == 1 and num_pi_minus == 1 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> B9 + 2p + n + pi_minus:
                    number_c12_b9_2p_n_piminus = number_c12_b9_2p_n_piminus + 0

                elif num_p == 0 and num_n == 3 and num_pi_minus == 0 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> B9 + 3n + pi_plus:
                    number_c12_b9_3n_piplus = number_c12_b9_3n_piplus + 0

                elif num_p == 1 and num_n == 2 and num_pi_minus == 2 and num_pi_plus == 2:
                    # interaction channel: nu + C12 -> B9 + p + 2n + 2*pi_minus + 2*pi_plus:
                    number_c12_b9_p_2n_2piminus_2piplus = number_c12_b9_p_2n_2piminus_2piplus + 0

                elif num_p == 2 and num_n == 1 and num_pi_minus == 2 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> B9 + 2p + n + 2*pi_minus + pi_plus:
                    number_c12_b9_2p_n_2piminus_piplus = number_c12_b9_2p_n_2piminus_piplus + 0

                else:
                    number_c12_b9_other = number_c12_b9_other + 1
                    # print("new interaction channel with nu + C12 -> B9: channel ID = {0:.0f}"
                    # .format(channel_id[index]))


            elif isotope_pdg[index] == 1000040090:
                # for Be9:
                num_p, num_n, num_pi_minus, num_pi_plus = get_number_of_particles_of_channelid(channel_id[index])

                if num_p == 2 and num_n == 1 and num_pi_minus == 0 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> Be9 + 2*p + n:
                    number_c12_be9_2p_n = number_c12_be9_2p_n + 1

                elif num_p == 1 and num_n == 2 and num_pi_minus == 0 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> Be9 + p + 2*n + pi_plus:
                    number_c12_be9_p_2n_piplus = number_c12_be9_p_2n_piplus + 1

                elif num_p == 3 and num_n == 0 and num_pi_minus == 1 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> Be9 + 3p + pi_minus:
                    number_c12_be9_3p_piminus = number_c12_be9_3p_piminus + 1

                elif num_p == 1 and num_n == 2 and num_pi_minus == 1 and num_pi_plus == 2:
                    # interaction channel: nu + C12 -> Be9 + p + 2*n + pi_minus + 2*pi_plus:
                    number_c12_be9_p_2n_piminus_2piplus = number_c12_be9_p_2n_piminus_2piplus + 1

                elif num_p == 2 and num_n == 1 and num_pi_minus == 1 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> Be9 + 2*p + n + pi_minus + pi_plus:
                    number_c12_be9_2p_n_piminus_piplus = number_c12_be9_2p_n_piminus_piplus + 1

                elif num_p == 2 and num_n == 1 and num_pi_minus == 3 and num_pi_plus == 3:
                    # interaction channel: nu + C12 -> Be9 + 2*p + n + 3*pi_minus + 3*pi_plus:
                    number_c12_be9_2p_n_3piminus_3piplus = number_c12_be9_2p_n_3piminus_3piplus + 1

                elif num_p == 2 and num_n == 1 and num_pi_minus == 2 and num_pi_plus == 2:
                    # interaction channel: nu + C12 -> Be9 + 2*p + n + 2*pi_minus + 2*pi_plus:
                    number_c12_be9_2p_n_2piminus_2piplus = number_c12_be9_2p_n_2piminus_2piplus + 1

                elif num_p == 0 and num_n == 3 and num_pi_minus == 0 and num_pi_plus == 2:
                    # interaction channel: nu + C12 -> Be9 + 3*n + 2*pi_plus:
                    number_c12_be9_3n_2piplus = number_c12_be9_3n_2piplus + 1

                elif num_p == 3 and num_n == 0 and num_pi_minus == 2 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> Be9 + 3*p + 2*pi_minus + pi_plus:
                    number_c12_be9_3p_2piminus_piplus = number_c12_be9_3p_2piminus_piplus + 1

                else:
                    number_c12_be9_other = number_c12_be9_other + 1
                    # print("new interaction channel with nu + C12 -> Be9: channel ID = {0:.0f}"
                    #       .format(channel_id[index]))


            elif isotope_pdg[index] == 1000040080:
                # for Be8:
                num_p, num_n, num_pi_minus, num_pi_plus = get_number_of_particles_of_channelid(channel_id[index])

                if num_p == 2 and num_n == 2 and num_pi_minus == 0 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> Be8 + 2*p + 2*n:
                    number_c12_be8_2p_2n = number_c12_be8_2p_2n + 1

                elif num_p == 3 and num_n == 1 and num_pi_minus == 1 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> Be8 + 3*p + n + pi_minus:
                    number_c12_be8_3p_n_piminus = number_c12_be8_3p_n_piminus + 1

                elif num_p == 1 and num_n == 3 and num_pi_minus == 0 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> Be8 + p + 3*n + pi_plus:
                    number_c12_be8_p_3n_piplus = number_c12_be8_p_3n_piplus + 1

                elif num_p == 2 and num_n == 2 and num_pi_minus == 2 and num_pi_plus == 2:
                    # interaction channel: nu + C12 -> Be8 + 2*p + 2*n + 2*pi_minus + 2*pi_plus:
                    number_c12_be8_2p_2n_2piminus_2piplus = number_c12_be8_2p_2n_2piminus_2piplus + 1

                elif num_p == 0 and num_n == 4 and num_pi_minus == 0 and num_pi_plus == 2:
                    # interaction channel: nu + C12 -> Be8 + 4*n + 2*pi_plus:
                    number_c12_be8_4n_2piplus = number_c12_be8_4n_2piplus + 1

                elif num_p == 2 and num_n == 2 and num_pi_minus == 1 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> Be8 + 2*p + 2*n + pi_minus + pi_plus:
                    number_c12_be8_2p_2n_piminus_piplus = number_c12_be8_2p_2n_piminus_piplus + 1

                elif num_p == 3 and num_n == 1 and num_pi_minus == 2 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> Be8 + 3*p + n + 2*pi_minus + pi_plus:
                    number_c12_be8_3p_n_2piminus_piplus = number_c12_be8_3p_n_2piminus_piplus + 1

                elif num_p == 4 and num_n == 0 and num_pi_minus == 2 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> Be8 + 4*p + 2*pi_minus:
                    number_c12_be8_4p_2piminus = number_c12_be8_4p_2piminus + 1

                else:
                    number_c12_be8_other = number_c12_be8_other + 1
                    # print("new interaction channel with nu + C12 -> Be8: channel ID = {0:.0f}"
                    #       .format(channel_id[index]))


            elif isotope_pdg[index] == 1000060090:
                # for C9:
                num_p, num_n, num_pi_minus, num_pi_plus = get_number_of_particles_of_channelid(channel_id[index])

                if num_p == 1 and num_n == 2 and num_pi_minus == 1 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> C9 + p + 2*n + pi_minus:
                    number_c12_c9_p_2n_piminus = number_c12_c9_p_2n_piminus + 1

                elif num_p == 0 and num_n == 3 and num_pi_minus == 0 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> C9 + 3*n:
                    number_c12_c9_3n = number_c12_c9_3n + 1

                elif num_p == 2 and num_n == 1 and num_pi_minus == 2 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> C9 + 2*p + n + 2*pi_minus:
                    number_c12_c9_2p_n_2piminus = number_c12_c9_2p_n_2piminus + 1

                elif num_p == 0 and num_n == 3 and num_pi_minus == 2 and num_pi_plus == 2:
                    # interaction channel: nu + C12 -> C9 + 3*n + 2*pi_minus + 2*pi_plus:
                    number_c12_c9_3n_2piminus_2piplus = number_c12_c9_3n_2piminus_2piplus + 1

                else:
                    number_c12_c9_other = number_c12_c9_other + 1
                    # print("new interaction channel with nu + C12 -> C9: channel ID = {0:.0f}"
                    # .format(channel_id[index]))


            elif isotope_pdg[index] == 1000040070:
                # for Be7:
                num_p, num_n, num_pi_minus, num_pi_plus = get_number_of_particles_of_channelid(channel_id[index])

                if num_p == 2 and num_n == 3 and num_pi_minus == 0 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> Be7 + 2*p + 3*n:
                    number_c12_be7_2p_3n = number_c12_be7_2p_3n + 1

                elif num_p == 1 and num_n == 4 and num_pi_minus == 0 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> Be7 + p + 4*n + pi_plus:
                    number_c12_be7_p_4n_piplus = number_c12_be7_p_4n_piplus + 1

                elif num_p == 2 and num_n == 3 and num_pi_minus == 2 and num_pi_plus == 2:
                    # interaction channel: nu + C12 -> Be7 + 2*p + 3*n + 2*pi_minus + 2*pi_plus:
                    number_c12_be7_2p_3n_2piminus_2piplus = number_c12_be7_2p_3n_2piminus_2piplus + 1

                elif num_p == 3 and num_n == 2 and num_pi_minus == 1 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> Be7 + 3*p + 2*n + pi_minus:
                    number_c12_be7_3p_2n_piminus = number_c12_be7_3p_2n_piminus + 1

                elif num_p == 4 and num_n == 1 and num_pi_minus == 2 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> Be7 + 4*p + n + 2*pi_minus:
                    number_c12_be7_4p_n_2piminus = number_c12_be7_4p_n_2piminus + 1

                elif num_p == 3 and num_n == 2 and num_pi_minus == 2 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> Be7 + 3*p + 2*n + 2*pi_minus + pi_plus:
                    number_c12_be7_3p_2n_2piminus_piplus = number_c12_be7_3p_2n_2piminus_piplus + 1

                else:
                    number_c12_be7_other = number_c12_be7_other + 1
                    # print("new interaction channel with nu + C12 -> Be7: channel ID = {0:.0f}"
                    #       .format(channel_id[index]))


            elif isotope_pdg[index] == 1000030060:
                # for Li6:
                num_p, num_n, num_pi_minus, num_pi_plus = get_number_of_particles_of_channelid(channel_id[index])

                if num_p == 3 and num_n == 3 and num_pi_minus == 0 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> Li6 + 3*p + 3*n:
                    number_c12_li6_3p_3n = number_c12_li6_3p_3n + 1

                elif num_p == 2 and num_n == 4 and num_pi_minus == 0 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> Li6 + 2*p + 4*n + pi_plus:
                    number_c12_li6_2p_4n_piplus = number_c12_li6_2p_4n_piplus + 1

                elif num_p == 5 and num_n == 1 and num_pi_minus == 2 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> Li6 + 5*p + n + 2*pi_minus:
                    number_c12_li6_5p_n_2piminus = number_c12_li6_5p_n_2piminus + 1

                elif num_p == 2 and num_n == 4 and num_pi_minus == 1 and num_pi_plus == 2:
                    # interaction channel: nu + C12 -> Li6 + 2*p + 4*n + pi_minus + 2*pi_plus:
                    number_c12_li6_2p_4n_piminus_2piplus = number_c12_li6_2p_4n_piminus_2piplus + 1

                elif num_p == 4 and num_n == 2 and num_pi_minus == 1 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> Li6 + 4*p + 2*n + pi_minus:
                    number_c12_li6_4p_2n_piminus = number_c12_li6_4p_2n_piminus + 1

                elif num_p == 3 and num_n == 3 and num_pi_minus == 1 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> Li6 + 3*p + 3*n + pi_minus + pi_plus:
                    number_c12_li6_3p_3n_piminus_piplus = number_c12_li6_3p_3n_piminus_piplus + 1

                else:
                    number_c12_li6_other = number_c12_li6_other + 1
                    # print("new interaction channel with nu + C12 -> Li6: channel ID = {0:.0f}"
                    #       .format(channel_id[index]))


            elif isotope_pdg[index] == 1000030080:
                # for Li8:
                num_p, num_n, num_pi_minus, num_pi_plus = get_number_of_particles_of_channelid(channel_id[index])

                if num_p == 3 and num_n == 1 and num_pi_minus == 0 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> Li8 + 3*p + n:
                    number_c12_li8_3p_n = number_c12_li8_3p_n + 1

                elif num_p == 4 and num_n == 0 and num_pi_minus == 1 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> Li8 + 4*p + pi_minus:
                    number_c12_li8_4p_piminus = number_c12_li8_4p_piminus + 1

                elif num_p == 4 and num_n == 0 and num_pi_minus == 2 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> Li8 + 4*p + 2*pi_minus + pi_plus:
                    number_c12_li8_4p_2piminus_piplus = number_c12_li8_4p_2piminus_piplus + 1

                elif num_p == 2 and num_n == 2 and num_pi_minus == 0 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> Li8 + 2*p + 2*n + pi_plus:
                    number_c12_li8_2p_2n_piplus = number_c12_li8_2p_2n_piplus + 1

                elif num_p == 3 and num_n == 1 and num_pi_minus == 1 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> Li8 + 3*p + n + pi_minus + pi_plus:
                    number_c12_li8_3p_n_piminus_piplus = number_c12_li8_3p_n_piminus_piplus + 1

                else:
                    number_c12_li8_other = number_c12_li8_other + 1
                    # print("new interaction channel with nu + C12 -> Li8: channel ID = {0:.0f}"
                    #       .format(channel_id[index]))


            elif isotope_pdg[index] == 1000030070:
                # for Li7:
                num_p, num_n, num_pi_minus, num_pi_plus = get_number_of_particles_of_channelid(channel_id[index])

                if num_p == 2 and num_n == 3 and num_pi_minus == 0 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> Li7 + 2*p + 3*n + pi_plus:
                    number_c12_li7_2p_3n_piplus = number_c12_li7_2p_3n_piplus + 1

                elif num_p == 4 and num_n == 1 and num_pi_minus == 1 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> Li7 + 4*p + n + pi_minus:
                    number_c12_li7_4p_n_piminus = number_c12_li7_4p_n_piminus + 1

                elif num_p == 3 and num_n == 2 and num_pi_minus == 0 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> Li7 + 3*p + 2*n:
                    number_c12_li7_3p_2n = number_c12_li7_3p_2n + 1

                elif num_p == 3 and num_n == 2 and num_pi_minus == 1 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> Li7 + 3*p + 2*n + pi_minus + pi_plus:
                    number_c12_li7_3p_2n_piminus_piplus = number_c12_li7_3p_2n_piminus_piplus + 1

                elif num_p == 4 and num_n == 1 and num_pi_minus == 2 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> Li7 + 4*p + n + 2*pi_minus + pi_plus:
                    number_c12_li7_4p_n_2piminus_piplus = number_c12_li7_4p_n_2piminus_piplus + 1

                elif num_p == 2 and num_n == 3 and num_pi_minus == 1 and num_pi_plus == 2:
                    # interaction channel: nu + C12 -> Li7 + 2*p + 3*n + pi_minus + 2*pi_plus:
                    number_c12_li7_2p_3n_piminus_2piplus = number_c12_li7_2p_3n_piminus_2piplus + 1

                else:
                    number_c12_li7_other = number_c12_li7_other + 1
                    # print("new interaction channel with nu + C12 -> Li7: channel ID = {0:.0f}"
                    #       .format(channel_id[index]))


            elif isotope_pdg[index] == 1000050080:
                # for B8:
                num_p, num_n, num_pi_minus, num_pi_plus = get_number_of_particles_of_channelid(channel_id[index])

                if num_p == 1 and num_n == 3 and num_pi_minus == 0 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> B8 + p + 3*n:
                    number_c12_b8_p_3n = number_c12_b8_p_3n + 1

                elif num_p == 1 and num_n == 3 and num_pi_minus == 1 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> B8 + p + 3*n + pi_minus + pi_plus:
                    number_c12_b8_p_3n_piminus_piplus = number_c12_b8_p_3n_piminus_piplus + 1

                elif num_p == 2 and num_n == 2 and num_pi_minus == 2 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> B8 + 2*p + 2*n + 2*pi_minus + pi_plus:
                    number_c12_b8_2p_2n_2piminus_piplus = number_c12_b8_2p_2n_2piminus_piplus + 1

                elif num_p == 2 and num_n == 2 and num_pi_minus == 1 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> B8 + 2*p + 2*n + pi_minus:
                    number_c12_b8_2p_2n_piminus = number_c12_b8_2p_2n_piminus + 1

                elif num_p == 0 and num_n == 4 and num_pi_minus == 0 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> B8 + 4*n + pi_plus:
                    number_c12_b8_4n_piplus = number_c12_b8_4n_piplus + 1

                else:
                    number_c12_b8_other = number_c12_b8_other + 1
                    # print("new interaction channel with nu + C12 -> B8: channel ID = {0:.0f}"
                    #       .format(channel_id[index]))


            elif isotope_pdg[index] == 1000030090:
                # for Li9:
                num_p, num_n, num_pi_minus, num_pi_plus = get_number_of_particles_of_channelid(channel_id[index])

                if num_p == 2 and num_n == 1 and num_pi_minus == 0 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> Li9 + 2*p + n + pi_plus:
                    number_c12_li9_2p_n_piplus = number_c12_li9_2p_n_piplus + 1

                elif num_p == 3 and num_n == 0 and num_pi_minus == 0 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> Li9 + 3*p:
                    number_c12_li9_3p = number_c12_li9_3p + 1

                elif num_p == 3 and num_n == 0 and num_pi_minus == 1 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> Li9 + 3*p + pi_minus + pi_plus:
                    number_c12_li9_3p_piminus_piplus = number_c12_li9_3p_piminus_piplus + 1

                elif num_p == 2 and num_n == 1 and num_pi_minus == 1 and num_pi_plus == 2:
                    # interaction channel: nu + C12 -> Li9 + 2*p + n + pi_minus + 2*pi_plus:
                    number_c12_li9_2p_n_piminus_2piplus = number_c12_li9_2p_n_piminus_2piplus + 1

                elif num_p == 1 and num_n == 2 and num_pi_minus == 1 and num_pi_plus == 3:
                    # interaction channel: nu + C12 -> Li9 + p + 2*n + pi_minus + 3*pi_plus:
                    number_c12_li9_p_2n_piminus_3piplus = number_c12_li9_p_2n_piminus_3piplus + 1

                else:
                    number_c12_li9_other = number_c12_li9_other + 1
                    # print("new interaction channel with nu + C12 -> Li9: channel ID = {0:.0f}"
                    #       .format(channel_id[index]))


            else:
                print("other isotope than expected: {0:.0f}, corresponding channel ID = {1:.0f}"
                      .format(isotope_pdg[index], channel_id[index]))
                # lala = 1

        else:
            # print("other target than C12: target PDG = {0:.0f}, channel ID = {1:.0f}, isotope PDG = {2:.0f}"
                  # .format(target_pdg[index], channel_id[index], isotope_pdg[index]))
            number_no_c12 = number_no_c12 + 1


    """ calculate the fraction of the different NC interaction channels in PERCENT (float): """
    # B11:
    frac_c12_b11_p = float(number_c12_b11_p) / float(number_entries) * 100
    frac_c12_b11_n_piplus = float(number_c12_b11_n_piplus) / float(number_entries) * 100
    frac_c12_b11_n_piminus_2piplus = float(number_c12_b11_n_piminus_2piplus) / float(number_entries) * 100
    frac_c12_b11_p_piminus_piplus = float(number_c12_b11_p_piminus_piplus) / float(number_entries) * 100
    frac_c12_b11_p_2piminus_2piplus = float(number_c12_b11_p_2piminus_2piplus) / float(number_entries) * 100
    frac_c12_b11_piplus = float(number_c12_b11_piplus) / float(number_entries) * 100
    frac_c12_b11_other = float(number_c12_b11_other) / float(number_entries) * 100

    # C11:
    frac_c12_c11_n = float(number_c12_c11_n) / float(number_entries) * 100
    frac_c12_c11_p_piminus = float(number_c12_c11_p_piminus) / float(number_entries) * 100
    frac_c12_c11_n_piminus_piplus = float(number_c12_c11_n_piminus_piplus) / float(number_entries) * 100
    frac_c12_c11_p_2piminus_piplus = float(number_c12_c11_p_2piminus_piplus) / float(number_entries) * 100
    frac_c12_c11_p_3piminus_2piplus = float(number_c12_c11_p_3piminus_2piplus) / float(number_entries) * 100
    frac_c12_c11_n_2piminus_2piplus = float(number_c12_c11_n_2piminus_2piplus) / float(number_entries) * 100
    frac_c12_c11_other = float(number_c12_c11_other) / float(number_entries) * 100

    # B10:
    frac_c12_b10_p_n = float(number_c12_b10_p_n) / float(number_entries) * 100
    frac_c12_b10_2p_piminus = float(number_c12_b10_2p_piminus) / float(number_entries) * 100
    frac_c12_b10_p_n_piminus_piplus = float(number_c12_b10_p_n_piminus_piplus) / float(number_entries) * 100
    frac_c12_b10_2n_piplus = float(number_c12_b10_2n_piplus) / float(number_entries) * 100
    frac_c12_b10_2n_piminus_2piplus = float(number_c12_b10_2n_piminus_2piplus) / float(number_entries) * 100
    frac_c12_b10_2p_2piminus_piplus = float(number_c12_b10_2p_2piminus_piplus) / float(number_entries) * 100
    frac_c12_b10_2p_3piminus_2piplus = float(number_c12_b10_2p_3piminus_2piplus) / float(number_entries) * 100
    frac_c12_b10_p_n_2piminus_2piplus = float(number_c12_b10_p_n_2piminus_2piplus) / float( number_entries) * 100
    frac_c12_b10_other = float(number_c12_b10_other) / float(number_entries) * 100

    # C10:
    frac_c12_c10_2n = float(number_c12_c10_2n) / float(number_entries) * 100
    frac_c12_c10_p_n_piminus = float(number_c12_c10_p_n_piminus) / float(number_entries) * 100
    frac_c12_c10_p_n_2piminus_piplus = float(number_c12_c10_p_n_2piminus_piplus) / float(number_entries) * 100
    frac_c12_c10_2n_piminus_piplus = float(number_c12_c10_2n_piminus_piplus) / float(number_entries) * 100
    frac_c12_c10_2p_2piminus = float(number_c12_c10_2p_2piminus) / float(number_entries) * 100
    frac_c12_c10_other = float(number_c12_c10_other) / float(number_entries) * 100

    # Be10:
    frac_c12_be10_2p = float(number_c12_be10_2p) / float(number_entries) * 100
    frac_c12_be10_p_n_piplus = float(number_c12_be10_p_n_piplus) / float(number_entries) * 100
    frac_c12_be10_p_n_piminus_2piplus = float(number_c12_be10_p_n_piminus_2piplus) / float(number_entries) * 100
    frac_c12_be10_2p_piminus_piplus = float(number_c12_be10_2p_piminus_piplus) / float(number_entries) * 100
    frac_c12_be10_2n_2piplus = float(number_c12_be10_2n_2piplus) / float(number_entries) * 100
    frac_c12_be10_p_n_2piminus_3piplus = float(number_c12_be10_p_n_2piminus_3piplus) / float(number_entries) * 100
    frac_c12_be10_2p_2piminus_2piplus = float(number_c12_be10_2p_2piminus_2piplus) / float(number_entries) * 100
    frac_c12_be10_2p_3piminus_3piplus = float(number_c12_be10_2p_3piminus_3piplus) / float(number_entries) * 100
    frac_c12_be10_other = float(number_c12_be10_other) / float(number_entries) * 100

    # B9:
    frac_c12_b9_p_2n = float(number_c12_b9_p_2n) / float(number_entries) * 100
    frac_c12_b9_p_2n_piminus_piplus = float(number_c12_b9_p_2n_piminus_piplus) / float(number_entries) * 100
    frac_c12_b9_2p_n_3piminus_2piplus = float(number_c12_b9_2p_n_3piminus_2piplus) / float(number_entries) * 100
    frac_c12_b9_2p_n_piminus = float(number_c12_b9_2p_n_piminus) / float(number_entries) * 100
    frac_c12_b9_3n_piplus = float(number_c12_b9_3n_piplus) / float(number_entries) * 100
    frac_c12_b9_p_2n_2piminus_2piplus = float(number_c12_b9_p_2n_2piminus_2piplus) / float(number_entries) * 100
    frac_c12_b9_2p_n_2piminus_piplus = float(number_c12_b9_2p_n_2piminus_piplus) / float(number_entries) * 100
    frac_c12_b9_other = float(number_c12_b9_other) / float(number_entries) * 100

    # Be9:
    frac_c12_be9_2p_n = float(number_c12_be9_2p_n) / float(number_entries) * 100
    frac_c12_be9_p_2n_piplus = float(number_c12_be9_p_2n_piplus) / float(number_entries) * 100
    frac_c12_be9_3p_piminus = float(number_c12_be9_3p_piminus) / float(number_entries) * 100
    frac_c12_be9_p_2n_piminus_2piplus = float(number_c12_be9_p_2n_piminus_2piplus) / float(number_entries) * 100
    frac_c12_be9_2p_n_piminus_piplus = float(number_c12_be9_2p_n_piminus_piplus) / float(number_entries) * 100
    frac_c12_be9_2p_n_3piminus_3piplus = float(number_c12_be9_2p_n_3piminus_3piplus) / float(number_entries) * 100
    frac_c12_be9_2p_n_2piminus_2piplus = float(number_c12_be9_2p_n_2piminus_2piplus) / float(number_entries) * 100
    frac_c12_be9_3n_2piplus = float(number_c12_be9_3n_2piplus) / float(number_entries) * 100
    frac_c12_be9_3p_2piminus_piplus = float(number_c12_be9_3p_2piminus_piplus) / float(number_entries) * 100
    frac_c12_be9_other = float(number_c12_be9_other) / float(number_entries) * 100

    # Be8:
    frac_c12_be8_2p_2n = float(number_c12_be8_2p_2n) / float(number_entries) * 100
    frac_c12_be8_3p_n_piminus = float(number_c12_be8_3p_n_piminus) / float(number_entries) * 100
    frac_c12_be8_p_3n_piplus = float(number_c12_be8_p_3n_piplus) / float(number_entries) * 100
    frac_c12_be8_2p_2n_2piminus_2piplus = float(number_c12_be8_2p_2n_2piminus_2piplus) / float(number_entries) * 100
    frac_c12_be8_4n_2piplus = float(number_c12_be8_4n_2piplus) / float(number_entries) * 100
    frac_c12_be8_2p_2n_piminus_piplus = float(number_c12_be8_2p_2n_piminus_piplus) / float(number_entries) * 100
    frac_c12_be8_3p_n_2piminus_piplus = float(number_c12_be8_3p_n_2piminus_piplus) / float(number_entries) * 100
    frac_c12_be8_4p_2piminus = float(number_c12_be8_4p_2piminus) / float(number_entries) * 100
    frac_c12_be8_other = float(number_c12_be8_other) / float(number_entries) * 100

    # C9:
    frac_c12_c9_p_2n_piminus = float(number_c12_c9_p_2n_piminus) / float(number_entries) * 100
    frac_c12_c9_3n = float(number_c12_c9_3n) / float(number_entries) * 100
    frac_c12_c9_2p_n_2piminus = float(number_c12_c9_2p_n_2piminus) / float(number_entries) * 100
    frac_c12_c9_3n_2piminus_2piplus = float(number_c12_c9_3n_2piminus_2piplus) / float(number_entries) * 100
    frac_c12_c9_other = float(number_c12_c9_other) / float(number_entries) * 100

    # Be7:
    frac_c12_be7_2p_3n = float(number_c12_be7_2p_3n) / float(number_entries) * 100
    frac_c12_be7_p_4n_piplus = float(number_c12_be7_p_4n_piplus) / float(number_entries) * 100
    frac_c12_be7_2p_3n_2piminus_2piplus = float(number_c12_be7_2p_3n_2piminus_2piplus) / float(number_entries) * 100
    frac_c12_be7_3p_2n_piminus = float(number_c12_be7_3p_2n_piminus) / float(number_entries) * 100
    frac_c12_be7_4p_n_2piminus = float(number_c12_be7_4p_n_2piminus) / float(number_entries) * 100
    frac_c12_be7_3p_2n_2piminus_piplus = float(number_c12_be7_3p_2n_2piminus_piplus) / float(number_entries) * 100
    frac_c12_be7_other = float(number_c12_be7_other) / float(number_entries) * 100

    # Li6:
    frac_c12_li6_3p_3n = float(number_c12_li6_3p_3n) / float(number_entries) * 100
    frac_c12_li6_2p_4n_piplus = float(number_c12_li6_2p_4n_piplus) / float(number_entries) * 100
    frac_c12_li6_5p_n_2piminus = float(number_c12_li6_5p_n_2piminus) / float(number_entries) * 100
    frac_c12_li6_2p_4n_piminus_2piplus = float(number_c12_li6_2p_4n_piminus_2piplus) / float(number_entries) * 100
    frac_c12_li6_4p_2n_piminus = float(number_c12_li6_4p_2n_piminus) / float(number_entries) * 100
    frac_c12_li6_3p_3n_piminus_piplus = float(number_c12_li6_3p_3n_piminus_piplus) / float(number_entries) * 100
    frac_c12_li6_other = float(number_c12_li6_other) / float(number_entries) * 100

    # Li8:
    frac_c12_li8_3p_n = float(number_c12_li8_3p_n) / float(number_entries) * 100
    frac_c12_li8_4p_piminus = float(number_c12_li8_4p_piminus) / float(number_entries) * 100
    frac_c12_li8_4p_2piminus_piplus = float(number_c12_li8_4p_2piminus_piplus) / float(number_entries) * 100
    frac_c12_li8_2p_2n_piplus = float(number_c12_li8_2p_2n_piplus) / float(number_entries) * 100
    frac_c12_li8_3p_n_piminus_piplus = float(number_c12_li8_3p_n_piminus_piplus) / float(number_entries) * 100
    frac_c12_li8_other = float(number_c12_li8_other) / float(number_entries) * 100

    # Li7:
    frac_c12_li7_2p_3n_piplus = float(number_c12_li7_2p_3n_piplus) / float(number_entries) * 100
    frac_c12_li7_4p_n_piminus = float(number_c12_li7_4p_n_piminus) / float(number_entries) * 100
    frac_c12_li7_3p_2n = float(number_c12_li7_3p_2n) / float(number_entries) * 100
    frac_c12_li7_3p_2n_piminus_piplus = float(number_c12_li7_3p_2n_piminus_piplus) / float(number_entries) * 100
    frac_c12_li7_4p_n_2piminus_piplus = float(number_c12_li7_4p_n_2piminus_piplus) / float(number_entries) * 100
    frac_c12_li7_2p_3n_piminus_2piplus = float(number_c12_li7_2p_3n_piminus_2piplus) / float(number_entries) * 100
    frac_c12_li7_other = float(number_c12_li7_other) / float(number_entries) * 100

    # B8:
    frac_c12_b8_p_3n = float(number_c12_b8_p_3n) / float(number_entries) * 100
    frac_c12_b8_p_3n_piminus_piplus = float(number_c12_b8_p_3n_piminus_piplus) / float(number_entries) * 100
    frac_c12_b8_2p_2n_2piminus_piplus = float(number_c12_b8_2p_2n_2piminus_piplus) / float(number_entries) * 100
    frac_c12_b8_2p_2n_piminus = float(number_c12_b8_2p_2n_piminus) / float(number_entries) * 100
    frac_c12_b8_4n_piplus = float(number_c12_b8_4n_piplus) / float(number_entries) * 100
    frac_c12_b8_other = float(number_c12_b8_other) / float(number_entries) * 100

    # Li9:
    frac_c12_li9_2p_n_piplus = float(number_c12_li9_2p_n_piplus) / float(number_entries) * 100
    frac_c12_li9_3p = float(number_c12_li9_3p) / float(number_entries) * 100
    frac_c12_li9_3p_piminus_piplus = float(number_c12_li9_3p_piminus_piplus) / float(number_entries) * 100
    frac_c12_li9_2p_n_piminus_2piplus = float(number_c12_li9_2p_n_piminus_2piplus) / float(number_entries) * 100
    frac_c12_li9_p_2n_piminus_3piplus = float(number_c12_li9_p_2n_piminus_3piplus) / float(number_entries) * 100
    frac_c12_li9_other = float(number_c12_li9_other) / float(number_entries) * 100


    frac_no_c12 = float(number_no_c12) / float(number_entries) * 100


    return frac_c12_b11_p, frac_c12_b11_n_piplus, frac_c12_b11_n_piminus_2piplus, frac_c12_b11_p_piminus_piplus, \
           frac_c12_b11_p_2piminus_2piplus, frac_c12_b11_piplus, frac_c12_b11_other, \
           frac_c12_c11_n, frac_c12_c11_p_piminus, frac_c12_c11_n_piminus_piplus, frac_c12_c11_p_2piminus_piplus, \
           frac_c12_c11_p_3piminus_2piplus, frac_c12_c11_n_2piminus_2piplus, frac_c12_c11_other, \
           frac_c12_b10_p_n, frac_c12_b10_2p_piminus, frac_c12_b10_p_n_piminus_piplus, frac_c12_b10_2n_piplus, \
           frac_c12_b10_2n_piminus_2piplus, frac_c12_b10_2p_2piminus_piplus, frac_c12_b10_2p_3piminus_2piplus, \
           frac_c12_b10_p_n_2piminus_2piplus, frac_c12_b10_other, \
           frac_c12_c10_2n, frac_c12_c10_p_n_piminus, frac_c12_c10_p_n_2piminus_piplus, frac_c12_c10_2n_piminus_piplus,\
           frac_c12_c10_2p_2piminus, frac_c12_c10_other, \
           frac_c12_be10_2p, frac_c12_be10_p_n_piplus, frac_c12_be10_p_n_piminus_2piplus, \
           frac_c12_be10_2p_piminus_piplus, frac_c12_be10_2n_2piplus, frac_c12_be10_p_n_2piminus_3piplus, \
           frac_c12_be10_2p_2piminus_2piplus, frac_c12_be10_2p_3piminus_3piplus, frac_c12_be10_other, \
           frac_c12_b9_p_2n, frac_c12_b9_p_2n_piminus_piplus, frac_c12_b9_2p_n_3piminus_2piplus, \
           frac_c12_b9_2p_n_piminus, frac_c12_b9_3n_piplus, frac_c12_b9_p_2n_2piminus_2piplus, \
           frac_c12_b9_2p_n_2piminus_piplus, frac_c12_b9_other, \
           frac_c12_be9_2p_n, frac_c12_be9_p_2n_piplus, frac_c12_be9_3p_piminus, frac_c12_be9_p_2n_piminus_2piplus, \
           frac_c12_be9_2p_n_piminus_piplus, frac_c12_be9_2p_n_3piminus_3piplus, frac_c12_be9_2p_n_2piminus_2piplus, \
           frac_c12_be9_3n_2piplus, frac_c12_be9_3p_2piminus_piplus, frac_c12_be9_other, \
           frac_c12_be8_2p_2n, frac_c12_be8_3p_n_piminus, frac_c12_be8_p_3n_piplus, \
           frac_c12_be8_2p_2n_2piminus_2piplus, frac_c12_be8_4n_2piplus, frac_c12_be8_2p_2n_piminus_piplus, \
           frac_c12_be8_3p_n_2piminus_piplus, frac_c12_be8_4p_2piminus, frac_c12_be8_other, \
           frac_c12_c9_p_2n_piminus, frac_c12_c9_3n, frac_c12_c9_2p_n_2piminus, frac_c12_c9_3n_2piminus_2piplus, \
           frac_c12_c9_other, \
           frac_c12_be7_2p_3n, frac_c12_be7_p_4n_piplus, frac_c12_be7_2p_3n_2piminus_2piplus, \
           frac_c12_be7_3p_2n_piminus, frac_c12_be7_4p_n_2piminus, frac_c12_be7_3p_2n_2piminus_piplus, \
           frac_c12_be7_other, \
           frac_c12_li6_3p_3n, frac_c12_li6_2p_4n_piplus, frac_c12_li6_5p_n_2piminus, \
           frac_c12_li6_2p_4n_piminus_2piplus, frac_c12_li6_4p_2n_piminus, frac_c12_li6_3p_3n_piminus_piplus, \
           frac_c12_li6_other, \
           frac_c12_li8_3p_n, frac_c12_li8_4p_piminus, frac_c12_li8_4p_2piminus_piplus, frac_c12_li8_2p_2n_piplus, \
           frac_c12_li8_3p_n_piminus_piplus, frac_c12_li8_other, \
           frac_c12_li7_2p_3n_piplus, frac_c12_li7_4p_n_piminus, frac_c12_li7_3p_2n, frac_c12_li7_3p_2n_piminus_piplus,\
           frac_c12_li7_4p_n_2piminus_piplus, frac_c12_li7_2p_3n_piminus_2piplus, frac_c12_li7_other, \
           frac_c12_b8_p_3n, frac_c12_b8_p_3n_piminus_piplus, frac_c12_b8_2p_2n_2piminus_piplus, \
           frac_c12_b8_2p_2n_piminus, frac_c12_b8_4n_piplus, frac_c12_b8_other, \
           frac_c12_li9_2p_n_piplus, frac_c12_li9_3p, frac_c12_li9_3p_piminus_piplus, \
           frac_c12_li9_2p_n_piminus_2piplus, frac_c12_li9_p_2n_piminus_3piplus, frac_c12_li9_other, \
           number_no_c12, frac_no_c12
