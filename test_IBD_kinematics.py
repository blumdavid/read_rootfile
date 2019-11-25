""" script to check the kinematics of IBD based on paper 'angular distribution of neutron inverse beta decay' by
    Vogel and Beacom (this kinematics is also used in IBD generator of JUNO offline):

    With this script you can also produce hep-evt files for one neutrino energy, where the positron momenta and the
    neutron momenta is stored for different angles of theta and phi.

    format of hepevt file (momentum and mass in GeV):

    number particles in event
    1   PDGID 0 0 px py pz mass
    1   PDGID 0 0 px py pz mass

    for example:
    2
    1   -11 0 0 px py pz 0.000511
    1   2112 0 0 px py pz 0.939570

"""

import datetime
import sys
from NC_background_functions import energy_resolution
import numpy as np
from matplotlib import pyplot as plt


def dsigma_dcostheta(cosinetheta, energy_pos_0, energy_pos_1, velocity_pos_1, momentum_pos_0, momentum_pos_1, sig_0,
                     gamma_fac, f_factor, g_factor, mass_proton):
    """
    differential cross-section at first order of 1/M (equ. 12)

    :param cosinetheta: cosine of angle theta between neutrino and positron direction
    :param energy_pos_0: positron energy at zeroth order of 1/M in MeV (equ. 6)
    :param energy_pos_1: positron energy at first order of 1/M in MeV (equ. 11)
    :param velocity_pos_1: positron velocity at first order of 1/M in MeV (velocity = momentum/energy)
    :param momentum_pos_0: positron momentum at zeroth order of 1/M in MeV (momentum = sqrt(energy^2 - m_pos^2))
    :param momentum_pos_1: positron momentum at first order of 1/M in MeV (momentum = sqrt(energy^2 - m_pos^2))
    :param sig_0: normalizing constant in 1/MeV^2, including the energy-independent inner radiative corrections (equ. 8)
    :param gamma_fac: gamma factor (equ. 13)
    :param f_factor: vector coupling constant
    :param g_factor: axial vector coupling constant
    :param mass_proton: proton mass in MeV
    :return:
    """
    # calculate first line of equation 12:
    first_line = (sig_0 / 2.0 * ((f_factor**2 + 3 * g_factor**2) + (f_factor**2 - g_factor**2) * velocity_pos_1
                                 * cosinetheta) * energy_pos_1 * momentum_pos_1)
    # calculate second line of equation 12:
    second_line = sig_0 / 2.0 * gamma_fac / mass_proton * energy_pos_0 * momentum_pos_0

    # calculate differential cross-section (equ. 12):
    cross_section = first_line - second_line

    return cross_section


def gamma_function(cosinetheta, energy_pos_0, velocity_pos_0, f_factor, f2_factor, g_factor, delta, mass_positron):
    """
    gamma function (equ. 13)

    :param cosinetheta: cosine of angle theta between neutrino and positron direction
    :param energy_pos_0: positron energy at zeroth order of 1/M in MeV (equ. 6)
    :param velocity_pos_0: positron velocity at zeroth order of 1/M in MeV (velocity = momentum/energy)
    :param f_factor: vector coupling constant
    :param f2_factor: anomalous nucleon iso-vector magnetic moment
    :param g_factor: axial vector coupling constant
    :param delta: mass_neutron - mass_proton in MeV
    :param mass_positron: positron mass in MeV
    :return:
    """
    first_line = (2 * (f_factor + f2_factor) * g_factor *
                  ((2 * energy_pos_0 + delta) * (1 - velocity_pos_0 * cosinetheta) - mass_positron ** 2 / energy_pos_0))

    second_line = (f_factor ** 2 + g_factor ** 2) * (delta * (1 + velocity_pos_0 * cosinetheta) + mass_positron ** 2 /
                                                     energy_pos_0)

    third_line = (f_factor ** 2 + 3 * g_factor ** 2) * ((energy_pos_0 + delta) * (1 - cosinetheta / velocity_pos_0)
                                                        - delta)

    fourth_line = (f_factor ** 2 - g_factor ** 2) * ((energy_pos_0 + delta) * (1 - cosinetheta / velocity_pos_0) -
                                                     delta) * velocity_pos_0 * cosinetheta

    gamma_fac = first_line + second_line + third_line + fourth_line

    return gamma_fac


def energy_positron_1(cosinetheta, energy_pos_0, energy_neutrino, velocity_pos_0, mass_proton, delta, mass_positron):
    """
    positron energy depending on scattering angle theta at first order in 1/M (equ. 11)

    :param cosinetheta: cosine of angle theta between neutrino and positron direction
    :param energy_pos_0: positron energy at zeroth order of 1/M in MeV (equ. 6)
    :param energy_neutrino: neutrino energy in MeV
    :param velocity_pos_0: positron velocity at zeroth order of 1/M in MeV (velocity = momentum/energy)
    :param mass_proton: proton mass in MeV
    :param delta: mass_neutron - mass_proton in MeV
    :param mass_positron: positron mass in MeV
    :return:
    """
    # calculate y_square in MeV^2:
    y_square = (delta**2 - mass_positron**2) / 2

    # calculate positron energy in MeV at first order of 1/M (equ. 11):
    energy_pos_1 = (energy_pos_0 * (1 - energy_neutrino / mass_proton * (1 - velocity_pos_0 * cosinetheta))
                    - y_square / mass_proton)

    return energy_pos_1


# Flag, if hepevt files are created:
CREATE_HEPFILES = False

# define neutrino energy in MeV:
E_nu = 60.0

# set the number of events, that should be stored in one hepevt file:
number_events_per_file = 100
# set the number of files that should be created:
number_of_files = 100
# set the path, where the hepevt files should be stored:
output_path = "/home/astro/blum/PhD/work/MeVDM_JUNO/gen_spectrum_v2/DM_signal_detsim/hepevt_files/"

# set angle in degree:
theta_min = 0.0
theta_max = 180.0
theta_interval = 0.5
# set axis for theta:
theta_axis = np.arange(theta_min, theta_max+theta_interval, theta_interval)
# get number of entries in theta_axis:
number_events = len(theta_axis)

# PDG ID of positron:
pdg_positron = -11
# PDG ID of neutron:
pdg_neutron = 2112
# mass of positron in MeV (float constant):
MASS_POSITRON = 0.51099892
# mass of proton in MeV (float constant):
MASS_PROTON = 938.27203
# mass of neutron in MeV (float constant):
MASS_NEUTRON = 939.56536
# difference MASS_NEUTRON - MASS_PROTON in MeV (float):
DELTA = MASS_NEUTRON - MASS_PROTON
# fermi constant in 1/MeV^2:
G_f = 1.16637 * 10**(-11)
# cosine of theta_c:
cos_theta_c = 0.974
# delta_R_inner:
Delta_R_inner = 0.024
# vector coupling constant:
f = 1.0
# axial vector coupling constant:
g = 1.26
# f2_factor: anomalous nucleon iso-vector magnetic moment:
f_2 = 3.706

# calculate positron energy in MeV at zeroth order of 1/M (equ. 6):
E_pos_0 = E_nu - DELTA
# calculate positron momentum in MeV:
p_pos_0 = np.sqrt(E_pos_0**2 - MASS_POSITRON**2)
# calculate positron velocity in MeV:
v_pos_0 = p_pos_0 / E_pos_0
# calculate sigma_0 (equ. 8):
sigma_0 = G_f**2 * cos_theta_c**2 / np.pi * (1 + Delta_R_inner)

# preallocate array, where differential cross-section is stored in cm^2:
array_DsigDcosTh = []

# generate neutrino events:
for theta in theta_axis:

    # calculate cosine of theta
    cosTheta = np.cos(np.deg2rad(theta))

    """ calculate positron energy (equ. 11): """
    E_pos_1 = energy_positron_1(cosTheta, E_pos_0, E_nu, v_pos_0, MASS_PROTON, DELTA, MASS_POSITRON)

    """ calculate IBD cross-section (equ. 12): """
    # positron momentum O(1/M) in MeV:
    p_pos_1 = np.sqrt(E_pos_1**2 - MASS_POSITRON**2)
    # positron velocity O(1/M) in MeV:
    v_pos_1 = p_pos_1 / E_pos_1

    # calculate gamma factor (equ. 13):
    gamma_factor = gamma_function(cosTheta, E_pos_0, v_pos_0, f, f_2, g, DELTA, MASS_POSITRON)

    # differential cross-section at first order of 1/M (equ.12):
    dSigmaDcosTheta = dsigma_dcostheta(cosTheta, E_pos_0, E_pos_1, v_pos_1, p_pos_0, p_pos_1,
                                       sigma_0, gamma_factor, f, g, MASS_PROTON)

    if dSigmaDcosTheta < 0.0:
        dSigmaDcosTheta = 0.0

    # append dSigmaDcosTheta to array:
    array_DsigDcosTh.append(dSigmaDcosTheta)

# normalize array_DsigDcosTh to 1 to get probability function:
array_DsigDcosTh = array_DsigDcosTh[:-1] / np.sum(array_DsigDcosTh[:-1])
# set the number of theta values, that should be generated:
number = 100000
# generate random values of theta with the differential cross-section as probability function in degree (array):
theta_pos_random = np.random.choice(theta_axis[:-1], p=array_DsigDcosTh, size=number)

# calculate cosine of these randomly generated theta values (array):
cosTheta_pos_random = np.cos(np.deg2rad(theta_pos_random))

# calculate positron energy in MeV with this randomly generated theta values (array):
E_pos_1_random = energy_positron_1(cosTheta_pos_random, E_pos_0, E_nu, v_pos_0, MASS_PROTON, DELTA, MASS_POSITRON)

# calculate positron momentum in MeV with this randomly generated theta values (array):
p_pos_1_random = np.sqrt(E_pos_1_random**2 - MASS_POSITRON**2)

# generate angle phi between neutrino and positron with uniformly distributed random number between 0 und 360 degree:
phi_pos_random = np.random.uniform(0.0, 360.0, size=number)

# calculate neutron momentum in MeV (array):
p_neutron_random = E_nu - p_pos_1_random
# set angle theta of the neutron in degree (neutrino is assumed as vector (E_nu, 0, 0) in spherical coordinates).
# Therefore theta_neutron = -theta_positron (proton is at rest) (array):
theta_neutron_random = - theta_pos_random
# set angle phi of the neutron in degree (array):
phi_neutron_random = - phi_pos_random

# convert total momentum p_pos_1_random to px, py, pz in MeV (array):
px_pos_random = p_pos_1_random * np.sin(np.deg2rad(theta_pos_random)) * np.cos(np.deg2rad(phi_pos_random))
py_pos_random = p_pos_1_random * np.sin(np.deg2rad(theta_pos_random)) * np.sin(np.deg2rad(phi_pos_random))
pz_pos_random = p_pos_1_random * np.cos(np.deg2rad(theta_pos_random))

# convert total momentum p_neutron_random to px, py, pz in MeV (array):
px_neutron_random = p_neutron_random * np.sin(np.deg2rad(theta_neutron_random)) * np.cos(np.deg2rad(phi_neutron_random))
py_neutron_random = p_neutron_random * np.sin(np.deg2rad(theta_neutron_random)) * np.sin(np.deg2rad(phi_neutron_random))
pz_neutron_random = p_neutron_random * np.cos(np.deg2rad(theta_neutron_random))

# momenta must be in GeV for hepevt file input:
px_pos_random = px_pos_random / 1000.0
py_pos_random = py_pos_random / 1000.0
pz_pos_random = pz_pos_random / 1000.0
px_neutron_random = px_neutron_random / 1000.0
py_neutron_random = py_neutron_random / 1000.0
pz_neutron_random = pz_neutron_random / 1000.0

# masses must be also in GeV for hepevt file input:
mass_pos_GeV = MASS_POSITRON / 1000.0
mass_neutron_GeV = MASS_NEUTRON / 1000.0

""" create hepevt files, where p_pos_random and p_neutron_random are stored: """
if CREATE_HEPFILES:
    # loop over the number of files:
    for filenumber in range(number_of_files):

        # preallocate hepevt array:
        hepevt_file = []

        # loop over the number of events per file:
        for event in range(number_events_per_file):

            # calculate index corresponding to arrays px_pos_random:
            index = filenumber*number_events_per_file + event

            """ append information to hepevt_file: """
            # 2 particles in the event:
            hepevt_file.append("2")
            # append positron information:
            hepevt_file.append("1\t{0} 0 0 {1} {2} {3} {4}".format(pdg_positron, px_pos_random[index], py_pos_random[index],
                                                                   pz_pos_random[index], mass_pos_GeV))
            # append neutron information to file:
            hepevt_file.append("1\t{0} 0 0 {1} {2} {3} {4}".format(pdg_neutron, px_neutron_random[index],
                                                                   py_neutron_random[index],
                                                                   pz_neutron_random[index], mass_neutron_GeV))

        # open file:
        MyFile = open(output_path + "DM{2:.0f}MeV_hepevt_{1:d}events_file{0:d}.txt"
                      .format(filenumber, number_events_per_file, E_nu), 'w')

        for element in hepevt_file:
            print >>MyFile, element
        MyFile.close()

# print("E_pos min")
# print(min(E_pos_1_random))
# print("E_pos max")
# print(max(E_pos_1_random))
# print("E_pos mean")
# print(np.mean(E_pos_1_random))

# calculate visible energy in MeV (array):
E_visible = E_pos_1_random + MASS_POSITRON

# get sigma from energy resolution in MeV (array):
sigma = energy_resolution(E_visible)

# smear E_visible with sigma (array):
E_visible_smeared = np.random.normal(E_visible, sigma)

# build histogram of array_E_pos:
h1 = plt.figure(1, figsize=(15, 8))
First_bin = 10.0
Last_bin = 100.0
Bin_width = 0.5
Bins = np.arange(First_bin, Last_bin+Bin_width, Bin_width)
E_visible_hist, bin_edges, patches = plt.hist(E_visible, bins=Bins, histtype="step", align='mid', color="r",
                                              linewidth=1.5,
                                              label="w/o energy resolution (entries = {0:d})".format(len(E_visible)))
E_visible_hist_smeared, bin_edges, patches = plt.hist(E_visible_smeared, bins=Bins, histtype="step", align='mid',
                                                      color="b", linewidth=1.5, label="w/ energy resolution "
                                                                                      "(entries = {0:d})"
                                                      .format(len(E_visible_smeared)))
plt.xlabel("visible energy of positron in MeV")
plt.ylabel("events")
plt.title("DM mass = {0:.0f} MeV".format(E_nu))
plt.legend()
plt.grid()
# plt.savefig()
# plt.close()

h2 = plt.figure(2, figsize=(15, 8))
hist_theta_random, bin_edges = np.histogram(theta_pos_random, bins=theta_axis)
hist_theta_random = hist_theta_random / float(np.sum(hist_theta_random))
plt.step(theta_axis[:-1], hist_theta_random, 'r', label='random numbers')
plt.plot(theta_axis[:-1], array_DsigDcosTh, label="entries = {0:d}".format(len(array_DsigDcosTh)))
plt.xlabel("theta in degree")
plt.ylabel("dSigma / dcosTheta")
# plt.title("DM mass = {0:.0f} MeV".format(E_nu))
plt.legend()
plt.grid()
# plt.savefig()
# plt.close()

h3 = plt.figure(3, figsize=(15, 8))
Bins2 = np.arange(0.0, Last_bin+Bin_width, Bin_width)
plt.hist(p_neutron_random, bins=Bins2, histtype="step", align='mid', color="r", linewidth=1.5,
         label="entries = {0:d}".format(len(p_neutron_random)))
plt.xlabel("kinetic energy of neutron in MeV")
plt.ylabel("events")
plt.title("DM mass = {0:.0f} MeV".format(E_nu))
plt.legend()
plt.grid()
# plt.savefig()
# plt.close()

# compare E_visible from calculation with the result of check_DMsignal_simulation.py of simulation:
# load result of simulation:
Qedep_simu = np.loadtxt("/home/astro/blum/PhD/work/MeVDM_JUNO/gen_spectrum_v2/DM_signal_detsim/"
                        "DMsignal_{0:.0f}MeV_R0mTo17m_Qedep.txt".format(E_nu))
# calculate number of analyzed events:
number_events_simu = np.sum(Qedep_simu)

# preallocate array where smeared Qedep are stored:
Qedep_simu_smeared = np.loadtxt("/home/astro/blum/PhD/work/MeVDM_JUNO/gen_spectrum_v2/DM_signal_detsim/"
                                "DMsignal_{0:.0f}MeV_R0mTo17m_Qedep_smeared.txt".format(E_nu))

# normalize spectra to 1:
E_vis_norm = E_visible_hist / np.sum(E_visible_hist)
E_vis_smeared_norm = E_visible_hist_smeared / np.sum(E_visible_hist_smeared)
Qedep_simu_norm = Qedep_simu / number_events_simu
Qedep_simu_smeared_norm = Qedep_simu_smeared / number_events_simu

# # shift Qedep_simu_norm in a way, that the maximum values of E_vis_norm and Qedep_simu_norm are equal
# # (correct error of conversion from nPE to MeV):
# # get index corresponding to last entry in E_vis_norm that differs from zero:
# index_max_E_vis_norm = np.max(np.nonzero(E_vis_norm))
# # get index corresponding to last entry in Qedep_simu_norm that differs from zero:
# index_max_E_vis_simu_norm = np.max(np.nonzero(Qedep_simu_norm))
# # calculate difference of indices:
# index_diff = index_max_E_vis_simu_norm - index_max_E_vis_norm
# # delete first index_diff entries of Qedep_simu_norm:
# Qedep_simu_norm = np.delete(Qedep_simu_norm, np.arange(index_diff))
# # append index_diff values of zero to Qedep_simu_norm:
# Qedep_simu_norm = np.append(Qedep_simu_norm, np.zeros(index_diff))

h4 = plt.figure(4, figsize=(15, 8))
plt.step(Bins[:-1], E_vis_norm, color="r", linewidth=1.5, where='post',
         label="E_vis of positron from calculation (entries = {0:.0f})".format(np.sum(E_visible_hist)))
plt.step(Bins[:-1], Qedep_simu_norm, color="b", linewidth=1.5, where='post',
         label="Qedep of event from simulation (entries = {0:.0f})".format(number_events_simu))
plt.xlabel("energy in MeV")
plt.ylabel("events")
plt.title("Normalized spectra for DM mass = {0:.0f} MeV".format(E_nu))
plt.legend()
plt.grid()
# plt.savefig()
# plt.close()

h5 = plt.figure(5, figsize=(15, 8))
plt.step(Bins[:-1], E_vis_smeared_norm, color="r", linewidth=1.5, where='post',
         label="smeared E_vis of positron from calculation (entries = {0:.0f})".format(np.sum(E_visible_hist_smeared)))
plt.step(Bins[:-1], Qedep_simu_smeared_norm, color="b", linewidth=1.5, where='post',
         label="smeared Qedep of event from simulation (entries = {0:.0f})".format(number_events_simu))
plt.xlabel("energy in MeV")
plt.ylabel("events")
plt.title("Normalized spectra for DM mass = {0:.0f} MeV".format(E_nu))
plt.legend()
plt.grid()
# plt.savefig()
# plt.close()

plt.show()
















