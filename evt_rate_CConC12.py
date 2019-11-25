""" script to calculate the theoretical event rate in JUNO for the CC interaction of atmospheric electron-antineutrino
    on C12 into positron and excited state of B12.

    This B12 de-excites then into B11 and a neutron.

    Event rate for the channel: nu_e_bar + C12 -> positron + B11 + neutron

    Equation to calculate event rate:

    dN/dt = A * flux * cross-section * 4*pi * P(nu_e_bar + C12 -> positron + B11 + neutron)

    script is based on function event_rate() of script NC_background_functions.py

"""
from NC_background_functions import number_c12_atoms, read_xml_xsec, convolution_new
import numpy as np
from matplotlib import pyplot as plt


# radius cut of fiducial volume in mm:
radius_cut = 17700.0

# minimum neutrino energy in MeV:
E_nu_min = 10.0
# maximum neutrino energy in MeV:
E_nu_max = 100.0
# energy interval in MeV:
interval_energy = 0.1
# define the energy range in MeV, that should be integrated (array of float):
energy_neutrino = np.arange(E_nu_min, E_nu_max + interval_energy, interval_energy)
# mass of proton in MeV:
mass_proton = 938.27203
# mass of neutron in MeV:
mass_neutron = 939.56536
# define delta = mass_neutron - mass_proton in MeV:
delta = mass_neutron - mass_proton
# define mass of positron in MeV:
mass_positron = 0.51099892
# energy interval of positron energy in MeV:
interval_energy_1 = 0.5
# define positron energy in MeV:
energy = np.arange(E_nu_min, E_nu_max+interval_energy_1, interval_energy_1)

""" Results of the FLUKA simulation (from the paper of Battistoni2005 'The atmospheric neutrino fluxes below 
100 MeV: The FLUKA results'): """
# Neutrino energy in MeV from table 3 from paper 1-s2.0-S0927650505000526-main (np.array of float):
energy_fluka = np.array([0, 13, 15, 17, 19, 21, 24, 27, 30, 33, 38, 42, 47, 53, 60, 67, 75, 84, 94, 106, 119, 133,
                         150, 168, 188, 211, 237, 266, 299, 335, 376, 422, 473, 531, 596, 668, 750, 841, 944])

# differential flux from FLUKA in energy for no oscillation for electron-antineutrinos for solar average at the site
# of Super-Kamiokande, in (MeV**(-1) * cm**(-2) * s**(-1)) (np.array of float).
# Assumption: for energy = 0 MeV, the flux is also 0!
flux_nuebar_fluka = 10 ** (-4) * np.array([0, 63.7, 69.7, 79.5, 84.2, 89.4, 95.0, 99.3, 103., 104., 101., 96.1,
                                           83.5, 65.9, 60.0, 56.4, 51.4, 46.3, 43.0, 37.2, 32.9, 28.8, 24.9, 21.3,
                                           18.3, 15.4, 12.9, 10.6, 8.80, 7.13, 5.75, 4.60, 3.68, 2.88, 2.28,
                                           1.87, 1.37, 1.06, 0.800])

# differential flux from FLUKA in energy for no oscillation for electron-neutrinos for solar average at the site
# of Super-Kamiokande, in (MeV**(-1) * cm**(-2) * s**(-1)) (np.array of float).
# Assumption: for energy = 0 MeV, the flux is also 0!
flux_nue_fluka = 10 ** (-4) * np.array([0, 69.6, 74.6, 79.7, 87.4, 94.2, 101., 103., 109., 108., 107., 101., 88.5,
                                        69.6, 64.4, 59.3, 54.3, 49.7, 45.1, 40.6, 35.8, 31.7, 27.3, 23.9, 20.4,
                                        17.0, 14.5, 12.0, 9.96, 8.11, 6.62, 5.27, 4.23, 3.37, 2.66, 2.09, 1.62,
                                        1.24, 0.950])

# differential flux from FLUKA in energy for no oscillation for muon-neutrinos for solar average at the site
# of Super-Kamiokande, in (MeV**(-1) * cm**(-2) * s**(-1)) (np.array of float).
# Assumption: for energy = 0 MeV, the flux is also 0!
flux_numu_fluka = 10 ** (-4) * np.array([0, 114., 124., 138., 146., 155., 159., 164., 181., 174., 179., 178., 176.,
                                         153., 131., 123., 114., 107., 96.3, 84.2, 72.7, 63.5, 55.2, 47.7, 41.2,
                                         34.4, 28.4, 23.6, 19.6, 15.8, 12.8, 10.3, 8.20, 6.49, 5.15, 3.98, 3.13,
                                         2.41, 1.82])

# differential flux from FLUKA in energy for no oscillation for muon-antineutrinos for solar average at the site of
# Super-K, in (MeV**(-1) * cm**(-2) * s**(-1)) (np.array of float).
# Assumption: for energy = 0 MeV, the flux is also 0!
flux_numubar_fluka = 10 ** (-4) * np.array([0, 116., 128., 136., 150., 158., 162., 170., 196., 177., 182., 183.,
                                            181., 155., 132., 123., 112., 101., 92.1, 82.2, 72.5, 64.0, 55.6,
                                            47.6, 40.8, 34.1, 28.6, 23.5, 19.3, 15.7, 12.6, 10.2, 8.15, 6.48,
                                            5.02, 3.94, 3.03, 2.33, 1.79])


""" Results of the HONDA simulation (based on the paper of Honda2015: 'Atmospheric neutrino flux calculation using
the NRLMSISE-00 atmospheric model'), (for solar maximum (HONDA_juno-ally-01-01-solmax.d)): """
# Neutrino energy in MeV from the table from file HONDA_juno-ally-01-01-solmax.d (np.array of float):
energy_honda = 10 ** 3 * np.array([1.0000E-01, 1.1220E-01, 1.2589E-01, 1.4125E-01, 1.5849E-01, 1.7783E-01,
                                   1.9953E-01, 2.2387E-01, 2.5119E-01, 2.8184E-01, 3.1623E-01, 3.5481E-01,
                                   3.9811E-01, 4.4668E-01, 5.0119E-01, 5.6234E-01, 6.3096E-01, 7.0795E-01,
                                   7.9433E-01, 8.9125E-01, 1.0000E+00, 1.1220E+00, 1.2589E+00, 1.4125E+00,
                                   1.5849E+00, 1.7783E+00, 1.9953E+00, 2.2387E+00, 2.5119E+00, 2.8184E+00,
                                   3.1623E+00, 3.5481E+00, 3.9811E+00, 4.4668E+00, 5.0119E+00, 5.6234E+00,
                                   6.3096E+00, 7.0795E+00, 7.9433E+00, 8.9125E+00, 1.0000E+01])

# all-direction averaged flux for no oscillation for electron-antineutrinos for solar maximum at the site of JUNO
# (WITHOUT mountain over the detector), in (MeV**(-1) * cm**(-2) * s**(-1) * sr**(-1)) (np.array of float):
flux_nuebar_max_honda = 10 ** (-7) * np.array([2.7733E+03, 2.4332E+03, 2.1124E+03, 1.8187E+03, 1.5545E+03,
                                               1.3190E+03, 1.1105E+03, 9.2820E+02, 7.7040E+02, 6.3403E+02,
                                               5.1790E+02, 4.1997E+02, 3.3811E+02, 2.7054E+02, 2.1539E+02,
                                               1.7049E+02, 1.3418E+02, 1.0499E+02, 8.1651E+01, 6.3166E+01,
                                               4.8654E+01, 3.7230E+01, 2.8329E+01, 2.1428E+01, 1.6121E+01,
                                               1.2064E+01, 8.9697E+00, 6.6258E+00, 4.8598E+00, 3.5435E+00,
                                               2.5685E+00, 1.8478E+00, 1.3252E+00, 9.4491E-01, 6.6836E-01,
                                               4.7226E-01, 3.3159E-01, 2.3192E-01, 1.6107E-01, 1.1131E-01,
                                               7.7823E-02])

# all-direction averaged flux for no oscillation for electron-neutrinos for solar maximum at the site of JUNO
# (WITHOUT mountain over the detector), in (MeV**(-1) * cm**(-2) * s**(-1) * sr**(-1)) (np.array of float):
flux_nue_max_honda = 10 ** (-7) * np.array([2.7743E+03, 2.4562E+03, 2.1530E+03, 1.8723E+03, 1.6181E+03, 1.3882E+03,
                                            1.1819E+03, 9.9837E+02, 8.3703E+02, 6.9614E+02, 5.7337E+02, 4.6903E+02,
                                            3.8128E+02, 3.0772E+02, 2.4680E+02, 1.9673E+02, 1.5600E+02, 1.2295E+02,
                                            9.6275E+01, 7.4975E+01, 5.8069E+01, 4.4733E+01, 3.4262E+01, 2.6103E+01,
                                            1.9770E+01, 1.4881E+01, 1.1137E+01, 8.2775E+00, 6.1088E+00, 4.4822E+00,
                                            3.2629E+00, 2.3653E+00, 1.7104E+00, 1.2266E+00, 8.7045E-01, 6.1557E-01,
                                            4.3368E-01, 3.0448E-01, 2.1286E-01, 1.4843E-01, 1.0281E-01])

# all-direction averaged flux for no oscillation for muon-neutrinos for solar maximum at the site of JUNO
# (WITHOUT mountain over the detector), in (MeV**(-1) * cm**(-2) * s**(-1) * sr^(-1)) (np.array of float):
flux_numu_max_honda = 10 ** (-7) * np.array([5.7913E+03, 5.0884E+03, 4.4520E+03, 3.8714E+03, 3.3388E+03, 2.8520E+03,
                                             2.4128E+03, 2.0226E+03, 1.6807E+03, 1.3858E+03, 1.1351E+03, 9.2472E+02,
                                             7.4912E+02, 6.0324E+02, 4.8323E+02, 3.8514E+02, 3.0543E+02, 2.4122E+02,
                                             1.8959E+02, 1.4845E+02, 1.1579E+02, 8.9940E+01, 6.9618E+01, 5.3647E+01,
                                             4.1114E+01, 3.1343E+01, 2.3751E+01, 1.7914E+01, 1.3453E+01, 1.0049E+01,
                                             7.4735E+00, 5.5296E+00, 4.0719E+00, 2.9889E+00, 2.1817E+00, 1.5909E+00,
                                             1.1558E+00, 8.3657E-01, 6.0575E-01, 4.3508E-01, 3.1237E-01])

# all-direction averaged flux for no oscillation for muon-antineutrinos for solar maximum at the site of JUNO
# (WITHOUT mountain over the detector), in (MeV**(-1) * cm**(-2) * s**(-1) * sr^(-1)) (np.array of float):
flux_numubar_max_honda = 10 ** (-7) * np.array([5.8966E+03, 5.1676E+03, 4.5104E+03, 3.9127E+03, 3.3665E+03,
                                                2.8701E+03, 2.4238E+03, 2.0277E+03, 1.6821E+03, 1.3857E+03,
                                                1.1333E+03, 9.2144E+02, 7.4476E+02, 5.9875E+02, 4.7865E+02,
                                                3.8024E+02, 3.0060E+02, 2.3645E+02, 1.8519E+02, 1.4444E+02,
                                                1.1204E+02, 8.6529E+01, 6.6529E+01, 5.0910E+01, 3.8731E+01,
                                                2.9299E+01, 2.2048E+01, 1.6504E+01, 1.2291E+01, 9.1084E+00,
                                                6.7266E+00, 4.9403E+00, 3.6136E+00, 2.6356E+00, 1.9115E+00,
                                                1.3828E+00, 9.9752E-01, 7.1482E-01, 5.1189E-01, 3.6743E-01,
                                                2.6256E-01])

# to compensate unit 1/sr from HONDA fluxes, multiply fluxes with 4*pi:
flux_nuebar_max_honda = flux_nuebar_max_honda * 4*np.pi
flux_nue_max_honda = flux_nue_max_honda * 4*np.pi
flux_numu_max_honda = flux_numu_max_honda * 4*np.pi
flux_numubar_max_honda = flux_numubar_max_honda * 4*np.pi

# Extrapolate the HONDA flux to the energies of the FLUKA simulation from 0 MeV to 100 MeV
# for electron-antineutrinos:
# Assumption:
# 1. the shape of the FLUKA flux as function of energy do NOT depend on the location
#     -> the shape of the flux at Super-K can also be used at JUNO site
#
# 2. the absolute value of the FLUKA flux at Super-K should be normalized to the location of JUNO
#     ->  therefore get the normalization factor by comparing the HONDA flux and the FLUKA flux in the energy
#         range from 100 MeV to 10 GeV
# define the energy-array, in which the normalization will be calculated (neutrino energy in MeV)
# (np.array of float):
energy_norm = np.arange(min(energy_honda), max(energy_fluka) + interval_energy, interval_energy)

# Interpolate the flux of FLUKA to get the differential flux in the energy range from 100 MeV to 950 MeV,
# in 1/(MeV * cm**2 * s) (np.array of float):
flux_nuebar_fluka_interpolated = np.interp(energy_norm, energy_fluka, flux_nuebar_fluka)
flux_nue_fluka_interpolated = np.interp(energy_norm, energy_fluka, flux_nue_fluka)
flux_numubar_fluka_interpolated = np.interp(energy_norm, energy_fluka, flux_numubar_fluka)
flux_numu_fluka_interpolated = np.interp(energy_norm, energy_fluka, flux_numu_fluka)

# Interpolate the flux of HONDA to get the differential flux in the energy range from 100 MeV to 950 MeV,
# in 1/(MeV * cm**2 * s) (np.array of float):
flux_nuebar_honda_interpolated = np.interp(energy_norm, energy_honda, flux_nuebar_max_honda)
flux_nue_honda_interpolated = np.interp(energy_norm, energy_honda, flux_nue_max_honda)
flux_numubar_honda_interpolated = np.interp(energy_norm, energy_honda, flux_numubar_max_honda)
flux_numu_honda_interpolated = np.interp(energy_norm, energy_honda, flux_numu_max_honda)

# Calculate the integral of the FLUKA flux in the energy range given by energy_norm (float):
integral_nuebar_fluka = np.trapz(flux_nuebar_fluka_interpolated, energy_norm)
integral_nue_fluka = np.trapz(flux_nue_fluka_interpolated, energy_norm)
integral_numubar_fluka = np.trapz(flux_numubar_fluka_interpolated, energy_norm)
integral_numu_fluka = np.trapz(flux_numu_fluka_interpolated, energy_norm)

# Calculate the integral of the HONDA flux in the energy range given by energy_norm (float):
integral_nuebar_honda = np.trapz(flux_nuebar_honda_interpolated, energy_norm)
integral_nue_honda = np.trapz(flux_nue_honda_interpolated, energy_norm)
integral_numubar_honda = np.trapz(flux_numubar_honda_interpolated, energy_norm)
integral_numu_honda = np.trapz(flux_numu_honda_interpolated, energy_norm)

# Interpolate the part of the FLUKA flux in the energy range from 0 MeV to 99.9 MeV, in 1/(MeV*s*cm**2)
# (np.array of float):
flux_nuebar_fluka_interesting = np.interp(np.arange(E_nu_min, E_nu_max+interval_energy, interval_energy),
                                          energy_fluka, flux_nuebar_fluka)
flux_nue_fluka_interesting = np.interp(np.arange(E_nu_min, E_nu_max+interval_energy, interval_energy),
                                       energy_fluka, flux_nue_fluka)
flux_numubar_fluka_interesting = np.interp(np.arange(E_nu_min, E_nu_max+interval_energy, interval_energy),
                                           energy_fluka, flux_numubar_fluka)
flux_numu_fluka_interesting = np.interp(np.arange(E_nu_min, E_nu_max+interval_energy, interval_energy),
                                        energy_fluka, flux_numu_fluka)

# Normalize flux_nuebar_fluka_interesting at Super-K to the electron-antineutrino flux at JUNO,
# in 1/(MeV * s * cm**2) (np.array of float):
flux_nuebar_fluka_norm = flux_nuebar_fluka_interesting * integral_nuebar_honda / integral_nuebar_fluka
flux_nue_fluka_norm = flux_nue_fluka_interesting * integral_nue_honda / integral_nue_fluka
flux_numubar_fluka_norm = flux_numubar_fluka_interesting * integral_numubar_honda / integral_numubar_fluka
flux_numu_fluka_norm = flux_numu_fluka_interesting * integral_numu_honda / integral_numu_fluka

# flux at JUNO site is given by normalized FLUKA flux in 1/(MeV * s * cm**2) (np.array of float):
flux_nuebar_juno = flux_nuebar_fluka_norm * 0.67 + flux_numubar_fluka_norm * 0.17
flux_nue_juno = flux_nue_fluka_norm * 0.67 + flux_numu_fluka_norm * 0.17

""" event rate of CC interaction of neutrinos with C12 from GENIE: 
"""
# path, where cross-sections are saved:
path_xsec_GENIE = "/home/astro/blum/juno/GENIE/genie_xsec_2.12.0_eventrate/genie_xsec/v2_12_0/NULL/" \
                  "DefaultPlusMECWithNC/data/"

""" nu_e_bar + C12 -> ... :"""
# define path, where cross-sections are saved (string):
path_xsec_nuebar_GENIE = path_xsec_GENIE + "gxspl_FNALsmall_nuebar_C12_CC.xml"

# calculate total cross-section with function 'read_xml_xsec()' (total cross-section in cm**2 for neutrino energies
# from 0 MeV to 10000 MeV, array of float):
xsec_nuebar_GENIE = read_xml_xsec(path_xsec_nuebar_GENIE, interval_energy)

# get index that corresponds to E_nu_min:
index_E_min = int((E_nu_min - 0.0) / interval_energy)
index_E_max = int((E_nu_max - 0.0) / interval_energy) + 1

# take only cross-section in energy range from E_nu_min to E_nu_max:
xsec_nuebar_GENIE = xsec_nuebar_GENIE[index_E_min:index_E_max]

""" Calculate the number of C12 atoms in the LS for fiducial volume with radius = radius_cut: """
# number of C12 atoms in fiducial volume (radius in meter) (float):
number_c12 = number_c12_atoms(radius_cut/1000.0)
print("number of C12 = {0}".format(number_c12))

# event rate as function of energy in events/(MeV * s) (array of float) (equ. 1 in AtmNeuBkgStudies_DocDB3884.pdf):
evt_rate_per_energy_GENIE = number_c12 * flux_nuebar_juno * xsec_nuebar_GENIE

# integrate evt_rate_per_energy over the energy to get the total event rate in events/s (float):
evt_rate_GENIE = np.trapz(evt_rate_per_energy_GENIE, energy_neutrino)
# print("event rate GENIE (nu_e_bar + C12 -> positron + ...) = {0} events/s".format(evt_rate_GENIE))

# event rate in events/10yr:
time = 10.0 * 3.156 * 10 ** 7
evt_rate_GENIE = evt_rate_GENIE * time

""" nu_e + C12 -> ... :"""
# define path, where cross-sections are saved (string):
path_xsec_nue_GENIE = path_xsec_GENIE + "gxspl_FNALsmall_nue_C12_CC.xml"

# calculate total cross-section with function 'read_xml_xsec()' (total cross-section in cm**2 for neutrino energies
# from 0 MeV to 10000 MeV, array of float):
xsec_nue_GENIE = read_xml_xsec(path_xsec_nue_GENIE, interval_energy)

# take only cross-section in energy range from E_nu_min to E_nu_max:
xsec_nue_GENIE = xsec_nue_GENIE[index_E_min:index_E_max]

# event rate as function of energy in events/(MeV * s) (array of float) (equ. 1 in AtmNeuBkgStudies_DocDB3884.pdf):
evt_rate_per_energy_nue_GENIE = number_c12 * flux_nue_juno * xsec_nue_GENIE

# integrate evt_rate_per_energy over the energy to get the total event rate in events/s (float):
evt_rate_nue_GENIE = np.trapz(evt_rate_per_energy_nue_GENIE, energy_neutrino)
# print("event rate GENIE (nu_e + C12 -> electron + ...) = {0} events/s".format(evt_rate_nue_GENIE))

# event rate in events/10yr:
evt_rate_nue_GENIE = evt_rate_nue_GENIE * time


""" event rate of CC interaction of neutrinos with C12 from paper of Fukugita, Kohyama, Kubodera from 1988: 
    "Neutrino reaction cross-sections on C12 Target" (1-s2.0-0370269388905138_IBDonC_CrossSection.pdf)
"""
# neutrino energy in MeV (add energy of 15 MeV, where cross-section is set to 0):
energy_nue_paper_1988 = np.array([17.7, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0, 35.0, 40.0, 45.0, 50.0, 60.0, 70.0,
                                  80.0, 90.0, 100.0])

""" nu_e + C12 -> electron + N12 (ground state): """
# cross-section of nu_e + C12 -> electron + N12 in cm**2:
xsec_nue_paper_1988 = 10 ** (-42) * np.array([0.0, 0.036, 0.287, 0.772, 1.49, 2.44, 3.62, 5.03, 9.47, 15.1, 21.8, 29.2,
                                              45.2, 60.8, 74.2, 84.2, 90.6])
# interpolate cross-section with energy_neutrino:
xsec_nue_paper_1988 = np.interp(energy_neutrino, energy_nue_paper_1988, xsec_nue_paper_1988)

# event rate as function of energy in events/(MeV * s) (array of float) (equ. 1 in AtmNeuBkgStudies_DocDB3884.pdf):
evt_rate_per_energy_nue_paper_1988 = number_c12 * flux_nue_juno * xsec_nue_paper_1988
# integrate evt_rate_per_energy over the energy to get the total event rate in events/s (float):
evt_rate_nue_paper_1988 = np.trapz(evt_rate_per_energy_nue_paper_1988, energy_neutrino)
# print("event rate paper (nu_e + C12 -> electron + N12) = {0} events/s".format(evt_rate_nue_paper_1988))
# event rate in events/10yr:
evt_rate_nue_paper_1988 = evt_rate_nue_paper_1988 * time

""" nu_e_bar + C12 -> positron + B12 (ground state): """
energy_nuebar_paper_1988 = np.array([15.7, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0, 35.0, 40.0, 45.0, 50.0, 60.0,
                                     70.0, 80.0, 90.0, 100.0])
# cross-section of nu_e_bar + C12 -> positron + B12 in cm**2:
xsec_nuebar_paper_1988 = 10 ** (-42) * np.array([0.0, 0.086, 0.327, 0.711, 1.23, 1.87, 2.62, 3.48, 4.42, 7.10, 10.1,
                                                 13.2, 16.4, 22.2, 27.0, 30.5, 32.8, 34.2])
# interpolate cross-section with energy_neutrino:
xsec_nuebar_paper_1988 = np.interp(energy_neutrino, energy_nuebar_paper_1988, xsec_nuebar_paper_1988)

# event rate as function of energy in events/(MeV * s) (array of float) (equ. 1 in AtmNeuBkgStudies_DocDB3884.pdf):
evt_rate_per_energy_nuebar_paper_1988 = number_c12 * flux_nuebar_juno * xsec_nuebar_paper_1988
# integrate evt_rate_per_energy over the energy to get the total event rate in events/s (float):
evt_rate_nuebar_paper_1988 = np.trapz(evt_rate_per_energy_nuebar_paper_1988, energy_neutrino)
# print("event rate paper (nu_e + C12 -> positron + B12) = {0} events/s".format(evt_rate_nuebar_paper_1988))
# event rate in events/10yr:
evt_rate_nuebar_paper_1988 = evt_rate_nuebar_paper_1988 * time

""" event rate of CC interaction of neutrinos with C12 from paper of Kolbe, Langanke, Vogel "Weak reactions on C12 
    within the continuum random phase approximation with partial occupancies" (C12_KolbeLangankeVogel_1999_9903022.pdf):
"""
# neutrino energy in MeV
energy_nue_paper_1999 = np.array([19.8, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0])

""" nu_e + C12 -> electron + N12 (ground state): """
# cross-section of nu_e + C12 -> electron + N12 (g.s.) in cm**2:
xsec_nue_1999_gs = 10 ** (-42) * np.array([0.0, 0.284, 4.9, 14.6, 27.9, 43.1, 58.2, 71.6, 82.5, 90.3])
# interpolate cross-section with energy_neutrino:
xsec_nue_1999_gs = np.interp(energy_neutrino, energy_nue_paper_1999, xsec_nue_1999_gs)

# event rate as function of energy in events/(MeV * s) (array of float) (equ. 1 in AtmNeuBkgStudies_DocDB3884.pdf):
evt_rate_per_energy_nue_1999_gs = number_c12 * flux_nue_juno * xsec_nue_1999_gs
# integrate evt_rate_per_energy over the energy to get the total event rate in events/s (float):
evt_rate_nue_1999_gs = np.trapz(evt_rate_per_energy_nue_1999_gs, energy_neutrino)
# print("event rate paper (nu_e + C12 -> electron + N12) = {0} events/s".format(evt_rate_nue_paper_1999))
# event rate in events/10yr:
evt_rate_nue_1999_gs = evt_rate_nue_1999_gs * time

""" nu_e + C12 -> electron + N12 (total cross-section): """
# total cross-section of nu_e + C12 -> electron + N12  in cm**2:
xsec_nue_1999 = 10 ** (-42) * np.array([0.0, 0.285, 5.65, 22.3, 59.9, 131.0, 247.0, 421.0, 663.0, 979.0])
# interpolate cross-section with energy_neutrino:
xsec_nue_1999 = np.interp(energy_neutrino, energy_nue_paper_1999, xsec_nue_1999)

# event rate as function of energy in events/(MeV * s) (array of float) (equ. 1 in AtmNeuBkgStudies_DocDB3884.pdf):
evt_rate_per_energy_nue_1999 = number_c12 * flux_nue_juno * xsec_nue_1999
# integrate evt_rate_per_energy over the energy to get the total event rate in events/s (float):
evt_rate_nue_1999 = np.trapz(evt_rate_per_energy_nue_1999, energy_neutrino)
# print("event rate paper (nu_e + C12 -> electron + N12) = {0} events/s".format(evt_rate_nue_paper_1999))
# event rate in events/10yr:
evt_rate_nue_1999 = evt_rate_nue_1999 * time

""" nu_e + C12 -> electron + N12* (e.s.): """
# cross-section of nu_e_ + C12 -> electron + N12* (e.s.) in cm**2 (excited state = total - ground state):
xsec_nue_1999_es = xsec_nue_1999 - xsec_nue_1999_gs

# event rate as function of energy in events/(MeV * s) (array of float) (equ. 1 in AtmNeuBkgStudies_DocDB3884.pdf):
evt_rate_per_energy_nue_1999_es = number_c12 * flux_nue_juno * xsec_nue_1999_es
# integrate evt_rate_per_energy over the energy to get the total event rate in events/s (float):
evt_rate_nue_1999_es = np.trapz(evt_rate_per_energy_nue_1999_es, energy_neutrino)
# print("event rate paper (nu_e + C12 -> electron + N12) = {0} events/s".format(evt_rate_nue_paper_1999))
# event rate in events/10yr:
evt_rate_nue_1999_es = evt_rate_nue_1999_es * time

""" nu_e_bar + C12 -> positron + B12 (ground state): """
energy_nuebar_1999 = np.array([18.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0,
                               140.0, 150.0])
# cross-section of nu_e_bar + C12 -> positron + B12 in cm**2:
xsec_nuebar_1999_gs_data = 10 ** (-42) * np.array([0.0, 0.79, 5.01, 11.7, 19.5, 27.2, 34.1, 39.7, 43.9, 46.9, 49.0,
                                                   50.4, 51.4, 52.2, 52.9])
# interpolate cross-section with energy_neutrino:
xsec_nuebar_1999_gs = np.interp(energy_neutrino, energy_nuebar_1999, xsec_nuebar_1999_gs_data)

# event rate as function of energy in events/(MeV * s) (array of float) (equ. 1 in AtmNeuBkgStudies_DocDB3884.pdf):
evt_rate_per_energy_nuebar_1999_gs = number_c12 * flux_nuebar_juno * xsec_nuebar_1999_gs
# integrate evt_rate_per_energy over the energy to get the total event rate in events/s (float):
evt_rate_nuebar_1999_gs = np.trapz(evt_rate_per_energy_nuebar_1999_gs, energy_neutrino)
# print("event rate paper (nu_e + C12 -> positron + B12) = {0} events/s".format(evt_rate_nuebar_paper_1988))
# event rate in events/10yr:
evt_rate_nuebar_1999_gs = evt_rate_nuebar_1999_gs * time

""" nu_e_bar + C12 -> positron + B12: """
# cross-section of nu_e_bar + C12 -> positron + B12 in cm**2:
xsec_nuebar_1999_data = 10 ** (-42) * np.array([0.0, 0.8, 6.05, 18.4, 41.4, 78.0, 130.0, 199.0, 284.0, 386.0, 501.0,
                                                629.0, 768.0, 917.0, 1080.0])
# interpolate cross-section with energy_neutrino:
xsec_nuebar_1999 = np.interp(energy_neutrino, energy_nuebar_1999, xsec_nuebar_1999_data)

# event rate as function of energy in events/(MeV * s) (array of float) (equ. 1 in AtmNeuBkgStudies_DocDB3884.pdf):
evt_rate_per_energy_nuebar_1999 = number_c12 * flux_nuebar_juno * xsec_nuebar_1999
# integrate evt_rate_per_energy over the energy to get the total event rate in events/s (float):
evt_rate_nuebar_1999 = np.trapz(evt_rate_per_energy_nuebar_1999, energy_neutrino)
# print("event rate paper (nu_e + C12 -> positron + B12) = {0} events/s".format(evt_rate_nuebar_paper_1988))
# event rate in events/10yr:
evt_rate_nuebar_1999 = evt_rate_nuebar_1999 * time

""" nu_e_bar + C12 -> positron + B12 (e.s.): """
# cross-section of nu_e_bar + C12 -> positron + B12* (e.s.) in cm**2 (excited state = total - ground state):
xsec_nuebar_1999_es = xsec_nuebar_1999 - xsec_nuebar_1999_gs

# event rate as function of energy in events/(MeV * s) (array of float) (equ. 1 in AtmNeuBkgStudies_DocDB3884.pdf):
evt_rate_per_energy_nuebar_1999_es = number_c12 * flux_nuebar_juno * xsec_nuebar_1999_es
# integrate evt_rate_per_energy over the energy to get the total event rate in events/s (float):
evt_rate_nuebar_1999_es = np.trapz(evt_rate_per_energy_nuebar_1999_es, energy_neutrino)
# print("event rate paper (nu_e + C12 -> positron + B12) = {0} events/s".format(evt_rate_nuebar_paper_1988))
# event rate in events/10yr:
evt_rate_nuebar_1999_es = evt_rate_nuebar_1999_es * time

""" event rate of CC interaction nu_e_bar + proton -> positron + neutron: """
# calculate number of free protons for fiducial volume defined by radius_cut:
number_proton = 1.45E+33 * radius_cut**3 / 17700.0**3
print("number of free protons = {0}".format(number_proton))

# calculate IBD cross section in energy window define by energy_neutrino
# (equation (25) from paper 0302005_IBDcrosssection):
# positron energy defined as energy_neutrino - delta in MeV (np.array of float64):
energy_positron = energy_neutrino - delta
# positron momentum defined as sqrt(energy_positron**2-mass_positron**2) in MeV (np.array of float64):
momentum_positron = np.sqrt(energy_positron ** 2 - mass_positron ** 2)
# IBD cross-section in cm**2 (np.array of float):
sigma = (10 ** (-43) * momentum_positron * energy_positron *
         energy_neutrino ** (-0.07056 + 0.02018 * np.log(energy_neutrino) - 0.001953 * np.log(energy_neutrino) ** 3))

# calculate event rate of atmospheric CC interaction of nu_e_bar on proton in events/(MeV * s):
evt_rate_per_energy_proton = number_proton * flux_nuebar_juno * sigma
# integrate evt_rate_per_energy over the energy to get the total event rate in events/s (float):
evt_rate_proton = np.trapz(evt_rate_per_energy_proton, energy_neutrino)
# print("event rate (nu_e_bar + proton -> positron + neutron) = {0} events/s".format(evt_rate_proton))
# event rate in events/10yr:
evt_rate_proton = evt_rate_proton * time

""" event rate of CC interaction of neutrinos with C12 from paper of Yoshida "NEUTRINO-NUCLEUS REACTION CROSS 
    SECTIONS FOR LIGHT ELEMENT SYNTHESIS IN SUPERNOVA EXPLOSIONS" from 2008 (Yoshida_2008_Reference2_of_Kim_2009.pdf):
"""
# information about cross-section from Yoshida paper in xsec_C12_Yoshida_2008.ods:
# neutrino energy in MeV:
energy_Yoshida = np.arange(1, 150+1, 1)

""" nu_e_bar + C12 -> positron + neutron + ...: """
# cross-section of nu_e_bar + C12 -> positron + neutron + ... (all channels, where 1 neutron is produced) in cm**2:
xsec_nuebar_Yoshida_data = 10 ** (-42) * np.array([0.00E+00, 0.00E+00, 0.00E+00, 0.00E+00, 0.00E+00, 0.00E+00, 0.00E+00,
                                                   0.00E+00, 0.00E+00, 0.00E+00, 0.00E+00, 0.00E+00, 0.00E+00, 0.00E+00,
                                                   0.00E+00, 0.00E+00, 0.00E+00, 0.00E+00, 1.17E-03, 5.77E-03, 1.41E-02,
                                                   3.00E-02, 5.66E-02, 9.92E-02, 1.59E-01, 2.40E-01, 3.44E-01, 4.75E-01,
                                                   6.36E-01, 8.32E-01, 1.07E+00, 1.34E+00, 1.66E+00, 2.02E+00, 2.43E+00,
                                                   2.90E+00, 3.43E+00, 4.01E+00, 4.66E+00, 5.37E+00, 6.15E+00, 7.00E+00,
                                                   7.93E+00, 8.93E+00, 1.00E+01, 1.12E+01, 1.24E+01, 1.37E+01, 1.51E+01,
                                                   1.66E+01, 1.82E+01, 1.99E+01, 2.17E+01, 2.35E+01, 2.55E+01, 2.75E+01,
                                                   2.96E+01, 3.18E+01, 3.41E+01, 3.65E+01, 3.89E+01, 4.15E+01, 4.41E+01,
                                                   4.68E+01, 4.96E+01, 5.25E+01, 5.54E+01, 5.85E+01, 6.15E+01, 6.47E+01,
                                                   6.79E+01, 7.12E+01, 7.46E+01, 7.80E+01, 8.14E+01, 8.50E+01, 8.85E+01,
                                                   9.21E+01, 9.58E+01, 9.95E+01, 1.03E+02, 1.07E+02, 1.11E+02, 1.15E+02,
                                                   1.19E+02, 1.22E+02, 1.26E+02, 1.30E+02, 1.34E+02, 1.38E+02, 1.42E+02,
                                                   1.46E+02, 1.50E+02, 1.54E+02, 1.58E+02, 1.62E+02, 1.66E+02, 1.70E+02,
                                                   1.74E+02, 1.78E+02, 1.82E+02, 1.86E+02, 1.90E+02, 1.94E+02, 1.98E+02,
                                                   2.02E+02, 2.06E+02, 2.09E+02, 2.13E+02, 2.17E+02, 2.21E+02, 2.25E+02,
                                                   2.28E+02, 2.32E+02, 2.36E+02, 2.40E+02, 2.43E+02, 2.47E+02, 2.51E+02,
                                                   2.54E+02, 2.58E+02, 2.61E+02, 2.65E+02, 2.69E+02, 2.72E+02, 2.76E+02,
                                                   2.79E+02, 2.82E+02, 2.86E+02, 2.89E+02, 2.93E+02, 2.96E+02, 2.99E+02,
                                                   3.03E+02, 3.06E+02, 3.09E+02, 3.12E+02, 3.16E+02, 3.19E+02, 3.22E+02,
                                                   3.25E+02, 3.28E+02, 3.32E+02, 3.35E+02, 3.38E+02, 3.41E+02, 3.44E+02,
                                                   3.47E+02, 3.50E+02, 3.53E+02])
# interpolate cross-section with energy_neutrino:
xsec_nuebar_Yoshida = np.interp(energy_neutrino, energy_Yoshida, xsec_nuebar_Yoshida_data)

# event rate as function of energy in events/(MeV * s) (array of float) (equ. 1 in AtmNeuBkgStudies_DocDB3884.pdf):
evt_rate_per_energy_nuebar_Yoshida = number_c12 * flux_nuebar_juno * xsec_nuebar_Yoshida
# integrate evt_rate_per_energy over the energy to get the total event rate in events/s (float):
evt_rate_nuebar_Yoshida = np.trapz(evt_rate_per_energy_nuebar_Yoshida, energy_neutrino)
# event rate in events/10yr:
evt_rate_nuebar_Yoshida = evt_rate_nuebar_Yoshida * time

""" nu_e + C12 -> electron + neutron + ...: """
# cross-section of nu_e + C12 -> electron + neutron + ... (all channels, where 1 neutron is produced) in cm**2:
xsec_nue_Yoshida_data = 10 ** (-42) * np.array([0.00E+00, 0.00E+00, 0.00E+00, 0.00E+00, 0.00E+00, 0.00E+00, 0.00E+00,
                                                0.00E+00, 0.00E+00, 0.00E+00, 0.00E+00, 0.00E+00, 0.00E+00, 0.00E+00,
                                                0.00E+00, 0.00E+00, 0.00E+00, 0.00E+00, 0.00E+00, 0.00E+00, 0.00E+00,
                                                0.00E+00, 0.00E+00, 0.00E+00, 0.00E+00, 0.00E+00, 0.00E+00, 0.00E+00,
                                                0.00E+00, 0.00E+00, 0.00E+00, 7.26E-16, 5.30E-10, 1.96E-06, 2.19E-05,
                                                9.12E-05, 2.32E-04, 4.67E-04, 8.20E-04, 1.33E-03, 2.02E-03, 2.94E-03,
                                                4.11E-03, 5.56E-03, 7.34E-03, 9.48E-03, 1.20E-02, 1.51E-02, 1.87E-02,
                                                2.28E-02, 2.76E-02, 3.31E-02, 3.93E-02, 4.62E-02, 5.40E-02, 6.26E-02,
                                                7.20E-02, 8.23E-02, 9.35E-02, 1.06E-01, 1.19E-01, 1.33E-01, 1.48E-01,
                                                1.64E-01, 1.81E-01, 1.99E-01, 2.18E-01, 2.38E-01, 2.59E-01, 2.81E-01,
                                                3.04E-01, 3.28E-01, 3.53E-01, 3.80E-01, 4.07E-01, 4.35E-01, 4.65E-01,
                                                4.95E-01, 5.27E-01, 5.60E-01, 5.93E-01, 6.28E-01, 6.64E-01, 7.00E-01,
                                                7.38E-01, 7.77E-01, 8.16E-01, 8.57E-01, 8.98E-01, 9.41E-01, 9.84E-01,
                                                1.03E+00, 1.07E+00, 1.12E+00, 1.17E+00, 1.21E+00, 1.26E+00, 1.31E+00,
                                                1.36E+00, 1.41E+00, 1.46E+00, 1.51E+00, 1.57E+00, 1.62E+00, 1.67E+00,
                                                1.72E+00, 1.78E+00, 1.83E+00, 1.89E+00, 1.94E+00, 2.00E+00, 2.05E+00,
                                                2.11E+00, 2.17E+00, 2.22E+00, 2.28E+00, 2.34E+00, 2.39E+00, 2.45E+00,
                                                2.51E+00, 2.57E+00, 2.62E+00, 2.68E+00, 2.74E+00, 2.79E+00, 2.85E+00,
                                                2.91E+00, 2.96E+00, 3.02E+00, 3.07E+00, 3.13E+00, 3.18E+00, 3.24E+00,
                                                3.29E+00, 3.35E+00, 3.40E+00, 3.46E+00, 3.51E+00, 3.56E+00, 3.62E+00,
                                                3.67E+00, 3.72E+00, 3.77E+00, 3.82E+00, 3.87E+00, 3.92E+00, 3.97E+00,
                                                4.02E+00, 4.07E+00, 4.12E+00])
# interpolate cross-section with energy_neutrino:
xsec_nue_Yoshida = np.interp(energy_neutrino, energy_Yoshida, xsec_nue_Yoshida_data)

# event rate as function of energy in events/(MeV * s) (array of float) (equ. 1 in AtmNeuBkgStudies_DocDB3884.pdf):
evt_rate_per_energy_nue_Yoshida = number_c12 * flux_nue_juno * xsec_nue_Yoshida
# integrate evt_rate_per_energy over the energy to get the total event rate in events/s (float):
evt_rate_nue_Yoshida = np.trapz(evt_rate_per_energy_nue_Yoshida, energy_neutrino)
# event rate in events/10yr:
evt_rate_nue_Yoshida = evt_rate_nue_Yoshida * time

print("\nevent rate (nu_e_bar + proton -> positron + neutron) = {0:.5f} events/(10 yr)\n".format(evt_rate_proton))
print("GENIE: ")
print("event rate GENIE (nu_e_bar + C12 -> positron + ...) = {0:.5f} events/(10 yr)".format(evt_rate_GENIE))
print("event rate GENIE (nu_e + C12 -> electron + ...) = {0:.5f} events/(10 yr)\n".format(evt_rate_nue_GENIE))
print("B12:")
print("event rate 1999 (nu_e_bar + C12 -> positron + B12 = {0:.5f} events/(10 yr)".format(evt_rate_nuebar_1999))
print("event rate paper 1988 (nu_e_bar + C12 -> positron + B12 (g.s.)) = {0:.5f} events/(10 yr)"
      .format(evt_rate_nuebar_paper_1988))
print("event rate 1999 (nu_e_bar + C12 -> positron + B12 (g.s.) = {0:.5f} events/(10 yr)"
      .format(evt_rate_nuebar_1999_gs))
print("event rate 1999 (nu_e_bar + C12 -> positron + B12 (e.s.) = {0:.5f} events/(10 yr)\n"
      .format(evt_rate_nuebar_1999_es))
print("N12:")
print("event rate 1999 (nu_e + C12 -> electron + N12 = {0:.5f} events/(10 yr)".format(evt_rate_nue_1999))
print("event rate paper 1988 (nu_e + C12 -> electron + N12 (g.s.)) = {0:.5f} events/(10 yr)"
      .format(evt_rate_nue_paper_1988))
print("event rate 1999 (nu_e + C12 -> electron + N12 (g.s.) = {0:.5f} events/(10 yr)".format(evt_rate_nue_1999_gs))
print("event rate 1999 (nu_e + C12 -> electron + N12 (e.s.) = {0:.5f} events/(10 yr)\n".format(evt_rate_nue_1999_es))
print("event rate Yoshida (nu_e_bar + C12 -> positron + neutron + ... = {0:.5f} events/(10 yr)"
      .format(evt_rate_nuebar_Yoshida))
print("event rate Yoshida (nu_e + C12 -> electron + neutron + ... = {0:.5f} events/(10 yr)"
      .format(evt_rate_nue_Yoshida))

# Display cross-section and event rate:
h1 = plt.figure(1, figsize=(15, 8))
# plt.semilogy(energy_neutrino, xsec_nuebar_GENIE, "r-", label="GENIE: nu_e_bar + C12 -> positron + ...")
# plt.semilogy(energy_neutrino, xsec_nue_GENIE, "b-", label="GENIE: nu_e + C12 -> electron + ...")
plt.semilogy(energy_neutrino, sigma, "g-", label="nu_e_bar + p -> positron + n")
# plt.semilogy(energy_neutrino, xsec_nuebar_paper_1988, "r--", label="paper 1988: nu_e_bar + C12 -> positron + "
#                                                                    "B12 (g.s.)")
# plt.semilogy(energy_neutrino, xsec_nue_paper_1988, "b--", label="paper 1988: nu_e + C12 -> electron + N12 (g.s.)")
# plt.semilogy(energy_neutrino, xsec_nuebar_1999_gs, "r:", label="paper 1999: nu_e_bar + C12 -> positron + "
#                                                                "B12 (g.s.)")
# plt.semilogy(energy_neutrino, xsec_nue_1999_gs, "b:", label="paper 1999: nu_e + C12 -> electron + N12 (g.s.)")
plt.semilogy(energy_neutrino, xsec_nuebar_1999, "r-.", label="paper 1999: nu_e_bar + C12 -> positron + B12")
plt.semilogy(energy_neutrino, xsec_nue_1999, "b-.", label="paper 1999: nu_e + C12 -> electron + N12")
plt.semilogy(energy_neutrino, xsec_nuebar_1999_es, "r--", label="paper 1999: nu_e_bar + C12 -> positron + "
                                                                "B12 (e.s.)")
plt.semilogy(energy_neutrino, xsec_nue_1999_es, "b--", label="paper 1999: nu_e + C12 -> electron + N12 (e.s.)")
plt.semilogy(energy_neutrino, xsec_nuebar_Yoshida, "r-", label="Yoshida: nu_e_bar + C12 -> positron + neutron + ...")
plt.semilogy(energy_neutrino, xsec_nue_Yoshida, "b-", label="Yoshida: nu_e + C12 -> electron + neutron + ...")
plt.ylim(ymin=1E-44)
plt.xlabel("neutrino energy in MeV")
plt.ylabel("cross-section in cm^2")
plt.title("cross-section of CC interaction")
plt.grid()
plt.legend()

h2 = plt.figure(2, figsize=(15, 8))
# plt.semilogy(energy_neutrino, evt_rate_per_energy_GENIE, "r-", label="GENIE: nu_e_bar + C12 -> positron + ...")
# plt.semilogy(energy_neutrino, evt_rate_per_energy_nue_GENIE, "b-", label="GENIE: nu_e + C12 -> electron + ...")
plt.plot(energy_neutrino, evt_rate_per_energy_proton, "g-",
         label="nu_e_bar + p -> positron + n,\nN = {0:.2f} evts/(10 yr * 20 kton)".format(evt_rate_proton))
# plt.semilogy(energy_neutrino, evt_rate_per_energy_nuebar_paper_1988, "r--",
#              label="paper 1988: nu_e_bar + C12 -> positron + B12 (g.s.)")
# plt.semilogy(energy_neutrino, evt_rate_per_energy_nue_paper_1988, "b--",
#              label="paper 1988: nu_e + C12 -> electron + N12 (g.s.)")
# plt.semilogy(energy_neutrino, evt_rate_per_energy_nuebar_1999_gs, "r:",
#              label="paper 1999: nu_e_bar + C12 -> positron + B12 (g.s.)")
# plt.semilogy(energy_neutrino, evt_rate_per_energy_nue_1999_gs, "b:",
#              label="paper 1999: nu_e + C12 -> electron + N12 (g.s.)")
# plt.semilogy(energy_neutrino, evt_rate_per_energy_nuebar_1999, "r-.",
#              label="paper 1999: nu_e_bar + C12 -> positron + B12")
# plt.semilogy(energy_neutrino, evt_rate_per_energy_nue_1999, "b-.", label="paper 1999: nu_e + C12 -> electron + N12")
plt.plot(energy_neutrino, evt_rate_per_energy_nuebar_1999_es, "r--",
         label="paper 1999: nu_e_bar + C12 -> positron + B12 (e.s.),\nN = {0:.2f} evts/(10 yr * 20 kton)"
         .format(evt_rate_nuebar_1999_es))
plt.semilogy(energy_neutrino, evt_rate_per_energy_nue_1999_es, "b--",
             label="paper 1999: nu_e + C12 -> electron + N12 (e.s.)")
plt.semilogy(energy_neutrino, evt_rate_per_energy_nuebar_Yoshida, "r-",
             label="Yoshida: nu_e_bar + C12 -> positron + neutron + ...")
plt.semilogy(energy_neutrino, evt_rate_per_energy_nue_Yoshida, "b-",
             label="Yoshida: nu_e + C12 -> electron + neutron + ... (e.s.)")
plt.ylim(ymin=1E-12)
plt.xlabel("neutrino energy in MeV")
plt.ylabel("event rate in event/(MeV * s)")
plt.title("event rate of CC interaction")
plt.grid()
plt.legend()

h3 = plt.figure(3, figsize=(15, 8))
plt.plot(energy_neutrino, flux_nuebar_juno, label="flux nu_e_bar")
plt.plot(energy_neutrino, flux_nue_juno, label="flux nu_e")
plt.xlabel("neutrino energy in MeV")
plt.ylabel("neutrino flux in 1/(MeV * s * cm**2)")
plt.title("neutrino flux at Juno site")
plt.grid()
plt.legend()

""" Get evt_rate as function of positron energy: """
# take a large neutrino energy interval:
E_nu_min1 = 10.0
E_nu_min2 = 17.3
E_nu_max1 = 150.0
energy_neutrino1 = np.arange(E_nu_min1, E_nu_max1+interval_energy, interval_energy)
energy_neutrino2 = np.arange(E_nu_min2, E_nu_max1+interval_energy, interval_energy)

# Interpolate the part of the FLUKA flux in the energy range from 14 MeV to 150 MeV, in 1/(MeV*s*cm**2)
# (np.array of float):
flux_nuebar_fluka_interesting1 = np.interp(np.arange(E_nu_min1, E_nu_max1+interval_energy, interval_energy),
                                           energy_fluka, flux_nuebar_fluka)
flux_nue_fluka_interesting1 = np.interp(np.arange(E_nu_min1, E_nu_max1+interval_energy, interval_energy),
                                        energy_fluka, flux_nue_fluka)
flux_numubar_fluka_interesting1 = np.interp(np.arange(E_nu_min1, E_nu_max1+interval_energy, interval_energy),
                                            energy_fluka, flux_numubar_fluka)
flux_numu_fluka_interesting1 = np.interp(np.arange(E_nu_min1, E_nu_max1+interval_energy, interval_energy),
                                         energy_fluka, flux_numu_fluka)
flux_nuebar_fluka_interesting2 = np.interp(np.arange(E_nu_min2, E_nu_max1+interval_energy, interval_energy),
                                           energy_fluka, flux_nuebar_fluka)
flux_nue_fluka_interesting2 = np.interp(np.arange(E_nu_min2, E_nu_max1+interval_energy, interval_energy),
                                        energy_fluka, flux_nue_fluka)
flux_numubar_fluka_interesting2 = np.interp(np.arange(E_nu_min2, E_nu_max1+interval_energy, interval_energy),
                                            energy_fluka, flux_numubar_fluka)
flux_numu_fluka_interesting2 = np.interp(np.arange(E_nu_min2, E_nu_max1+interval_energy, interval_energy),
                                         energy_fluka, flux_numu_fluka)

# Normalize flux_nuebar_fluka_interesting at Super-K to the electron-antineutrino flux at JUNO,
# in 1/(MeV * s * cm**2) (np.array of float):
flux_nuebar_fluka_norm1 = flux_nuebar_fluka_interesting1 * integral_nuebar_honda / integral_nuebar_fluka
flux_nue_fluka_norm1 = flux_nue_fluka_interesting1 * integral_nue_honda / integral_nue_fluka
flux_numubar_fluka_norm1 = flux_numubar_fluka_interesting1 * integral_numubar_honda / integral_numubar_fluka
flux_numu_fluka_norm1 = flux_numu_fluka_interesting1 * integral_numu_honda / integral_numu_fluka
flux_nuebar_fluka_norm2 = flux_nuebar_fluka_interesting2 * integral_nuebar_honda / integral_nuebar_fluka
flux_nue_fluka_norm2 = flux_nue_fluka_interesting2 * integral_nue_honda / integral_nue_fluka
flux_numubar_fluka_norm2 = flux_numubar_fluka_interesting2 * integral_numubar_honda / integral_numubar_fluka
flux_numu_fluka_norm2 = flux_numu_fluka_interesting2 * integral_numu_honda / integral_numu_fluka

# flux at JUNO site is given by normalized FLUKA flux in 1/(MeV * s * cm**2) (consider oscillation) (np.array of float):
flux_nuebar_juno1 = flux_nuebar_fluka_norm1 * 0.67 + flux_numubar_fluka_norm1 * 0.17
flux_nue_juno1 = flux_nue_fluka_norm1 * 0.67 + flux_numu_fluka_norm1 * 0.17
flux_nuebar_juno2 = flux_nuebar_fluka_norm2 * 0.67 + flux_numubar_fluka_norm2 * 0.17
flux_nue_juno2 = flux_nue_fluka_norm2 * 0.67 + flux_numu_fluka_norm2 * 0.17

# nu_e_bar + proton -> positron + neutron:
# calculate IBD cross section in energy window define by energy_neutrino1
# (equation (25) from paper 0302005_IBDcrosssection):
# positron energy defined as energy_neutrino - delta in MeV (np.array of float64):
energy_positron1 = energy_neutrino1 - delta
# positron momentum defined as sqrt(energy_positron**2-mass_positron**2) in MeV (np.array of float64):
momentum_positron1 = np.sqrt(energy_positron1 ** 2 - mass_positron ** 2)
# IBD cross-section in cm**2 (np.array of float):
sigma1 = (10 ** (-43) * momentum_positron1 * energy_positron1 *
          energy_neutrino1 ** (-0.07056 + 0.02018 * np.log(energy_neutrino1) -
                               0.001953 * np.log(energy_neutrino1) ** 3))

# calculate event rate of atmospheric CC interaction of nu_e_bar on proton in events/(MeV * s):
evt_rate_per_energy_proton1 = number_proton * flux_nuebar_juno1 * sigma1
# convolve spectrum with energy resolution and interaction kinematics to get spectrum as function of positron energy:
spectrum_proton = convolution_new(energy_neutrino1, energy, 0.5, evt_rate_per_energy_proton1,
                                  mass_proton, mass_neutron, mass_positron, "proton")
# integrate spectrum_proton over the energy to get the total event rate in events/s (float):
evt_rate_proton1 = np.trapz(spectrum_proton, energy)
# event rate in events/10yr:
evt_rate_proton1 = evt_rate_proton1 * time

# nu_e_bar + C12 -> positron + neutron + ...:
# interpolate cross-section with energy_neutrino1:
xsec_nuebar_Yoshida_1 = np.interp(energy_neutrino2, energy_Yoshida, xsec_nuebar_Yoshida_data)
# event rate as function of energy in events/(MeV * s) (array of float) (equ. 1 in AtmNeuBkgStudies_DocDB3884.pdf):
evt_rate_per_energy_nuebar_Yoshida_1 = number_c12 * flux_nuebar_juno2 * xsec_nuebar_Yoshida_1
# convolve spectrum with energy resolution and interaction kinematics to get spectrum as function of positron energy:
spectrum_Yoshida = convolution_new(energy_neutrino2, energy, 0.5, evt_rate_per_energy_nuebar_Yoshida_1, mass_proton,
                                   mass_neutron, mass_positron, "n+B11")
# integrate evt_rate_per_energy over the energy to get the total event rate in events/s (float):
evt_rate_Yoshida = np.trapz(spectrum_Yoshida, energy)
# event rate in events/10yr:
evt_rate_Yoshida = evt_rate_Yoshida * time

h4 = plt.figure(4, figsize=(15, 8))
plt.plot(energy, spectrum_proton, "g-",
         label="nu_e_bar + p -> positron + n,\nN = {0:.2f} evts/(10 yr * 20 kton)".format(evt_rate_proton1))
plt.plot(energy, spectrum_Yoshida, "r-",
         label="Yoshida: nu_e_bar + C12 -> positron + neutron + ...,\nN = {0:.2f} evts/(10 yr * 20 kton)"
         .format(evt_rate_Yoshida))
plt.plot(energy, spectrum_proton+spectrum_Yoshida, "b-", label="atmospheric CC background")
plt.ylim(ymin=1E-12)
plt.xlabel("positron energy in MeV")
plt.ylabel("event rate in event/(MeV * s)")
plt.title("event rate of CC interaction")
plt.grid()
plt.legend()

plt.show()

















