""" script to calculate the number of background events of atmospheric muon neutrinos (and muon antineutrinos)
    interacting via charged current interaction with JUNO LS.

    All channels, where muons are generated, has to be considered:
        - nu_mu_bar + proton -> neutron + mu_plus
        - nu_mu_bar + C12 -> B12* + mu_plus
        - nu_mu + C12 -> N12* + mu_minus

    Neutrino fluxes are taken from Honda and FLUKA simulations.

    Cross-sections are taken from GENIE generator.

"""
import numpy as np
from matplotlib import pyplot as plt
import NC_background_functions

# neutrino energy in MeV:
E_neutrino_interval = 0.1
# INFO-me: take only neutrino energies between 100 and 207 MeV, because of kinematics of nu_mu_bar + p -> n + mu_plus
# minimum neutrino energy (threshold of reaction) is 107 MeV (= delta + mass_muon = 1.3 MeV + 105.66 MeV).
# maximum neutrino energy to get 100 MeV kinetic energy of muon is 207 MeV (= delta + mass_muon + 100 MeV = 207 MeV).
E_neutrino = np.arange(105, 223+E_neutrino_interval, E_neutrino_interval)
E_visible = np.arange(10, 100+0.5+0.5, 0.5)

# exposure time in seconds:
time_seconds = 10 * 3.156 * 10 ** 7
# fiducial volume cut in mm:
radius_cut = 16000
# mass of proton in MeV (reference PDG 2016) (float constant):
MASS_PROTON = 938.27203
# mass of neutron in MeV (reference PDG 2016) (float constant):
MASS_NEUTRON = 939.56536
# mass of muon in MeV (reference PDG 2019):
MASS_MUON = 105.658

# number of C12 in JUNO LS for specific radius cut:
number_C12 = NC_background_functions.number_c12_atoms(radius_cut/1000)

# number of free protons in total JUNO LS (17.7 m and 20 ktons):
N_free_protons_total = 1.45 * 10 ** 33
# number of free protons in JUNO LS for specific radius cut:
number_free_protons = N_free_protons_total * radius_cut**3 / 17700**3

""" Results of the HONDA simulation (based on the paper of Honda2015: 'Atmospheric neutrino flux calculation using
the NRLMSISE-00 atmospheric model'): """
# Neutrino energy in MeV from the table from file HONDA_juno-ally-01-01-solmin.d (is equal to neutrino energy
# in HONDA_juno-ally-01-01-solmax.d) (np.array of float):
energy_honda = 10 ** 3 * np.array([1.0000E-01, 1.1220E-01, 1.2589E-01, 1.4125E-01, 1.5849E-01, 1.7783E-01,
                                   1.9953E-01, 2.2387E-01])

""" for solar minimum (HONDA_juno-ally-01-01-solmin.d): """
# all-direction averaged flux for no oscillation for electron-neutrinos for solar minimum at the site of JUNO
# (WITHOUT mountain over the detector), in (MeV**(-1) * cm**(-2) * s**(-1)) (np.array of float):
# INFO-me: Solid angle (Raumwinkel) of the whole spherical angle  is = 4*pi sr! -> factor 4*pi must be correct!
flux_nue_min_honda = 10 ** (-7) * 4 * np.pi * np.array([2.9804E+03, 2.6365E+03, 2.3087E+03, 2.0054E+03, 1.7305E+03,
                                                        1.4821E+03, 1.2596E+03, 1.0619E+03])

# all-direction averaged flux for no oscillation for electron-antineutrinos for solar minimum at the site of JUNO
# (WITHOUT mountain over the detector), in (MeV**(-1) * cm**(-2) * s**(-1)) (np.array of float):
# INFO-me: Solid angle (Raumwinkel) of the whole spherical angle  is = 4*pi sr! -> factor 4*pi must be correct!
flux_nuebar_min_honda = 10 ** (-7) * 4 * np.pi * np.array([2.9367E+03, 2.5746E+03, 2.2332E+03, 1.9206E+03,
                                                           1.6395E+03, 1.3891E+03, 1.1679E+03, 9.7454E+02])

# all-direction averaged flux for no oscillation for muon-neutrinos for solar minimum at the site of JUNO
# (WITHOUT mountain over the detector), in (MeV**(-1) * cm**(-2) * s**(-1)) (np.array of float):
flux_numu_min_honda = 10 ** (-7) * 4 * np.pi * np.array([6.1827E+03, 5.4276E+03, 4.7443E+03, 4.1210E+03, 3.5492E+03,
                                                         3.0267E+03, 2.5556E+03, 2.1375E+03])

# all-direction averaged flux for no oscillation for muon-antineutrinos for solar minimum at the site of JUNO
# (WITHOUT mountain over the detector), in (MeV**(-1) * cm**(-2) * s**(-1)) (np.array of float):
flux_numubar_min_honda = 10 ** (-7) * 4 * np.pi * np.array([6.2903E+03, 5.5084E+03, 4.8032E+03, 4.1620E+03,
                                                            3.5763E+03, 3.0444E+03, 2.5663E+03, 2.1426E+03])

""" for solar maximum (HONDA_juno-ally-01-01-solmax.d): """
# all-direction averaged flux for no oscillation for electron-neutrinos for solar maximum at the site of JUNO
# (WITHOUT mountain over the detector), in (MeV**(-1) * cm**(-2) * s**(-1)) (np.array of float):
flux_nue_max_honda = 10 ** (-7) * 4 * np.pi * np.array([2.7743E+03, 2.4562E+03, 2.1530E+03, 1.8723E+03, 1.6181E+03,
                                                        1.3882E+03, 1.1819E+03, 9.9837E+02])

# all-direction averaged flux for no oscillation for electron-antineutrinos for solar maximum at the site of JUNO
# (WITHOUT mountain over the detector), in (MeV**(-1) * cm**(-2) * s**(-1)) (np.array of float):
flux_nuebar_max_honda = 10 ** (-7) * 4 * np.pi * np.array([2.7733E+03, 2.4332E+03, 2.1124E+03, 1.8187E+03,
                                                           1.5545E+03, 1.3190E+03, 1.1105E+03, 9.2820E+02])

# all-direction averaged flux for no oscillation for muon-neutrinos for solar maximum at the site of JUNO
# (WITHOUT mountain over the detector), in (MeV**(-1) * cm**(-2) * s**(-1)) (np.array of float):
flux_numu_max_honda = 10 ** (-7) * 4 * np.pi * np.array([5.7913E+03, 5.0884E+03, 4.4520E+03, 3.8714E+03, 3.3388E+03,
                                                         2.8520E+03, 2.4128E+03, 2.0226E+03])

# all-direction averaged flux for no oscillation for muon-antineutrinos for solar maximum at the site of JUNO
# (WITHOUT mountain over the detector), in (MeV**(-1) * cm**(-2) * s**(-1)) (np.array of float):
flux_numubar_max_honda = 10 ** (-7) * 4 * np.pi * np.array([5.8966E+03, 5.1676E+03, 4.5104E+03, 3.9127E+03,
                                                            3.3665E+03, 2.8701E+03, 2.4238E+03, 2.0277E+03])

# all-direction averaged flux for no oscillation for electron-neutrinos for solar AVERAGE at the site of JUNO
# (WITHOUT mountain over the detector), in (MeV**(-1) * cm**(-2) * s**(-1)) (np.array of float):
flux_atmo_nue_honda = (flux_nue_min_honda + flux_nue_max_honda) / 2

# all-direction averaged flux for no oscillation for electron-antineutrinos for solar AVERAGE at the site of JUNO
# (WITHOUT mountain over the detector), in (MeV**(-1) * cm**(-2) * s**(-1)) (np.array of float):
flux_atmo_nuebar_honda = (flux_nuebar_min_honda + flux_nuebar_max_honda) / 2

# all-direction averaged flux for no oscillation for muon-neutrinos for solar AVERAGE at the site of JUNO
# (WITHOUT mountain over the detector), in (MeV**(-1) * cm**(-2) * s**(-1)) (np.array of float):
flux_atmo_numu_honda = (flux_numu_min_honda + flux_numu_max_honda) / 2

# all-direction averaged flux for no oscillation for muon-antineutrinos for solar AVERAGE at the site of JUNO
# (WITHOUT mountain over the detector), in (MeV**(-1) * cm**(-2) * s**(-1)) (np.array of float):
flux_atmo_numubar_honda = (flux_numubar_min_honda + flux_numubar_max_honda) / 2

""" Interpolate the HONDA flux to the energies given by E_neutrino: """
flux_nue_juno = np.interp(E_neutrino, energy_honda, flux_atmo_nue_honda)
flux_nuebar_juno = np.interp(E_neutrino, energy_honda, flux_atmo_nuebar_honda)
flux_numu_juno = np.interp(E_neutrino, energy_honda, flux_atmo_numu_honda)
flux_numubar_juno = np.interp(E_neutrino, energy_honda, flux_atmo_numubar_honda)

""" Taking account neutrino oscillation from the Appendix B of the paper of Fogli et al. from 2004 with title
"Three-generation flavor transitions and decays of supernova relic neutrinos" (like in ccatmospheric_background_v2):
"""
# oscillation parameters from paper:
# survival probability nu_mu -> nu_mu:
p_numu_numu = 0.41
# probability of nu_e -> nu_mu (p_nue_nue + p_nue_numu + p_nue_nutau = 0.67 + p_nue_numu + 0 = 1
# -> p_nue_numu = 1 - 0.67):
p_nue_numu = 0.33

# total muon-neutrino flux in the INTERESTING part (10 to 100 MeV) of FLUKA simulation normalized to
# JUNO site (HONDA) in 1/(MeV * cm**2 * s), (np.array of float):
flux_total_ccatmospheric_nu_mu = p_nue_numu * flux_nue_juno + p_numu_numu * flux_numu_juno

# total muon-antineutrino flux in the INTERESTING part (10 to 100 MeV) of FLUKA simulation normalized to
# JUNO site (HONDA) in 1/(MeV * cm**2 * s), (np.array of float):
flux_total_ccatmospheric_nu_mu_bar = p_nue_numu * flux_nuebar_juno + p_numu_numu * flux_numubar_juno

""" get cross-section of nu_mu_bar + proton -> neutron + mu_plus, nu_mu + C12 -> N12 + mu_minus 
    and nu_mu_bar + C12 -> B12 + mu_plus from xml files from GENIE: """
# path, where cross-sections are saved:
path_xsec = "/home/astro/blum/juno/GENIE/genie_xsec_2.12.0_eventrate/genie_xsec/v2_12_0/NULL/DefaultPlusMECWithNC/data/"

# nu_mu_bar + proton -> neutron + mu_plus:
path_numubar_p = path_xsec + "gxspl-FNALsmall_numubar_proton_CC.xml"
# calculate total cross-section with function 'read_xml_xsec()' (total cross-section in cm**2, array of float):
# INFO-me: cross.section is calculated for energies from 0 MeV to 10 GeV:
xsec_numubar_p = NC_background_functions.read_xml_xsec(path_numubar_p, E_neutrino_interval)
# take cross-section only for E_neutrino:
xsec_numubar_p = np.interp(E_neutrino, np.arange(0, 10000+E_neutrino_interval, E_neutrino_interval), xsec_numubar_p)

# nu_mu_bar + C12 -> B12 + mu_plus:
path_numubar_C12 = path_xsec + "gxspl-FNALsmall_numubar_C12_CC.xml"
# calculate total cross-section with function 'read_xml_xsec()' (total cross-section in cm**2, array of float):
# INFO-me: cross.section is calculated for energies from 0 MeV to 10 GeV:
xsec_numubar_C12 = NC_background_functions.read_xml_xsec(path_numubar_C12, E_neutrino_interval)
# take cross-section only for E_neutrino:
xsec_numubar_C12 = np.interp(E_neutrino, np.arange(0, 10000+E_neutrino_interval, E_neutrino_interval), xsec_numubar_C12)

# nu_mu + C12 -> N12 + mu_minus:
path_numu_C12 = path_xsec + "gxspl-FNALsmall_numu_C12_CC.xml"
# calculate total cross-section with function 'read_xml_xsec()' (total cross-section in cm**2, array of float):
# INFO-me: cross.section is calculated for energies from 0 MeV to 10 GeV:
xsec_numu_C12 = NC_background_functions.read_xml_xsec(path_numu_C12, E_neutrino_interval)
# take cross-section only for E_neutrino:
xsec_numu_C12 = np.interp(E_neutrino, np.arange(0, 10000+E_neutrino_interval, E_neutrino_interval), xsec_numu_C12)

""" calculate number of events: """
# number of randomly generated values:
random_number = 1000000

# numubar + p -> neutron + mu_plus:
# spectrum (number of events as function of energy) in 1/MeV:
N_numubar_p_E = number_free_protons * flux_total_ccatmospheric_nu_mu_bar * xsec_numubar_p * time_seconds
# spectrum in 1/bin (not 1/MeV):
N_numubar_p_E = N_numubar_p_E * E_neutrino_interval
# normalize spectrum to 1:
N_numubar_p_E_norm = N_numubar_p_E / np.sum(N_numubar_p_E)
# generate 10000 random values of E_neutrino from theoretical spectrum (array):
array_e_neutrino_1 = np.random.choice(E_neutrino, p=N_numubar_p_E_norm, size=random_number)
# preallocate array, where visible energies are stored:
array_Evis_1 = []
# loop over entries in array_e_neutrino:
for index in range(len(array_e_neutrino_1)):
    # interaction channel (nu_mu_bar + p -> neutron + mu_plus;
    # muon energy in MeV (very simple kinematics (nu_mu_bar + p -> neutron + mu_plus;
    # E_mu = E_nu + m_p - m_n - m_muon, kinetic energy of neutron is neglected.
    Evis = array_e_neutrino_1[index] + MASS_PROTON - MASS_NEUTRON - MASS_MUON
    # check, if vis_neutrino is positive:
    if Evis <= 0:
        continue
    # append vis_neutrino to array_vis_neutrino:
    array_Evis_1.append(Evis)
# build histogram of visible spectrum in 1/bin from array_vis_neutrino in whole visible energy range:
Spectrum_numubar_p, bin_edges = np.histogram(array_Evis_1, bins=E_visible)
# normalize spectrum with number of events and number of random numbers:
Spectrum_numubar_p = Spectrum_numubar_p * np.sum(N_numubar_p_E) / random_number
# number of events:
N_numubar_p = np.sum(Spectrum_numubar_p)

# numubar + C12 -> B12 + mu_plus:
# spectrum (number of events as function of energy) in 1/MeV:
N_numubar_C12_E = number_C12 * flux_total_ccatmospheric_nu_mu_bar * xsec_numubar_C12 * time_seconds
# spectrum in 1/bin (not 1/MeV):
N_numubar_C12_E = N_numubar_C12_E * E_neutrino_interval
# normalize spectrum to 1:
N_numubar_C12_E_norm = N_numubar_C12_E / np.sum(N_numubar_C12_E)
# generate 10000 random values of E_neutrino from theoretical spectrum (array):
array_e_neutrino_2 = np.random.choice(E_neutrino, p=N_numubar_C12_E_norm, size=random_number)
# preallocate array, where visible energies are stored:
array_Evis_2 = []
# loop over entries in array_e_neutrino:
for index in range(len(array_e_neutrino_2)):
    # interaction channel (nu_mu_bar + C12 -> B12 + mu_plus;
    # muon energy in MeV (very simple kinematics (nu_mu_bar + C12 -> B12 + mu_plus;
    # E_mu = E_nu + m_C12 - m_B12 - m_muon, kinetic energy of B12 is neglected.
    Evis = array_e_neutrino_2[index] + 12.0*931.494 - 12.014*931.494 - MASS_MUON
    # check, if vis_neutrino is positive:
    if Evis <= 0:
        continue
    # append vis_neutrino to array_vis_neutrino:
    array_Evis_2.append(Evis)
# build histogram of visible spectrum in 1/bin from array_vis_neutrino in whole visible energy range:
Spectrum_numubar_C12, bin_edges = np.histogram(array_Evis_2, bins=E_visible)
# normalize spectrum with number of events and number of random numbers:
Spectrum_numubar_C12 = Spectrum_numubar_C12 * np.sum(N_numubar_C12_E) / random_number
# number of events:
N_numubar_C12 = np.sum(Spectrum_numubar_C12)

# numu + C12 -> N12 + mu_minus:
# spectrum (number of events as function of energy) in 1/MeV:
N_numu_C12_E = number_C12 * flux_total_ccatmospheric_nu_mu * xsec_numu_C12 * time_seconds
# spectrum in 1/bin (not 1/MeV):
N_numu_C12_E = N_numu_C12_E * E_neutrino_interval
# normalize spectrum to 1:
N_numu_C12_E_norm = N_numu_C12_E / np.sum(N_numu_C12_E)
# generate 10000 random values of E_neutrino from theoretical spectrum (array):
array_e_neutrino_3 = np.random.choice(E_neutrino, p=N_numu_C12_E_norm, size=random_number)
# preallocate array, where visible energies are stored:
array_Evis_3 = []
# loop over entries in array_e_neutrino:
for index in range(len(array_e_neutrino_3)):
    # interaction channel (nu_mu + C12 -> N12 + mu_minus;
    # muon energy in MeV (very simple kinematics (nu_mu + C12 -> N12 + mu_minus;
    # E_mu = E_nu + m_C12 - m_N12 - m_muon, kinetic energy of N12 is neglected.
    Evis = array_e_neutrino_3[index] + 12.0*931.494 - 12.018*931.494 - MASS_MUON

    # check, if vis_neutrino is positive:
    if Evis <= 0:
        continue

    # append vis_neutrino to array_vis_neutrino:
    array_Evis_3.append(Evis)
# build histogram of visible spectrum in 1/bin from array_vis_neutrino in whole visible energy range:
Spectrum_numu_C12, bin_edges = np.histogram(array_Evis_3, bins=E_visible)
# normalize spectrum with number of events and number of random numbers:
Spectrum_numu_C12 = Spectrum_numu_C12 * np.sum(N_numu_C12_E) / random_number
# number of events:
N_numu_C12 = np.sum(Spectrum_numu_C12)

# total number of events as function of energy (in 1/MeV):
Spectrum_total = Spectrum_numubar_p + Spectrum_numubar_C12 + Spectrum_numu_C12
# total number of events:
N_total = N_numubar_p + N_numubar_C12 + N_numu_C12

""" consider cut efficiency from Giulio (eff for muCC = 6.2 %): """
cut_efficiency = 0.062
Spectrum_numubar_p_eff = Spectrum_numubar_p * cut_efficiency
N_numubar_p_eff = N_numubar_p * cut_efficiency
Spectrum_numubar_C12_eff = Spectrum_numubar_C12 * cut_efficiency
N_numubar_C12_eff = N_numubar_C12 * cut_efficiency
Spectrum_numu_C12_eff = Spectrum_numu_C12 * cut_efficiency
N_numu_C12_eff = N_numu_C12 * cut_efficiency
Spectrum_total_eff = Spectrum_total * cut_efficiency
N_total_eff = N_total * cut_efficiency

h1 = plt.figure(1)
plt.step(E_visible[:-1], Spectrum_total, "k-",
         label="total spectrum ({0:.1f} events)".format(N_total))
plt.step(E_visible[:-1], Spectrum_numubar_p, "b--",
         label="$\\bar{\\nu}_{\\mu}$ + p $\\rightarrow$ n + $\\mu^+$\n" + "({0:.1f} events)".format(N_numubar_p))
plt.step(E_visible[:-1], Spectrum_numubar_C12, "r--",
         label="$\\bar{\\nu}_{\\mu} + ^{12}C \\rightarrow ^{12}B* + \\mu^+$\n" +
               "({0:.1f} events)".format(N_numubar_C12))
plt.step(E_visible[:-1], Spectrum_numu_C12, "g--",
         label="$\\nu_{\\mu} + ^{12}C \\rightarrow ^{12}N* + \\mu^-$\n" + "({0:.1f} events)".format(N_numu_C12))
plt.xlabel("visible energy in MeV")
plt.ylabel("events / bin")
plt.title("Atmospheric $\\nu_{\\mu}$/$\\bar{\\nu}_{\\mu}$ CC background spectrum in JUNO after 10 years\n")
plt.grid()
plt.legend()

h2 = plt.figure(2)
plt.step(E_visible[:-1], Spectrum_total_eff, "k-",
         label="total spectrum ({0:.1f} events)".format(N_total_eff))
plt.step(E_visible[:-1], Spectrum_numubar_p_eff, "b--",
         label="$\\bar{\\nu}_{\\mu}$ + p $\\rightarrow$ n + $\\mu^+$\n" + "({0:.1f} events)".format(N_numubar_p_eff))
plt.step(E_visible[:-1], Spectrum_numubar_C12_eff, "r--",
         label="$\\bar{\\nu}_{\\mu} + ^{12}C \\rightarrow ^{12}B* + \\mu^+$\n" +
               "({0:.1f} events)".format(N_numubar_C12_eff))
plt.step(E_visible[:-1], Spectrum_numu_C12_eff, "g--",
         label="$\\nu_{\\mu} + ^{12}C \\rightarrow ^{12}N* + \\mu^-$\n" + "({0:.1f} events)".format(N_numu_C12_eff))
plt.xlabel("visible energy in MeV")
plt.ylabel("events / bin")
plt.title("Atmospheric $\\nu_{\\mu}$/$\\bar{\\nu}_{\\mu}$ CC background spectrum in JUNO after 10 years\n"
          "($\\nu_{\\mu}$ CC cut efficiency = 6.2 %)")
plt.grid()
plt.legend()
plt.show()
