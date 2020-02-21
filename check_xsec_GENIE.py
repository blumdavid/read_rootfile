""" script to check the used cross-section from GENIE event generator: """
import NC_background_functions
import numpy as np
from matplotlib import pyplot as plt

# path, where GENIE cross-sections are saved:
path_xsec = "/home/astro/blum/juno/GENIE/genie_xsec_2.12.0_eventrate/genie_xsec/v2_12_0/NULL/DefaultPlusMECWithNC/data/"
# set energy interval in MeV:
interval_energy = 10
# set energy range in MeV (do NOT change this, it is hardcoded in function read_xml_xsec()):
energy = np.arange(0, 10000+interval_energy, interval_energy)

""" Neutral Current interaction cross-sections of neutrinos with C12 for each neutrino flavour: """
# NC interaction nu_e + C12 -> nu_e + ...:
# define path, where cross-sections are saved (string):
path_xsec_NC_nue_C12 = path_xsec + "gxspl-FNALsmall_nue.xml"
# calculate total cross-section with function 'read_xml_xsec()' (total cross-section in cm**2, array of float):
xsec_NC_nue_C12 = NC_background_functions.read_xml_xsec(path_xsec_NC_nue_C12, interval_energy)

# NC interaction nu_e_bar + C12 -> nu_e_bar + ...:
# define path, where cross-sections are saved (string):
path_xsec_NC_nuebar_C12 = path_xsec + "gxspl-FNALsmall_nuebar.xml"
# calculate total cross-section with function 'read_xml_xsec()' (total cross-section in cm**2, array of float):
xsec_NC_nuebar_C12 = NC_background_functions.read_xml_xsec(path_xsec_NC_nuebar_C12, interval_energy)

# NC interaction nu_mu + C12 -> nu_mu + ...:
# define path, where cross-sections are saved (string):
path_xsec_NC_numu_C12 = path_xsec + "gxspl-FNALsmall_numu.xml"
# calculate total cross-section with function 'read_xml_xsec()' (total cross-section in cm**2, array of float):
xsec_NC_numu_C12 = NC_background_functions.read_xml_xsec(path_xsec_NC_numu_C12, interval_energy)

# NC interaction nu_mu_bar + C12 -> nu_mu_bar + ...:
# define path, where cross-sections are saved (string):
path_xsec_NC_numubar_C12 = path_xsec + "gxspl-FNALsmall_numubar.xml"
# calculate total cross-section with function 'read_xml_xsec()' (total cross-section in cm**2, array of float):
xsec_NC_numubar_C12 = NC_background_functions.read_xml_xsec(path_xsec_NC_numubar_C12, interval_energy)

plt.figure(1, figsize=(11, 6))
plt.semilogy(energy, xsec_NC_nue_C12, "b-", label="NC: nu_e + C12 -> nu_e + ...")
plt.semilogy(energy, xsec_NC_nuebar_C12, "b--", label="NC: nu_e_bar + C12 -> nu_e_bar + ...")
plt.semilogy(energy, xsec_NC_numu_C12, "r-", label="NC: nu_mu + C12 -> nu_mu + ...")
plt.semilogy(energy, xsec_NC_numubar_C12, "r--", label="NC: nu_mu_bar + C12 -> nu_mu_bar + ...")
plt.xlabel("neutrino energy in MeV")
plt.ylabel("cross-section in $cm^2$")
plt.grid()
plt.legend()
plt.show()
