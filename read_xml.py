""" script to read xml file and store values in arrays: """

import xml.etree.ElementTree as ET
import numpy as np
from matplotlib import pyplot as plt

# xml file:
filename = '/home/astro/blum/juno/GENIE/genie_xsec/v2_12_0/NULL/DefaultPlusMECWithNC/data/gxspl-FNALsmall_numu.xml'
# read xml file:
tree = ET.parse(filename)
# get root from xml file, root = genie_xsec_spline_list:
root = tree.getroot()

# energy interval (bin-width) in MeV:
interval_energy = 0.1
# preallocate energy-array in MeV:
energy = np.arange(0, 10000+interval_energy, interval_energy)
# preallocate cross-section array in cm**2:
cross_section = np.zeros(len(energy))

# preallocate cross-section array in cm**2 for QES:
cross_section_qes = np.zeros(len(energy))
# preallocate cross-section array in cm**2 for DIS:
cross_section_dis = np.zeros(len(energy))
# preallocate cross-section array in cm**2 for RES:
cross_section_res = np.zeros(len(energy))
# preallocate cross-section array in cm**2 for COH:
cross_section_coh = np.zeros(len(energy))

# loop over root, child = spline
for child in root:
    # get attrib of child (dict):
    spline = child.attrib

    # get the name attribute from spline-dictionary
    # (e.g. spline_name = 'genie::AhrensNCELPXSec/Default/nu:14;tgt:1000060120;N:2112;proc:Weak[NC],QES;', type=string)
    spline_name = spline.get('name')
    print(spline_name)

    # INFO-me: for E=0MeV, set xsec=0cm**2 -> if not, interp sets all values below 10 MeV to value of xsec(10MeV)!
    # preallocate energy-array for this knot (MeV):
    e_spline = np.array([0])
    # preallocate cross-section-array for this knot (in cm**2):
    xsec_spline = np.array([0])

    # loop over child, subchild = knot, subchild is list with 2 values, first value = E, second value = xsec:
    for subchild in child:
        # get the value of E in GeV (string)
        E = subchild[0].text
        # convert to float and MeV:
        E = float(E) * 1000

        # get the value of cross-section in natural units (xsec) (string):
        xsec = subchild[1].text
        # convert to float and cm**2 (natural unit: 1/GeV**2 entspricht 3.89391289*10**(-28) cm**2):
        xsec = float(xsec) * 3.89391289 * 10**(-28)

        # append E to energy array:
        e_spline = np.append(e_spline, E)
        # append xsec to cross-section-array:
        xsec_spline = np.append(xsec_spline, xsec)

    # interpolate the cross-section of this spline to array "energy" to have the same bin-width and energy range for
    # all cross-section (also same bin-width and energy range like the fluxes):
    xsec_interp = np.interp(energy, e_spline, xsec_spline)

    # check the cross-section for different interactions (QES, DIS, RES, COH):
    if spline_name[:13] == "genie::Ahrens":
        # only QES interactions:
        cross_section_qes = cross_section_qes + xsec_interp

    elif spline_name[:13] == "genie::QPMDIS":
        # only DIS interactions:
        cross_section_dis = cross_section_dis + xsec_interp

    elif spline_name[:20] == "genie::ReinSehgalRES":
        # only RES interactions:
        cross_section_res = cross_section_res + xsec_interp

    else:
        # only COH interactions:
        cross_section_coh = cross_section_coh + xsec_interp

    # add cross-section to the total cross-section in cm**2
    cross_section = cross_section + xsec_interp

plt.plot(energy, cross_section_qes, 'b--', label='QES')
plt.plot(energy, cross_section_dis, 'r--', label='DIS')
plt.plot(energy, cross_section_res, 'g--', label='RES')
plt.plot(energy, cross_section_coh, '--', label='COH')
plt.plot(energy, cross_section, 'k', label='total')
plt.xlabel("neutrino energy in MeV")
plt.ylabel("cross-section of $\\nu_\\mu$ on $^{12}$C in cm$^2$")
plt.title("Muon-neutrino cross-section on $^{12}$C from GENIE")
plt.xlim(xmin=0.0)
plt.ylim(ymin=0.0)
plt.legend()
plt.grid()
plt.show()



