""" script to check the kinematics of IBD based on paper 'Precise quasielastic neutrino/nucleon cross section'
    by Strumia and Vissani

"""
import datetime
import sys
from scipy.stats import rv_discrete
from NC_background_functions import energy_resolution
import numpy as np
from matplotlib import pyplot as plt

E_nu = 100.0

# mass of positron in MeV (float constant):
MASS_POSITRON = 0.51099892
# mass of proton in MeV (float constant):
MASS_PROTON = 938.27203
# mass of neutron in MeV (float constant):
MASS_NEUTRON = 939.56536
# difference MASS_NEUTRON - MASS_PROTON in MeV (float):
DELTA = MASS_NEUTRON - MASS_PROTON

delta = (MASS_NEUTRON**2 - MASS_PROTON**2 - MASS_POSITRON**2) / (2 * MASS_PROTON)

s = 2 * MASS_PROTON * E_nu + MASS_PROTON**2

E_nu_CM = (s - MASS_PROTON**2) / (2 * np.sqrt(s))
E_pos_CM = (s - MASS_NEUTRON**2 + MASS_POSITRON**2) / (2 * np.sqrt(s))
p_pos_CM = np.sqrt((s - (MASS_NEUTRON - MASS_POSITRON)**2) * (s - (MASS_NEUTRON + MASS_POSITRON)**2)) / (2 * np.sqrt(s))

E_1 = E_nu - delta - E_nu_CM * (E_pos_CM + p_pos_CM) / MASS_PROTON

E_2 = E_nu - delta - E_nu_CM * (E_pos_CM - p_pos_CM) / MASS_PROTON

E_mean = (E_1 + E_2) / 2

print("E_1")
print(E_1)
print("E_2")
print(E_2)
print("E_mean")
print(E_mean)

















