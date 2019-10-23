""" scritp to check the kinematics of IBD based on paper 'angular distribution of neutron inverse beta decay' by
    Vogel and Beacom (this kinematics is also used in IBD generator of JUNO offline): """

import datetime
import sys
from NC_background_functions import energy_resolution
import numpy as np
from matplotlib import pyplot as plt

# define neutrino energy in MeV:
E_nu = 40.0
# number of neutrino events:
number_events = 100000

# set angle in degree:
theta_min = 0.0
theta_max = 180.0

# mass of positron in MeV (float constant):
MASS_POSITRON = 0.51099892
# mass of proton in MeV (float constant):
MASS_PROTON = 938.27203
# mass of neutron in MeV (float constant):
MASS_NEUTRON = 939.56536
# difference MASS_NEUTRON - MASS_PROTON in MeV (float):
DELTA = MASS_NEUTRON - MASS_PROTON

# calculate positron energy in MeV at zeroth order of 1/M (equ. 6):
E_pos_0 = E_nu - DELTA
# calculate positron momentum in MeV:
p_pos_0 = np.sqrt(E_pos_0**2 - MASS_POSITRON**2)
# calculate positron velocity in MeV:
v_pos_0 = p_pos_0 / E_pos_0

# preallocate array, where positron energies are stored:
array_E_pos = []
array_E_pos_smeared = []
array_theta = []

# generate neutrino events:
for index in range(number_events):

    # generate angle theta between positron and neutron from uniform distribution:
    theta = np.random.uniform(theta_min, theta_max)

    array_theta.append(theta)

    # calculate positron energy in MeV at first order of 1/M (equ. 11):
    E_pos_1 = (E_pos_0 * (1 - E_nu / MASS_PROTON * (1 - v_pos_0 * np.cos(np.deg2rad(theta)))) -
               (DELTA**2 - MASS_POSITRON**2) / (2 * MASS_PROTON))

    # calculate visible energy of positron in MeV:
    E_visible = E_pos_1 + MASS_POSITRON

    array_E_pos.append(E_visible)

    # get sigma from energy resolution in MeV:
    sigma = energy_resolution(E_visible)

    # smear E_visible with sigma:
    E_visible_smeared = np.random.normal(E_visible, sigma)

    # append to array:
    array_E_pos_smeared.append(E_visible_smeared)

# build histogram of array_E_pos:
h1 = plt.figure(1, figsize=(15, 8))
First_bin = 0.0
Last_bin = 100.0
Bin_width = 0.5
Bins = np.arange(First_bin, Last_bin+Bin_width, Bin_width)
plt.hist(array_E_pos, bins=Bins, histtype="step", align='mid', color="r", linewidth=1.5,
         label="w/o energy resolution (entries = {0:d})".format(number_events))
plt.hist(array_E_pos_smeared, bins=Bins, histtype="step", align='mid', color="b", linewidth=1.5,
         label="w/ energy resolution (entries = {0:d})".format(number_events))
plt.xlabel("visible energy of positron in MeV")
plt.ylabel("events")
plt.title("DM mass = {0:.0f} MeV".format(E_nu))
plt.legend()
plt.grid()
# plt.savefig()
# plt.close()

h2 = plt.figure(2, figsize=(15, 8))
First_bin = theta_min
Last_bin = theta_max
Bin_width = 0.5
Bins = np.arange(First_bin, Last_bin+Bin_width, Bin_width)
plt.hist(array_theta, bins=Bins, histtype="step", align='mid', color="r", linewidth=1.5,
         label="entries = {0:d}".format(number_events))
plt.xlabel("theta in degree")
plt.ylabel("events")
# plt.title("DM mass = {0:.0f} MeV".format(E_nu))
plt.legend()
plt.grid()
# plt.savefig()
# plt.close()

h3 = plt.figure(3, figsize=(15, 8))
Bin_width = 0.1
Bins_theta = np.arange(theta_min, theta_max+Bin_width, Bin_width)
Bins_E = np.arange(33.0, 41.0+Bin_width, Bin_width)
plt.hist2d(array_theta, array_E_pos_smeared, [Bins_theta, Bins_E])
plt.xlabel("theta in degree")
plt.ylabel("visible energy in MeV")
# plt.title("DM mass = {0:.0f} MeV".format(E_nu))
plt.grid()
# plt.savefig()
# plt.close()



plt.show()
















