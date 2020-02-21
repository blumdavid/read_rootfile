""" script to create hepevt files (as input for JUNO detsim) containing IBD signals.

    IBD signal:
    - positron with initial momentum between 10 MeV and 100 MeV
    - neutron with initial momentum corresponding to the positron energy

    format of hepevt file (momentum and mass in GeV):

    number particles in event
    1   PDGID 0 0 px py pz mass
    1   PDGID 0 0 px py pz mass

    for example:
    2
    1   -11 0 0 px py pz 0.000511
    1   2112 0 0 px py pz 0.939570



"""
import numpy as np
from matplotlib import pyplot as plt


def calculate_neutron_momentum(mom_positron, theta, mass_p, delta, mass_pos):
    """
    function to calculate the neutron momentum as function of the positron momentum and the angle theta

    from paper "Angular distribution of neutron inverse beta decay" of Beacom and Vogel (PhysRevD.60.053003.pdf),
    page 7, equation 29.

    In the case of large positron energies (E_pos > 10 MeV): momentum of positron = energy of positron
    (mass can be neglected) and therefore velocity_positron = momentum / energy = 1.

    For example: mom_pos(E_pos = 10 MeV) = 9.987 MeV,   mom_pos(E_pos = 100 MeV) = 99.9987 MeV

    With this approximation equation 29 of the paper becomes function below.

    :param mom_positron: momentum of positron in GeV
    :param theta: angle theta in degree
    :param mass_p: mass of proton in GeV
    :param delta: mass neutron - mass proton in GeV
    :param mass_pos: mass of positron in GeV
    :return: momentum of neutron in GeV
    """
    mom_neutron = ((mom_positron - delta) * mom_positron / mass_p * (1 - np.cos(np.deg2rad(theta))) +
                   (delta**2 - mass_pos**2)/(2 * mass_p))

    return mom_neutron


# mass of positron in GeV (float constant):
MASS_POSITRON = 0.51099892/1000.0
# mass of proton in GeV (float constant):
MASS_PROTON = 938.27203/1000.0
# mass of neutron in GeV (float constant):
MASS_NEUTRON = 939.56536/1000.0
# difference MASS_NEUTRON - MASS_PROTON in GeV (float):
DELTA = MASS_NEUTRON - MASS_PROTON

# set minimum and maximum momentum of positron in GeV:
E_min = 0.010
E_max = 0.100

# set path, where the output files (hepevt files should be saved):
output_path = "/home/astro/blum/juno/IBD_events/IBD_hepevt_files/"

# set the number of events, that should be stored in one hepevt file:
number_events_per_file = 100

# set the number of files that should be created:
number_of_files = 200

# parameters for the hepevt file:
# PDG ID of positron:
pdg_positron = str(-11)
# mass of positron in GeV:
mass_positron = str(0.000511)
# PDG ID of neutron:
pdg_neutron = str(2112)
# mass of neutron in GeV:
mass_neutron = str(0.93957)

# preallocate array to check positron momentum:
array_pos = []
array_pos_test = []
array_neutron = []

# loop over the number of files:
for filenumber in range(number_of_files):

    # preallocate hepevt array:
    hepevt_file = []

    # loop over the number of events per file:
    for event in range(number_events_per_file):

        # generate total momentum of positron in GeV:
        momentum_positron = np.random.uniform(E_min, E_max)

        # calculate square of the total momentum in GeV**2:
        momentum_positron_square = momentum_positron**2

        # generate square of x-momentum of positron in GeV**2:
        momentum_positron_x_square = np.random.uniform(0.0, momentum_positron_square)

        # generate square of y-momentum of positron in GeV**2:
        momentum_positron_y_square = np.random.uniform(0.0, (momentum_positron_square - momentum_positron_x_square))

        # calculate square of z-momentum of positron in GeV**2:
        momentum_positron_z_square = momentum_positron_square - momentum_positron_x_square - momentum_positron_y_square

        # calculate the momentum in x, y, z direction in GeV (consider that it can also be negative):
        # generate random 0 or 1:
        sign_x = np.random.randint(2)
        if sign_x == 0:
            momentum_positron_x = np.sqrt(momentum_positron_x_square)
        else:
            momentum_positron_x = -np.sqrt(momentum_positron_x_square)

        sign_y = np.random.randint(2)
        if sign_y == 0:
            momentum_positron_y = np.sqrt(momentum_positron_y_square)
        else:
            momentum_positron_y = -np.sqrt(momentum_positron_y_square)

        sign_z = np.random.randint(2)
        if sign_z == 0:
            momentum_positron_z = np.sqrt(momentum_positron_z_square)
        else:
            momentum_positron_z = -np.sqrt(momentum_positron_z_square)

        # calculate momentum of the generated x,y,z momenta in GeV as cross-check:
        momentum_positron_crosscheck = np.sqrt(momentum_positron_x**2 + momentum_positron_y**2 + momentum_positron_z**2)

        # generate angle theta in degree:
        theta_degree = np.random.uniform(0.0, 180.0)

        # calculate the momentum of neutron in GeV:
        momentum_neutron = calculate_neutron_momentum(momentum_positron, theta_degree, MASS_PROTON, DELTA,
                                                      MASS_POSITRON)

        # calculate square of neutron momentum in GeV**2:
        momentum_neutron_square = momentum_neutron**2

        # generate square of x-momentum of neutron in GeV**2:
        momentum_neutron_x_square = np.random.uniform(0.0, momentum_neutron_square)

        # generate square of y-momentum of neutron in GeV**2:
        momentum_neutron_y_square = np.random.uniform(0.0, (momentum_neutron_square - momentum_neutron_x_square))

        # calculate square of z-momentum of neutron in GeV**2:
        momentum_neutron_z_square = momentum_neutron_square - momentum_neutron_x_square - momentum_neutron_y_square

        # calculate the momentum in x, y, z, direction in GeV (consider that it can also be negative):
        # generate random 0 or 1:
        sign_x = np.random.randint(2)
        if sign_x == 0:
            momentum_neutron_x = np.sqrt(momentum_neutron_x_square)
        else:
            momentum_neutron_x = -np.sqrt(momentum_neutron_x_square)

        sign_y = np.random.randint(2)
        if sign_y == 0:
            momentum_neutron_y = np.sqrt(momentum_neutron_y_square)
        else:
            momentum_neutron_y = -np.sqrt(momentum_neutron_y_square)

        sign_z = np.random.randint(2)
        if sign_z == 0:
            momentum_neutron_z = np.sqrt(momentum_neutron_z_square)
        else:
            momentum_neutron_z = -np.sqrt(momentum_neutron_z_square)

        """ fill array as cross-check: """
        array_pos.append(momentum_positron)
        array_pos_test.append(momentum_positron_crosscheck)
        array_neutron.append(momentum_neutron)

        """ append information to hepevt_file: """
        # 2 particles in the event:
        hepevt_file.append("2")
        # append positron information:
        hepevt_file.append("1\t{0} 0 0 {1} {2} {3} {4}".format(pdg_positron, momentum_positron_x, momentum_positron_y,
                                                               momentum_positron_z, mass_positron))
        # append neutron information to file:
        hepevt_file.append("1\t{0} 0 0 {1} {2} {3} {4}".format(pdg_neutron, momentum_neutron_x, momentum_neutron_y,
                                                               momentum_neutron_z, mass_neutron))

    # open file:
    MyFile = open(output_path + "IBD_hepevt_{1:d}events_file{0:d}.txt".format(filenumber, number_events_per_file), 'w')

    for element in hepevt_file:
        print >>MyFile, element
    MyFile.close()

bins = np.arange(0.005, 0.120, 0.001)

h1 = plt.figure(1)
plt.hist(array_pos, bins, labeL='positron momentum in GeV')
plt.legend()


h2 = plt.figure(2)
plt.hist(array_pos_test, bins, labeL='cross-check positron momentum in GeV')
plt.legend()


h3 = plt.figure(3)
plt.hist(array_neutron, np.arange(0.0, 0.100, 0.001), labeL='neutron momentum in GeV')
plt.legend()

plt.show()











