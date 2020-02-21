""" Script to check conversion from nPE to MeV of neutrons and protons, respectively, which were simulated
    with tut_detsim.py of JUNO offline version J18v1r1-pre1.

    Results of this script are used to convert neutron/proton/positron with specific energy in MeV to number of PE
    in the JUNO detector.

    With this conversion the cut on the energy of a possible prompt signal can be made in the PE-regime and efficiency
    of this cut can be calculated.

    More information: info_conversion_proton_neutron.odt (/home/astro/blum/juno/atmoNC/data_NC/conversion_nPE_MeV/)
"""
import datetime
import NC_background_functions
import numpy as np
from matplotlib import pyplot as plt
from decimal import Decimal
from matplotlib.colors import LogNorm

""" define gaussian function: """


def gaussian(x, a, b, c):
    return a * np.exp(- (x-b)**2 / (2*c**2))


def get_info_from_file(start, stop, filename, num_entries, radius_cut):
    """

    :param start: number of first file
    :param stop: number of last file
    :param filename: path and name of the file
    :param num_entries: number of entries per file
    :param radius_cut: radius, that define the volume cut, in mm
    :return:
    """
    """ preallocate arrays: """
    # number of PE of each event:
    number_pe = np.array([])
    # initial total momentum of each event in MeV:
    momentum_init = np.array([])
    # deposit energy in each event in MeV:
    edep = np.array([])
    # quenched deposit energy in each event in MeV:
    qedep = np.array([])

    # loop over files of proton = 10 MeV:
    for num in range(start, stop+1, 1):
        # path to file:
        input_file = filename + "_{0:d}.root".format(num)

        # get number of PE per event (array of int), hit-times of the last event in ns (array of float),
        # initial momentum per event in MeV (array of float), deposit energy per event in MeV and quenched deposit
        # energy per event in MeV:
        num_pe, momentum, e, qe = NC_background_functions.conversion_npe_mev(input_file, num_entries, radius_cut)

        # append arrays to array:
        number_pe = np.append(number_pe, num_pe)
        momentum_init = np.append(momentum_init, momentum)
        edep = np.append(edep, e)
        qedep = np.append(qedep, qe)

    return number_pe, momentum_init, edep, qedep


def save_array_to_file(arr, out_path, file_name, number_events):
    """
    function to save an array (either number_pe or qedep) to txt file to save time, because you must read the file only
    once
    :param arr: array that should be saved (array of float)
    :param out_path: path, where the txt file should be saved (string)
    :param file_name: file name of the txt file (string)
    :param number_events: number of events in the array/root-file
    :return:
    """
    np.savetxt(out_path + file_name + ".txt", arr, fmt='%1.5f',
               header="{0} of {1:d} events analyzed with function get_info_from_file() in script\n "
                      "check_conversion_npe_mev.py (number of photo-electron per event OR \n"
                      "quenched deposited energy/ visible energy per event in MeV).\n"
                      "({2})\n"
                      "(volume cut (R <= {3:d} mm) applied on initial position):"
               .format(file_name, number_events, now, r_cut))

    return


# get the date and time, when the script was run:
date = datetime.datetime.now()
now = date.strftime("%Y-%m-%d %H:%M")

# set the path of the input files:
input_path = "/local/scratch1/pipc51/astro/blum/conversion_nPE_MeV/"
input_proton = input_path + "proton_output/"
input_neutron = input_path + "neutron_output/"

# set path, where results should be saved:
output_path = "/home/astro/blum/juno/atmoNC/data_NC/conversion_nPE_MeV/"

# set the number of the first file and number of the last file that should be read:
start_number = 0
stop_number_p = 99
stop_number_n = 99
# number of entries in the input files:
Number_entries_input = 10
# total number of events:
number_events_p = (stop_number_p - start_number + 1) * Number_entries_input
number_events_n = (stop_number_n - start_number + 1) * Number_entries_input

# set the radius for the volume cut in mm:
r_cut = 16000

# set maximum visible energy for plots and fit in MeV:
max_evis = 120.0

# Set boolean variables:
PLOT_INITENERGY = False
READ_P_10MEV = True
READ_N_10MEV = True
READ_P_100MEV = True
READ_N_100MEV = True
READ_N_300MEV = True
READ_N_500MEV = True
READ_N_500MEV_2 = True
READ_N_500MEV_3 = True
READ_N_500MEV_4 = True
READ_N_500MEV_5 = True
READ_N_500MEV_6 = True
READ_N_500MEV_7 = True
READ_N_500MEV_8 = True
READ_P_1GEV = True
READ_N_1GEV = True

""" 10 MeV proton: """
if READ_P_10MEV:
    print("\nstart reading 10 MeV proton files...")
    # file name:
    file_p_10MeV = input_proton + "user_proton_10_MeV"
    # read info of all files of 10 MeV protons:
    number_pe_p_10MeV, momentum_init_p_10MeV, edep_p_10MeV, qedep_p_10MeV = \
        get_info_from_file(start_number, stop_number_p, file_p_10MeV, Number_entries_input, r_cut)
    # save number of pe to txt file:
    save_array_to_file(number_pe_p_10MeV, output_path, "number_pe_p_10MeV", number_events_p)
    # save qedep to txt file:
    save_array_to_file(qedep_p_10MeV, output_path, "qedep_p_10MeV", number_events_p)
else:
    # load number of pe and qedep array from txt file:
    number_pe_p_10MeV = np.loadtxt(output_path + "number_pe_p_10MeV.txt")
    qedep_p_10MeV = np.loadtxt(output_path + "qedep_p_10MeV.txt")

""" 10 MeV neutron: """
if READ_N_10MEV:
    print("\nstart reading 10 MeV neutron files...")
    # file name:
    file_n_10MeV = input_neutron + "user_neutron_10_MeV"
    # read info of all files of 10 MeV neutrons:
    number_pe_n_10MeV, momentum_init_n_10MeV, edep_n_10MeV, qedep_n_10MeV = \
        get_info_from_file(start_number, stop_number_n, file_n_10MeV, Number_entries_input, r_cut)
    # save number of pe to txt file:
    save_array_to_file(number_pe_n_10MeV, output_path, "number_pe_n_10MeV", number_events_n)
    # save qedep to txt file:
    save_array_to_file(qedep_n_10MeV, output_path, "qedep_n_10MeV", number_events_n)
else:
    # load number of pe and qedep array from txt file:
    number_pe_n_10MeV = np.loadtxt(output_path + "number_pe_n_10MeV.txt")
    qedep_n_10MeV = np.loadtxt(output_path + "qedep_n_10MeV.txt")

""" 100 MeV proton: """
if READ_P_100MEV:
    print("\nstart reading 100 MeV proton files...")
    # file name:
    file_p_100MeV = input_proton + "user_proton_100_MeV"
    # read info of all files of 100 MeV protons:
    number_pe_p_100MeV, momentum_init_p_100MeV, edep_p_100MeV, qedep_p_100MeV = \
        get_info_from_file(start_number, stop_number_p, file_p_100MeV, Number_entries_input, r_cut)
    # save number of pe to txt file:
    save_array_to_file(number_pe_p_100MeV, output_path, "number_pe_p_100MeV", number_events_p)
    # save qedep to txt file:
    save_array_to_file(qedep_p_100MeV, output_path, "qedep_p_100MeV", number_events_p)
else:
    # load number of pe and qedep array from txt file:
    number_pe_p_100MeV = np.loadtxt(output_path + "number_pe_p_100MeV.txt")
    qedep_p_100MeV = np.loadtxt(output_path + "qedep_p_100MeV.txt")

""" 100 MeV neutron: """
if READ_N_100MEV:
    print("\nstart reading 100 MeV neutron files...")
    # file name:
    file_n_100MeV = input_neutron + "user_neutron_100_MeV"
    # read info of all files of 100 MeV neutrons:
    number_pe_n_100MeV, momentum_init_n_100MeV, edep_n_100MeV, qedep_n_100MeV = \
        get_info_from_file(start_number, stop_number_n, file_n_100MeV, Number_entries_input, r_cut)
    # save number of pe to txt file:
    save_array_to_file(number_pe_n_100MeV, output_path, "number_pe_n_100MeV", number_events_n)
    # save qedep to txt file:
    save_array_to_file(qedep_n_100MeV, output_path, "qedep_n_100MeV", number_events_n)
else:
    # load number of pe and qedep array from txt file:
    number_pe_n_100MeV = np.loadtxt(output_path + "number_pe_n_100MeV.txt")
    qedep_n_100MeV = np.loadtxt(output_path + "qedep_n_100MeV.txt")

""" 300 MeV neutron: """
if READ_N_300MEV:
    print("\nstart reading 300 MeV neutron files...")
    # file name:
    file_n_300MeV = input_neutron + "user_neutron_300_MeV"
    # read info of all files of 300 MeV neutrons:
    number_pe_n_300MeV, momentum_init_n_300MeV, edep_n_300MeV, qedep_n_300MeV = \
        get_info_from_file(start_number, stop_number_n, file_n_300MeV, Number_entries_input, r_cut)
    # save number of pe to txt file:
    save_array_to_file(number_pe_n_300MeV, output_path, "number_pe_n_300MeV", number_events_n)
    # save qedep to txt file:
    save_array_to_file(qedep_n_300MeV, output_path, "qedep_n_300MeV", number_events_n)
else:
    # load number of pe and qedep array from txt file:
    number_pe_n_300MeV = np.loadtxt(output_path + "number_pe_n_300MeV.txt")
    qedep_n_300MeV = np.loadtxt(output_path + "qedep_n_300MeV.txt")

""" 500 MeV neutron (user_neutron_500_MeV_0.root to user_neutron_500_MeV_99.root): """
if READ_N_500MEV:
    print("\nstart reading 500 MeV neutron files...")
    # file name:
    file_n_500MeV = input_neutron + "user_neutron_500_MeV"
    # read info of all files of 500 MeV neutrons:
    number_pe_n_500MeV, momentum_init_n_500MeV, edep_n_500MeV, qedep_n_500MeV = \
        get_info_from_file(start_number, stop_number_n, file_n_500MeV, Number_entries_input, r_cut)
    # save number of pe to txt file:
    save_array_to_file(number_pe_n_500MeV, output_path, "number_pe_n_500MeV", number_events_n)
    # save qedep to txt file:
    save_array_to_file(qedep_n_500MeV, output_path, "qedep_n_500MeV", number_events_n)
else:
    # load number of pe and qedep array from txt file:
    number_pe_n_500MeV = np.loadtxt(output_path + "number_pe_n_500MeV.txt")
    qedep_n_500MeV = np.loadtxt(output_path + "qedep_n_500MeV.txt")

""" 500 MeV neutron (user_neutron_500_MeV_100.root to user_neutron_500_MeV_199.root): """
if READ_N_500MEV_2:
    print("\nstart reading 500 MeV neutron files 2...")
    # file name:
    file_n_500MeV_2 = input_neutron + "user_neutron_500_MeV"
    # read info of all files of 500 MeV neutrons:
    number_pe_n_500MeV_2, momentum_init_n_500MeV_2, edep_n_500MeV_2, qedep_n_500MeV_2 = \
        get_info_from_file(start_number+100, stop_number_n+100, file_n_500MeV_2, Number_entries_input, r_cut)
    # save number of pe to txt file:
    save_array_to_file(number_pe_n_500MeV_2, output_path, "number_pe_n_500MeV_2", number_events_n)
    # save qedep to txt file:
    save_array_to_file(qedep_n_500MeV_2, output_path, "qedep_n_500MeV_2", number_events_n)
else:
    # load number of pe and qedep array from txt file:
    number_pe_n_500MeV_2 = np.loadtxt(output_path + "number_pe_n_500MeV_2.txt")
    qedep_n_500MeV_2 = np.loadtxt(output_path + "qedep_n_500MeV_2.txt")

""" 500 MeV neutron (user_neutron_500_MeV_200.root to user_neutron_500_MeV_299.root): """
if READ_N_500MEV_3:
    print("\nstart reading 500 MeV neutron files 3...")
    # file name:
    file_n_500MeV_3 = input_neutron + "user_neutron_500_MeV"
    # read info of all files of 500 MeV neutrons:
    number_pe_n_500MeV_3, momentum_init_n_500MeV_3, edep_n_500MeV_3, qedep_n_500MeV_3 = \
        get_info_from_file(start_number+200, stop_number_n+200, file_n_500MeV_3, Number_entries_input, r_cut)
    # save number of pe to txt file:
    save_array_to_file(number_pe_n_500MeV_3, output_path, "number_pe_n_500MeV_3", number_events_n)
    # save qedep to txt file:
    save_array_to_file(qedep_n_500MeV_3, output_path, "qedep_n_500MeV_3", number_events_n)
else:
    # load number of pe and qedep array from txt file:
    number_pe_n_500MeV_3 = np.loadtxt(output_path + "number_pe_n_500MeV_3.txt")
    qedep_n_500MeV_3 = np.loadtxt(output_path + "qedep_n_500MeV_3.txt")

""" 500 MeV neutron (user_neutron_500_MeV_300.root to user_neutron_500_MeV_399.root): """
if READ_N_500MEV_4:
    print("\nstart reading 500 MeV neutron files 4...")
    # file name:
    file_n_500MeV_4 = input_neutron + "user_neutron_500_MeV"
    # read info of all files of 500 MeV neutrons:
    number_pe_n_500MeV_4, momentum_init_n_500MeV_4, edep_n_500MeV_4, qedep_n_500MeV_4 = \
        get_info_from_file(start_number+300, stop_number_n+300, file_n_500MeV_4, Number_entries_input, r_cut)
    # save number of pe to txt file:
    save_array_to_file(number_pe_n_500MeV_4, output_path, "number_pe_n_500MeV_4", number_events_n)
    # save qedep to txt file:
    save_array_to_file(qedep_n_500MeV_4, output_path, "qedep_n_500MeV_4", number_events_n)
else:
    # load number of pe and qedep array from txt file:
    number_pe_n_500MeV_4 = np.loadtxt(output_path + "number_pe_n_500MeV_4.txt")
    qedep_n_500MeV_4 = np.loadtxt(output_path + "qedep_n_500MeV_4.txt")

""" 500 MeV neutron (user_neutron_500_MeV_400.root to user_neutron_500_MeV_499.root): """
if READ_N_500MEV_5:
    print("\nstart reading 500 MeV neutron files 5...")
    # file name:
    file_n_500MeV_5 = input_neutron + "user_neutron_500_MeV"
    # read info of all files of 500 MeV neutrons:
    number_pe_n_500MeV_5, momentum_init_n_500MeV_5, edep_n_500MeV_5, qedep_n_500MeV_5 = \
        get_info_from_file(start_number+400, stop_number_n+400, file_n_500MeV_5, Number_entries_input, r_cut)
    # save number of pe to txt file:
    save_array_to_file(number_pe_n_500MeV_5, output_path, "number_pe_n_500MeV_5", number_events_n)
    # save qedep to txt file:
    save_array_to_file(qedep_n_500MeV_5, output_path, "qedep_n_500MeV_5", number_events_n)
else:
    # load number of pe and qedep array from txt file:
    number_pe_n_500MeV_5 = np.loadtxt(output_path + "number_pe_n_500MeV_5.txt")
    qedep_n_500MeV_5 = np.loadtxt(output_path + "qedep_n_500MeV_5.txt")

""" 500 MeV neutron (user_neutron_500_MeV_500.root to user_neutron_500_MeV_599.root): """
if READ_N_500MEV_6:
    print("\nstart reading 500 MeV neutron files 6...")
    # file name:
    file_n_500MeV_6 = input_neutron + "user_neutron_500_MeV"
    # read info of all files of 500 MeV neutrons:
    number_pe_n_500MeV_6, momentum_init_n_500MeV_6, edep_n_500MeV_6, qedep_n_500MeV_6 = \
        get_info_from_file(start_number+500, stop_number_n+500, file_n_500MeV_6, Number_entries_input, r_cut)
    # save number of pe to txt file:
    save_array_to_file(number_pe_n_500MeV_6, output_path, "number_pe_n_500MeV_6", number_events_n)
    # save qedep to txt file:
    save_array_to_file(qedep_n_500MeV_6, output_path, "qedep_n_500MeV_6", number_events_n)
else:
    # load number of pe and qedep array from txt file:
    number_pe_n_500MeV_6 = np.loadtxt(output_path + "number_pe_n_500MeV_6.txt")
    qedep_n_500MeV_6 = np.loadtxt(output_path + "qedep_n_500MeV_6.txt")

""" 500 MeV neutron (user_neutron_500_MeV_600.root to user_neutron_500_MeV_699.root): """
if READ_N_500MEV_7:
    print("\nstart reading 500 MeV neutron files 7...")
    # file name:
    file_n_500MeV_7 = input_neutron + "user_neutron_500_MeV"
    # read info of all files of 500 MeV neutrons:
    number_pe_n_500MeV_7, momentum_init_n_500MeV_7, edep_n_500MeV_7, qedep_n_500MeV_7 = \
        get_info_from_file(start_number+600, stop_number_n+600, file_n_500MeV_7, Number_entries_input, r_cut)
    # save number of pe to txt file:
    save_array_to_file(number_pe_n_500MeV_7, output_path, "number_pe_n_500MeV_7", number_events_n)
    # save qedep to txt file:
    save_array_to_file(qedep_n_500MeV_7, output_path, "qedep_n_500MeV_7", number_events_n)
else:
    # load number of pe and qedep array from txt file:
    number_pe_n_500MeV_7 = np.loadtxt(output_path + "number_pe_n_500MeV_7.txt")
    qedep_n_500MeV_7 = np.loadtxt(output_path + "qedep_n_500MeV_7.txt")

# """ 500 MeV neutron (user_neutron_500_MeV_700.root to user_neutron_500_MeV_799.root): """
# if READ_N_500MEV_8:
#     print("\nstart reading 500 MeV neutron files 8...")
#     # file name:
#     file_n_500MeV_8 = input_neutron + "user_neutron_500_MeV"
#     # read info of all files of 500 MeV neutrons:
#     number_pe_n_500MeV_8, momentum_init_n_500MeV_8, edep_n_500MeV_8, qedep_n_500MeV_8 = \
#         get_info_from_file(start_number+700, stop_number_n+700, file_n_500MeV_8, Number_entries_input, r_cut)
#     # save number of pe to txt file:
#     save_array_to_file(number_pe_n_500MeV_8, output_path, "number_pe_n_500MeV_8", number_events_n)
#     # save qedep to txt file:
#     save_array_to_file(qedep_n_500MeV_8, output_path, "qedep_n_500MeV_8", number_events_n)
# else:
#     # load number of pe and qedep array from txt file:
#     number_pe_n_500MeV_8 = np.loadtxt(output_path + "number_pe_n_500MeV_8.txt")
#     qedep_n_500MeV_8 = np.loadtxt(output_path + "qedep_n_500MeV_8.txt")

""" 1 GeV proton: """
if READ_P_1GEV:
    print("\nstart reading 1 GeV proton files...")
    # file name:
    file_p_1GeV = input_proton + "user_proton_1000_MeV"
    # read info of all files of 1 GeV protons:
    number_pe_p_1GeV, momentum_init_p_1GeV, edep_p_1GeV, qedep_p_1GeV = \
        get_info_from_file(start_number, stop_number_p, file_p_1GeV, Number_entries_input, r_cut)
    # save number of pe to txt file:
    save_array_to_file(number_pe_p_1GeV, output_path, "number_pe_p_1GeV", number_events_p)
    # save qedep to txt file:
    save_array_to_file(qedep_p_1GeV, output_path, "qedep_p_1GeV", number_events_p)
else:
    # load number of pe and qedep array from txt file:
    number_pe_p_1GeV = np.loadtxt(output_path + "number_pe_p_1GeV.txt")
    qedep_p_1GeV = np.loadtxt(output_path + "qedep_p_1GeV.txt")

""" 1 GeV neutron: """
if READ_N_1GEV:
    print("\nstart reading 1 GeV neutron files...")
    # file name:
    file_n_1GeV = input_neutron + "user_neutron_1000_MeV"
    # read info of all files of 1 GeV neutrons:
    number_pe_n_1GeV, momentum_init_n_1GeV, edep_n_1GeV, qedep_n_1GeV = \
        get_info_from_file(start_number, stop_number_n, file_n_1GeV, Number_entries_input, r_cut)
    # save number of pe to txt file:
    save_array_to_file(number_pe_n_1GeV, output_path, "number_pe_n_1GeV", number_events_n)
    # save qedep to txt file:
    save_array_to_file(qedep_n_1GeV, output_path, "qedep_n_1GeV", number_events_n)
else:
    # load number of pe and qedep array from txt file:
    number_pe_n_1GeV = np.loadtxt(output_path + "number_pe_n_1GeV.txt")
    qedep_n_1GeV = np.loadtxt(output_path + "qedep_n_1GeV.txt")

""" linear fit to qedep vs. nPE diagram: """
# build one array for qedep:
qedep_total = np.concatenate((qedep_p_10MeV, qedep_n_10MeV, qedep_p_100MeV, qedep_n_100MeV, qedep_n_300MeV,
                              qedep_n_500MeV, qedep_n_500MeV_2, qedep_n_500MeV_3, qedep_n_500MeV_4, qedep_n_500MeV_5,
                              qedep_n_500MeV_6, qedep_n_500MeV_7, qedep_p_1GeV, qedep_n_1GeV))
# build one array for number of p.e.:
number_pe_total = np.concatenate((number_pe_p_10MeV, number_pe_n_10MeV, number_pe_p_100MeV, number_pe_n_100MeV,
                                  number_pe_n_300MeV, number_pe_n_500MeV, number_pe_n_500MeV_2, number_pe_n_500MeV_3,
                                  number_pe_n_500MeV_4, number_pe_n_500MeV_5, number_pe_n_500MeV_6,
                                  number_pe_n_500MeV_7, number_pe_p_1GeV,
                                  number_pe_n_1GeV))

""" take only values for qedep below max_evis: """
# preallocate arrays:
qedep_total_interesting = np.array([])
number_pe_total_interesting = np.array([])
# loop over qedep_total:
for index in range(len(qedep_total)):
    if qedep_total[index] <= max_evis:
        qedep_total_interesting = np.append(qedep_total_interesting, qedep_total[index])
        number_pe_total_interesting = np.append(number_pe_total_interesting, number_pe_total[index])

""" do linear fit """
# do linear fit with np.linalg.lstsq:
# The model is y = a * x; x = number_pe_total_interesting, y = qedep_total_interesting
# x needs to be a column vector instead of a 1D vector for this, however.
number_pe_total_interesting_columnvector = number_pe_total_interesting[:, np.newaxis]
# first value of output is slope of linear fit (fir_result is array with one entry):
fit_result = np.linalg.lstsq(number_pe_total_interesting_columnvector, qedep_total_interesting, rcond=None)[0]
# take first entry of fit_result:
fit_result = fit_result[0]
# set x axis for linear fit:
fit_x_axis = np.arange(0, max(number_pe_total_interesting), 100)
# set y axis for linear fit:
fit_y_axis = fit_result * fit_x_axis

print("n_500MeV = {0:d}".format(len(number_pe_n_500MeV) + len(number_pe_n_500MeV_2) +
                                len(number_pe_n_500MeV_3) + len(number_pe_n_500MeV_4) + len(number_pe_n_500MeV_5) +
                                len(number_pe_n_500MeV_6) + len(number_pe_n_500MeV_7)))

""" plot Qedep as function of nPE for all energies: """
h1 = plt.figure(1, figsize=(15, 8))
num_proton = len(number_pe_p_10MeV) + len(number_pe_p_100MeV) + len(number_pe_p_1GeV)
num_neutron = len(number_pe_n_10MeV) + len(number_pe_n_100MeV) + len(number_pe_n_300MeV) + len(number_pe_n_1GeV) + \
              len(number_pe_n_500MeV) + len(number_pe_n_500MeV_2) + len(number_pe_n_500MeV_3) + \
              len(number_pe_n_500MeV_4) + len(number_pe_n_500MeV_5) + len(number_pe_n_500MeV_6) + \
              len(number_pe_n_500MeV_7)
plt.plot(number_pe_p_10MeV, qedep_p_10MeV, "rx", label="proton ({0:d} entries)".format(num_proton))
plt.plot(number_pe_n_10MeV, qedep_n_10MeV, "bx", label="neutron ({0:d} entries)".format(num_neutron))
plt.plot(number_pe_p_100MeV, qedep_p_100MeV, "rx")
plt.plot(number_pe_n_100MeV, qedep_n_100MeV, "bx")
plt.plot(number_pe_n_300MeV, qedep_n_300MeV, "bx")
plt.plot(number_pe_n_500MeV, qedep_n_500MeV, "bx")
plt.plot(number_pe_n_500MeV_2, qedep_n_500MeV_2, "bx")
plt.plot(number_pe_n_500MeV_3, qedep_n_500MeV_3, "bx")
plt.plot(number_pe_n_500MeV_4, qedep_n_500MeV_4, "bx")
plt.plot(number_pe_n_500MeV_5, qedep_n_500MeV_5, "bx")
plt.plot(number_pe_n_500MeV_6, qedep_n_500MeV_6, "bx")
plt.plot(number_pe_n_500MeV_7, qedep_n_500MeV_7, "bx")
# plt.plot(number_pe_n_500MeV_8, qedep_n_500MeV_8, "bx")
plt.plot(number_pe_p_1GeV, qedep_p_1GeV, "rx")
plt.plot(number_pe_n_1GeV, qedep_n_1GeV, "bx")
plt.xlabel("number of p.e.")
plt.ylabel("visible energy in JUNO detector (in MeV)")
plt.title("Visible energy vs. number of p.e.")
plt.legend()
plt.grid()
plt.savefig(output_path + "qedep_vs_nPE_all_energies.png")

""" plot Qedep as function of nPE for qedep <= max_evis: """
h3 = plt.figure(3, figsize=(15, 8))
plt.plot(number_pe_total_interesting, qedep_total_interesting, "rx",
         label="{0:d} entries".format(len(number_pe_total_interesting)))
plt.xlabel("number of p.e.")
plt.ylabel("visible energy in JUNO detector (in MeV)")
plt.title("Visible energy vs. number of p.e.")
plt.legend()
plt.grid()
plt.savefig(output_path + "qedep_vs_nPE_interesting.png")

""" plot Qedep as function of nPE with fit for qedep <= max_evis: """
h4 = plt.figure(4, figsize=(15, 8))
plt.plot(number_pe_total_interesting, qedep_total_interesting, "rx",
         label="{0:d} entries".format(len(number_pe_total_interesting)))
plt.plot(fit_x_axis, fit_y_axis, "b", label="linear fit: f(x) = {0:.3E} * x"
         .format(fit_result))
plt.xlabel("number of p.e.")
plt.ylabel("visible energy in JUNO detector (in MeV)")
plt.title("Visible energy vs. number of p.e.\n(with linear fit)")
plt.legend()
plt.grid()
plt.savefig(output_path + "fit_qedep_vs_nPE_interesting.png")

""" display Qedep as function of nPE in 2D histogram for qedep <= max_evis: """
h5 = plt.figure(5, figsize=(15, 8))
bins_edges_nPE = np.arange(0, max(number_pe_total_interesting), 2000)
bins_edges_Qedep = np.arange(0, max_evis+2, 2)
plt.hist2d(number_pe_total_interesting, qedep_total_interesting, [bins_edges_nPE, bins_edges_Qedep], norm=LogNorm(),
           cmap="rainbow")
plt.xlabel("number of p.e.")
plt.ylabel("visible energy in JUNO detector (in MeV)")
plt.title("Visible energy vs. number of p.e.")
plt.colorbar()
plt.legend()
plt.grid()
plt.savefig(output_path + "hist2d_Qedep_vs_nPE_interesting.png")

""" display Qedep as function of nPE in 2D histogram for qedep <= max_evis with fit: """
h6 = plt.figure(6, figsize=(15, 8))
bins_edges_nPE = np.arange(0, max(number_pe_total_interesting), 2000)
bins_edges_Qedep = np.arange(0, max_evis+2, 2)
plt.hist2d(number_pe_total_interesting, qedep_total_interesting, [bins_edges_nPE, bins_edges_Qedep], norm=LogNorm(),
           cmap="rainbow")
plt.plot(fit_x_axis, fit_y_axis, "k", label="{1:d} entries\nlinear fit: f(x) = {0:.3E} * x"
         .format(fit_result, len(number_pe_total_interesting)))
plt.xlabel("number of p.e.")
plt.ylabel("visible energy in JUNO detector (in MeV)")
plt.title("Visible energy vs. number of p.e.\nwith linear fit")
plt.colorbar()
plt.legend()
plt.grid()
plt.savefig(output_path + "hist2d_Qedep_vs_nPE_interesting_fit.png")

# plot initial energy of proton/neutron:
if PLOT_INITENERGY and READ_P_10MEV and READ_N_10MEV and READ_P_100MEV and READ_N_100MEV and READ_N_300MEV and \
        READ_P_1GEV and READ_N_1GEV:
    h3 = plt.figure(3, figsize=(15, 8))
    bin_width = 0.1
    Bins = np.arange(9.5, 1000.5, bin_width)
    plt.hist(momentum_init_p_10MeV, bins=Bins, color='r', align='mid',
             label="{0:d} protons".format(number_events_p) + " with $E_{kin}$ = 10 MeV")
    plt.hist(momentum_init_n_10MeV, bins=Bins, color='b', align='mid',
             label="{0:d} neutrons".format(number_events_n) + " with $E_{kin}$ = 10 MeV")
    plt.hist(momentum_init_p_100MeV, bins=Bins, color='r', linestyle="--", align='mid',
             label="{0:d} protons".format(number_events_p) + " with $E_{kin}$ = 100 MeV")
    plt.hist(momentum_init_n_100MeV, bins=Bins, color='b', linestyle="--", align='mid',
             label="{0:d} neutrons".format(number_events_n) + " with $E_{kin}$ = 100 MeV")
    plt.hist(momentum_init_n_300MeV, bins=Bins, color='b', linestyle="-.", align='mid',
             label="{0:d} neutrons".format(number_events_n) + " with $E_{kin}$ = 300 MeV")
    plt.hist(momentum_init_p_1GeV, bins=Bins, color='r', linestyle=":", align='mid',
             label="{0:d} protons".format(number_events_p) + " with $E_{kin}$ = 1 GeV")
    plt.hist(momentum_init_n_1GeV, bins=Bins, color='b', linestyle=":", align='mid',
             label="{0:d} neutrons".format(number_events_n) + " with $E_{kin}$ = 1 GeV")
    plt.xlabel("initial kinetic energy in MeV", fontsize=13)
    plt.ylabel("entries per bin (bin-width = {0:0.3f} MeV)".format(bin_width), fontsize=13)
    plt.title("Initial neutron/proton energy", fontsize=18)
    plt.legend()
    plt.grid()
    plt.savefig(output_path + "init_energy.png")

""" Efficiency of the conversion fit: """
# the prompt energy cut is defined by min_ecut and max_ecut in MeV:
min_ecut = 10.0
max_ecut = 100.0
# total number of simulated events:
number_entries = len(number_pe_total)
# preallocate number of 'real' entries inside energy window:
number_entries_real = 0
# preallocate number of 'calculated' entries inside energy window:
number_entries_calculated = 0
# preallocate number of events that are counted too less:
number_too_less = 0
# preallocate number of events that are counted too much:
number_too_much = 0

# calculate Qedep for each entry in number_pe_total with the function of the linear fit:
qedep_calculated = fit_result * number_pe_total

# loop over qedep_total (same like looping over qedep_calculated, since len(qedep_total) == len(qedep_calculated)).
# Therefore check lengths before:
if len(qedep_total) != len(qedep_calculated):
    print("--------------------ERROR: len(qedep_total) != len(qedep_calculated)")

for index in range(len(qedep_total)):
    # preallocate indices to check difference between real and calculated data:
    index_real = 0
    index_calc = 0

    # get the number of entries from the simulated data, where min_ecut <= Qedep <= max_ecut:
    # check min_ecut <= Qedep <= max_ecut:
    if min_ecut <= qedep_total[index] <= max_ecut:
        # entry inside energy window:
        number_entries_real += 1
        # get index:
        index_real = index

    # get the number of entries from the calculated data, where min_ecut <= qedep_calc <= max_ecut:
    # check min_ecut <= Qedep <= max_ecut:
    if min_ecut <= qedep_calculated[index] <= max_ecut:
        # entry inside energy window:
        number_entries_calculated += 1
        # get index:
        index_calc = index

    # check entry, where there is an entry in real data, but not in calculated data (events are counted too much):
    if index_real != 0 and index_calc == 0:
        number_too_less += 1
        print("\nindex_real = {0:d}, index_calc = {1:d}".format(index_real, index_calc))
        print("qedep_total = {0:.2f} MeV".format(qedep_total[index]))
        print("qedep_calculated = {0:.2f} MeV".format(qedep_calculated[index]))

    # check entry, where there is an entry in calculated data, but not in real data (events are counted too less):
    if index_real == 0 and index_calc != 0:
        number_too_much += 1
        print("\nindex_real = {0:d}, index_calc = {1:d}".format(index_real, index_calc))
        print("qedep_total = {0:.2f} MeV".format(qedep_total[index]))
        print("qedep_calculated = {0:.2f} MeV".format(qedep_calculated[index]))

# calculate the efficiency of the prompt energy cut (describes the 'error', when using the conversion from nPE to Qedep)
# in percent:
efficiency_prompt_energy_cut = 100 + float(number_too_less - number_too_much) / float(number_entries_calculated) * 100

print("total number of simulated events = {0:d}\n".format(number_entries))
print("number of 'real' entries from simulated data with {0:.1f} MeV <= Qedep_real <= {1:.1f} MeV: {2:d}\n"
      .format(min_ecut, max_ecut, number_entries_real))
print("number of entries calculated with linear fit with {0:.1f} MeV <= Qedep_calc <= {1:.1f} MeV: {2:d}\n"
      .format(min_ecut, max_ecut, number_entries_calculated))
print("number of events counted too less = {0:d}".format(number_too_less))
print("number of events counted too much = {0:d}".format(number_too_much))
print("efficiency of prompt energy cut = {0:.4f} % (1 + (number_too_less - number_too_much) / number_calculated)"
      .format(efficiency_prompt_energy_cut))



