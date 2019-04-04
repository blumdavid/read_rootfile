""" script to check, if there can be a delayed signal without a neutron in the events of user_atmoNC_.root:

    To get a delayed signal, there should be at least one neutron in the event. But to check this, around 1000 events
    with no neutron from user_atmoNC_.root must be analyzed and checked for a delayed signal.


"""
import datetime
import NC_background_functions
import numpy as np
from matplotlib import pyplot as plt

# get the date and time, when the script was run:
date = datetime.datetime.now()
now = date.strftime("%Y-%m-%d %H:%M")

# set the path of the input files (filename must be 'user_atmoNC_{}.root'):
input_path = "/local/scratch1/pipc51/astro/blum/detsim_output_data/"

# set path, where results should be saved:
output_path = "/home/astro/blum/juno/atmoNC/data_NC/output_neutron_cut/"

# set the number of the first file and number of the last file that should be read:
start_number = 21
stop_number = 24
# number of entries in the input files:
Number_entries_input = 100

# Set SAVE_HITTIME and SAVE_TXT flag:
SAVE_HITTIME = True
SAVE_TXT = True

# Set minimum and maximum hittime in ns (define time window, where delayed signal should lie):
min_hittime = 10000
max_hittime = 1000000
# Set threshold of number of PE per bin for possible delayed signal (bin-width = 5 ns):
threshold = 30
# set threshold2 of number of PEs per bin (signal peak is summed as long as nPE is above threshold2):
threshold2 = 1
# Set bin-width of hittime histograms in ns:
binwidth = 5.0
# min and max number of PE for delayed energy cut (from check_delayed_energy.py):
min_PE_delayed = 2647.80
max_PE_delayed = 3524.22

# preallocate variables:
# number of total events, that are analyzed:
number_events = 0
# number of events with at least one neutron:
number_neutron = 0
# number of events without neutron:
number_no_neutron = 0
# number of events without neutron and without possible delayed signal:
number_no_delayed = 0
# number of events without neutron, but with delayed signal (agree with time and delayed energy cut):
number_delayed = 0
# number of events without neutron, but with possible delayed signal (agree with time but NOT with energy cut):
number_possible_delayed = 0
# number of events without neutron, but with possible second delayed signal (agree only with time cut) after one
# delayed or possible delayed cut:
number_possible_second_delayed = 0

# loop over files:
for file_number in range(start_number, stop_number+1, 1):
    # path to file:
    input_file = "user_atmoNC_{0:d}.root".format(file_number)
    print("Start reading {0} ...".format(input_file))

    # analyze file with function check_neutron_cut():
    num_events, num_neutron, num_no_neutron, num_no_delayed, num_delayed, num_pos_delayed, num_pos_second_delayed = \
        NC_background_functions.check_neutron_cut(input_path, file_number, output_path, min_hittime, max_hittime,
                                                  threshold, threshold2, binwidth, min_PE_delayed, max_PE_delayed,
                                                  Number_entries_input, SAVE_HITTIME)

    # add variables:
    number_events = number_events + num_events
    number_neutron = number_neutron + num_neutron
    number_no_neutron = number_no_neutron + num_no_neutron
    number_no_delayed = number_no_delayed + num_no_delayed
    number_delayed = number_delayed + num_delayed
    number_possible_delayed = number_possible_delayed + num_pos_delayed
    number_possible_second_delayed = number_possible_second_delayed + num_pos_second_delayed

print("\nnumber_events = {0}".format(number_events))
print("\nnumber_neutron = {0}".format(number_neutron))
print("\nnumber_no_neutron = {0}".format(number_no_neutron))
print("\nnumber_no_delayed = {0}".format(number_no_delayed))
print("\nnumber_delayed = {0}".format(number_delayed))
print("\nnumber_possible_delayed = {0}".format(number_possible_delayed))
print("\nnumber_possible_second_delayed = {0}".format(number_possible_second_delayed))

# calculate neutron cut efficiency in percent (I lost number_delayed possible delayed signals when I only look
# at event with at least 1 neutron):
efficiency = float(number_delayed) / float(number_events) * 100

if SAVE_TXT:
    # save numbers from above to txt file:
    np.savetxt(output_path + "result_neutron_cut_atmoNC_{0:d}_to_{1:d}.txt".format(start_number, stop_number),
               np.array([number_events, number_neutron, number_no_neutron, number_no_delayed, number_delayed,
                         number_possible_delayed, number_possible_second_delayed, efficiency]),
               fmt='%.3f',
               header="Results of script check_neutron_cut.py (only 20inch PMTs) ({0}):\n"
                      "Input-files: user_atmoNC_{1:d}.root to user_atmoNC_{2:d}.root;\n"
                      "Time window: min_hittime = {3:d} ns, max_hittime = {4:d} ns;\n"
                      "Threshold (nPE per bin, bin-width = {6:0.1f} ns) = {5:d};\n"
                      "Threshold2 (nPE per bin, define integration of pulse) = {7:d};\n"
                      "Delayed energy cut: min. PE = {8:0.2f}, max. PE = {9:0.2f}.\n"
                      "Results:\n"
                      "\n"
                      "Total number of analyzed events;\n"
                      "Number of events with at least 1 neutron;\n"
                      "Number of events without a neutron;\n"
                      "Number of events without neutron and without possible delayed signal;\n"
                      "Number of events without neutron, but with delayed signal (agrees with time and delayed energy "
                      "cut);\n"
                      "Number of events without neutron, but with possible delayed signal (agrees with time but NOT "
                      "with energy cut);\n"
                      "Number of events without neutron, but with possible second delayed signal (agree only with time "
                      "cut) after one delayed or possible delayed cut;\n"
                      "Cut efficiency in percent (number_delayed/number_events):"
               .format(now, start_number, stop_number, min_hittime, max_hittime, threshold, binwidth, threshold2,
                       min_PE_delayed, max_PE_delayed))
