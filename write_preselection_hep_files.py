""" script to write hep-evt files, which contains the event that pass the preselection of script
    preselection_detsim_user.py.

    1.  The old hep-evt file gen_NC_onlyC12_250000evts_seed1.txt are read (from output_generator/)
    2.  The files evtID_preselected_{}.txt (from /output_preselection/), where the event ID's of the preselected events
        are saved, are read.
    3.  The info about the event of the old hep-evt file for the event ID, that passes preselection is added to a new
        hep-evt file.
    4.  This new hep-evt file is then saved and can be snipped into files with 100 events each with script
        snip_hepevt_file.py

"""


def write_new_hep_file(out_file, input_arr):
    output_file = open(out_file, 'w')

    for line in input_arr:

        output_file.write(line)

    output_file.close()

    return


# input path (string), where old hep-evt files are stored:
# input_path_hep = "/home/astro/blum/juno/atmoNC/data_NC/output_generator/"
input_path_hep = "/home/astro/blum/juno/atmoNC/data_NC/output_generator/folder_NC_onlyC12_250000evts_seed1_1000/"

# input file name of old hep-evt file (the whole file with 250000 events is read, not the snippets):
# input_name_hep = input_path_hep + "gen_NC_onlyC12_250000evts_seed1.txt"
input_name_hep = input_path_hep + "out_gen_NC_onlyC12_1000evts_seed1_0.txt"

# input path, where evtID_preselected_{}.txt files are stored:
input_path = "/home/astro/blum/juno/atmoNC/data_NC/output_preselection/"

# input file name of preselected evtID's:
input_name = input_path + "evtID_preselected_"

# index of first and last evtID_preselected_{}.txt file:
first_index = 0
last_index = 0

# number of events, that are simulated with old hep-files out_gen_NC_onlyC12_1000evts_seed1.txt:
number_events = 1000
# total number of events:
number_total_events = (last_index + 1) * number_events

# output path (string), where new hep-evt files should be stored:
output_path_hep = input_path_hep

# output file name of new hep-evt file:
output_name = output_path_hep + "NC_onlyC12_250000evts_seed1_preselect_test3_0.txt"

""" read old hep file """
# open file:
hep_file = open(input_name_hep, "r")
# preallocate array, where hep events are stored (array):
hep_file_array = []
# loop over lines of hep-file and append lines to array:
for index1, line1 in enumerate(hep_file):
    hep_file_array.append(line1)

# close file:
hep_file.close()

""" loop over evtID_preselected files and add evtID to evtID_arr: """
# preallocate array, where preselected evtID's are stored:
evtID_arr = []

for index in range(first_index, last_index+1, 1):

    """ read evtID_preselected file: """
    # open file:
    file_evtID = open(input_name + str(index) + ".txt", "r")
    # loop over lines of evtID file and check content:
    for index1, line1 in enumerate(file_evtID):

        if line1[0] != '#':
            # check if first character of line1 is not "#" (this is to exclude the header).
            # remove the last character of line1 (all lines are of type '2\n', therefore remove '\n'
            # (\n is one character))
            line1 = line1[:-1]
            # change string-type to integer:
            line1 = int(line1, base=10)
            # add index*1000 to current evt ID to have the evtID's in increment number and not always from 0:
            line1 = line1 + index*number_events
            # append line1 to array
            evtID_arr.append(line1)

    # close file:
    file_evtID.close()

""" create list, where one array per event is saved: """
# preallocated list where info of preselected events is stored:
hep_file_list_new = []
# start index of line in hep_file_arr:
start_line_index = 0
# start index of evtID_arr:
start_index_evtID = 0

# loop over total number of events:
for index in range(0, number_total_events, 1):

    # preallocate array, where information about one event is stored (array):
    hep_file_arr = []
    # flag for second short line:
    flag_short_line = 0

    # loop over lines of hep-file and append lines to array, then add this array to the list:
    for index1, line1 in enumerate(hep_file_array[start_line_index:]):

        if len(line1) < 11:
            # short line: increment flag_short_line:
            flag_short_line += 1

        if flag_short_line < 2:
            # append line to array:
            hep_file_arr.append(line1)
        else:
            # second short line is reached:
            # add index1 to start_line_index to define the new start index:
            start_line_index = start_line_index + index1
            break

    # now hep_file_arr contains all information about the event specified by index.

    flag_event_preselected = False

    # loop over evtID_arr to check if this event is preselected:
    for index2 in range(start_index_evtID, len(evtID_arr), 1):

        if evtID_arr[index2] == index:
            # event passes preselection!
            # append hep_file_arr to hep_file_list_new:
            hep_file_list_new.append(hep_file_arr)

            flag_event_preselected = True

            # increment start_index_evtID by 1:
            start_index_evtID = start_index_evtID + 1
            break
        else:
            continue

    if not flag_event_preselected:
        # get number of particles in the event (first parameter in array):
        number_particles_in_event = hep_file_arr[0]
        # remove \n from string:
        number_particles_in_event = number_particles_in_event[:-1]
        # change string to integer:
        number_particles_in_event = int(number_particles_in_event, base=10)

        # create new "dummy" hep_array:
        hep_array = [hep_file_arr[0]]

        # loop over number of particles in the event:
        for index in range(number_particles_in_event):
            # append 'dummy' B9 isotope:
            hep_array.append("1\t1000050090 0 0 0.000000 0.000000 0.000000 8.393500\n")

        hep_file_list_new.append(hep_array)


# print length of hep_file_list_new:
print(len(hep_file_list_new))

# hep_file_list_new is a list of lists. To be able to write it to a new file, this list of lists must be flattened to a
# single list:
hep_file_flatlist_new = []
for sublist in hep_file_list_new:
    for item in sublist:
        hep_file_flatlist_new.append(item)

# write flatten list to txt file:
write_new_hep_file(output_name, hep_file_flatlist_new)






















# write_new_hep_file(output_name, hep_file_arr_new)
