""" script to snip a large hepevt-file into smaller hepevt-files:

    e.g. snip the hepevt file from DSNB-NC generator (gen_NC_onlyC12_250000evts_seed1.txt), that contains 250000 NC
    events, into 250 files with 1000 events each, which then can be used as input for JUNO detector simulation.

"""


def write_new_hep_file(out_file, num, input_arr):
    output_file = open(out_file + str(num) + '.txt', 'w')

    for line in input_arr:

        output_file.write(line)

    output_file.close()

    return


# input path (string):
input_path = "/home/astro/blum/juno/atmoNC/data_NC/output_generator/"

# input file name of the hep file (string):
input_file = input_path + "gen_NC_onlyC12_250000evts_seed1.txt"
# input_file = input_path + "gen_test.txt"

# output path (string):
output_path = input_path + "folder_NC_onlyC12_250000evts_seed1/"

# output file name (string):
output_name = output_path + "out_gen_NC_onlyC12_1000evts_seed1_"


""" read hep file """
# open file:
hep_file = open(input_file, "r")
# preallocate array, where hep events are stored (array):
hep_file_arr = []
# loop over lines of hep-file and append lines to array:
for index1, line1 in enumerate(hep_file):
    hep_file_arr.append(line1)

# close file:
hep_file.close()


""" write snippets of hep file to output file: """
# preallocate number of events per output file:
num_evts_out = 1000
# preallocate number of output files:
number_files = 250
# set start index of the line in hep_file_arr for the new output file:
index_line = 0

# loop over output files:
for index2 in range(number_files):

    # preallocate array, that is written to output file:
    hep_arr_output = []

    # preallocate number of events:
    num_events = 0

    # preallocate index:
    index_line2 = 0

    # loop over hep_file_arr:
    for index3, line2 in enumerate(hep_file_arr[index_line:]):

        # index of the line in hep_file_arr:
        index_line2 = index3

        if len(line2) < 11:
            # this line contains the number of particles in one event:
            num_events = num_events + 1

        if num_events == num_evts_out+1:
            break

        # append line to output array:
        hep_arr_output.append(line2)

    # write hep_arr_output to file:
    write_new_hep_file(output_name, index2, hep_arr_output)

    # add index_line2:
    index_line = index_line + index_line2

    print(index2)

















