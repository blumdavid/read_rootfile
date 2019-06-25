""" script to check the number of entries in root-files (to check, if simulation worked correctly): """
import ROOT

# path to the root files:
input_path = "/local/scratch1/pipc51/astro/blum/detsim_output_data/"
# number of entries in each file:
number_entries = 100
# index of first file:
start_file = 900
# index of last file:
stop_file = 999

# loop over every file:
for index in range(start_file, stop_file+1):
    # file name:
    file_name = "user_atmoNC_{0:d}.root".format(index)

    # input name:
    input_name = input_path + file_name

    # load ROOT file:
    rfile = ROOT.TFile(input_name)
    # get the "evt"-TTree from the TFile:
    rtree_evt = rfile.Get("evt")
    # get the number of events in the geninfo Tree:
    number_events = rtree_evt.GetEntries()

    if number_events == number_entries:
        continue
    else:
        print("\nnumber of events ({0:d}) != {1:d}".format(number_events, number_entries))
        print("failed file: {0}".format(file_name))
