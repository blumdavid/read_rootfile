""" Script to read the 'original' ROOT file from Julia's GENIE simulation and convert it to a ROOT-file, which can be
    read from the DSNB-NC.exe generator of the JUNO offline software.

    The ROOT-file, which is generated with this script can be used as input for the DSNB-NC.exe generator.

"""



# import ROOT
import datetime
# import glob
import NC_background_functions
import numpy as np
from matplotlib import pyplot as plt


# set the path of the inputs:
input_path = "/home/astro/blum/juno/atmoNC/data_Julia/"

# file name of the input file:
input_name = input_path + "gntp.101.gst.root"

# set the path, where the outputs are saved:
output_path = "/home/astro/blum/juno/atmoNC/data_NC/"

NC_background_functions.convert_genie_file_for_generator(input_name, output_path)


