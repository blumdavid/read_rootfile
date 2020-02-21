""" compare file and event number of file{}_evt{}_prompt_signal.txt with evtID_IBDlike_{}.txt and Evis_file{}.txt:

    The problem is, that in script OLD_IBDlike_events_channels_v1.py 9001 IBD-like events are anaylzed, in script
    OLD_atmoNC_spectrum_v1.py only 8946 events.

"""

import ROOT
import numpy as np
from matplotlib import pyplot as plt
import NC_background_functions
import datetime
import os
import re
import sys


def get_numbers_from_filename(filename):
    """
    function to get number as integer out of a string. For example: filename='file235' -> num = 235 of type integer

    :param filename: string of one part of the filename 'file{}_evt{}_prompt_signal.txt'
    :return:
    """
    # get the number out of filename and convert it into integer:
    num = int(re.search(r'\d+', filename).group(0))

    return num


# path, where the txt files with the number of pe of prompt signal are saved:
input_path = "/home/astro/blum/juno/atmoNC/data_NC/output_detsim/"

# set the file number of the first file to be analyzed:
first_file = 0
# set the file number of the last file to be analyzed:
last_file = 999

# number of IBD-like events from evtID_IBDlike_{}.txt:
number_from_evtID_IBDlike = 0

# number of IBD-like events from Evis_file{}.txt:
number_from_Evis_file = 0

# number of IBD-like events from file{}_evt{}_prompt_signal.txt:
number_from_file_evt_prompt_signal = 0

""" read file{}_evt{}_prompt_signal.txt files: """
# preallocate array, where the event number is stored (event number from 0 to 99999):
evt_number_from_file_evt_prompt_signal = []

# loop over all files in folder input_path, that start with 'file' and end with '_prompt_signal.txt' (the name
# of these files contains the file number and evtID):
for file_ibdlike in os.listdir(input_path):

    if file_ibdlike.startswith("file") and file_ibdlike.endswith("_prompt_signal.txt"):

        # increment number_from_file_evt_prompt_signal:
        number_from_file_evt_prompt_signal += 1

        # split string file_ibdlike into two parts: 1. part: 'file{}', 2. part: 'vt{}_prompt_signal.txt'
        x = file_ibdlike.split("_e")

        # x[0] is string 'file{}':
        file_string = x[0]
        # x[1] is string 'vt{}_prompt_signal.txt':
        event_string = x[1]

        # get file_number of file_string:
        file_number = get_numbers_from_filename(file_string)
        # get evtID of event_string:
        evtID = get_numbers_from_filename(event_string)

        # calculate the event_num with file_number, evtID:
        event_num = file_number * 100 + evtID

        # append event_num to event_number array:
        evt_number_from_file_evt_prompt_signal.append(event_num)

# sort evt_number_from_file_evt_prompt_signal. Therefore sort array in ascending:
evt_number_from_file_evt_prompt_signal.sort()

""" read evtID_IBDlike_{}.txt files: """
evt_number_from_evtID_IBDlike = []

for index in range(first_file, last_file+1, 1):
    # file name:
    input_file = input_path + "Evis_file{0:d}.txt".format(index)
    # read this file:
    Evis_arr = np.loadtxt(input_file)
    # increment number_from_Evis_file:
    number_from_Evis_file += len(Evis_arr)

    # file name, where evtID of IBD-like events are saved:
    input_file_evtID = input_path + "evtID_IBDlike_{0:d}.txt".format(index)
    # read this file:
    evtID_arr = np.loadtxt(input_file_evtID)
    # increment number_from_evtID_IBDlike:
    number_from_evtID_IBDlike += len(evtID_arr)

    # loop over evtID, calculate the total event number (with filenumber and evtID_arr[]) and append total event
    # number to array:
    for index1 in range(len(evtID_arr)):

        evt_num = index * 100 + evtID_arr[index1]

        evt_number_from_evtID_IBDlike.append(evt_num)

""" compare evt_number_from_file_evt_prompt_signal and evt_number_from_evtID_IBDlike: """
print("number of IBD-like events from evtID_IBDlike_.txt = {0:d}".format(number_from_evtID_IBDlike))
print("number of IBD-like events from Evis_file.txt = {0:d}".format(number_from_Evis_file))
print("number of IBD-like events from file_evt_prompt_signal.txt = {0:d}\n".format(number_from_file_evt_prompt_signal))

print("len(evt_number_from_evtID_IBDlike) = {0:d}".format(len(evt_number_from_evtID_IBDlike)))
print("len(evt_number_from_file_evt_prompt_signal) = {0:d}\n".format(len(evt_number_from_file_evt_prompt_signal)))

# loop over entries of evt_number_from_file_evt_prompt_signal:
for index2 in range(len(evt_number_from_file_evt_prompt_signal)):

    # set flag, that event number is missing in evt_number_from_evtID_IBDlike:
    flag_event_number_missing = True

    # loop over entries of evt_number_from_evtID_IBDlike and check if event number are equal:
    for index3 in range(len(evt_number_from_evtID_IBDlike)):

        if evt_number_from_evtID_IBDlike[index3] == evt_number_from_file_evt_prompt_signal[index2]:
            # set flag:
            flag_event_number_missing = False
            break

    # print event number, if it is missing in evt_number_from_evtID_IBDlike:
    if flag_event_number_missing:
        print("event number missing in evt_number_from_evtID_IBDlike = {0:d}"
              .format(evt_number_from_file_evt_prompt_signal[index2]))

# loop over entries of evt_number_from_evtID_IBDlike:
for index4 in range(len(evt_number_from_evtID_IBDlike)):

    # set flag, that event number is missing in evt_number_from_file_evt_prompt_signal:
    flag_event_number_miss = True

    # loop over entries of evt_number_from_file_evt_prompt_signal and check if event number are equal:
    for index5 in range(len(evt_number_from_file_evt_prompt_signal)):

        if evt_number_from_file_evt_prompt_signal[index5] == evt_number_from_evtID_IBDlike[index4]:
            # set flag:
            flag_event_number_miss = False
            break

    # print event number, if it is missing in evt_number_from_file_evt_prompt_signal:
    if flag_event_number_miss:
        print("event number missing in evt_number_from_file_evt_prompt_signal = {0:d}"
              .format(evt_number_from_evtID_IBDlike[index4]))
