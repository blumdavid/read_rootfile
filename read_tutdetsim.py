""" Script to read the sample_detsim_user.root file from JUNO offline detector simulation (tut_detsim.py) and to
    calculate the visible spectrum of events that mimic IBD events:

    - file sample_detsim_user.root is generated with the JUNO detector simulation:
        python $TUTORIALROOT/share/tut_detsim.py
    - the hepevt file from the DSNB-NC.exe generator is used as input for the detector simulation

    This ROOT file is read and analyzed to get the spectrum of visible energy of the NC QEL events that mimic an IBD
    signal in the JUNO detector.

    Criteria for IBD signal (from JUNO PhysicsReport page 39):
        1.  fiducial volume cut: r < 17 m
        2.  prompt energy cut: 10 MeV < E_prompt < 105 MeV
        3.  delayed energy cut: 1.9 MeV < E_delayed < 2.5 MeV (neutron is captured on Carbon and releases gamma of
            around 2.2 MeV)
        4.  time cut: time interval between prompt and delayed signal: 600 ns < delta_T < 1.0 ms
            (Julia has taken 600 ns, in PhysicsReport there is no further specification)
        5.  distance cut between prompt and delayed signal: R_(prompt-delayed) < 1.5 m
        6.  neutron multiplicity cut: only 1 delayed signal (neutron capture with 1.9 MeV to 2.5 MeV) in the time window
            (600 ns to 1.0 ms)


        6. Muon veto criteria: -> Do I have to consider this??????????????


"""

# import ROOT
import datetime
# import glob
import NC_background_functions
import numpy as np
from matplotlib import pyplot as plt

# get the date and time, when the script was run:
date = datetime.datetime.now()
now = date.strftime("%Y-%m-%d %H:%M")

# set SAVE_TXT, defines if txt files are saved:
SAVE_TXT = True

""" set the different cut properties: """
# fiducial volume cut in mm:
R_cut = 17000
# prompt energy cut in MeV:
E_prompt_min = 10.0
E_prompt_max = 105.0
# delayed energy cut in MeV:
E_delayed_min = 1.9
E_delayed_max = 2.5
# time cut between prompt and delayed signal in ns:
time_cut_min = 600
time_cut_max = 1000000
# distance cut between prompt and delayed signal in mm:
Distance_cut = 1500
# time in ns, where two prompt signals cannot be separated anymore (very conservative limit 1000 ns = 1 microsecond):
# TODO-me: what is the time resolution where you can separate 2 prompt signals?
time_resolution = 10

""" set the number of the first file and number of the last file that should be read: """
start_number = 22
stop_number = 99

# set the path of the inputs:
input_path = "/home/astro/blum/juno/atmoNC/data_NC/detsim_output_data/"

# preallocate array, where the visible energy in MeV is stored:
E_visible = np.array([])
# preallocate number of simulated events:
number_events = 0
# preallocate number of IBD-like events:
number_IBDevts = 0

# loop over the files that are read:
for index in range(start_number, stop_number+1):

    # file name of the input file:
    input_name = input_path + "user_atmoNC_{0:d}.root".format(index)
    print("------------------------------------------------------------------------------------")
    print(input_name)


    # get the visible energy of the prompt signal from events that mimic IBD signals (E_vis in MeV) (np.array):
    num_evts, evt_ID_IBD, E_vis = NC_background_functions.read_sample_detsim_user(input_name, R_cut, E_prompt_min,
                                                                                  E_prompt_max, E_delayed_min,
                                                                                  E_delayed_max, time_cut_min,
                                                                                  time_cut_max, Distance_cut,
                                                                                  time_resolution)

    # add number_evts to number_events:
    number_events = number_events + num_evts

    # calculate number of IBD-like events:
    num_IBD = len(E_vis)
    # add num_IBD to number_IBDevts:
    number_IBDevts = number_IBDevts + num_IBD

    # append E_vis to E_visible:
    E_visible = np.append(E_visible, E_vis)

    print(evt_ID_IBD)
    print(E_vis)
    print(num_evts)
    print(num_IBD)


# number of simulated events:
print("Total number of simulated events = {0:d}".format(number_events))
# number of IBD-like events:
print("Number of IBD-like events = {0:d}".format(number_IBDevts))
