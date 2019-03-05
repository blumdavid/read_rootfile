""" Script to read sample_detsim_user.root files (user_atmoNC.root), to make preselection of possible IBD-like
    events (in user_atmoNC.root) and to compare preselection with possible IBD-like events from Electronics Simulation
    (user_elecsim_atmoNC.root).

    Only user_atmoNC_0.root to user_atmoNC_9.root is read (1000 events). The files atmoNC_0.root to atmoNC_9.root are
    used as input for elecsim (tut_det2elec.py). Output of elecsim is saved in files user_elecsim_atmoNC_0.root to
    user_elecsim_atmoNC_9.root (also 1000 events).

    IBD preselection is then checked by comparing results from detsim with results from elecsim.

    Preselction criteria for possible IBD-like events:
    Criteria for IBD signal (from JUNO PhysicsReport page 39):
        1.  fiducial volume cut: r < 17 m (InitX, InitY, InitZ, ExitX, ExitY, ExitZ from root file)
        2.  prompt energy cut: 10 MeV < E_prompt < 105 MeV (edep from root file)
        3.  delayed energy cut: 1.9 MeV < E_delayed < 2.5 MeV (neutron is captured on Carbon and releases gamma of
            around 2.2 MeV) (edep from roo file)
"""

import datetime
import NC_background_functions
import numpy as np
from matplotlib import pyplot as plt

# get the date and time, when the script was run:
date = datetime.datetime.now()
now = date.strftime("%Y-%m-%d %H:%M")

""" set the different cut properties: """
# fiducial volume cut in mm:
R_cut_mm = 17000
# prompt energy cut in MeV:
E_prompt_min = 10.0
E_prompt_max = 105.0
# delayed energy cut in MeV:
E_delayed_min = 1.9
E_delayed_max = 2.5
# time cut between prompt and delayed signal in ns:
time_cut_min = 600
time_cut_max = 1000000

""" set the number of the first file and number of the last file that should be read: """
start_number = 0
stop_number = 699
# number of entries in the input files:
Number_entries_input = 100
# set the path of the inputs:
input_path = "/local/scratch1/pipc51/astro/blum/detsim_output_data/"

""" parameters to check """
# event ID of preselected events (array of float):
event_id_preselected = np.array([])
# number of events, that pass preselection (float):
number_preselected = 0
# number of events, that are rejected (float):
number_rejected = 0
# preallocate number of simulated events:
number_events = 0

# loop over the files that are read:
for index in range(start_number, stop_number+1):

    # file name of the input file:
    input_name = input_path + "user_atmoNC_{0:d}.root".format(index)
    print("------------------------------------------------------------------------------------")
    print(input_name)

    # preselection of detsim events:
    num_events, evt_id_preselected, num_preselected, num_rejected = \
        NC_background_functions.preselect_sample_detsim_user(input_name, R_cut_mm, E_prompt_min, E_prompt_max,
                                                             E_delayed_min, E_delayed_max, time_cut_max,
                                                             Number_entries_input)

    # add numbers to parameters:
    number_events = number_events + num_events
    event_id_preselected = np.append(event_id_preselected, evt_id_preselected)
    number_preselected = number_preselected + num_preselected
    number_rejected = number_rejected + num_rejected

    print("event ID of preselected events:")
    print(evt_id_preselected)

print("\n-----------------------------------------------------------")
print("results from user_atmoNC_{0:d}.root to user_atmoNC_{1:d}.root".format(start_number, stop_number))
print("number of total events = {0:d}".format(number_events))
print("number of preselected events = {0:d}".format(number_preselected))
print("number of rejected events = {0:d}".format(number_rejected))









