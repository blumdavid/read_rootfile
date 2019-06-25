""" Script to read sample_detsim_user.root files (user_atmoNC_noopt_{}.root) and to make preselection of possible
    IBD-like events (in user_atmoNC_noopt_{}.root).

    These pre-selected events should then be simulated again, but now with optical properties, to be able to get the
    energy of the prompt signal and to be able to get the hittime of all photons (important for later PSD).

    Files that should be read: user_atmoNC_noopt_0.root to user_atmoNC_noopt_249.root
    Each file contains 1000 events. Therefore, 250000 events are analyzed.

    Preselection criteria for possible IBD-like events:
    Criteria for IBD signal (from JUNO PhysicsReport page 39):
        1.  fiducial volume cut: r < 16 m (because of fast neutron background reduction)
        2.  prompt energy cut -> is done later after preselection
        3.  delayed energy cut: only 1 neutron is captured in total 1 ms time-window
                                (neutron-multiplicity cut with nCapture tree)
        4.  time cut between prompt and delayed signal: deltaT < 1 ms
        5.  Distance cut between prompt and delayed signal: deltaR < 1.5 m
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
R_cut_mm = 16000
# minimum of total deposit energy of one event in MeV:
dep_energy_min = 10.0
# maximum of total deposit energy of one event in MeV:
dep_energy_max = 200.0
# time cut between prompt and delayed signal in ns:
time_cut_min = 500
time_cut_max = 1000000
# distance cut between prompt and delayed signal in mm:
dist_cut_mm = 1500

""" set the number of the first file and number of the last file that should be read: """
start_number = 0
stop_number = 999
# number of entries in the input files:
Number_entries_input = 100
# set the path of the inputs:
input_path = "/local/scratch1/pipc51/astro/blum/detsim_output_data/"
# set the path of the output:
output_path = "/home/astro/blum/juno/atmoNC/data_NC/output_preselection/preselection_detsim/"

""" parameters to check """
# event ID of preselected events (array of float):
event_id_preselected = np.array([])
# total deposit energy of the event:
e_dep_total = np.array([])
# number of events, that pass preselection (float):
number_preselected = 0
# number of events, that are rejected (float):
number_rejected = 0
# preallocate number of simulated events:
number_events = 0
# number of events, that pass volume cut:
number_vol_pass = 0
# number of events, that are rejected by volume cut:
number_vol_reject = 0
# number of events that are rejected by volume cut of init_geninfo position:
number_vol_reject_init_geninfo = 0
# number of events that are rejected by volume cut of exit_geninfo position:
# number_vol_reject_exit_geninfo = 0
# number of events that are rejected by volume cut of edep_prmtrkdep position:
# number_vol_reject_edep_prmtrkdep = 0
# number of events that are rejected by volume cut of start_ncap position:
number_vol_reject_start_ncap = 0
# number of events that are rejected by volume cut of stop_ncap position:
number_vol_reject_stop_ncap = 0
# number of events, that pass energy cut:
number_e_pass = 0
# number of events, that are rejected by minimum energy cut:
number_mine_reject = 0
# number of events with total deposit energy above dep_energy_max:
number_maxe_reject = 0
# number of events, that pass neutron-multiplicity cut:
number_nmult_pass = 0
# number of events, that are rejected by neutron-multiplicity cut:
number_nmult_reject = 0
# number of events without nCapture:
number_without_ncap = 0
# number of events, that pass time cut:
number_time_pass = 0
# number of events, that are rejected by time cut:
number_time_reject = 0
# number of events, that pass distance cut:
number_dist_pass = 0
# number of events, that are rejected by distance cut:
number_dist_reject = 0

# loop over the files that are read:
for index in range(start_number, stop_number+1):

    # file name of the input file:
    input_name = input_path + "user_atmoNC_{0:d}.root".format(index)
    print("------------------------------------------------------------------------------------")
    print(input_name)

    # preselection of detsim events:
    (num_events, evt_id_preselected, edep_total, num_preselected, num_rejected, num_vol_pass, num_vol_reject,
     num_vol_reject_init_geninfo, num_vol_reject_start_ncap, num_vol_reject_stop_ncap,
     num_e_pass, num_mine_reject, num_maxe_reject, num_nmult_pass, num_nmult_reject, num_without_ncap, num_time_pass,
     num_time_reject, num_dist_pass, num_dist_reject) = \
        NC_background_functions.preselect_sample_detsim_user(input_name, R_cut_mm, dep_energy_min, dep_energy_max,
                                                             time_cut_min, time_cut_max, dist_cut_mm,
                                                             Number_entries_input)

    # add numbers to parameters:
    number_events += num_events
    event_id_preselected = np.append(event_id_preselected, evt_id_preselected)
    e_dep_total = np.append(e_dep_total, edep_total)
    number_preselected += num_preselected
    number_rejected += num_rejected
    number_vol_pass += num_vol_pass
    number_vol_reject += num_vol_reject
    number_vol_reject_init_geninfo += num_vol_reject_init_geninfo
    # number_vol_reject_exit_geninfo += num_vol_reject_exit_geninfo
    # number_vol_reject_edep_prmtrkdep += num_vol_reject_edep_prmtrkdep
    number_vol_reject_start_ncap += num_vol_reject_start_ncap
    number_vol_reject_stop_ncap += num_vol_reject_stop_ncap
    number_e_pass += num_e_pass
    number_mine_reject += num_mine_reject
    number_maxe_reject += num_maxe_reject
    number_nmult_pass += num_nmult_pass
    number_nmult_reject += num_nmult_reject
    number_without_ncap += num_without_ncap
    number_time_pass += num_time_pass
    number_time_reject += num_time_reject
    number_dist_pass += num_dist_pass
    number_dist_reject += num_dist_reject

    """ save event ID of preselected events to txt file together with information about cut parameters: """
    np.savetxt(output_path + "evtID_preselected_{0:d}.txt".format(index), evt_id_preselected, fmt='%i',
               header="event ID's of NC events, that are preselected with script 'preselection_detsim_user.py' ({0})\n"
                      "{1:d} events analyzed from file: user_atmoNC_{2:d}.root\n"
                      "{3:d} events pass preselection, {4:d} events rejected by preselection.\n"
                      "Cut Parameters:\n"
                      "radius cut: {5:d} mm (initial pos. from geninfo-tree, start pos. and stop pos. of n-Capture),\n"
                      "total deposit energy cut: {6:0.1f} MeV <= edep <= {10:0.1f} MeV (edep of total event from "
                      "evt-tree),\n"
                      "time cut between prompt (initial time from geninfo-tree) and delayed signal (capture time from "
                      "nCapture-tree and kin. energy of gamma between 2.2 and 2.25 MeV): min = {7:0.1f} ns, "
                      "max = {8:0.1f} ns,\n"
                      "neutron multiplicity cut: only 1 neutron capture on Hydrogen in time window from above,\n"
                      "distance cut (initial pos. from geninfo-tree to stop pos. of nCapture, only for events that pass"
                      " neutron cut): {9:d} mm:"
               .format(now, num_events, index, num_preselected, num_rejected, R_cut_mm, dep_energy_min, time_cut_min,
                       time_cut_max, dist_cut_mm, dep_energy_max))

    # print("event ID of preselected events:")
    # print(evt_id_preselected)
    # print("total deposited energy per event:")
    # print(edep_total)

print("\n-----------------------------------------------------------")
print("results from user_atmoNC_{0:d}.root to user_atmoNC_{1:d}.root".format(start_number, stop_number))
print("number of total events = {0:d}".format(number_events))
print("number of preselected events = {0:d}".format(number_preselected))
print("number of rejected events = {0:d}".format(number_rejected))
print("------------")
print("number of events, that pass volume cut = {0:d}".format(number_vol_pass))
print("number of events, that are rejected by volume cut = {0:d}".format(number_vol_reject))
print("init_geninfo rejection = {0:d}".format(number_vol_reject_init_geninfo))
# print("exit_geninfo rejection = {0:d}".format(number_vol_reject_exit_geninfo))
# print("edep_prmtrkdep rejection = {0:d}".format(number_vol_reject_edep_prmtrkdep))
print("start_ncap rejection = {0:d}".format(number_vol_reject_start_ncap))
print("stop_ncap rejection = {0:d}".format(number_vol_reject_stop_ncap))
print("------------")
print("number of events, that pass energy cut = {0:d}".format(number_e_pass))
print("number of events, that are rejected by minimum energy cut = {0:d}".format(number_mine_reject))
print("number of events, with total deposit energy above {0:0.1f} MeV that pass preselection = {1:d}"
      .format(dep_energy_max, number_maxe_reject))
print("------------")
print("number of events without nCapture = {0:d}".format(number_without_ncap))
print("number of events, that pass neutron-multiplicity cut = {0:d}".format(number_nmult_pass))
print("number of events, that are rejected by neutron-multiplicity cut = {0:d}".format(number_nmult_reject))
print("------------")
print("number of events, that pass time cut = {0:d}".format(number_time_pass))
print("number of events, that are rejected by time cut = {0:d}".format(number_time_reject))
print("------------")
print("number of events, that pass distance cut (and n-mult. cut) = {0:d}".format(number_dist_pass))
print("number of events, that are rejected by distance cut (but pass n-mult. cut) = {0:d}".format(number_dist_reject))
