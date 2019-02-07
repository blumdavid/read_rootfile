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


        6. Muon veto criteria (83% efficiency by yellow book): -> Do I have to consider this??????????????


"""

# import ROOT
import datetime
import NC_background_functions
import numpy as np
from matplotlib import pyplot as plt

# get the date and time, when the script was run:
date = datetime.datetime.now()
now = date.strftime("%Y-%m-%d %H:%M")

# set SAVE_TXT, defines if txt files are saved:
SAVE_TXT = True

""" set properties like used in the simulation of reactor, DSNB, CCatmo background and DM signal): 
"""
# minimum of visible energy in MeV (float):
E_vis_min = 10.0
# maximum of visible energy in MeV (float):
E_vis_max = 100.0
# interval (bin-width) of visible energy in MeV (float):
interval_E_vis = 0.1
# total exposure time in years (float):
time_exposure_years = 10
# total exposure time in seconds (float):
time_exposure_sec = time_exposure_years * 3.156 * 10 ** 7

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
# distance cut between prompt and delayed signal in mm:
Distance_cut = 1500
# time in ns, where two prompt signals cannot be separated anymore (very conservative limit 1000 ns = 1 microsecond):
# TODO-me: what is the time resolution where you can separate 2 prompt signals?
time_resolution = 10

""" set the number of the first file and number of the last file that should be read: """
start_number = 0
stop_number = 499
# number of entries in the input files:
Number_entries_input = 100
# set the path of the inputs:
input_path = "/local/scratch1/pipc51/astro/blum/detsim_output_data/"
# set boolean, that define, if spectrum of simulated IBD-like events is shown or not:
SHOW_SIMU_SPECTRUM = True

""" set parameters, which define the calculation of the event rate: """
# interval-width (bin-width) of the energy array (energy from 0 MeV to 10 GeV) in MeV (float):
interval_energy = 0.1
# radius of the fiducial volume cut in meter (float):
R_cut_meter = R_cut_mm / 1000
# set booleans, that define, which plots are shown or saved (boolean):
PLOT_FLUX = True
SHOW_FLUXPLOT = False
SAVE_FLUXPLOT = True
PLOT_EVT_RATE = True
SHOW_EVT_RATE = False
SAVE_EVT_RATE = True


""" Spectrum of IBD-like atmospheric NC neutrino background events: """
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
    num_evts, evt_ID_IBD, E_vis = NC_background_functions.read_sample_detsim_user(input_name, R_cut_mm, E_prompt_min,
                                                                                  E_prompt_max, E_delayed_min,
                                                                                  E_delayed_max, time_cut_min,
                                                                                  time_cut_max, Distance_cut,
                                                                                  time_resolution, Number_entries_input)

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


""" get spectrum of simulated IBD-like NC events: """
# Array of visible energy with finer binning in MeV (array of float):
array_visible_energy = np.arange(E_vis_min, E_vis_max+interval_E_vis, interval_E_vis)
# put E_visible in histogram to get the number of events as function of the visible energy (spectrum_NC: values of the
# histogram bins (IBD-like NC events per MeV, NOT normalized to 10 years); edges_E_vis: edges of the bins
# (length nbins + 1), visible energy in MeV):
# TODO-me: is spectrum_NC really 1/MeV?
h3 = plt.figure(3, figsize=(15, 8))
spectrum_NC, edges_E_vis, patches = plt.hist(E_visible, np.append(array_visible_energy,
                                                                  array_visible_energy[-1]+interval_E_vis),
                                             histtype="step", align="mid", color="blue",
                                             label="number of events = {0:.2f}".format(number_IBDevts))
plt.xlabel("Visible energy in MeV")
plt.ylabel("IBD-like events per MeV")
plt.title("Spectrum of IBD-like atmospheric NC neutrino background events (not normalized)")
plt.xlim(xmin=E_vis_min, xmax=E_vis_max)
plt.ylim(ymin=0.0)
plt.legend()
if SHOW_SIMU_SPECTRUM:
    plt.show()
else:
    plt.close(h3)


""" Event rate calculation: """
# calculate the theoretical event rate in events/sec for neutrino energies from 0 MeV to 10 GeV (float)
# (event_rate = A * (flux_nue*xsec_nue + flux_nuebar*xsec_nuebar + flux_numu*xsec_numu + flux_numubar*xsec_numubar)):
event_rate = NC_background_functions.event_rate(interval_energy, R_cut_meter, PLOT_FLUX, SHOW_FLUXPLOT, SAVE_FLUXPLOT,
                                                PLOT_EVT_RATE, SHOW_EVT_RATE, SAVE_EVT_RATE)

# This event rate (event_rate) corresponds to the number of simulated events (number_events) per time "time" ("time" is
# unknown). The time "time, in which 'number_events' NC interactions proceed (stattfinden), is given by
# number_events/event_rate, so number of simulated NC events / event rate.
# time in seconds (float):
time = number_events / event_rate

# number of NC interactions during time exposure (normally 10 years) (float):
number_NC_interaction_exposure = event_rate * time_exposure_sec

# number of IBD-like events of the NC interactions during time exposure (float):
number_IBD_exposure = number_NC_interaction_exposure * number_IBDevts / number_events
print("number of IBD-like events after {0:d} years = {1:.2f}".format(time_exposure_years, number_IBD_exposure))

# Normalize the spectrum of IBD-like events (spectrum_NC). Normalized spectrum of IBD-like NC events in evts/MeV
# (array of float):
spectrum_NC_norm = number_IBD_exposure / number_IBDevts * spectrum_NC

h4 = plt.figure(4, figsize=(15, 8))
plt.plot(array_visible_energy, spectrum_NC_norm, "b", drawstyle="steps-mid", label="number of events = {0:.2f}"
         .format(number_IBD_exposure))
# plt.plot(edges_E_vis[:-1], spectrum_NC_norm, "r--", drawstyle="steps-mid")
plt.xlabel("Visible energy in MeV")
plt.ylabel("IBD-like events per MeV")
plt.title("Spectrum of IBD-like atmospheric NC neutrino background events in JUNO after {0:d} years"
          .format(time_exposure_years))
plt.xlim(xmin=E_vis_min, xmax=E_vis_max)
plt.ylim(ymin=0.0)
plt.legend()
plt.grid()
plt.show()




