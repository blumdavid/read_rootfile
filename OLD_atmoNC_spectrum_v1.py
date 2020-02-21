""" Script to read the energy of the IBD-like events of atmospheric NC neutrino background and
    calculate the spectrum of atmospheric NC neutrino background as function of the visible energy.

    1.  Read txt files where the energy of the IBD-like events are saved.
        Read also txt files where the event ID of the IBD-like events are saved.

    2.  Put all visible energies into histogram to get the spectrum of atmospheric NC neutrino
        background as function of the visible energy

    3.  Consider the efficiencies of all cuts applied on the NC events by displaying the spectrum

    4.  Consider the event rate of NC interactions on C12 inside the detector (calculated with cross-sections and
        neutrino fluxes) and calculate the 'real' spectrum of atmospheric NC background, JUNO will measure after 10
        years of data taking

    5.  Consider the PSD efficiency calculated with OLD_pulse_shape_analysis_v1.py and calculate the spectrum of atmospheric
        NC background, JUNO will measure after 10 years of data taking, after Pulse Shape Analysis

"""
import datetime
import sys
import NC_background_functions
import numpy as np
from matplotlib import pyplot as plt

# get the date and time, when the script was run:
date = datetime.datetime.now()
now = date.strftime("%Y-%m-%d %H:%M")

# path, where the txt files with the number of pe of prompt signal are saved:
input_path = "/home/astro/blum/juno/atmoNC/data_NC/output_detsim_v1/"

# path, where the output files are saved:
output_path = "/home/astro/blum/juno/atmoNC/data_NC/output_detsim_v1/results_atmoNC/"

# set the file number of the first file to be analyzed:
first_file = 0
# set the file number of the last file to be analyzed:
last_file = 999
# number of simulated NC events per user_atmoNC_{}.root file:
number_events_in_rootfile = 100

""" prompt signal energy window: """
# minimal visible energy of prompt signal in MeV:
min_energy = 10.0
# maximal visible energy of prompt signal in MeV:
max_energy = 100.0
# bin-width of visible energy for histogram in MeV (must be the same like for DM signal, Reactor, CCatmo and DSNB;
# 100 keV = 0.1 MeV):
bin_width_energy = 0.5

""" Pulse Shape discrimination: """
# efficiency in percent, of how many NC events are cut away by PSD:
efficiency_PSD = 95.0

""" cut efficiencies (see detection_efficiency_NCevents.odt): """
efficiency_volume_cut = 100.074
efficiency_prompt_energy_cut = 98.873
efficiency_time_cut = 100.050
efficiency_delayed_energy_cut = 100.050
efficiency_neutron_multiplicity_cut = 99.947
efficiency_distance_cut = 101.723

error_efficiency_volume_cut = 0.370
error_efficiency_prompt_energy_cut = 2.006
error_efficiency_time_cut = 0.744
error_efficiency_delayed_energy_cut = 0.746
error_efficiency_neutron_multiplicity_cut = 1.028
error_efficiency_distance_cut = 1.051

# TODO-me: muon-veto?????????
efficiency_muon_veto = 100.00
error_efficiency_muon_veto = 0.00


# radius cut in m:
r_cut = 16.0
# time exposure in years:
time_in_years = 10
# time exposure in seconds:
time_seconds = 10 * 3.156 * 10 ** 7
# set booleans, that define, which plots are shown or saved (boolean):
PLOT_FLUX = False
SHOW_FLUXPLOT = True
SAVE_FLUXPLOT = False
PLOT_EVT_RATE = False
SHOW_EVT_RATE = True
SAVE_EVT_RATE = False

# number of IBD-like events from simulation (cut efficiency is not yet considered):
number_IBDlike_events_simu_without_eff = 0
number_cross_check = 0

# total number of simulated NC events:
number_NC_events_simu = (last_file + 1 - first_file) * number_events_in_rootfile

# preallocate array, where Evis is saved:
Evis_array = []

# loop over the files with the visible energy:
for index in range(first_file, last_file+1, 1):
    # print("read file {0:d}.....".format(index))

    # file name:
    input_file = input_path + "Evis_file{0:d}.txt".format(index)
    # read this file:
    Evis_arr = np.loadtxt(input_file)

    # file name, where evtID of IBD-like events are saved:
    input_file_evtID = input_path + "evtID_IBDlike_{0:d}.txt".format(index)
    # read this file:
    evtID = np.loadtxt(input_file_evtID)

    # number of IBD-like events:
    number_IBDlike_events_simu_without_eff += len(Evis_arr)
    # number of IBD-like events from evtID_IBDlike file as cross-check:
    number_cross_check += len(evtID)

    # loop over all entries in Evis_arr:
    for index1 in range(len(Evis_arr)):

        # append Evis_arr[index1] to array:
        Evis_array.append(Evis_arr[index1])

# set bin-edges of e_vis histogram in MeV:
bins_evis = np.arange(min_energy, max_energy + 2*bin_width_energy, bin_width_energy)
# build histogram:
Evis_histo_without_eff, bin_edges_evis = np.histogram(Evis_array, bins_evis)

print("number of NC events from simulation = {0:d}".format(number_NC_events_simu))
print("number of IBD-like events from simulation (without cut efficiency) = {0:d}"
      .format(number_IBDlike_events_simu_without_eff))
print("number of IBD-like events (cross-check) = {0:d}".format(number_cross_check))

""" consider cut efficiencies: """
# consider the cut efficiencies for volume, prompt energy, time, delayed energy, multiplicity, distance cut and muon
# veto cut in percent:
cut_efficiency = efficiency_volume_cut/100.0 * efficiency_prompt_energy_cut/100.0 * efficiency_time_cut/100.0 * \
                 efficiency_delayed_energy_cut/100.0 * efficiency_neutron_multiplicity_cut/100.0 \
                 * efficiency_distance_cut/100.0 * efficiency_muon_veto/100.0 * 100.0
# calculate statistical error of cut_efficiency with Gaussian error propagation in percent:
error_cut_efficiency = (efficiency_prompt_energy_cut/100.0 * efficiency_time_cut/100.0 *
                        efficiency_delayed_energy_cut/100.0 * efficiency_neutron_multiplicity_cut/100.0 *
                        efficiency_distance_cut/100.0 * efficiency_muon_veto/100.0 * error_efficiency_volume_cut +
                        efficiency_volume_cut/100.0 * efficiency_time_cut/100.0 * efficiency_delayed_energy_cut/100.0 *
                        efficiency_neutron_multiplicity_cut/100.0 * efficiency_distance_cut/100.0 *
                        efficiency_muon_veto/100.0 * error_efficiency_prompt_energy_cut +
                        efficiency_volume_cut/100.0 * efficiency_prompt_energy_cut/100.0 *
                        efficiency_delayed_energy_cut/100.0 * efficiency_neutron_multiplicity_cut/100.0 *
                        efficiency_distance_cut/100.0 * efficiency_muon_veto/100.0 * error_efficiency_time_cut +
                        efficiency_volume_cut/100.0 * efficiency_prompt_energy_cut/100.0 *
                        efficiency_time_cut/100.0 * efficiency_neutron_multiplicity_cut/100.0 *
                        efficiency_distance_cut/100.0 * efficiency_muon_veto/100.0 *
                        error_efficiency_delayed_energy_cut +
                        efficiency_volume_cut/100.0 * efficiency_prompt_energy_cut/100.0 *
                        efficiency_time_cut/100.0 * efficiency_delayed_energy_cut/100.0 * efficiency_distance_cut/100.0
                        * efficiency_muon_veto/100.0 * error_efficiency_neutron_multiplicity_cut +
                        efficiency_volume_cut/100.0 * efficiency_prompt_energy_cut/100.0 *
                        efficiency_time_cut/100.0 * efficiency_delayed_energy_cut/100.0 *
                        efficiency_neutron_multiplicity_cut/100.0 * efficiency_muon_veto/100.0 *
                        error_efficiency_distance_cut +
                        efficiency_volume_cut/100.0 * efficiency_prompt_energy_cut/100.0 *
                        efficiency_time_cut/100.0 * efficiency_delayed_energy_cut/100.0 *
                        efficiency_neutron_multiplicity_cut/100.0 * efficiency_distance_cut/100.0 *
                        error_efficiency_muon_veto)

# spectrum of all simulated IBD-like events (cut efficiencies are considered):
Evis_histo = Evis_histo_without_eff * cut_efficiency/100.0
# number of simulated IBD-like events (cut efficiencies are considered):
number_IBDlike_events_simu = number_IBDlike_events_simu_without_eff * cut_efficiency/100.0
print("number of IBD-like events from simulation (with cut efficiency) = {0:.2f}"
      .format(number_IBDlike_events_simu))

""" Event rate calculation: """
# calculate the theoretical event rate in NC events/sec in JUNO for neutrino energies from 0 MeV to 10 GeV (float)
# (event_rate = A * (flux_nue*xsec_nue + flux_nuebar*xsec_nuebar + flux_numu*xsec_numu + flux_numubar*xsec_numubar)):
event_rate = NC_background_functions.event_rate(bin_width_energy, r_cut, output_path, PLOT_FLUX, SHOW_FLUXPLOT,
                                                SAVE_FLUXPLOT, PLOT_EVT_RATE, SHOW_EVT_RATE, SAVE_EVT_RATE)

# number of NC events in JUNO after 10 years:
number_NC_events_JUNO = event_rate * time_seconds

# number of IBD-like events in JUNO after 10 years (cut efficiencies are considered):
number_IBDlike_events_JUNO = int(number_NC_events_JUNO * number_IBDlike_events_simu / number_NC_events_simu)

# normalize the spectrum of IBD-like events to the spectrum, JUNO will measure after 10 years (cut efficiencies are
# considered):
Evis_histo_JUNO = float(number_IBDlike_events_JUNO) / float(number_IBDlike_events_simu) * Evis_histo

""" display simulated spectrum: """
h1 = plt.figure(1, figsize=(15, 8))
plt.plot(bins_evis[:-1], Evis_histo, drawstyle="steps", linestyle="-", color="orange",
         label="atmospheric NC background\n(number of events = {0:.1f},\ncut efficiency = {1:.2f} % $\\pm$ {2:.2f} %)"
         .format(number_IBDlike_events_simu, cut_efficiency, error_cut_efficiency))
plt.xlim(xmin=min_energy, xmax=max_energy)
plt.ylim(ymin=0.0)
plt.xlabel("visible energy of prompt signal in MeV")
plt.ylabel("number of IBD-like events per bin (bin-width = {0:.2f} MeV)".format(bin_width_energy))
plt.title("Simulated spectrum of atmospheric NC neutrino events with IBD-like signature")
plt.legend()
plt.grid()
plt.savefig(output_path + "atmoNC_spectrum_simulated_bins{0:.0f}keV.png".format(bin_width_energy*1000))
# plt.show()
plt.close()

""" display visible spectrum in JUNO after 10 years: """
h2 = plt.figure(2, figsize=(15, 8))
plt.plot(bins_evis[:-1], Evis_histo_JUNO, drawstyle="steps", linestyle="-", color="orange",
         label="atmospheric NC background\n(number of events = {0:.1f},\ncut efficiency = {1:.2f} % $\\pm$ {2:.2f} %)"
         .format(number_IBDlike_events_JUNO, cut_efficiency, error_cut_efficiency))
plt.xlim(xmin=min_energy, xmax=max_energy)
plt.ylim(ymin=0.0)
plt.xlabel("visible energy of prompt signal in MeV")
plt.ylabel("number of IBD-like events per bin (bin-width = {0:.2f} MeV)".format(bin_width_energy))
plt.title("Expected spectrum of atmospheric NC neutrino events with IBD-like signature in JUNO after {0:.0f} years"
          .format(time_in_years))
plt.legend()
plt.grid()
plt.savefig(output_path + "atmoNC_spectrum_JUNO_bins{0:.0f}keV.png".format(bin_width_energy*1000))
# plt.show()
plt.close()

""" display visible spectrum in JUNO after 10 years after Pulse Shape Discrimination: """
# calculate spectrum after PSD:
Evis_histo_JUNO_PSD = Evis_histo_JUNO * (100.0-efficiency_PSD)/100.0
# number of events after PSD:
number_IBDlike_events_JUNO_PSD = number_IBDlike_events_JUNO * (100.0-efficiency_PSD)/100.0

h3 = plt.figure(3, figsize=(15, 8))
plt.plot(bins_evis[:-1], Evis_histo_JUNO_PSD, drawstyle="steps", linestyle="-", color="orange",
         label="atmospheric NC background after PSD\n(number of events = {0:0.1f})"
         .format(number_IBDlike_events_JUNO_PSD))
plt.xlim(xmin=min_energy, xmax=max_energy)
plt.ylim(ymin=0.0)
plt.xlabel("visible energy of prompt signal in MeV")
plt.ylabel("number of IBD-like events per bin (bin-width = {0:.2f} MeV)".format(bin_width_energy))
plt.title("Expected spectrum of atmospheric NC neutrino events with IBD-like signature in JUNO after {0:.0f} years\n"
          "after Pulse Shape Discrimination".format(time_in_years) + " (PSD efficiency $\\epsilon_{NC}$ = " +
          "{0:0.1f} %)"
          .format(efficiency_PSD))
plt.legend()
plt.grid()
plt.savefig(output_path + "atmoNC_spectrum_JUNO_afterPSD{1:.0f}_bins{0:.0f}keV.png".format(bin_width_energy*1000,
                                                                                           efficiency_PSD))
# plt.show()
plt.close()

""" save e_vis_array_JUNO to txt file (txt file must have same shape like files in folder 
    /home/astro/blum/PhD/work/MeVDM_JUNO/gen_spectrum_v2/).
    Save the array before Pulse Shape Discrimination. PSD efficiencies for real IBD signals and NC events is then 
    applied afterwards before calculating the Limits: """
# save Evis_histo_JUNO to txt-spectrum-file and information about simulation in txt-info-file:
print("... save data of spectrum to file...")
np.savetxt(output_path + 'NCatmo_onlyC12_bin{0:.0f}keV.txt'
           .format(bin_width_energy * 1000), Evis_histo_JUNO, fmt='%1.5e',
           header='Spectrum in IBD-like events/bin of atmospheric NC background events that mimic IBD signals '
                  '(calculated with OLD_atmoNC_spectrum_v1.py, {0}):'
                  '\n{3:d} NC events are simulated with JUNO detector software (tut_detsim.py).'
                  '\nNumber of IBD-like NC events in JUNO = {1:.2f}, binning of E_visible = {2:.3f} MeV,'
                  '\nNC interactions of nu_e, nu_e_bar, nu_mu and nu_mu_bar with C12 of liquid scintillator are '
                  'simulated with GENIE.'
                  '\nDeexcitation of residual isotopes are simulated with modified DSNB-NC.exe generator.'
                  '\nThen the final products are simulated with JUNO detector simulation and cuts are applied to get'
                  '\nthe number of NC events that mimic an IBD signal:'
                  '\ntotal cut efficiency = {4:.2f} %,'
                  '\nstatistical error of total cut efficiency = {11:.2f} %'
                  '\nefficiency volume cut = {5:.2f} %,'
                  '\nefficiency prompt energy cut = {6:.2f} %,'
                  '\nefficiency time cut = {7:.2f} %,'
                  '\nefficiency delayed energy cut = {8:.2f} %,'
                  '\nefficiency neutron multiplicity cut = {9:.2f} %,'
                  '\nefficiency distance cut = {10:.2f} %,'
           .format(now, number_IBDlike_events_JUNO, bin_width_energy, number_NC_events_simu, cut_efficiency,
                   efficiency_volume_cut, efficiency_prompt_energy_cut, efficiency_time_cut,
                   efficiency_delayed_energy_cut, efficiency_neutron_multiplicity_cut, efficiency_distance_cut,
                   error_cut_efficiency))
np.savetxt(output_path + 'NCatmo_info_onlyC12_bin{0:.0f}keV.txt'
           .format(bin_width_energy * 1000),
           np.array([min_energy, max_energy, bin_width_energy, time_in_years, r_cut, number_NC_events_simu,
                     number_IBDlike_events_JUNO, event_rate]),
           fmt='%1.9e',
           header='Information to simulation NCatmo_onlyC12_bin{0:.0f}keV.txt (analyzed files: user_atmoNC_{1:d}.root '
                  'to user_atmoNC_{2:d}.root):\n'
                  'values below: E_visible[0] in MeV, E_visible[-1] in MeV, interval_E_visible in MeV,'
                  '\nexposure time t_years in years, applied volume cut for radius in meter,'
                  '\nnumber of simulated NC events, number of IBD-like NC events in spectrum JUNO will measure, '
                  '\ntheoretical NC event rate in JUNO detector in NC events/sec,'
           .format(bin_width_energy * 1000, first_file, last_file))

