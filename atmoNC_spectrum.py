""" Script to read the energy of the prompt signals of preselected events of atmospheric NC neutrino background and
    calculate the spectrum of atmospheric NC neutrino background as function of the visible energy.

    1.  Read txt files where the energy of the prompt signals of the preselected events are saved.
        Read also txt files where the energy of the delayed signals of preselected events are saved.

    2.  check, if delayed energy passes delayed energy cut as a cross-check of the preselection (cuts used in script
        prompt_signal_preselected_evts.py is 2300 nPE < E_delayed < 3600 nPE)

    2.  Convert the energy of the prompt signal from number of pe to visible energy in the detector in MeV
        This conversion is based on the analysis in script check_conversion_npe_mev.py, where protons and neutrons with
        different kinetic energy are simulated and the visible energy as function of number of PE was plotted. With a
        linear fit, you get the conversion function E_vis(nPE) = 0.0007872*nPE in MeV.
        This also define the energy window of the prompt signal:
        nPE(10 MeV) = 12703.3 nPE
        nPE(100 MeV) = 127032.5 nPE

    3.  Put all these calculated visible energies into histogram to get the spectrum of atmospheric NC neutrino
        background as function of the visible energy

    4.  Consider the event rate of NC interactions on C12 inside the detector (calculated with cross-sections and
        neutrino fluxes) and calculate the 'real' spectrum of atmospheric NC background, JUNO will measure after 10
        years of data taking
"""
import datetime
import sys
import NC_background_functions
import numpy as np
from matplotlib import pyplot as plt


def conversion_npe_to_evis(number_photo_electron):
    """
    Function to calculate the visible energy in MeV for a given number of p.e. for the prompt signal.
    This function is the result of linear fit from script check_conversion_npe_mev.py.

    :param number_photo_electron: number of photo-electrons of the prompt signal
    :return: quenched deposited energy (visible energy in MeV)
    """
    # TODO-me: parameters of the fit has to be checked!!!!!!!!!
    # first fit parameter (slope) in MeV/nPE:
    parameter_a = 0.0007872

    energy = parameter_a * number_photo_electron

    return energy, parameter_a


# get the date and time, when the script was run:
date = datetime.datetime.now()
now = date.strftime("%Y-%m-%d %H:%M")

# path, where the txt files with the number of pe of prompt signal are saved:
input_path = "/home/astro/blum/juno/atmoNC/data_NC/output_detsim/"
# path, where corresponding evtID is saved:
input_path_evtID = "/home/astro/blum/juno/atmoNC/data_NC/output_preselection/preselection_detsim/"

# path, where the output files are saved:
output_path = "/home/astro/blum/juno/atmoNC/data_NC/output_detsim/results_atmoNC/"

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
bin_width_energy = 0.4
# preallocate number of events that are rejected by prompt energy cut:
number_rejected_prompt_cut = 0
# preallocate number of events where nPE of prompt signal is below min_energy:
number_rejected_prompt_cut_min = 0
# preallocate number of events where nPE of prompt signal is above max_energy:
number_rejected_prompt_cut_max = 0

""" delayed energy cut parameters: """
# minimum energy of delayed signal in nPE:
# TODO-me: check the values for delayed energy cut!!!!!!!!
min_nPE_delayed = 2600
# maximum energy of delayed signal in nPE:
max_nPE_delayed = 3600
# preallocate number of events that are rejected by delayed energy cut:
number_rejected_del_cut = 0
# preallocate number of events with nPE=0 for delayed signal:
number_nPE0_delayed = 0
# preallocate number of events where nPE of delayed signal is below min_PE_delayed:
number_rejected_del_cut_min = 0
# preallocate number of events where nPE of delayed signal is above max_PE_delayed:
number_rejected_del_cut_max = 0

""" Pulse Shape discrimination: """
# efficiency in percent, of how many NC events are cut away by PSD:
efficiency_PSD = 95.0

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

# number of preselected events:
number_preselected = 0

# total number of simulated NC events:
number_NC_events_simu = (last_file + 1 - first_file) * number_events_in_rootfile

# preallocate histogram, where visible energy is saved:
# set bin-edges of e_vis histogram in MeV:
bins_evis = np.arange(min_energy, max_energy + 2*bin_width_energy, bin_width_energy)
# preallocate empty array to build default e_vis-histogram:
e_vis_empty = np.array([])
# build default e_vis histogram:
e_vis_array, bin_edges_evis = np.histogram(e_vis_empty, bins_evis)

# loop over the files with the number of pe:
for index in range(first_file, last_file+1, 1):
    # print("read file {0:d}.....".format(index))

    # file name:
    input_file_delayed = input_path + "number_pe_delayed_file{0:d}.txt".format(index)
    # read this file:
    number_pe_delayed = np.loadtxt(input_file_delayed)

    # file name:
    input_file = input_path + "number_pe_file{0:d}.txt".format(index)
    # read this file:
    number_pe_file = np.loadtxt(input_file)

    # file name, where evtID of preselected events are saved:
    input_file_evtID = input_path_evtID + "evtID_preselected_{0:d}.txt".format(index)
    # read this file:
    evtID_preselected = np.loadtxt(input_file_evtID)

    # number of preselected events:
    number_preselected = number_preselected + len(number_pe_file)

    # loop over all entries in number_pe_file:
    for index1 in range(len(number_pe_file)):
        # convert number_pe to E_vis:
        e_vis, slope_of_linear_fit = conversion_npe_to_evis(number_pe_file[index1])

        # check, if energy is in the correct time window:
        if min_energy <= e_vis <= max_energy:
            # add e_vis to default evis histogram:
            e_vis_array += np.histogram(e_vis, bins_evis)[0]

            # get evtID of this event:
            evtID = evtID_preselected[index1]

        else:
            # event is rejected by prompt energy cut:
            number_rejected_prompt_cut += 1

            if e_vis < min_energy:
                # nPE below min_energy:
                number_rejected_prompt_cut_min += 1
            elif e_vis > max_energy:
                # nPE above max_energy:
                number_rejected_prompt_cut_max += 1

        # check delayed energy cut:
        # TODO-me: Apply also delayed energy cut -> Until now, all preselected events pass this cut!
        if number_pe_delayed[index1] < min_nPE_delayed or number_pe_delayed[index1] > max_nPE_delayed:
            # event is rejected by delayed energy cut:
            number_rejected_del_cut += 1

        if number_pe_delayed[index1] == 0:
            # delayed signal lies in prompt signal and is not calculated correctly:
            number_nPE0_delayed += 1
            # print("nPE=0: file = {0:d}, evt = {1:.0f}".format(index, evtID_preselected[index1]))

        elif 0 < number_pe_delayed[index1] < min_nPE_delayed:
            # nPE below min_nPE_delayed:
            number_rejected_del_cut_min += 1
            # print("nPE<2300: file = {0:d}, evt = {1:.0f}".format(index, evtID_preselected[index1]))

        elif number_pe_delayed[index1] > max_nPE_delayed:
            # nPE above max_nPE_delayed:
            number_rejected_del_cut_max += 1

# number of NC events in e_vis_array (from simulation):
number_IBDlike_events_simu = np.sum(e_vis_array)

print("number of NC events from simulation = {0:d}".format(number_NC_events_simu))
print("number of preselected events from simulation = {0:d}".format(number_preselected))
print("number of IBD-like from simulation = {0:d}\n".format(number_IBDlike_events_simu))
print("number of events (of preselect. evts), that would be rejected by delayed energy cut = {0:d}"
      .format(number_rejected_del_cut))
print("number of events (of preselect. evts) with nPE=0 of delayed signal = {0:d}"
      .format(number_nPE0_delayed))
print("number of events (of preselect. evts) with nPE < min_PE_delayed = {0:d}"
      .format(number_rejected_del_cut_min))
print("number of events (of preselect. evts) with nPE > max_PE_delayed = {0:d}"
      .format(number_rejected_del_cut_max))
print("number of events (of preselect. evts) rejected by prompt energy cut = {0:d}".format(number_rejected_prompt_cut))
print("number of events (of preselect. evts) with nPE of prompt signal < min_energy = {0:d}"
      .format(number_rejected_prompt_cut_min))
print("number of events (of preselect. evts) with nPE of prompt signal > max_energy = {0:d}"
      .format(number_rejected_prompt_cut_max))

""" Event rate calculation: """
# calculate the theoretical event rate in NC events/sec in JUNO for neutrino energies from 0 MeV to 10 GeV (float)
# (event_rate = A * (flux_nue*xsec_nue + flux_nuebar*xsec_nuebar + flux_numu*xsec_numu + flux_numubar*xsec_numubar)):
event_rate = NC_background_functions.event_rate(bin_width_energy, r_cut, output_path, PLOT_FLUX, SHOW_FLUXPLOT,
                                                SAVE_FLUXPLOT, PLOT_EVT_RATE, SHOW_EVT_RATE, SAVE_EVT_RATE)

# number of NC events in JUNO after 10 years:
number_NC_events_JUNO = event_rate * time_seconds

# number of IBD-like events in JUNO after 10 years:
number_IBDlike_events_JUNO = int(number_NC_events_JUNO * number_IBDlike_events_simu / number_NC_events_simu)

# normalize the spectrum of simulated events to the spectrum, JUNO will measure after 10 years:
e_vis_array_JUNO = float(number_IBDlike_events_JUNO) / float(number_IBDlike_events_simu) * e_vis_array

""" display simulated spectrum: """
h1 = plt.figure(1, figsize=(15, 8))
plt.plot(bins_evis[:-1], e_vis_array, drawstyle="steps", linestyle="-", color="orange",
         label="atmospheric NC background\n(number of events = {0:d})".format(number_IBDlike_events_simu))
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
plt.plot(bins_evis[:-1], e_vis_array_JUNO, drawstyle="steps", linestyle="-", color="orange",
         label="atmospheric NC background\n(number of events = {0:d})".format(number_IBDlike_events_JUNO))
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
e_vis_array_JUNO_PSD = e_vis_array_JUNO * (100.0-efficiency_PSD)/100.0
# number of events after PSD:
number_IBDlike_events_JUNO_PSD = number_IBDlike_events_JUNO * (100.0-efficiency_PSD)/100.0

h3 = plt.figure(3, figsize=(15, 8))
plt.plot(bins_evis[:-1], e_vis_array_JUNO_PSD, drawstyle="steps", linestyle="-", color="orange",
         label="atmospheric NC background after PSD\n(number of events = {0:0.1f})"
         .format(number_IBDlike_events_JUNO_PSD))
plt.xlim(xmin=min_energy, xmax=max_energy)
plt.ylim(ymin=0.0)
plt.xlabel("visible energy of prompt signal in MeV")
plt.ylabel("number of IBD-like events per bin (bin-width = {0:.2f} MeV)".format(bin_width_energy))
plt.title("Expected spectrum of atmospheric NC neutrino events with IBD-like signature in JUNO after {0:.0f} years\n"
          "after Pulse Shape Discrimination".format(time_in_years) + " (cut efficiency $\\epsilon_{NC}$ = " +
          "{0:0.1f} %)"
          .format(efficiency_PSD))
plt.legend()
plt.grid()
plt.savefig(output_path + "atmoNC_spectrum_JUNO_afterPSD_bins{0:.0f}keV.png".format(bin_width_energy*1000))
# plt.show()
plt.close()

""" save e_vis_array_JUNO to txt file (txt file must have same shape like files in folder 
    /home/astro/blum/PhD/work/MeVDM_JUNO/gen_spectrum_v2/).
    Save the array before Pulse Shape Discrimination. PSD efficiencies for real IBD signals and NC events is then 
    applied afterwards before calculating the Limits: """
# save e_vis_array_JUNO to txt-spectrum-file and information about simulation in txt-info-file:
print("... save data of spectrum to file...")
np.savetxt(output_path + 'NCatmo_onlyC12_bin{0:.0f}keV.txt'
           .format(bin_width_energy * 1000), e_vis_array_JUNO, fmt='%1.5e',
           header='Spectrum in IBD-like events/bin of atmospheric NC background events that mimic IBD signals '
                  '(calculated with atmoNC_spectrum.py, {0}):'
                  '\n{3:d} NC events are simulated with JUNO detector software (tut_detsim.py).'
                  '\nNumber of IBD-like NC events from spectrum = {1:.5f}, binning of E_visible = {2:.3f} MeV,'
                  '\nNC interactions of nu_e, nu_e_bar, nu_mu and nu_mu_bar with C12 of liquid scintillator are '
                  'simulated with GENIE.'
                  '\nDeexcitation of residual isotopes are simulated with modified DSNB-NC.exe generator.'
                  '\nThen the final products are simulated with JUNO detector simulation and cuts are applied to get'
                  '\nthe number of NC events that mimic an IBD signal:'
           .format(now, number_IBDlike_events_JUNO, bin_width_energy, number_NC_events_simu))
np.savetxt(output_path + 'NCatmo_info_onlyC12_bin{0:.0f}keV.txt'
           .format(bin_width_energy * 1000),
           np.array([min_energy, max_energy, bin_width_energy, time_in_years, r_cut, number_NC_events_simu,
                     number_IBDlike_events_JUNO, event_rate, slope_of_linear_fit]),
           fmt='%1.9e',
           header='Information to simulation NCatmo_onlyC12_bin{0:.0f}keV.txt (analyzed files: user_atmoNC_{1:d}.root '
                  'to user_atmoNC_{2:d}.root):\n'
                  'values below: E_visible[0] in MeV, E_visible[-1] in MeV, interval_E_visible in MeV,'
                  '\nexposure time t_years in years, applied volume cut for radius in meter,'
                  '\nnumber of simulated NC events, number of IBD-like NC events in spectrum, '
                  '\ntheoretical NC event rate in JUNO detector in NC events/sec,'
                  '\nslope A of the linear conversion function E_vis[MeV] = A * nPE:'
           .format(bin_width_energy * 1000, first_file, last_file))


