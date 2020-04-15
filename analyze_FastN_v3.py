""" Script to analyze Fast Neutron Background in JUNO detector:

    Version 3:
    -   consider shape of energy spectrum of used neutrons -> must be flat distribution
    -   use different tail windows and TTR cut values for different energy regions (results from PSD analysis of
        NC and IBD events)

    In the first part of the script:
    -   Calculate the pulse shape of the prompt signal of fast neutrons.
    -   As simulated neutrons, the files user_neutron_10_MeV_{}.root, user_neutron_100_MeV_{}.root,
        user_neutron_300_MeV_{}.root and user_neutron_500_MeV_{}.root (from 'conversion_nPE_MeV') are used.
    -   Files user_neutron_10_MeV_0.root to user_neutron_10_MeV_99.root can not be used, because they deposit too less
        energy in the detector and therefore the hittime distribution has only a few ten entries (Qedep only from
        0.2 MeV to 5 MeV).
        For files user_neutron_100_MeV_0.root to user_neutron_100_MeV_99.root, Qedep is from 1.2 MeV to 7.4 MeV.
        For files user_neutron_300_MeV_0.root to user_neutron_300_MeV_99.root, Qedep is from 1.5 MeV to 40 MeV.
        For files user_neutron_500_MeV_0.root to user_neutron_500_MeV_99.root, Qedep is from 0.0 MeV to 113 MeV.
    -   Procedure to get the hittime distribution with vertex reconstruction and time smearing of PMTs:

    1.  apply same cuts like on prompt signals of NC events:
        1.1 energy cut on prompt signal: only neutrons with energy from 10 MeV to 100 MeV (here either 10 or 100 MeV)
            are simulated -> energy cut is applied automatically
        1.2 volume cut: must be same like for NC events (take initial position of initial particle -> smear it with
            vertex resolution with function position_smearing())

    2.  calculate time of flight:
        2.1 for every photon, that hits a PMT (20inch and 3inch), take the PMT position (via PMT ID from file
            PMT_position.root) and calculate the time-of-flight distance with the reconstructed position from above.
        2.2 with the time-of-flight distance, calculate the time-of-flight of this photon from production to PMT by
            considering an effective speed of light in the LS.

    3.  consider TTS of PMTs:
        3.1 for every photon, that hits a PMT (20inch and 3inch), take the time resolution (sigma) of the PMT
            (via PMT ID either from file PmtData.root for the 20inch PMTs or set TTS = 5 ns for 3inch PMTs.)
        3.2 the TTS of the PMT is FWHM. Therefore calculate sigma from TTS (FWHM = 2*sqrt(2*ln(2)) * sigma).
        3.3 smear hittime of detsim with gaussian of sigma (time resolution) around the value of detsim hittime to get
        the smeared hittime

    4.  for every photon, calculate the 'real' hittime (= smeared hittime - time_of_flight) and store it in array

    5.  Do points 2. to 4. for every photon. Then you get the correct hittime of this event. Build histogram with
        correct hittimes and save histogram value in txt file and display histogram in png file


    In the second part of the script:
    -   use calculated pulse shapes of fast neutrons and



"""
import datetime
import ROOT
import sys
import NC_background_functions
import numpy as np
from matplotlib import pyplot as plt
import os

# get the date and time, when the script was run:
date = datetime.datetime.now()
now = date.strftime("%Y-%m-%d %H:%M")

""" set flag, if hittimes must be calculated or read from file: """
flag_read_hittime_from_file = True

""" define time window and bin width: """
# set time window of whole signal in ns:
min_time = -50
max_time = 2000
# set time in ns, where the prompt signal should be 0:
time_limit_prompt = 1000
# Set bin-width of hittime histogram in ns:
binwidth = 5.0

""" PSD parameter from analyze_PSD_cut_IBDlike_v2.py and atmoNC_spectrum_v3.py, respectively: """
# total NC PSD suppression in percent:
NC_supp_total = 99.12
# total IBD PSD suppression in percent:
IBD_supp_total = 11.60
# PSD suppression of NC and IBD events in %, tail window in ns and TTR cut value for energies between 10 MeV and 20 MeV:
PSD_NC_supp_10_20 = 98.93
PSD_IBD_supp_10_20 = 74.24
tail_start_10_20 = 25.0
tail_stop_10_20 = 600.0
TTR_cut_10_20 = 0.27355
# PSD suppression of NC and IBD events in %, tail window in ns and TTR cut value for energies between 20 MeV and 30 MeV:
PSD_NC_supp_20_30 = 98.25
PSD_IBD_supp_20_30 = 32.96
tail_start_20_30 = 325.0
tail_stop_20_30 = 600.0
TTR_cut_20_30 = 0.01125
# PSD suppression of NC and IBD events in %, tail window in ns and TTR cut value for energies between 30 MeV and 40 MeV:
PSD_NC_supp_30_40 = 98.63
PSD_IBD_supp_30_40 = 0.6
tail_start_30_40 = 325.0
tail_stop_30_40 = 600.0
TTR_cut_30_40 = 0.01137
# PSD suppression of NC and IBD events in %, tail window in ns and TTR cut value for energies
# between 40 MeV and 100 MeV:
PSD_NC_supp_40_100 = 99.96
PSD_IBD_supp_40_100 = 0.04
tail_start_40_100 = 325.0
tail_stop_40_100 = 800.0
TTR_cut_40_100 = 0.01489

""" set cut parameters from analysis of NC events and IBD events: """
# cut-radius in mm:
radius_cut = 16000
# prompt energy cut in MeV:
E_prompt_min = 10.0
E_prompt_max = 100.0
# time cut in ns:
T_cut_min = 1000.0
T_cut_max = 1000000.0
# multiplicity:
multiplicity = 1
# delayed energy cut in MeV:
E_delayed_min_MeV = 1.80
E_delayed_max_MeV = 2.55
# distance cut in mm:
distance_cut = 500.0
# delayed volume cut in mm:
delayed_volume = 17700.0
# PSD suppression to get the correct path to the folder with the TTR values of NC and IBD events:
suppression_for_folder = 99.0

# path, where root files of neutron simulation are saved:
input_path_neutron = "/local/scratch1/pipc51/astro/blum/conversion_nPE_MeV/neutron_output/"

# path, where TTR values of NC and IBD events, that pass all cuts, are stored:
input_path_TTR = ("/home/astro/blum/juno/atmoNC/data_NC/output_detsim_v2/DCR_results_{0:.0f}mm_{1:.0f}MeVto{2:.0f}MeV_"
                  "{3:.0f}nsto{4:.0f}ms_mult{5}_{6:.0f}keVto{7:.0f}keV_dist{8:.0f}mm_R{9:.0f}mm_PSD{10:.0f}/"
                  "test_10to20_20to30_30to40_40to100_final/"
                  .format(radius_cut, E_prompt_min, E_prompt_max, T_cut_min, T_cut_max/1000000, multiplicity,
                          E_delayed_min_MeV*1000, E_delayed_max_MeV*1000, distance_cut, delayed_volume,
                          suppression_for_folder))

# path, where hittime distributions (png and txt) are saved:
input_path_hittimes = "/home/astro/blum/PhD/work/MeVDM_JUNO/fast_neutrons/"
# output_path_neutron = "/home/astro/blum/PhD/work/MeVDM_JUNO/fast_neutrons/"
output_path_neutron = "/home/astro/blum/PhD/work/MeVDM_JUNO/fast_neutrons/DCR/"

""" analyze neutron hittime distribution for kinetic energy of 10 MeV and 100 MeV: """
# first file of neutron simulation:
first_file_neutron = 0
# last file of neutron simulation (100 corresponds to user_neutron_300_MeV_0.root, 200 corresponds to
# user_neutron_500_MeV_0.root and 299 corresponds to user_neutron_500_MeV_99.root):
last_file_neutron = 1099
# number of events per file:
number_evts_per_file = 10
# total number of neutron events (factor 2 for 10 MeV and 100 MeV):
number_evts_total = (last_file_neutron - first_file_neutron + 1) * number_evts_per_file
# bin width of energy spectrum:
binwidth_energy = 1.0
# array corresponding to energy range:
energy_range = np.arange(E_prompt_min, E_prompt_max+binwidth_energy, binwidth_energy)
# preallocate number of events that are analyzed (pass volume cut):
number_analyzed_total = 0
number_analyzed_array = np.zeros(len(energy_range))
# array, where TTR values of neutron events are saved:
array_TTR = []
array_TTR_10_20 = []
array_TTR_20_30 = []
array_TTR_30_40 = []
array_TTR_40_100 = []
# array, where prompt energy of neutron events is saved in MeV:
array_E_prompt_neutron = []
array_E_10_20 = []
array_E_20_30 = []
array_E_30_40 = []
array_E_40_100 = []
# number of neutron events, that pass the PSD cut:
number_pass_PSD_total = 0
number_pass_PSD_array = np.zeros(len(energy_range))

""" load position of the PMTs and corresponding PMT ID from file PMT_position.root: """
file_PMT_position = "/home/astro/blum/juno/atmoNC/PMT_information/PMT_position.root"
# array with PMT ID and corresponding x, y, z position in mm:
pmtID_pos_file, x_pos_pmt, y_pos_pmt, z_pos_pmt = NC_background_functions.get_pmt_position(file_PMT_position)

""" load 'time resolution' in ns of the 20 inch PMTs and corresponding PMT ID from file PmtData.root: """
file_PMT_time = "/home/astro/blum/juno/atmoNC/PMT_information/PmtData_old.root"
# array with PMT ID and corresponding sigma in ns:
pmtID_time_file, sigma_time_20inch = NC_background_functions.get_20inchpmt_tts(file_PMT_time)
# set TTS (FWHM) of the 3inch PMTs in ns:
tts_3inch = 5.0
# calculate time resolution (sigma) for the 3inch PMTs in ns:
sigma_time_3inch = tts_3inch / (2 * np.sqrt(2 * np.log(2)))
# set effective speed of light in the liquid scintillator in mm/ns (see page 12 of
# 20200111_zli_VertexReconstruction_page20.pdf in folder /home/astro/blum/PhD/paper/reconstruction/).
# The effective refraction index in LS depends on the TTS of the PMT (Hamamatsu with TTS ~ 2.7 ns,
# NNVT with TTS ~ 18 ns).
# for Hamamatsu and 3inch PMTs (TTS ~ 2.7 ns and 5 ns) use n_eff = 1.544 (c/n_eff = 299792458 m / 1.544 s
# = 194166100 * 10**(-6) mm/ns ~ 194.17 mm/ns):
c_effective_smallTTS = 194.17
# for MCP PMTs (TTS ~ 18 ns) use n_eff = 1.578 (c/n_eff = 299792458 m / 1.578 s = 189982546 * 10**(-6) mm/ns ~
# 189.98 mm/ns):
c_effective_largeTTS = 189.98

# if pulse shapes are not calculated already, read ROOT files and save pulse shape to txt file:
if not flag_read_hittime_from_file:

    # loop over root files with neutron simulation:
    for index in range(first_file_neutron, last_file_neutron + 1, 1):
        # load user_neutron_{}.root file:
        if index < 100:
            rfile = ROOT.TFile(input_path_neutron + "user_neutron_10_MeV_{0:d}.root".format(index))
        elif 100 <= index < 200:
            rfile = ROOT.TFile(input_path_neutron + "user_neutron_100_MeV_{0:d}.root".format(index-100))
        elif 200 <= index < 300:
            rfile = ROOT.TFile(input_path_neutron + "user_neutron_300_MeV_{0:d}.root".format(index-200))
        elif 300 <= index < 1000:
            rfile = ROOT.TFile(input_path_neutron + "user_neutron_500_MeV_{0:d}.root".format(index-300))
        elif 1000 <= index <= 1099:
            rfile = ROOT.TFile(input_path_neutron + "user_neutron_1000_MeV_{0:d}.root".format(index-1000))
        else:
            print("--------------------------------- ERROR-------------------")

        print(rfile)

        # get the "evt"-TTree from the TFile:
        rtree_evt = rfile.Get("evt")
        # get geninfo tree from TFile:
        rtree_geninfo = rfile.Get("geninfo")
        # get prmtrkdep tree from TFile:
        rtree_prmtrkdep = rfile.Get("prmtrkdep")
        # get the number of events in the 'evt' Tree:
        num_events_evt = rtree_evt.GetEntries()
        # check number of events:
        if num_events_evt != number_evts_per_file:
            sys.exit("ERROR: number of events in root file ({0:d}) != {1:d}"
                     .format(num_events_evt, number_evts_per_file))

        # loop over the events:
        for event in range(num_events_evt):

            """ check volume cut: """
            # get current event in prmtrkdep tree:
            rtree_prmtrkdep.GetEntry(event)
            # get number of initial particles:
            n_init_part = int(rtree_prmtrkdep.GetBranch('nInitParticles').GetLeaf('nInitParticles').GetValue())

            if n_init_part != 1:
                # check if there is just one initial positron:
                sys.exit("ERROR: more than 1 initial particles in event {0:d}".format(event))

            # get quenched deposited energy of the initial particle in MeV:
            qedep_prmtrkdep = float(rtree_prmtrkdep.GetBranch("Qedep").GetLeaf("Qedep").GetValue())

            # do qedep cut (qedep must be between 10 and 100 MeV):
            if qedep_prmtrkdep < E_prompt_min:
                continue
            if qedep_prmtrkdep > E_prompt_max:
                continue

            # get current event in geninfo tree:
            rtree_geninfo.GetEntry(event)

            # get number of initial particles:
            n_init_particles = int(rtree_geninfo.GetBranch('nInitParticles').GetLeaf('nInitParticles').GetValue())

            if n_init_particles != 1:
                # check if there is just one initial positron:
                sys.exit("ERROR: more than 1 initial particles in event {0:d}".format(event))

            # get initial x, y, z position:
            x_init = float(rtree_geninfo.GetBranch('InitX').GetLeaf('InitX').GetValue())
            y_init = float(rtree_geninfo.GetBranch('InitY').GetLeaf('InitY').GetValue())
            z_init = float(rtree_geninfo.GetBranch('InitZ').GetLeaf('InitZ').GetValue())

            # do vertex reconstruction with function position_smearing():
            # Smear x,y and z position of the initial position (returns reconstructed position in mm):
            x_reconstructed = NC_background_functions.position_smearing(x_init, qedep_prmtrkdep)
            y_reconstructed = NC_background_functions.position_smearing(y_init, qedep_prmtrkdep)
            z_reconstructed = NC_background_functions.position_smearing(z_init, qedep_prmtrkdep)

            # calculate distance to detector center in mm:
            r_reconstructed = np.sqrt(x_reconstructed**2 + y_reconstructed**2 + z_reconstructed**2)

            # check if reconstructed position is within 17 m:
            if r_reconstructed >= 17000.0:
                continue
            else:
                number_analyzed += 1

            """ calculate the photon emission time distribution (time of flight correction with reconstructed position 
            and time smearing with TTS for each hit): """
            # get event of 'evt'-tree:
            rtree_evt.GetEntry(event)
            # get evtID of the tree and compare with event:
            evt_id = int(rtree_evt.GetBranch('evtID').GetLeaf('evtID').GetValue())
            if evt_id != event:
                sys.exit("ERROR: evtID of tree ({0:d}) != {1:d}".format(evt_id, event))

            # get number of photons of this event:
            n_photons = int(rtree_evt.GetBranch('nPhotons').GetLeaf('nPhotons').GetValue())

            # preallocate list, where corrected (real) hit-times are saved:
            hittime_array = []

            # loop over every photon in the event:
            for index1 in range(n_photons):

                # get number of pe per photon and check if it is equal to 1:
                npe = int(rtree_evt.GetBranch('nPE').GetLeaf('nPE').GetValue(index1))
                if npe != 1:
                    sys.exit("ERROR: more than one p.e. per photon in event {0:d}, file {1}".format(event, index))

                # get the pmtID of the hit PMT:
                pmtID = int(rtree_evt.GetBranch('pmtID').GetLeaf('pmtID').GetValue(index1))

                """ time of flight correction: """
                # get hittime of PMT from tree in ns:
                hittime = float(rtree_evt.GetBranch('hitTime').GetLeaf('hitTime').GetValue(index1))

                # get position of the PMT with specific pmtID (pmtID is ascending number from 0 to 17738
                # (17739 large PMTs) and from 300000 to 336571 (36572 small PMTs)).
                # For large PMTs -> For 20inch PMTs, the pmtID is equal to index of x,y,z_pos_pmt array.
                # For small PMTs -> For 3inch PMTs, the pmtID - (300000 - 17739) is equal to index of
                # x,y,z_pos_pmt array.
                # check if PMT is 20 inch or 3inch (pmtID < 50000 means 20inch PMT):
                if pmtID < 50000:
                    # 20inch PMT:
                    # get PMT position in mm from arrays:
                    x_pmt = x_pos_pmt[pmtID]
                    y_pmt = y_pos_pmt[pmtID]
                    z_pmt = z_pos_pmt[pmtID]

                else:
                    # 3inch PMT:
                    # calculate index of pos_pmt array that correspond to pmtID of 3inch PMTs (for example:
                    # first small PMT: 300000-282261 = 17739, last small PMT: 336571-282261 = 54310)
                    index_3inch = pmtID - 282261
                    # get PMT position in mm from arrays:
                    x_pmt = x_pos_pmt[index_3inch]
                    y_pmt = y_pos_pmt[index_3inch]
                    z_pmt = z_pos_pmt[index_3inch]

                # calculate distance between reconstructed position of event and position of PMT (in mm):
                distance_tof = np.sqrt((x_reconstructed - x_pmt)**2 + (y_reconstructed - y_pmt)**2 +
                                       (z_reconstructed - z_pmt)**2)

                """ time resolution of PMT: """
                # get time resolution of PMT with specific pmtID (pmtID is ascending number from 0 to 17738 (17739 large
                # PMTs)) -> For 20inch PMTs, the pmtID is equal to index of sigma_time_20inch array.
                # check if PMT is 20 inch or 3inch (pmtID < 50000 means 20inch PMT):
                if pmtID < 20000:
                    # 20inch PMT:
                    # get time resolution (sigma) of PMT in ns from array:
                    sigma_pmt_1 = sigma_time_20inch[pmtID]

                    if sigma_pmt_1 < 3:
                        # Hamamatsu PMT:
                        # sigma = TTS / (2*np.sqrt(2*np.log(2))) -> TTS = 7 ns as edge between Hamamatsu and MCP
                        # -> sigma = 3 ns: Hamamatsu if sigma < 3 ns, MCP if sigma > 3 ns:
                        # For Hamamatsu PMTs use sigma_t / TTS of old PmtData_old.root file:
                        sigma_pmt = sigma_pmt_1

                        # Calculate time of flight in ns for the small TTS:
                        time_of_flight = distance_tof / c_effective_smallTTS

                    else:
                        # MCP PMT:
                        # do NOT use sigma_t / TTS from old PmtData_old.root file, because there the TTS is
                        # around 12 ns.
                        # Use TTS of 18 ns and calculate sigma_pmt:
                        TTS_MCP = 18.0
                        sigma_pmt = TTS_MCP / (2 * np.sqrt(2 * np.log(2)))

                        # Calculate time of flight in ns for the large TTS:
                        time_of_flight = distance_tof / c_effective_largeTTS

                elif 20000 < pmtID < 40000:
                    # there are some PMTs with ID around 30000 (user_atmoNC_7.root, event=32: 30637, 30276, 30573,30561,
                    # 30377) -> PMTs with ID above 30000 are Water Pool PMTs!!
                    # go to next photon:
                    continue

                else:
                    # 3inch PMT:
                    sigma_pmt = sigma_time_3inch

                    # Calculate time of flight in ns for the small TTS:
                    time_of_flight = distance_tof / c_effective_smallTTS

                # consider time resolution of PMT by generating normal distributed random number with mu = hittime and
                # sigma = sigma_pmt (only the hittime at the PMT must be smeared, not the time-of-flight):
                hittime_tts = np.random.normal(hittime, sigma_pmt)

                """ calculate the 'real' hittime of the photon in ns: """
                hittime_real = hittime_tts - time_of_flight

                # take only hittimes that are within time window specified by min_time and max_time:
                if min_time <= hittime_real <= max_time:
                    # append real hittime to array:
                    hittime_array.append(hittime_real)

            # build histogram, where hittimes are saved:
            # set bin-edges of hittime histogram in ns:
            bins_hittime = np.arange(min_time, max_time + binwidth, binwidth)
            # build hittime histogram:
            npe_per_hittime, bin_edges_hittime = np.histogram(hittime_array, bins_hittime)

            # before saving the histogram of the hittimes check the reconstructed distance:
            if r_reconstructed < radius_cut:
                # event within 16 m:
                event_position = 16
            elif radius_cut <= r_reconstructed < 17000.0:
                # event between 16 m and 17 m:
                event_position = 17

            """ save time distribution/ pulse shape to file: """
            # save hittime distribution of the event to txt file:
            # build list, where 0th entry is start-hittime in ns, 1st entry is last-hittime in ns, 2nd entry is
            # binwidth in ns and the following entries are nPE of each hittime-bin of whole signal:
            npe_per_hittime_save = [x_reconstructed, y_reconstructed, z_reconstructed]
            npe_per_hittime_save.extend([min_time, max_time, binwidth])
            npe_per_hittime_save.extend(npe_per_hittime)
            np.savetxt(input_path_hittimes + "hittimes/file{0:d}_evt{1:d}_pulse_shape_R{2:d}.txt"
                       .format(index, event, event_position),
                       npe_per_hittime_save, fmt='%1.2f',
                       header="Pulse shape: Number of pe as function of the time "
                              "(time-of-flight correction and TTS smearing) of file user_neutron_{0:d}.root,"
                              "\nevent {1:d}, {2}:"
                              "\ntime window of pulse shape: from {3:.3f} ns to {4:.3f} ns with "
                              "bin-width = {5:0.3f} ns,"
                       .format(index, event, now, min_time, max_time, binwidth))

else:
    # read pulse shapes for each event from txt file:

    # loop over all files in folder input_path_neutron + "hittimes/", that start with 'file' and end with 'R16.txt'
    # (files where hittime distribution is saved, each file is equal to one event):
    for file_neutron in os.listdir(input_path_hittimes + "hittimes/"):

        # if file_neutron.startswith("file") and file_neutron.endswith("R{0:d}.txt".format(int(radius_cut/1000.0))):
        if file_neutron.startswith("file") and file_neutron.endswith("R{0:d}_DCR.txt".format(int(radius_cut/1000.0))):

            # get the file name:
            file_name_neutron = input_path_hittimes + "hittimes/" + file_neutron
            # read txt file:
            npe_from_file = np.loadtxt(file_name_neutron)

            # get reconstructed position in mm:
            x_reconstructed = npe_from_file[0]
            y_reconstructed = npe_from_file[1]
            z_reconstructed = npe_from_file[2]

            # get min_time, max_time and binwidth from txt file and compare it with the values set above:
            min_time_total_txt = npe_from_file[3]
            if min_time != min_time_total_txt:
                sys.exit("ERROR: min_time_total from file differ from the value set in script")
            max_time_total_txt = npe_from_file[4]
            if max_time != max_time_total_txt:
                sys.exit("ERROR: max_time_total from file differ from the value set in script")
            binwidth_txt = npe_from_file[5]
            if binwidth != binwidth_txt:
                sys.exit("ERROR: binwidth from file differ from the value set in script")

            # the rest of npe_from_file is the hittime distribution histogram in nPE per bin. Take only the
            # prompt signal defined by start_time and end_time:
            nPE_per_bin = npe_from_file[6:(int((time_limit_prompt + binwidth + np.abs(min_time)) / binwidth)+3)]

            # calculate the total number of PE of the prompt signal:
            nPE_total = np.sum(nPE_per_bin)

            # convert nPE_total to energy in MeV:
            E_prompt = NC_background_functions.conversion_npe_to_evis(nPE_total)
            # only take energies between 10 MeV and 100 MeV:
            if E_prompt < E_prompt_min or E_prompt > E_prompt_max:
                continue

            # set the time window corresponding to nPE_per_bin:
            time_window = np.arange(min_time, time_limit_prompt + binwidth, binwidth)

            # append E_prompt to array:
            array_E_prompt_neutron.append(E_prompt)
            # increment number of analyzed events:
            number_analyzed_total += 1

            # get index of number_analyzed_array depending on E_prompt (index 0 = 10-11 MeV, index 1 = 11-12 MeV, ...):
            index_energy = int(E_prompt) - 10
            # increment index index_energy of number_analyzed_array by 1:
            np.add.at(number_analyzed_array, [index_energy], 1)

            # calculate TTR value depending on the energy:
            if E_prompt <= 20.0:
                # get TTR and normalized pulse shape of neutron event:
                TTR_neutron, npe_norm_neutron = NC_background_functions.pulse_shape(time_window, nPE_per_bin,
                                                                                    tail_start_10_20, tail_stop_10_20)

                # append E_prompt to array:
                array_E_10_20.append(E_prompt)
                # append TTR value to array:
                array_TTR.append(TTR_neutron)
                array_TTR_10_20.append(TTR_neutron)

                # check if event passes PSD cut:
                if TTR_neutron <= TTR_cut_10_20:
                    # event passes PSD cut:
                    number_pass_PSD_total += 1
                    # increment index index_energy of number_pass_PSD_array by 1:
                    np.add.at(number_pass_PSD_array, [index_energy], 1)

            elif 20.0 < E_prompt <= 30.0:
                # get TTR and normalized pulse shape of neutron event:
                TTR_neutron, npe_norm_neutron = NC_background_functions.pulse_shape(time_window, nPE_per_bin,
                                                                                    tail_start_20_30, tail_stop_20_30)

                # append E_prompt to array:
                array_E_20_30.append(E_prompt)
                # append TTR value to array:
                array_TTR.append(TTR_neutron)
                array_TTR_20_30.append(TTR_neutron)

                # check if event passes PSD cut:
                if TTR_neutron <= TTR_cut_20_30:
                    # event passes PSD cut:
                    number_pass_PSD_total += 1
                    # increment index index_energy of number_pass_PSD_array by 1:
                    np.add.at(number_pass_PSD_array, [index_energy], 1)

            elif 30.0 < E_prompt <= 40.0:
                # get TTR and normalized pulse shape of neutron event:
                TTR_neutron, npe_norm_neutron = NC_background_functions.pulse_shape(time_window, nPE_per_bin,
                                                                                    tail_start_30_40, tail_stop_30_40)

                # append E_prompt to array:
                array_E_30_40.append(E_prompt)
                # append TTR value to array:
                array_TTR.append(TTR_neutron)
                array_TTR_30_40.append(TTR_neutron)

                # check if event passes PSD cut:
                if TTR_neutron <= TTR_cut_30_40:
                    # event passes PSD cut:
                    number_pass_PSD_total += 1
                    # increment index index_energy of number_pass_PSD_array by 1:
                    np.add.at(number_pass_PSD_array, [index_energy], 1)

            else:
                # energies between 40 MeV and 100 MeV
                # get TTR and normalized pulse shape of neutron event:
                TTR_neutron, npe_norm_neutron = NC_background_functions.pulse_shape(time_window, nPE_per_bin,
                                                                                    tail_start_40_100, tail_stop_40_100)

                # append E_prompt to array:
                array_E_40_100.append(E_prompt)
                # append TTR value to array:
                array_TTR.append(TTR_neutron)
                array_TTR_40_100.append(TTR_neutron)

                # check if event passes PSD cut:
                if TTR_neutron <= TTR_cut_40_100:
                    # event passes PSD cut:
                    number_pass_PSD_total += 1
                    # increment index index_energy of number_pass_PSD_array by 1:
                    np.add.at(number_pass_PSD_array, [index_energy], 1)

    """ Build histograms of E_prompt for the different energy regions: """
    binwidth_energy = 1.0
    N_energy, bins_energy = np.histogram(array_E_prompt_neutron, bins=np.arange(E_prompt_min, E_prompt_max+
                                                                                2*binwidth_energy, binwidth_energy))

    """ calculate weights for every energy bin: """
    # get mean value of N_energy (mean number of events per energy bin):
    N_energy_mean = np.mean(N_energy)
    print(N_energy_mean)
    # preallocate array, where all weights are stored (weight per energy bin, in this case index0 = 10-11 MeV, ...):
    weights_array = []
    # preallocate array, where number of events per energy bin are stored to get flat energy distribution:
    N_energy_flat = []
    # loop over N_energy and calculate weight for every bin:
    for index in range(len(N_energy)):
        if N_energy[index] == 0:
            # set weight = 1
            weights_array.append(1.0)
            N_energy_flat.append(0.0)
            continue
        # weight = N_mean per bin/ number events per bin):
        weight = N_energy_mean / N_energy[index]
        # append weight to array:
        weights_array.append(weight)
        # calculate weighted N_energy value:
        N_weighted = N_energy[index] * weight
        # append N_weighted to N_energy_flat:
        N_energy_flat.append(N_weighted)
    # convert weights_array to numpy array:
    weights_array = np.asarray(weights_array)

    """ display energy of prompt signal in histogram: """
    h1 = plt.figure(1, figsize=(11, 6))
    plt.step(bins_energy[:-1], N_energy, color="r", label="entries = {0:d}".format(len(array_E_prompt_neutron)))
    plt.step(bins_energy[:-1], N_energy_flat, color="b", label="weighted spectrum")
    plt.xlim(xmin=E_prompt_min, xmax=E_prompt_max)
    plt.xlabel("Visible energy of prompt signal in MeV")
    plt.ylabel("events per bin (bin width = {0:.1f} MeV)")
    plt.title("Energy spectrum of prompt signal of simulated neutron events")
    plt.legend()
    plt.grid()
    plt.savefig(output_path_neutron + "E_prompt_{0:.0f}MeVto{1:.0f}MeV_R{2:.0f}mm.png"
                .format(E_prompt_min, E_prompt_max, radius_cut))
    plt.close()

    """ calculate PSD efficiencies for fast neutron events: """
    # total efficiency without weight in percent:
    PSD_supp_wo_weight = 100.0 - float(number_pass_PSD_total) / float(number_analyzed_total) * 100.0

    # total efficiency (all events) with weight in percent:
    PSD_supp_w_weight = (100.0 - np.sum(number_pass_PSD_array * weights_array) /
                         np.sum(number_analyzed_array * weights_array) * 100.0)
    # efficiency between 10 and 20 MeV with weight in percent:
    PSD_supp_w_weight_10_20 = (100.0 - np.sum(number_pass_PSD_array[0:10] * weights_array[0:10]) /
                               np.sum(number_analyzed_array[0:10] * weights_array[0:10]) * 100.0)
    # efficiency between 20 and 30 MeV with weight in percent:
    PSD_supp_w_weight_20_30 = (100.0 - np.sum(number_pass_PSD_array[10:20] * weights_array[10:20]) /
                               np.sum(number_analyzed_array[10:20] * weights_array[10:20]) * 100.0)
    # efficiency between 30 and 40 MeV with weight in percent:
    PSD_supp_w_weight_30_40 = (100.0 - np.sum(number_pass_PSD_array[20:30] * weights_array[20:30]) /
                               np.sum(number_analyzed_array[20:30] * weights_array[20:30]) * 100.0)
    # efficiency between 40 and 100 MeV with weight in percent:
    PSD_supp_w_weight_40_100 = (100.0 - np.sum(number_pass_PSD_array[30:] * weights_array[30:]) /
                                np.sum(number_analyzed_array[30:] * weights_array[30:]) * 100.0)

    # print information about PSD efficiency of neutron events:
    print("total PSD suppression:")
    print("NC = {0:.2f} %".format(NC_supp_total))
    print("IBD = {0:.2f} %\n".format(IBD_supp_total))
    print("total number of analyzed fast neutron events = {0}".format(number_analyzed_total))
    print("total number of fast neutron events passed PSD cut = {0}".format(number_pass_PSD_total))
    print("Fast neutron without weighting = {0:.2f} %".format(PSD_supp_wo_weight))
    print("Fast neutron with weighting = {0:.2f} %\n".format(PSD_supp_w_weight))

    print("number of analyzed events (10 to 20 MeV) = {0}".format(np.sum(number_analyzed_array[0:10])))
    print("number events passed PSD (10 to 20 MeV) = {0}".format(np.sum(number_pass_PSD_array[0:10])))
    print("FN suppression with weighting between 10 and 20 MeV = {0:.3f} %\n".format(PSD_supp_w_weight_10_20))

    print("number of analyzed events (20 to 30 MeV) = {0}".format(np.sum(number_analyzed_array[10:20])))
    print("number events passed PSD (20 to 30 MeV) = {0}".format(np.sum(number_pass_PSD_array[10:20])))
    print("FN suppression with weighting between 20 and 30 MeV = {0:.3f} %\n".format(PSD_supp_w_weight_20_30))

    print("number of analyzed events (30 to 40 MeV) = {0}".format(np.sum(number_analyzed_array[20:30])))
    print("number events passed PSD (30 to 40 MeV) = {0}".format(np.sum(number_pass_PSD_array[20:30])))
    print("FN suppression with weighting between 30 and 40 MeV = {0:.3f} %\n".format(PSD_supp_w_weight_30_40))

    print("number of analyzed events (40 to 100 MeV) = {0}".format(np.sum(number_analyzed_array[30:])))
    print("number events passed PSD (40 to 100 MeV) = {0}".format(np.sum(number_pass_PSD_array[30:])))
    print("FN suppression with weighting between 40 and 100 MeV = {0:.3f} %".format(PSD_supp_w_weight_40_100))

    """ 2D histograms of Evis vs TTR values: """
    h2 = plt.figure(2)
    plt.hist2d(array_E_10_20, array_TTR_10_20, bins=100, cmap='rainbow')
    plt.hlines(TTR_cut_10_20, xmin=min(array_TTR_10_20), xmax=max(array_E_10_20),
               label="TTR cut value = {0:.5f}".format(TTR_cut_10_20))
    plt.ylim(ymin=0.0)
    plt.xlabel("visible energy of prompt signal in MeV")
    plt.ylabel("TTR value for tail window from {0:.0f} ns to {1:.0f}".format(tail_start_10_20, tail_stop_10_20))
    plt.title("TTR vs visible energy of simulated neutron events")
    plt.grid()
    plt.legend()
    plt.savefig(output_path_neutron + "TTR_vs_Evis_10_20_{0:.0f}ns_{1:.0f}ns_R{2:.0f}mm.png"
                .format(tail_start_10_20, tail_stop_10_20, radius_cut))
    plt.close()

    h3 = plt.figure(3)
    plt.hist2d(array_E_20_30, array_TTR_20_30, bins=100, cmap='rainbow')
    plt.hlines(TTR_cut_20_30, xmin=min(array_TTR_20_30), xmax=max(array_E_20_30),
               label="TTR cut value = {0:.5f}".format(TTR_cut_20_30))
    plt.ylim(ymin=0.0)
    plt.xlabel("visible energy of prompt signal in MeV")
    plt.ylabel("TTR value for tail window from {0:.0f} ns to {1:.0f}".format(tail_start_20_30, tail_stop_20_30))
    plt.title("TTR vs visible energy of simulated neutron events")
    plt.grid()
    plt.legend()
    plt.savefig(output_path_neutron + "TTR_vs_Evis_20_30_{0:.0f}ns_{1:.0f}ns_R{2:.0f}mm.png"
                .format(tail_start_20_30, tail_stop_20_30, radius_cut))
    plt.close()

    h4 = plt.figure(4)
    plt.hist2d(array_E_30_40, array_TTR_30_40, bins=100, cmap='rainbow')
    plt.hlines(TTR_cut_30_40, xmin=min(array_TTR_30_40), xmax=max(array_E_30_40),
               label="TTR cut value = {0:.5f}".format(TTR_cut_30_40))
    plt.ylim(ymin=0.0)
    plt.xlabel("visible energy of prompt signal in MeV")
    plt.ylabel("TTR value for tail window from {0:.0f} ns to {1:.0f}".format(tail_start_30_40, tail_stop_30_40))
    plt.title("TTR vs visible energy of simulated neutron events")
    plt.grid()
    plt.legend()
    plt.savefig(output_path_neutron + "TTR_vs_Evis_30_40_{0:.0f}ns_{1:.0f}ns_R{2:.0f}mm.png"
                .format(tail_start_30_40, tail_stop_30_40, radius_cut))
    plt.close()

    h5 = plt.figure(5)
    plt.hist2d(array_E_40_100, array_TTR_40_100, bins=100, cmap='rainbow')
    plt.hlines(TTR_cut_40_100, xmin=min(array_TTR_40_100), xmax=max(array_E_40_100),
               label="TTR cut value = {0:.5f}".format(TTR_cut_40_100))
    plt.ylim(ymin=0.0)
    plt.xlabel("visible energy of prompt signal in MeV")
    plt.ylabel("TTR value for tail window from {0:.0f} ns to {1:.0f}".format(tail_start_40_100, tail_stop_40_100))
    plt.title("TTR vs visible energy of simulated neutron events")
    plt.grid()
    plt.legend()
    plt.savefig(output_path_neutron + "TTR_vs_Evis_40_100_{0:.0f}ns_{1:.0f}ns_R{2:.0f}mm.png"
                .format(tail_start_40_100, tail_stop_40_100, radius_cut))
    plt.close()

    """ compare TTR values of neutron events with TTR values of NC and IBD events: """
    # 10 to 20 MeV:
    # load TTR values of NC events that pass all cuts (before PSD):
    array_TTR_NC = np.loadtxt(input_path_TTR + "TTR_IBDlike_{0:.0f}ns_to_{1:.0f}ns_10_20.txt"
                              .format(tail_start_10_20, tail_stop_10_20))
    array_TTR_IBD = np.loadtxt(input_path_TTR + "TTR_IBD_{0:.0f}ns_to_{1:.0f}ns_10_20.txt"
                               .format(tail_start_10_20, tail_stop_10_20))

    # display tot-values for positrons, NC and fast neutron events for the given configuration with efficiencies:
    h6 = plt.figure(6, figsize=(11, 6))
    First_bin = 0.0
    Last_bin = max(array_TTR_NC)
    Bin_width = (Last_bin-First_bin) / 200
    Bins = np.arange(First_bin, Last_bin+Bin_width, Bin_width)
    n_pos_1, bins_pos_1, patches_pos_1 = plt.hist(array_TTR_IBD, bins=Bins, histtype="step", align='mid',
                                                  color="r", linewidth=1.5,
                                                  label="prompt signal of IBD events after all cuts "
                                                        "(entries = {0:d})".format(len(array_TTR_IBD)))
    n_NC_1, bins_NC_1, patches_NC_1 = plt.hist(array_TTR_NC, bins=Bins, histtype="step", align='mid',
                                               color="b", linewidth=1.5,
                                               label="prompt signal of NC events that mimic IBD signal "
                                                     "(entries = {0:d})".format(len(array_TTR_NC)))
    n_n_1, bins_n_1, patches_n_1 = plt.hist(array_TTR_10_20, bins=Bins, histtype="step", align='mid',
                                            color="g", linewidth=1.5,
                                            label="prompt signal of neutrons representing fast neutron "
                                                  "events (entries = {0:d})"
                                            .format(len(array_TTR_10_20)))
    plt.vlines(TTR_cut_10_20, 0, max(n_NC_1)+max(n_NC_1)/10, colors="k", linestyles="--",
               label="$\\epsilon_{IBD}$ = "+"{0:0.2f} %\n".format(PSD_IBD_supp_10_20) +
                     "$\\epsilon_{NC}$ = "+"{0:0.2f} %\n".format(PSD_NC_supp_10_20) +
                     "$\\epsilon_{fastN}$ = "+"{0:0.2f} %\n".format(PSD_supp_w_weight_10_20) +
                     "tot value = {0:.5f}".format(TTR_cut_10_20))
    plt.xlabel("tail-to-total ratio")
    plt.ylabel("events")
    plt.title("Tail-to-total ratio of prompt signals of IBD, NC and fast neutron events\n" +
              "with energies between 10 MeV and 20 MeV (tail window {0:0.1f} ns to {1:0.1f} ns)"
              .format(tail_start_10_20, tail_stop_10_20))
    plt.legend()
    plt.grid()
    plt.savefig(output_path_neutron + "TTR_{0:.0f}nsto{1:.0f}ns_PosNCfastN_R{2:.0f}mm_10_20.png"
                .format(tail_start_10_20, tail_stop_10_20, radius_cut))
   # plt.close()

    # 20 to 30 MeV:
    # load TTR values of NC events that pass all cuts (before PSD):
    array_TTR_NC = np.loadtxt(input_path_TTR + "TTR_IBDlike_{0:.0f}ns_to_{1:.0f}ns_20_30.txt"
                              .format(tail_start_20_30, tail_stop_20_30))
    array_TTR_IBD = np.loadtxt(input_path_TTR + "TTR_IBD_{0:.0f}ns_to_{1:.0f}ns_20_30.txt"
                               .format(tail_start_20_30, tail_stop_20_30))

    # display tot-values for positrons, NC and fast neutron events for the given configuration with efficiencies:
    h7 = plt.figure(7, figsize=(11, 6))
    First_bin = 0.0
    Last_bin = 0.05
    Bin_width = (Last_bin - First_bin) / 200
    Bins = np.arange(First_bin, Last_bin + Bin_width, Bin_width)
    n_pos_1, bins_pos_1, patches_pos_1 = plt.hist(array_TTR_IBD, bins=Bins, histtype="step", align='mid',
                                                  color="r", linewidth=1.5,
                                                  label="prompt signal of IBD events after all cuts "
                                                        "(entries = {0:d})".format(len(array_TTR_IBD)))
    n_NC_1, bins_NC_1, patches_NC_1 = plt.hist(array_TTR_NC, bins=Bins, histtype="step", align='mid',
                                               color="b", linewidth=1.5,
                                               label="prompt signal of NC events that mimic IBD signal "
                                                     "(entries = {0:d})".format(len(array_TTR_NC)))
    n_n_1, bins_n_1, patches_n_1 = plt.hist(array_TTR_20_30, bins=Bins, histtype="step", align='mid',
                                            color="g", linewidth=1.5,
                                            label="prompt signal of neutrons representing fast neutron "
                                                  "events (entries = {0:d})"
                                            .format(len(array_TTR_20_30)))
    plt.vlines(TTR_cut_20_30, 0, max(n_pos_1) + max(n_pos_1) / 10, colors="k", linestyles="--",
               label="$\\epsilon_{IBD}$ = " + "{0:0.2f} %\n".format(PSD_IBD_supp_20_30) +
                     "$\\epsilon_{NC}$ = " + "{0:0.2f} %\n".format(PSD_NC_supp_20_30) +
                     "$\\epsilon_{fastN}$ = " + "{0:0.2f} %\n".format(PSD_supp_w_weight_20_30) +
                     "tot value = {0:.5f}".format(TTR_cut_20_30))
    plt.xlabel("tail-to-total ratio")
    plt.ylabel("events")
    plt.title("Tail-to-total ratio of prompt signals of IBD, NC and fast neutron events\n" +
              "with energies between 20 MeV and 30 MeV (tail window {0:0.1f} ns to {1:0.1f} ns)"
              .format(tail_start_20_30, tail_stop_20_30))
    plt.legend()
    plt.grid()
    plt.savefig(output_path_neutron + "TTR_{0:.0f}nsto{1:.0f}ns_PosNCfastN_R{2:.0f}mm_20_30.png"
                .format(tail_start_20_30, tail_stop_20_30, radius_cut))
 #   plt.close()

    # 30 to 40 MeV:
    # load TTR values of NC events that pass all cuts (before PSD):
    array_TTR_NC = np.loadtxt(input_path_TTR + "TTR_IBDlike_{0:.0f}ns_to_{1:.0f}ns_30_40.txt"
                              .format(tail_start_30_40, tail_stop_30_40))
    array_TTR_IBD = np.loadtxt(input_path_TTR + "TTR_IBD_{0:.0f}ns_to_{1:.0f}ns_30_40.txt"
                               .format(tail_start_30_40, tail_stop_30_40))

    # display tot-values for positrons, NC and fast neutron events for the given configuration with efficiencies:
    h8 = plt.figure(8, figsize=(11, 6))
    First_bin = 0.0
    Last_bin = 0.05
    Bin_width = (Last_bin - First_bin) / 200
    Bins = np.arange(First_bin, Last_bin + Bin_width, Bin_width)
    n_pos_1, bins_pos_1, patches_pos_1 = plt.hist(array_TTR_IBD, bins=Bins, histtype="step", align='mid',
                                                  color="r", linewidth=1.5,
                                                  label="prompt signal of IBD events after all cuts "
                                                        "(entries = {0:d})".format(len(array_TTR_IBD)))
    n_NC_1, bins_NC_1, patches_NC_1 = plt.hist(array_TTR_NC, bins=Bins, histtype="step", align='mid',
                                               color="b", linewidth=1.5,
                                               label="prompt signal of NC events that mimic IBD signal "
                                                     "(entries = {0:d})".format(len(array_TTR_NC)))
    n_n_1, bins_n_1, patches_n_1 = plt.hist(array_TTR_30_40, bins=Bins, histtype="step", align='mid',
                                            color="g", linewidth=1.5,
                                            label="prompt signal of neutrons representing fast neutron "
                                                  "events (entries = {0:d})"
                                            .format(len(array_TTR_30_40)))
    plt.vlines(TTR_cut_30_40, 0, max(n_pos_1) + max(n_pos_1) / 10, colors="k", linestyles="--",
               label="$\\epsilon_{IBD}$ = " + "{0:0.2f} %\n".format(PSD_IBD_supp_30_40) +
                     "$\\epsilon_{NC}$ = " + "{0:0.2f} %\n".format(PSD_NC_supp_30_40) +
                     "$\\epsilon_{fastN}$ = " + "{0:0.2f} %\n".format(PSD_supp_w_weight_30_40) +
                     "tot value = {0:.5f}".format(TTR_cut_30_40))
    plt.xlabel("tail-to-total ratio")
    plt.ylabel("events")
    plt.title("Tail-to-total ratio of prompt signals of IBD, NC and fast neutron events\n" +
              "with energies between 30 MeV and 40 MeV (tail window {0:0.1f} ns to {1:0.1f} ns)"
              .format(tail_start_30_40, tail_stop_30_40))
    plt.legend()
    plt.grid()
    plt.savefig(output_path_neutron + "TTR_{0:.0f}nsto{1:.0f}ns_PosNCfastN_R{2:.0f}mm_30_40.png"
                .format(tail_start_30_40, tail_stop_30_40, radius_cut))
 #   plt.close()

    # 40 to 100 MeV:
    # load TTR values of NC events that pass all cuts (before PSD):
    array_TTR_NC = np.loadtxt(input_path_TTR + "TTR_IBDlike_{0:.0f}ns_to_{1:.0f}ns_40_100.txt"
                              .format(tail_start_40_100, tail_stop_40_100))
    array_TTR_IBD = np.loadtxt(input_path_TTR + "TTR_IBD_{0:.0f}ns_to_{1:.0f}ns_40_100.txt"
                               .format(tail_start_40_100, tail_stop_40_100))

    # display tot-values for positrons, NC and fast neutron events for the given configuration with efficiencies:
    h9 = plt.figure(9, figsize=(11, 6))
    First_bin = 0.0
    Last_bin = 0.05
    Bin_width = (Last_bin - First_bin) / 200
    Bins = np.arange(First_bin, Last_bin + Bin_width, Bin_width)
    n_pos_1, bins_pos_1, patches_pos_1 = plt.hist(array_TTR_IBD, bins=Bins, histtype="step", align='mid',
                                                  color="r", linewidth=1.5,
                                                  label="prompt signal of IBD events after all cuts "
                                                        "(entries = {0:d})".format(len(array_TTR_IBD)))
    n_NC_1, bins_NC_1, patches_NC_1 = plt.hist(array_TTR_NC, bins=Bins, histtype="step", align='mid',
                                               color="b", linewidth=1.5,
                                               label="prompt signal of NC events that mimic IBD signal "
                                                     "(entries = {0:d})".format(len(array_TTR_NC)))
    n_n_1, bins_n_1, patches_n_1 = plt.hist(array_TTR_40_100, bins=Bins, histtype="step", align='mid',
                                            color="g", linewidth=1.5,
                                            label="prompt signal of neutrons representing fast neutron "
                                                  "events (entries = {0:d})"
                                            .format(len(array_TTR_40_100)))
    plt.vlines(TTR_cut_40_100, 0, max(n_pos_1) + max(n_pos_1) / 10, colors="k", linestyles="--",
               label="$\\epsilon_{IBD}$ = " + "{0:0.2f} %\n".format(PSD_IBD_supp_40_100) +
                     "$\\epsilon_{NC}$ = " + "{0:0.2f} %\n".format(PSD_NC_supp_40_100) +
                     "$\\epsilon_{fastN}$ = " + "{0:0.2f} %\n".format(PSD_supp_w_weight_40_100) +
                     "tot value = {0:.5f}".format(TTR_cut_40_100))
    plt.xlabel("tail-to-total ratio")
    plt.ylabel("events")
    plt.title("Tail-to-total ratio of prompt signals of IBD, NC and fast neutron events\n" +
              "with energies between 40 MeV and 100 MeV (tail window {0:0.1f} ns to {1:0.1f} ns)"
              .format(tail_start_40_100, tail_stop_40_100))
    plt.legend()
    plt.grid()
    plt.savefig(output_path_neutron + "TTR_{0:.0f}nsto{1:.0f}ns_PosNCfastN_R{2:.0f}mm_40_100.png"
                .format(tail_start_40_100, tail_stop_40_100, radius_cut))
  #  plt.close()
    plt.show()






