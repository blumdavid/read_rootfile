""" Script to analyze the prompt energy cut and the delayed energy cut (time cut, delayed energy cut, multiplicity cut
    and distance cut) and their efficiencies (either of real IBD events (positron + neutron) or NC events):

    Version 3:

    -   the prompt energy and delayed energy is converted to MeV, but NOT smeared with energy resolution, because
    the energy resolution is already considered in the number of PE

    -   the delayed energy cut is done on MeV and not on number of PE

    Version 2 (v2) means, that:
    -   prompt energy cut is calculated independently of the other cuts
    -   filenumber, evtID and Q_smeared of each event passing the prompt energy cut is stored in file
    -   all values of Q_smeared are stored in a histogram with the prompt energy cut parameters
        E_min and E_max

    -   delayed cut is calculated independently of other cuts
    -   all values of begin_pulse and end_pulse are stored in a histogram with the time cut parameters t_min and t_max
    -   all values of Q_smeared of delayed signal are stored in a histogram with the delayed energy cut parameters
        PE_min and PE_max of delayed signal
    -   all values of N_PEsmeared are stored in a histogram with the neutron multiplicity parameters
    -   all values of distance_reco are stored in a histogram with the distance cut parameters R_distance
    -   filenumber and evtID of each event passing the whole delayed cut is stored in file

    -   save prompt signal of pulse shape of all events in png and txt file

    1.  set the cut parameters (prompt_energy_min, prompt_energy_max, time_cut_min, time_cut_max, multiplicity,
        min_PE_delayed, max_PE_delayed, distance_cut)

    2.  read all NC events (user_atmoNC_0.root to user_atmoNC_999.root, each file contains 100 events.
        Therefore, 100000 events are analyzed.)
        or read all IBD events (user_IBD_hepevt_0.root to user_IBD_hepevt_199.root, each file
        contains 100 events.)

    3.  calculate corrected time distribution (pulse shape) of each event to be able to do prompt energy cut and
        delayed cut

    4.  analyze pulse shape:
        4.1.    take prompt signal and do prompt energy cut

        4.2.    take delayed signal and do delayed cut (time cut, multiplicity cut, delayed energy cut, distance cut,
        volume cut on delayed signal)

"""
import datetime
import ROOT
import sys
import NC_background_functions
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit


def gaus(x, a, x0, sigma):
    """
    gaussian fit function
    :param x: values on x axis
    :param a: normalization factor
    :param x0: mu
    :param sigma: standard deviation
    :return:
    """

    return a * np.exp(-(x-x0)**2 / (2 * sigma**2))


# get the date and time, when the script was run:
date = datetime.datetime.now()
now = date.strftime("%Y-%m-%d %H:%M")

""" set flag, if hittimes must be calculated or read from file: """
flag_read_hittime_from_file = True

""" file information: """
# first file to be read:
start_number = 0
# last file to be read:
stop_number = 999
# number of entries in the input files:
Number_entries_input = 100
# set string, that define, if positrons or NC events are analyzed:
event_type = "atmoNC"
# set integer, that defines the number of the analyzed files (hittime analysis needs much time, therefore separate
# analysis in different parts)
analysis_part = 100
# set the path of the input root files:
if event_type == "atmoNC":
    input_path = "/local/scratch1/pipc51/astro/blum/detsim_output_data/user_atmoNC_"
elif event_type == "IBD":
    input_path = "/local/scratch1/pipc51/astro/blum/IBD_hepevt/user_IBD_hepevt_"

# set the path of the output, where the txt file with visible energy of events that pass cuts is saved:
if event_type == "atmoNC":
    output_path = "/home/astro/blum/juno/atmoNC/data_NC/output_detsim_v2/"
elif event_type == "IBD":
    output_path = "/home/astro/blum/juno/IBD_events/"

""" define parameters depending on the cuts: """
# fiducial volume cut on delayed signal position in mm:
R_cut_delayed_mm = 17700.0
# minimum of prompt energy in MeV:
prompt_energy_min = 10.0
# maximum of prompt energy in MeV:
prompt_energy_max = 100.0
# set time window of whole signal in ns:
min_time_total = -50
max_time_total = 2000000
# time window of delayed signal in ns:
time_cut_min = 1000
time_cut_max = 1000000
# Set bin-width of hittime histogram in ns:
binwidth = 5.0
# set multiplicity cut (only one delayed signal in event):
multiplicity = 1
# min and max of energy in MeV for delayed energy cut:
min_PE_delayed_MeV = 1.80
max_PE_delayed_MeV = 2.55
# minimum number of PE that defines a delayed signal (threshold, nPE_delayed > min_PE -> delayed signal):
min_PE = 250
# Set threshold of number of PE per bin for possible delayed signal (bin-width = 5 ns):
threshold1_del = 50
# set threshold2 of number of PEs per bin (signal peak is summed as long as nPE is above threshold2):
threshold2_del = 0
# distance cut between prompt and delayed signal in mm:
distance_cut = 500

# set the path, where the numbers.txt files and histograms for the delayed cut (time, multiplicity, delayed E,
# distance and delayed volume cut) should be stored:
if event_type == "atmoNC":
    output_path_del = (output_path + "delayed_cut_{0:.0f}nsto{1:.0f}ms_mult{2:d}_{3:.0f}keVto{4:.0f}keV_dist{5:.0f}mm_"
                                     "R{6:.0f}mm/".format(time_cut_min, time_cut_max/1000000.0, multiplicity,
                                                          min_PE_delayed_MeV*1000, max_PE_delayed_MeV*1000,
                                                          distance_cut, R_cut_delayed_mm))
elif event_type == "IBD":
    output_path_del = (output_path + "delayed_cut_{0:.0f}nsto{1:.0f}ms_mult{2:d}_{3:.0f}keVto{4:.0f}keV_dist{5:.0f}mm_"
                                     "R{6:.0f}mm/".format(time_cut_min, time_cut_max/1000000.0, multiplicity,
                                                          min_PE_delayed_MeV*1000, max_PE_delayed_MeV*1000,
                                                          distance_cut, R_cut_delayed_mm))

# number of analyzed events (total number of events):
number_total_events = (stop_number + 1 - start_number) * Number_entries_input

""" preallocate variables: """
# number of events without initial particles:
number_without_particles = 0

""" prompt energy """
# number of events, where Q_converted pass prompt energy cut:
number_prompt_energy_pass_reco = 0
# number of events, where Q_converted is rejected by prompt energy cut:
number_prompt_energy_rejected_reco = 0
# number of events, where Q_converted is rejected by prompt energy cut (E < prompt_energy_min):
number_prompt_energy_rejected_min_reco = 0
# number of events, where Q_converted is rejected by prompt energy cut (E > prompt_energy_max):
number_prompt_energy_rejected_max_reco = 0
# array, where filenumber of events that pass the prompt energy cut (Q_converted, real) are stored:
array_filenumber_prompt_energy = []
# array, where corresponding evtID of events that pass the prompt energy cut (Q_converted, real) are stored:
array_evtID_prompt_energy = []
# array, where corresponding prompt energy of events that pass the prompt energy cut (Q_converted, real) are stored in
# MeV:
array_Evis_prompt_energy = []
# array, where all Q_converted values are stored to build a histogram:
array_Q_converted_prompt_energy = []
# efficiency of conversion of prompt energy from nPE to MeV (efficiency = N_before_conversion / N_after_conversion):
efficiency_conversion = 0.9879

""" delayed cut: """
# number of events that pass the delayed cut:
number_delayed_cut_pass = 0
# number of events that are rejected by delayed cut:
number_delayed_cut_rejected = 0
# array, where filenumber of events that pass delayed cut (real) are stored:
array_filenumber_delayed_cut = []
# array, where corresponding evtID of events that pass delayed cut (real) are stored:
array_evtID_delayed_cut = []
# array, where filenumber of events that pass delayed cut (ideal) are stored:
array_filenumber_delayed_cut_MCtruth = []
# array, where corresponding evtID of events that pass delayed cut (ideal) are stored:
array_evtID_delayed_cut_MCtruth = []

""" time """
# number of events, where begin/end time pass time cut:
number_time_pass_reco = 0
# number of events, where begin/end time is rejected by time cut:
number_time_rejected_reco = 0
# number of events, where NeutronCaptureT pass time cut:
number_time_pass_MC = 0
# number of events, where NeutronCaptureT is rejected by time cut:
number_time_rejected_MC = 0
# number of events, where begin/end time pass time cut, but NeutronCaptureT is rejected (events counted too much):
number_time_toomuch = 0
# number of events, where begin/end time is rejected by time cut, but NeutronCaptureT pass (events counted too less):
number_time_tooless = 0
# array, where all begin times are stored to build 2d histogram:
array_begin_times = []
# array, where all end times are stored to build 2d histogram:
array_end_times = []
# array, where all NeutronCaptureT are stored:
array_NeutronCaptureT = []

""" multiplicity """
# number of events that pass multiplicity cut (and time cut) (only 1 delayed signal (begin/end time) in delayed time
# window):
number_multiplicity_pass_reco = 0
# number of events that are rejected by multiplicity cut (more than 1 delayed signal (begin/end time) in delayed time
# window):
number_multiplicity_rejected_reco = 0
# number of events that would pass multiplicity cut (only 1 delayed signal (NeutronCaptureT) in delayed time window):
number_multiplicity_pass_MC = 0
# number of events that would be rejected by multiplicity cut (more than 1 delayed signal (NeutronCaptureT) in delayed
# time window):
number_multiplicity_rejected_MC = 0
# number of events, where only 1 delayed signal from begin/end,
# but 0 or more than 1 delayed signal from NeutronCaptureT:
number_multiplicity_toomuch = 0
# number of events, where more than 1 delayed signal from begin/end, but only 1 delayed signal from NeutronCaptureT:
number_multiplicity_tooless = 0
# array, where all numbers of delayed signals inside delayed time window are stored:
array_multiplicity = []

""" delayed energy """
# number of events, where Evis_converted pass delayed energy cut (and time cut and multiplicity cut):
number_delayed_energy_pass_reco = 0
# number of events, where Evis_converted is rejected by delayed energy cut:
number_delayed_energy_rejected_reco = 0
# array, where all Evis_converted values are stored:
array_Evis_delayed_energy = []

""" distance """
# number of events, where distance_reco pass distance cut (and time, multiplicity and delayed energy cut):
number_distance_cut_pass_reco = 0
# number of events, where distance_reco is rejected by distance cut:
number_distance_cut_rejected_reco = 0
# number of events, where distance_MC pass distance cut:
number_distance_cut_pass_MC = 0
# number of events, where distance_MC is rejected by distance cut:
number_distance_cut_rejected_MC = 0
# number of events, where reconstructed distance pass the cut, but "real" (MC truth) distance is rejected:
number_distance_cut_toomuch = 0
# number of events, where reconstructed distance is rejected by cut, but "real" distance pass:
number_distance_cut_tooless = 0
# array, where all distance_reco values are stored to build histogram:
array_distance_reco = []
# array, where filenumber of events that pass distance cut and all delayed cuts before (real) are stored:
array_filenumber_distance_cut = []
# array, where corresponding evtID of events that pass distance cut and all delayed cuts before (real) are stored:
array_evtID_distance_cut = []

""" volume cut on delayed position: """
# number of events that pass volume cut on delayed recon. position:
number_volume_pass_delayed_reco = 0
# number of events that are rejected by volume cut on delayed recon. position:
number_volume_rejected_delayed_reco = 0
# number of events that pass volume cut on delayed MC truth position:
number_volume_pass_delayed_MC = 0
# number of events that are rejected by volume cut on delayed MC position:
number_volume_rejected_delayed_MC = 0
# number of leak-in events (initial position outside, but recon. position inside fiducial volume):
number_volume_leak_in = 0
# number of leak-out events (initial position inside, but recon. position outside fiducial volume):
number_volume_leak_out = 0
# array, where delayed recon. position is stored:
array_volume_delayed_position = []

""" load position of the PMTs and corresponding PMT ID from file PMT_position.root: """
file_PMT_position = "/home/astro/blum/juno/atmoNC/PMT_information/PMT_position.root"
# array with PMT ID and corresponding x, y, z position in mm:
pmtID_pos_file, x_pos_pmt, y_pos_pmt, z_pos_pmt = NC_background_functions.get_pmt_position(file_PMT_position)

""" load 'time resolution' in ns of the 20 inch PMTs and corresponding PMT ID from file PmtData.root: """
# The simulation is done with version J18v1r1-Pre1. In this version 17738 20 inch PMTs are set in the detector with TTS
# for Hamamatsu ~ 2.7 ns and MCP ~ 12 ns (see PmtData_old.root).
# The PMT positions for this version are saved in file PMT_position.root.
# In the latest version J19v1r1-Pre4, the number of 20 inch PMTs and their TTS have been updated (Hamamatsu ~ 2.7 ns,
# MCP ~ 18 ns, TTS of MCP larger!). Therefore use TTS = 18 ns for all MCP PMTs!
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

# loop over the files that are read:
for filenumber in range(start_number, stop_number+1):

    # file name of the input file:
    input_name = input_path + "{0:d}.root".format(filenumber)

    print(input_name)

    # load user_atmoNC_index.root file:
    rfile = ROOT.TFile(input_name)
    # get the "evt"-TTree from the TFile:
    rtree_evt = rfile.Get("evt")
    # get the "geninfo"-TTree from the TFile:
    rtree_geninfo = rfile.Get("geninfo")
    # get the "prmtrkdep"-TTree from the TFile:
    rtree_prmtrkdep = rfile.Get("prmtrkdep")
    # get the "nCapture"-TTree from TFile:
    rtree_ncapture = rfile.Get("nCapture")

    # get number of events in evt tree:
    number_events_evt = rtree_evt.GetEntries()
    # get the number of events in the geninfo Tree:
    number_events_geninfo = rtree_geninfo.GetEntries()
    # get the number of events in the prmtrkdep tree:
    number_events_prmtrkdep = rtree_prmtrkdep.GetEntries()
    # get number of events in nCapture tree:
    number_events_ncapture = rtree_ncapture.GetEntries()

    if (number_events_evt == number_events_geninfo and number_events_evt == number_events_ncapture and
            number_events_evt == number_events_prmtrkdep):
        number_events = number_events_geninfo
    else:
        sys.exit("ERROR: number of events in the Trees are NOT equal!!")

    # check if number_events is equal to number_entries_input (if not, the detector simulation was incorrect!!):
    if number_events != Number_entries_input:
        sys.exit("ERROR: number of events {0:d} are not equal to {1:d} -> Detector Simulation not correct!"
                 .format(number_events, Number_entries_input))

    # loop over every event in the file:
    for event in range(number_events):

        """ preallocate variables corresponding to this event: """
        # array where corrected hittimes are stored (hittimes in ns):
        hittime_array = []
        # sum of Qedep of all initial particles (in MeV):
        Qedep_sum = 0
        # array where NeutronCaptureT of each NeutronN is stored (in ns):
        ncaptureT_array = np.array([])

        """ get reconstructed position of event """
        # get the current event in the TTree:
        rtree_geninfo.GetEntry(event)
        # get the value of the event ID:
        event_id = int(rtree_geninfo.GetBranch('evtID').GetLeaf('evtID').GetValue())

        # get the value of the number of initial particles:
        nInitParticles_geninfo = int(rtree_geninfo.GetBranch('nInitParticles').GetLeaf('nInitParticles').GetValue())

        # check if there are initial particles in this event:
        if nInitParticles_geninfo < 1:
            # there are NO particles in the event:
            number_without_particles += 1
            # go to next event
            continue

        # get initial position of the first initial particle:
        init_x = float(rtree_geninfo.GetBranch('InitX').GetLeaf('InitX').GetValue(0))
        init_y = float(rtree_geninfo.GetBranch('InitY').GetLeaf('InitY').GetValue(0))
        init_z = float(rtree_geninfo.GetBranch('InitZ').GetLeaf('InitZ').GetValue(0))

        # check, if initial positions of the other initial particles are equal to position of first initial:
        for index in range(1, nInitParticles_geninfo):
            # get initial x position:
            init_x_test = float(rtree_geninfo.GetBranch('InitX').GetLeaf('InitX').GetValue(index))
            # get initial y position:
            init_y_test = float(rtree_geninfo.GetBranch('InitY').GetLeaf('InitY').GetValue(index))
            # get initial z position:
            init_z_test = float(rtree_geninfo.GetBranch('InitZ').GetLeaf('InitZ').GetValue(index))

            if init_x != init_x_test or init_y != init_y_test or init_z != init_z_test:
                sys.exit("ERROR: initial positions are NOT equal (event = {0:d}, filenumber = {1:d})"
                         .format(event, filenumber))

        if not flag_read_hittime_from_file:
            # do not read hittime from file, but calculate it:

            # get current event in TTree:
            rtree_prmtrkdep.GetEntry(event)
            # get number of initial particles in prmtrkdep:
            nInitParticles_prmtrkdep = int(
                rtree_geninfo.GetBranch('nInitParticles').GetLeaf('nInitParticles').GetValue())

            # get quenched deposit energy of all initial particles in MeV:
            for index in range(nInitParticles_prmtrkdep):
                Qedep = float(rtree_prmtrkdep.GetBranch('Qedep').GetLeaf('Qedep').GetValue(index))
                Qedep_sum += Qedep

            # Calculate vertex smeared position.
            if Qedep_sum != 0:
                # Smear initial x,y and z position with function position_smearing(). (returns reconstructed position
                # in mm):
                x_reconstructed = NC_background_functions.position_smearing(init_x, Qedep_sum)
                y_reconstructed = NC_background_functions.position_smearing(init_y, Qedep_sum)
                z_reconstructed = NC_background_functions.position_smearing(init_z, Qedep_sum)

            else:
                # Qedep_sum = 0, use initial position:
                x_reconstructed = init_x
                y_reconstructed = init_y
                z_reconstructed = init_z

            # calculate reconstructed distance to detector center in mm:
            r_reconstructed = np.sqrt(x_reconstructed ** 2 + y_reconstructed ** 2 + z_reconstructed ** 2)

            """ calculate the corrected time distribution/pulse shape (time of flight correction with reconstructed 
            position and time smearing with TTS for each hit): """
            # get current event of evt-tree:
            rtree_evt.GetEntry(event)
            # get number of photons of this event:
            n_photons = int(rtree_evt.GetBranch('nPhotons').GetLeaf('nPhotons').GetValue())

            # loop over every photon in the event:
            for index in range(n_photons):

                # get nPE for this photon:
                n_pe = int(rtree_evt.GetBranch('nPE').GetLeaf('nPE').GetValue(index))
                # check, if photon produces only 1 PE:
                if n_pe != 1:
                    print("{1:d} PE for 1 photon in event {0:d} in file user_atmoNC_{2:d}.root"
                          .format(event_id, n_pe, index))

                # get the pmtID of the hit PMT:
                pmtID = int(rtree_evt.GetBranch('pmtID').GetLeaf('pmtID').GetValue(index))

                # get hittime of this photon:
                hit_time = float(rtree_evt.GetBranch('hitTime').GetLeaf('hitTime').GetValue(index))

                # get position of the PMT with specific pmtID (pmtID is ascending number from 0 to 17738
                # (17739 large PMTs) and from 300000 to 336571 (36572 small PMTs)).
                # For large PMTs -> For 20inch PMTs, the pmtID is equal to index of x,y,z_pos_pmt array.
                # For small PMTs -> For 3inch PMTs, the pmtID - (300000 - 17739) is equal to index of x,y,z_pos_pmt
                # array.
                # check if PMT is 20 inch or 3inch (pmtID < 20000 means 20inch PMT):
                if pmtID < 20000:
                    # 20inch PMT:
                    # get PMT position in mm from arrays:
                    x_pmt = x_pos_pmt[pmtID]
                    y_pmt = y_pos_pmt[pmtID]
                    z_pmt = z_pos_pmt[pmtID]

                elif 20000 < pmtID < 40000:
                    # there are some PMTs with ID around 30000 (user_atmoNC_7.root, event=32: 30637, 30276, 30573,
                    # 30561, 30377) -> PMTs with ID above 30000 are Water Pool PMTs!!
                    # go to next photon:
                    continue

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
                # check if PMT is 20 inch or 3inch (pmtID < 20000 means 20inch PMT):
                if pmtID < 20000:
                    # 20inch PMT:
                    # get time resolution (sigma) of PMT in ns from array:
                    sigma_pmt_1 = sigma_time_20inch[pmtID]
                    # check if PMT is Hamamatsu or MCP:
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

                # consider time resolution of PMT by generating normal distributed random number with mu = hit_time and
                # sigma = sigma_pmt (only the hit_time at the PMT must be smeared, not the time-of-flight):
                hittime_tts = np.random.normal(hit_time, sigma_pmt)

                # calculate the 'real' hittime of the photon in ns:
                hittime_real = hittime_tts - time_of_flight

                # append hittime to array:
                hittime_array.append(hittime_real)

            # hittime_array contains now the corrected hittimes of all PMTs in ns!
            # now all cuts can be done with hittime_array!

            # build histogram, where hittimes are saved:
            # set bin-edges of hittime histogram in ns in whole time window:
            bins_hittime = np.arange(min_time_total, max_time_total + 2 * binwidth, binwidth)
            # build hittime histogram:
            npe_per_hittime, bin_edges_hittime = np.histogram(hittime_array, bins_hittime)

        else:
            # read corrected hittime from txt file:
            npe_from_file = np.loadtxt(output_path + "hittimes/file{0:d}_evt{1:d}_pulse_shape.txt"
                                       .format(filenumber, event_id))

            # get reconstructed x, y and z position in mm:
            x_reconstructed = npe_from_file[0]
            y_reconstructed = npe_from_file[1]
            z_reconstructed = npe_from_file[2]

            # calculate reconstructed distance to detector center in mm:
            r_reconstructed = np.sqrt(x_reconstructed ** 2 + y_reconstructed ** 2 + z_reconstructed ** 2)

            # get min_time_total, max_time_total and binwidth from txt file and compare it with the values set above:
            min_time_total_txt = npe_from_file[3]
            if min_time_total != min_time_total_txt:
                sys.exit("ERROR: min_time_total from file differ from the value set in script")

            max_time_total_txt = npe_from_file[4]
            if max_time_total != max_time_total_txt:
                sys.exit("ERROR: max_time_total from file differ from the value set in script")

            binwidth_txt = npe_from_file[5]
            if binwidth != binwidth_txt:
                sys.exit("ERROR: binwidth from file differ from the value set in script")

            # set bin-edges of hittime histogram in ns in whole time window:
            bins_hittime = np.arange(min_time_total, max_time_total + 2 * binwidth, binwidth)

            # set bin_edges_hittime:
            bin_edges_hittime = bins_hittime

            # get hittime histogram from file:
            npe_per_hittime = npe_from_file[6:]

        """ analyze prompt signal: """
        # get index of bins_hittime corresponding to min_time (should be index = 0):
        index_min_hittime_prompt = 0

        # Where does prompt signal end?
        # get index of bins_hittime, where prompt time window ends/delayed time window starts:
        index_time_cut_min = int((time_cut_min + np.abs(min_time_total)) / binwidth)
        # check if npe_per_hittime (and the following two bins) are 0 for this index:
        if (npe_per_hittime[index_time_cut_min] == npe_per_hittime[index_time_cut_min + 1]
                == npe_per_hittime[index_time_cut_min + 2] == 0):
            # prompt signal already 0:
            index_max_hittime_prompt = index_time_cut_min
        else:
            # prompt signal not yet 0.
            # loop over npe_per_hittime from index_time_cut_min until npe_per_hittime (and the following two bins)
            # are 0:
            for index in range(index_time_cut_min, index_time_cut_min + 500):
                if npe_per_hittime[index] == npe_per_hittime[index + 1] == npe_per_hittime[index + 2] == 0:
                    index_max_hittime_prompt = index
                    break

        # calculate nPE as function of hittime only for prompt time window (last index should be included):
        npe_per_hittime_prompt = npe_per_hittime[index_min_hittime_prompt:index_max_hittime_prompt + 1]
        # bin edges of hittime histogram only for prompt time window:
        bins_hittime_prompt = bin_edges_hittime[index_min_hittime_prompt:index_max_hittime_prompt + 1]

        # save pulse shape of prompt signal to txt file (this is used in analyze_PSD_cut_v2.py for pulse shape
        # analysis):
        npe_per_hittime_save_prompt = [x_reconstructed, y_reconstructed, z_reconstructed]
        npe_per_hittime_save_prompt.extend([min_time_total, 1000.0+20.0, binwidth])
        npe_per_hittime_save_prompt.extend(npe_per_hittime[index_min_hittime_prompt:
                                                           (int((1000.0 + binwidth + np.abs(min_time_total))
                                                                / binwidth)+3)])
        np.savetxt(output_path + "hittimes/file{0:d}_evt{1:d}_prompt_signal.txt".format(filenumber, event),
                   npe_per_hittime_save_prompt, fmt='%1.2f',
                   header="Pulse shape of prompt signal: Number of pe as function of the time "
                          "(time-of-flight correction and TTS smearing) of file user_{6}_{0:d}.root,"
                          "\nevent {1:d}, {2}:"
                          "\ntime window of pulse shape: from {3:.3f} ns to {4:.3f} ns with bin-width = {5:0.3f} "
                          "ns,"
                   .format(filenumber, event, now, min_time_total, 1000.0+20.0, binwidth, event_type))

        # get the minimum and maximum time of the prompt signal time window in ns:
        min_time_prompt = bins_hittime_prompt[0]
        max_time_prompt = bins_hittime_prompt[-1]

        # sum up the values of npe_per_hittime_prompt to get the total number of pe of the prompt signal:
        number_pe_prompt = np.sum(npe_per_hittime_prompt)

        # convert the total number of pe to quenched deposited energy in MeV:
        Qedep_converted = NC_background_functions.conversion_npe_to_evis(number_pe_prompt)

        """ store Qedep_converted in array_Q_converted_prompt_energy: """
        array_Q_converted_prompt_energy.append(Qedep_converted)

        """ save time distribution/ pulse shape to file: """
        if not flag_read_hittime_from_file:
            # save hittime distribution of the event to txt file:
            # build list, where 0th entry is reconstructed x position, 1st entry is reconstructed y position, 2nd entry
            # is reconstructed z-position, 3rd entry is start-hittime in ns, 4th entry is last-hittime in ns, 5th
            # entry is binwidth in ns and the following entries are nPE of each hittime-bin of whole signal:
            npe_per_hittime_save = [x_reconstructed, y_reconstructed, z_reconstructed]
            npe_per_hittime_save.extend([min_time_total, max_time_total, binwidth])
            npe_per_hittime_save.extend(npe_per_hittime)
            np.savetxt(output_path + "hittimes/file{0:d}_evt{1:d}_pulse_shape.txt".format(filenumber, event_id),
                       npe_per_hittime_save, fmt='%1.2f',
                       header="Pulse shape: Number of pe as function of the time "
                              "(time-of-flight correction and TTS smearing) of file user_{0}_{1:d}.root,"
                              "\nevent {2:d}, {3}:"
                              "\ntime window of pulse shape: from {4:.3f} ns to {5:.3f} ns with bin-width = {6:0.3f} "
                              "ns,"
                       .format(event_type, filenumber, event_id, now, min_time_total, max_time_total, binwidth))

        """ do prompt energy cut: """
        if prompt_energy_min <= Qedep_converted <= prompt_energy_max:
            # converted energy pass the prompt energy cut:
            number_prompt_energy_pass_reco += 1
            # append filenumber, evtID and Qedep_converted to the arrays:
            array_filenumber_prompt_energy.append(filenumber)
            array_evtID_prompt_energy.append(event_id)
            array_Evis_prompt_energy.append(Qedep_converted)
        else:
            # converted energy is rejected by prompt energy cut:
            number_prompt_energy_rejected_reco += 1

            if Qedep_converted < prompt_energy_min:
                number_prompt_energy_rejected_min_reco += 1
            elif Qedep_converted > prompt_energy_max:
                number_prompt_energy_rejected_max_reco += 1

        """ analyze delayed signal: """
        # get index, where delayed time window starts:
        index_min_hittime_delayed = index_max_hittime_prompt+2
        # get time corresponding to index_min_hittime_delayed in ns:
        min_time_delayed = bins_hittime[index_min_hittime_delayed]

        # get npe_per_hittime and bins_hittime for the delayed time window (from index_max_hittime_prompt+2 to
        # max_time_total (end of array)):
        npe_per_hittime_delayed = npe_per_hittime[index_min_hittime_delayed:]
        bins_hittime_delayed = bins_hittime[index_min_hittime_delayed:-1]

        # preallocate number of delayed signals with correct energy in delayed time window in this event:
        number_delayed_signal = 0
        # preallocate first index of npe_per_hittime_delayed:
        index_first_del = 0
        # preallocate array, where number of pe of each delayed signal is stored:
        number_pe_delayed_array = np.array([])
        # preallocate array, where begin time of each signal is stored in ns:
        begin_pulse_array = np.array([])
        # preallocate array, where end time of each signal is stored in ns:
        end_pulse_array = np.array([])

        # analyze npe_per_hittime_delayed for possible delayed signals. As long as index has not reached the end of
        # npe_per_hittime_delayed, check event for possible delayed signals:
        while index_first_del < len(npe_per_hittime_delayed):
            # is_delayed_signal (=0 if no signal, =1 if there is delayed signal with min_PE < nPE,
            # index_first_del (index after delayed signal and start of next analysis),
            # num_pe_delayed (nPE of delayed signal), begin_pulse (time, where
            # signal begins in ns), end_pulse (time, where signal ends in ns):
            is_delayed_signal, index_first_del, num_pe_delayed, begin_pulse, end_pulse = \
                NC_background_functions.analyze_delayed_signal_v2(npe_per_hittime_delayed, bins_hittime_delayed,
                                                                  index_first_del, threshold1_del, threshold2_del,
                                                                  min_PE, event_id)

            number_delayed_signal += is_delayed_signal
            # if there was a delayed signal, append variables to arrays:
            if num_pe_delayed > 0:
                # append nPE of delayed signal:
                number_pe_delayed_array = np.append(number_pe_delayed_array, num_pe_delayed)

                # append time where signal begins:
                begin_pulse_array = np.append(begin_pulse_array, begin_pulse)
                if begin_pulse == 0:
                    print("**** WARNING: num_pe_delayed > 0, but begin_pulse = 0 (event = {0:d}, file = {1:d})"
                          .format(event_id, filenumber))

                # append time where signal ends:
                end_pulse_array = np.append(end_pulse_array, end_pulse)
                if end_pulse == 2000000:
                    print("**** WARNING: num_pe_delayed > 0, but end_pulse = 2 ms (event = {0:d}, file = {1:d})"
                          .format(event_id, filenumber))

        # delayed time window is analyzed!

        """ do delayed cut on real data (reconstructed) and on ideal (MC truth) data: """
        for cut_type in range(2):
            # cut_type = 0 means delayed cut on real (reconstructed) data.
            # cut_type = 1 means delayed cut on ideal (MC truth) data.

            if cut_type == 0:
                # do delayed cut on real data (including efficiencies and alpha of all single cuts)

                """ time cut: """
                # preallocate array, where indices of number_pe_delayed_array, begin_pulse_array, end_pulse_array,
                # that pass time cut are stored (important for further cuts):
                array_index_pass_time_cut = np.array([])
                # preallocate number of delayed signals inside delayed time window for NeutronCaptureT:
                n_delayed_signal_NeutronCaptureT = 0

                """ store begin_pulse_array and end_pulse array to array_begin_times or array_end_times: """
                array_begin_times.extend(begin_pulse_array)
                array_end_times.extend(end_pulse_array)

                # check efficiency of time cut compared to cut on NeutronCaptureT:
                # get current event in nCapture-tree:
                rtree_ncapture.GetEntry(event)
                # preallocate array, where NeutronCaptureT is stored:
                NeutronCaptureT_array = np.array([])
                # get number of neutron captures NeutronN:
                NeutronN = int(rtree_ncapture.GetBranch('NeutronN').GetLeaf('NeutronN').GetValue())
                # loop over NeutronN:
                for index in range(NeutronN):
                    # get neutron capture time in ns:
                    NeutronCaptureT = float(rtree_ncapture.GetBranch("NeutronCaptureT").GetLeaf("NeutronCaptureT")
                                            .GetValue(index))
                    # append to NeutronCaptureT array:
                    NeutronCaptureT_array = np.append(NeutronCaptureT_array, NeutronCaptureT)

                """ store NeutronCaptureT_array (NeutronCaptureT of one event) to array_NeutronCaptureT 
                (NeutronCaptureT of all events): """
                array_NeutronCaptureT.extend(NeutronCaptureT_array)

                # set flag to determine time cut efficiency:
                flag_pass_NeutronCaptureT = 0
                flag_rejected_NeutronCaptureT = 0
                flag_pass_begin_end_pulse = 0
                flag_rejected_begin_end_pulse = 0

                # check time cut corresponding to NeutronCaptureT:
                if len(NeutronCaptureT_array) == 0:
                    # no neutron capture in event -> no delayed signal:
                    flag_rejected_NeutronCaptureT = 1
                else:
                    # at least one neutron capture in event:
                    # loop over NeutronN and check if NeutronCaptureT pass time cut:
                    for index in range(len(NeutronCaptureT_array)):
                        if min_time_delayed <= NeutronCaptureT_array[index] <= time_cut_max:
                            # at least neutron capture inside delayed time window:
                            flag_pass_NeutronCaptureT = 1
                            # increment n_delayed_signal_NeutronCaptureT:
                            n_delayed_signal_NeutronCaptureT += 1

                    if flag_pass_NeutronCaptureT == 0:
                        flag_rejected_NeutronCaptureT = 1

                # Are there signals with begin_pulse and end_pulse inside delayed time window (min_time_delayed to
                # time_cut_max)?
                if len(number_pe_delayed_array) == 0:
                    # no signal in delayed time window (either there is no delayed signal at all or the delayed signal
                    # is before min_time_delayed and therefore inside prompt time window (overlap with prompt signal)).
                    # event is rejected by time cut:
                    flag_rejected_begin_end_pulse = 1

                else:
                    # there is at least one signal in delayed time window. check, if begin_pulse and end_pulse of one
                    # signal is inside time window (min_time_delayed to time_cut_max).
                    # loop over delayed signals:
                    for index in range(len(number_pe_delayed_array)):

                        if (min_time_delayed <= begin_pulse_array[index] <= time_cut_max and
                                min_time_delayed <= end_pulse_array[index] <= time_cut_max):
                            # signal passes time cut:
                            # stored index in array_index_pass_time_cut:
                            array_index_pass_time_cut = np.append(array_index_pass_time_cut, index)

                    # check array_index_pass_time_cut:
                    if len(array_index_pass_time_cut) == 0:
                        # NO delayed signal pass time cut:
                        flag_rejected_begin_end_pulse = 1
                    else:
                        # at least one delayed signal pass time cut:
                        flag_pass_begin_end_pulse = 1

                # check efficiency of time cut:
                if flag_pass_begin_end_pulse == 1 and flag_rejected_NeutronCaptureT == 1:
                    # event passes time cut with begin and end time, but would be rejected by NeutronCaptureT time cut
                    # -> event counted too much
                    number_time_toomuch += 1
                elif flag_rejected_begin_end_pulse == 1 and flag_pass_NeutronCaptureT == 1:
                    # event is rejected by begin and end time, but would pass NeutronCaptureT time cut
                    # -> event counted too less
                    number_time_tooless += 1

                # check time cut with NeutronCaptureT:
                if flag_pass_NeutronCaptureT == 1:
                    # event pass time cut:
                    number_time_pass_MC += 1
                elif flag_rejected_NeutronCaptureT == 1:
                    # event is rejected by time cut:
                    number_time_rejected_MC += 1

                # check time cut with begin time and end time:
                if flag_pass_begin_end_pulse == 1:
                    # event pass time cut (at least one delayed signal in correct time window):
                    number_time_pass_reco += 1
                elif flag_rejected_begin_end_pulse == 1:
                    # event is rejected by time cut:
                    number_time_rejected_reco += 1
                    number_delayed_cut_rejected += 1
                    # go to next event
                    continue

                """ multiplicity cut: """
                # event passes time cut -> at least one delayed signal with nPE > min_PE.
                # number of delayed signals in delayed time window:
                n_delayed_signal_timewindow = len(array_index_pass_time_cut)

                """ store n_delayed_signal_timewindow in array_multiplicity: """
                array_multiplicity.append(n_delayed_signal_timewindow)

                # do multiplicity cut:
                # check efficiency of multiplicity cut:
                if n_delayed_signal_timewindow == multiplicity and n_delayed_signal_NeutronCaptureT != multiplicity:
                    # events are counted too much:
                    number_multiplicity_toomuch += 1
                elif n_delayed_signal_timewindow != multiplicity and n_delayed_signal_NeutronCaptureT == multiplicity:
                    # events are counted too less:
                    number_multiplicity_tooless += 1

                if n_delayed_signal_NeutronCaptureT == multiplicity:
                    # only 1 delayed signal from NeutronCaptureT:
                    number_multiplicity_pass_MC += 1
                else:
                    # 0 or more than 1 delayed signal from NeutronCaptureT:
                    number_multiplicity_rejected_MC += 1

                if n_delayed_signal_timewindow == multiplicity:
                    # only 1 delayed signal from begin/end time -> event passes multiplicity cut:
                    number_multiplicity_pass_reco += 1
                else:
                    # more than 1 delayed signal from begin/end time -> event is rejected by multiplicity cut:
                    number_multiplicity_rejected_reco += 1
                    number_delayed_cut_rejected += 1
                    # go to next event:
                    continue

                """ delayed energy cut: """
                # after time and multiplicity cut, there is only 1 signal in delayed time window!
                # check if len(array_index_pass_time_cut) is really = 1:
                if len(array_index_pass_time_cut) == 1:
                    # set index for delayed energy analysis:
                    index_delayed_signal = int(array_index_pass_time_cut[0])
                else:
                    print("ERROR: 0 or more than 1 signal in delayed time window ----------------file {0:d}, evt {1:d}"
                          .format(filenumber, event_id))

                # get nPE of delayed signal corresponding to index_delayed_signal:
                nPE_del_MC = number_pe_delayed_array[index_delayed_signal]

                """ convert nPE to visible energy in MeV: """
                Evis_delayed = NC_background_functions.conversion_npe_to_evis(nPE_del_MC)

                """ store Evis_delayed in array_Evis_delayed_energy: """
                array_Evis_delayed_energy.append(Evis_delayed)

                if min_PE_delayed_MeV <= Evis_delayed <= max_PE_delayed_MeV:
                    # Evis_del_smeared pass cut:
                    number_delayed_energy_pass_reco += 1
                else:
                    # Evis_del_smeared is rejected:
                    number_delayed_energy_rejected_reco +=1
                    number_delayed_cut_rejected += 1
                    # go to next event:
                    continue

                """ distance cut: """
                # set Qedep of delayed signal to 2.2 MeV (neutron capture on Hydrogen):
                Qedep_capture = 2.2

                # get the start position of neutron capture in mm:
                x_ncapture = float(rtree_ncapture.GetBranch("NCStartX").GetLeaf("NCStartX")
                                   .GetValue(index_delayed_signal))
                y_ncapture = float(rtree_ncapture.GetBranch("NCStartY").GetLeaf("NCStartY")
                                   .GetValue(index_delayed_signal))
                z_ncapture = float(rtree_ncapture.GetBranch("NCStartZ").GetLeaf("NCStartZ")
                                   .GetValue(index_delayed_signal))

                # calculate distance of MC truth neutron capture position to detector center in mm:
                r_MC_ncapture = np.sqrt(x_ncapture**2 + y_ncapture**2 + z_ncapture**2)

                # do vertex reconstruction of neutron capture position with function position_smearing():
                x_reco_ncapture = NC_background_functions.position_smearing(x_ncapture, Qedep_capture)
                y_reco_ncapture = NC_background_functions.position_smearing(y_ncapture, Qedep_capture)
                z_reco_ncapture = NC_background_functions.position_smearing(z_ncapture, Qedep_capture)

                # calculate distance of reconstructed neutron capture position to detector center in mm:
                r_reconstructed_ncapture = np.sqrt(x_reco_ncapture**2 + y_reco_ncapture**2 + z_reco_ncapture**2)

                # calculate distance between initial position and MC truth nCapture position in mm:
                distance_MC = np.sqrt((init_x - x_ncapture)**2 + (init_y - y_ncapture)**2 + (init_z - z_ncapture)**2)

                # calculate distance between reconstructed initial position and reconstructed nCapture position in mm:
                distance_reco = np.sqrt((x_reconstructed - x_reco_ncapture)**2 + (y_reconstructed - y_reco_ncapture)**2
                                        + (z_reconstructed - z_reco_ncapture)**2)

                """ store distance_reco to array_distance_reco: """
                array_distance_reco.append(distance_reco)

                # check efficiency of distance cut:
                if distance_reco < distance_cut and distance_MC >= distance_cut:
                    # event is counted too much:
                    number_distance_cut_toomuch += 1
                elif distance_reco >= distance_cut and distance_MC < distance_cut:
                    # event is counted too less:
                    number_distance_cut_tooless += 1

                if distance_MC < distance_cut:
                    number_distance_cut_pass_MC += 1
                else:
                    number_distance_cut_rejected_MC += 1

                if distance_reco < distance_cut:
                    # event passes distance cut:
                    number_distance_cut_pass_reco += 1
                else:
                    # event is rejected by distance cut:
                    number_distance_cut_rejected_reco += 1
                    number_delayed_cut_rejected += 1
                    # go to next event:
                    continue

                """ event has passed time, multiplicity, delayed energy and distance cut: """
                # store filenumber and event ID of this event to the array:
                array_filenumber_distance_cut.append(filenumber)
                array_evtID_distance_cut.append(event_id)

                """ volume cut on delayed position: """
                # store r_reconstructed_ncapture in array_volume_delayed_position:
                array_volume_delayed_position.append(r_reconstructed_ncapture)

                # check efficiency:
                if r_reconstructed_ncapture < R_cut_delayed_mm and r_MC_ncapture >= R_cut_delayed_mm:
                    # event is counted too much:
                    number_volume_leak_in += 1
                elif r_reconstructed_ncapture >= R_cut_delayed_mm and r_MC_ncapture < R_cut_delayed_mm:
                    # event is counted too less:
                    number_volume_leak_out += 1

                if r_MC_ncapture < R_cut_delayed_mm:
                    number_volume_pass_delayed_MC += 1
                else:
                    number_volume_rejected_delayed_MC += 1

                if r_reconstructed_ncapture < R_cut_delayed_mm:
                    # event passes volume cut:
                    number_volume_pass_delayed_reco += 1
                else:
                    # event is rejected by volume cut:
                    number_volume_rejected_delayed_reco += 1
                    number_delayed_cut_rejected += 1
                    # go to next event:
                    continue

                """ event has passed complete delayed cut (time, multiplicity, delayed energy, distance and 
                    delayed volume cut): """
                # increment number_delayed_cut_pass:
                number_delayed_cut_pass += 1

                # store filenumber and event ID of this event to the array:
                array_filenumber_delayed_cut.append(filenumber)
                array_evtID_delayed_cut.append(event_id)

            elif cut_type == 1:
                # do same delayed cut on ideal (MC truth) data like above to get the filenumber and evtID of the event.

                """ time cut: """
                # preallocate array, where indices of NeutronCaptureT that pass time cut are stored
                # (important for further cuts):
                array_index_pass_time_cut_MC = np.array([])

                # preallocate number of delayed signals inside delayed time window for NeutronCaptureT:
                n_delayed_signal_NeutronCaptureT = 0

                # check efficiency of time cut of NeutronCaptureT:
                # get current event in nCapture-tree:
                rtree_ncapture.GetEntry(event)
                # preallocate array, where NeutronCaptureT is stored:
                NeutronCaptureT_array = np.array([])
                # get number of neutron captures NeutronN:
                NeutronN = int(rtree_ncapture.GetBranch('NeutronN').GetLeaf('NeutronN').GetValue())
                # loop over NeutronN:
                for index in range(NeutronN):
                    # get neutron capture time in ns:
                    NeutronCaptureT = float(rtree_ncapture.GetBranch("NeutronCaptureT").GetLeaf("NeutronCaptureT")
                                            .GetValue(index))
                    # append to NeutronCaptureT array:
                    NeutronCaptureT_array = np.append(NeutronCaptureT_array, NeutronCaptureT)

                # check time cut corresponding to NeutronCaptureT:
                if len(NeutronCaptureT_array) == 0:
                    # no neutron capture in event -> no delayed signal -> go to next event
                    continue
                else:
                    # at least one neutron capture in event:
                    # loop over NeutronN and check if NeutronCaptureT pass time cut:
                    for index in range(len(NeutronCaptureT_array)):
                        if min_time_delayed <= NeutronCaptureT_array[index] <= time_cut_max:
                            # at least neutron capture inside delayed time window:
                            # append index to array_index_pass_time_cut_MC:
                            array_index_pass_time_cut_MC = np.append(array_index_pass_time_cut_MC, index)
                            # increment n_delayed_signal_NeutronCaptureT:
                            n_delayed_signal_NeutronCaptureT += 1

                    if n_delayed_signal_NeutronCaptureT == 0:
                        # no Neutron Capture inside delayed time window:
                        continue

                """ multiplicity cut: """
                # event passes time cut -> at least one neutron capture in delayed time window:
                # do multiplicity cut:
                if n_delayed_signal_NeutronCaptureT != multiplicity:
                    # more than 1 delayed signal from NeutronCaptureT -> go to next event
                    continue

                """ delayed energy cut: """
                # after time and multiplicity cut, there is only 1 signal from NeutronCaptureT in the time window,
                # check if len(array_index_pass_time_cut_MC) is really = 1:
                if len(array_index_pass_time_cut_MC) == 1:
                    # set index for delayed energy analysis:
                    index_delayed_signal = int(array_index_pass_time_cut_MC[0])
                else:
                    print("ERROR: 0 or more than 1 signal in delayed time window (MC truth)-------file {0:d}, evt {1:d}"
                          .format(filenumber, event_id))

                # set flag, if there is neutron capture on H (if event passes delayed energy cut):
                flag_pass_delayed_energy_MC = 0

                # get NeutronTrkid of the neutron that is captured inside delayed time window:
                NeutronTrkid = int(rtree_ncapture.GetBranch('NeutronTrkid').GetLeaf('NeutronTrkid')
                                   .GetValue(index_delayed_signal))

                # get number of particles after neutron capture:
                n = int(rtree_ncapture.GetBranch('n').GetLeaf('n').GetValue())
                # loop over n to get kinetic energy (kine) of the gamma (pdgid==22) that corresponds to the captured
                # neutron specified by NeutronTrkid
                for index in range(n):
                    # get pdgid:
                    pdgid = int(rtree_ncapture.GetBranch('pdgid').GetLeaf('pdgid').GetValue(index))
                    # get trkid:
                    trkid = int(rtree_ncapture.GetBranch('trkid').GetLeaf('trkid').GetValue(index))
                    # get kinetic energy in MeV:
                    kine = float(rtree_ncapture.GetBranch('kine').GetLeaf('kine').GetValue(index))
                    if pdgid == 22 and trkid == NeutronTrkid and (2.18 < kine < 2.25):
                        # gamma with around 2.2 MeV is produced by neutron capture in delayed time window:
                        # -> event passes delayed energy cut:
                        flag_pass_delayed_energy_MC = 1
                        break

                # do delayed energy cut:
                if flag_pass_delayed_energy_MC == 0:
                    # no neutron capture on Hydrogen in delayed time window -> event is rejected by delayed energy cut
                    # -> go to next event:
                    continue

                """ distance cut: """
                # get the start position of neutron capture in mm:
                x_ncapture = float(rtree_ncapture.GetBranch("NCStartX").GetLeaf("NCStartX")
                                   .GetValue(index_delayed_signal))
                y_ncapture = float(rtree_ncapture.GetBranch("NCStartY").GetLeaf("NCStartY")
                                   .GetValue(index_delayed_signal))
                z_ncapture = float(rtree_ncapture.GetBranch("NCStartZ").GetLeaf("NCStartZ")
                                   .GetValue(index_delayed_signal))

                # calculate distance of MC truth neutron capture position to detector center in mm:
                r_MC_ncapture = np.sqrt(x_ncapture**2 + y_ncapture**2 + z_ncapture**2)

                # calculate distance between initial position and MC truth nCapture position in mm:
                distance_MC = np.sqrt((init_x - x_ncapture)**2 + (init_y - y_ncapture)**2 + (init_z - z_ncapture)**2)

                if distance_MC >= distance_cut:
                    # event is rejected by distance cut -> go to next event!
                    continue

                """ volume cut on delayed position: """
                if r_MC_ncapture >= R_cut_delayed_mm:
                    # event is rejected by delayed volume cut -> go to next event !
                    continue

                """ event has passed complete delayed cut (time, multiplicity, delayed energy, distance and 
                    delayed volume cut) on MC truth data:
                """
                # store filenumber and event ID of this event to the array for MC truth data:
                array_filenumber_delayed_cut_MCtruth.append(filenumber)
                array_evtID_delayed_cut_MCtruth.append(event_id)

""" SAVE INFORMATION ABOUT PROMPT ENERGY CUT: """

""" build histogram from array_Q_smeared_prompt_energy: """
h1 = plt.figure(1, figsize=(11, 6))
bin_width = 1.0
bins_prompt_energy = np.arange(0.0, 200+bin_width, bin_width)
values1, bins1, patches1 = plt.hist(array_Q_converted_prompt_energy, bins_prompt_energy, align='mid', histtype='step',
                                    linewidth='1.5',
                                    label='entries = {0:d},\n'
                                          'events passing cut = {1:d},\n'
                                          'events rejected by cut = {2:d}'
                                    .format(number_total_events, number_prompt_energy_pass_reco,
                                            number_prompt_energy_rejected_reco))
plt.vlines(prompt_energy_min, ymin=0.0, ymax=(max(values1)+max(values1)*0.1), colors='k', linestyles='solid',
           label='lower edge of prompt energy window = {0:.2f} MeV'.format(prompt_energy_min))
plt.vlines(prompt_energy_max, ymin=0.0, ymax=(max(values1)+max(values1)*0.1), colors='k', linestyles='dashed',
           label='upper edge of prompt energy window = {0:.2f} MeV'.format(prompt_energy_max))
plt.xlim(xmin=0.0, xmax=200)
plt.xlabel("visible energy of prompt signal in MeV")
plt.ylabel("number of events per bin (bin-width = {0:.1f} MeV)".format(bin_width))
plt.title("Visible energy of prompt signal of {0} events".format(event_type))
plt.legend()
plt.grid()
plt.savefig(output_path_del + "histo_prompt_energy_{0}_{1:.0f}MeV_to_{2:.0f}MeV_{3:d}.png"
            .format(event_type, prompt_energy_min, prompt_energy_max, analysis_part))
plt.close()

""" save array_filenumber_prompt_energy, array_evtID_prompt_energy and array_Evis_prompt_energy of Q_converted (real) to 
txt file: """
np.savetxt(output_path_del + 'filenumber_evtID_Evis_prompt_energy_cut_{0}_{1:.0f}MeV_to_{2:.0f}MeV_{3:d}.txt'
           .format(event_type, prompt_energy_min, prompt_energy_max, analysis_part),
           np.c_[array_filenumber_prompt_energy, array_evtID_prompt_energy, array_Evis_prompt_energy], fmt='%.3f',
           header='filenumber | evtID | E_vis in MeV fo events that pass prompt energy cut')

""" save different number of events of the prompt energy cut to txt file: """
np.savetxt(output_path_del + "numbers_prompt_energy_cut_{0}_{1:.0f}MeV_to_{2:.0f}MeV_{3:d}.txt"
           .format(event_type, prompt_energy_min, prompt_energy_max, analysis_part),
           np.array([number_total_events, number_without_particles, number_total_events - number_without_particles,
                     number_prompt_energy_pass_reco, number_prompt_energy_rejected_reco,
                     efficiency_conversion,
                     number_prompt_energy_rejected_min_reco, number_prompt_energy_rejected_max_reco]), fmt='%i',
           header='number of events from analyze_prompt_delayed_cut_v2.py ({0}):'
                  '\nnumber of events of the prompt energy cut.'
                  '\nanalyzed files: user_{1}_{4:d}.root to user_{1}_{5:d}.root'
                  '\n{1} events are analyzed for E_min = {2:.2f} MeV to E_max = {3:.2f} MeV.'
                  '\nValues below:'
                  '\nnumber of total events,'
                  '\nnumber of events without initial particles,'
                  '\nnumber of events before cut (these events were analyzed),'
                  '\nnumber of events that pass prompt energy cut on prompt Q_converted,'
                  '\nnumber of events that are rejected by prompt energy cut on prompt Q_converted,'
                  '\nefficiency of conversion from nPE to MeV (N_before_conversion / N_after_conversion)'
                  '\nnumber of events that are rejected by prompt energy cut (Q_converted < E_min),'
                  '\nnumber of events that are rejected by prompt energy cut (Q_converted > E_max):'
                  .format(now, event_type, prompt_energy_min, prompt_energy_max, start_number, stop_number))

""" SAVE INFORMATION ABOUT DELAYED CUT: """

""" build 2D histogram of the begin and end times of the delayed pulse: """
# h2 = plt.figure(2, figsize=(11, 6))
# bins_times = np.arange(0, max_time_total-max_time_total/4, 10000)
# plt.hist2d(array_begin_times, array_end_times, bins=bins_times, cmap=plt.cm.YlOrBr)
# plt.vlines(time_cut_min, ymin=bins_times[0], ymax=bins_times[-1], colors='k', linestyles='dashed')
# plt.vlines(time_cut_max, ymin=bins_times[0], ymax=bins_times[-1], colors='k', linestyles='solid')
# plt.hlines(time_cut_min, xmin=bins_times[0], xmax=bins_times[-1], colors='k', linestyles='dashed',
#            label='$t_{min}$'+' = {0:.1f} ns'.format(time_cut_min))
# plt.hlines(time_cut_max, xmin=bins_times[0], xmax=bins_times[-1], colors='k', linestyles='solid',
#            label='$t_{max}$'+' = {0:.1f} ms'.format(time_cut_max/1000000.0))
# plt.xlabel("begin time of the delayed pulse in ns")
# plt.ylabel("end time of the delayed pulse in ns")
# plt.title("Begin and end time of the delayed pulses of {0} events".format(event_type))
# plt.colorbar()
# plt.legend()
# plt.savefig(output_path + "histo_times_{0}_{1:.0f}ns_to_{2:.0f}ms_{3:d}.png"
#             .format(event_type, time_cut_min, time_cut_max/1000000.0, analysis_part))
# plt.show()
# plt.close()

""" build 1D histogram of the begin times of the delayed pulses with x axis in log scale: """
h7 = plt.figure(7, figsize=(11, 6))
bin_width = 500.0
bins_times = np.arange(0.0, max_time_total, bin_width)
values7, bins7, patches7 = plt.hist(array_begin_times, bins_times, align='mid', histtype='step',
                                    linewidth='1.5',
                                    label='entries = {0:d},\n'
                                          'events passing cut = {1:d},\n'
                                          'events rejected by cut = {2:d}'
                                    .format(len(array_begin_times), number_time_pass_reco,
                                            number_time_rejected_reco))
plt.vlines(time_cut_min, ymin=0.0, ymax=(max(values7)+max(values7)*0.1), colors='k', linestyles='solid',
           label='$t_{min}$'+' = {0:.1f} ns'.format(time_cut_min))
plt.vlines(time_cut_max, ymin=0.0, ymax=(max(values7)+max(values7)*0.1), colors='k', linestyles='dashed',
           label='$t_{max}$'+' = {0:.1f} ms'.format(time_cut_max/1000000.0))
plt.xlim(xmin=50.0, xmax=(max_time_total-bin_width))
plt.xscale('log')
plt.xlabel("begin time of delayed pulses in ns")
plt.ylabel("entries per bin (bin-width = {0:.0f} ns)".format(bin_width))
plt.title("Begin times of delayed pulses of {0} events".format(event_type))
plt.legend()
plt.grid()
plt.savefig(output_path_del + "histo_begin_times_{0}_{1:.0f}ns_to_{2:.0f}ms_log_{3:d}.png"
            .format(event_type, time_cut_min, time_cut_max/1000000.0, analysis_part))
plt.close()

""" build 1D histogram of the begin times of the delayed pulses: """
h9 = plt.figure(9, figsize=(11, 6))
bin_width_9 = 10.0
bins_times_9 = np.arange(0.0, 2000, bin_width_9)
values9, bins9, patches9 = plt.hist(np.asarray(array_begin_times)/1000.0, bins_times_9, align='mid', histtype='step',
                                    linewidth='1.5',
                                    label='entries = {0:d},\nmean = {1:.0f} $\\mu$s'
                                    .format(len(array_begin_times), np.mean(array_begin_times)/1000.0))
plt.xlim(xmin=0.0, xmax=2000)
plt.xlabel("begin time of delayed pulses in $\\mu$s")
plt.ylabel("entries per bin (bin-width = {0:.0f} $\\mu$s)".format(bin_width_9))
plt.title("Begin times of delayed pulses of {0} events".format(event_type))
plt.legend()
plt.grid()
plt.savefig(output_path_del + "histo_begin_times_{0}_{1:.0f}ns_to_{2:.0f}ms_{3:d}.png"
            .format(event_type, time_cut_min, time_cut_max/1000000.0, analysis_part))
plt.close()

""" build 1D histogram of the end times of the delayed pulses with x axis in log scale: """
h8 = plt.figure(8, figsize=(11, 6))
values8, bins8, patches8 = plt.hist(array_end_times, bins_times, align='mid', histtype='step',
                                    linewidth='1.5',
                                    label='entries = {0:d},\n'
                                          'events passing cut = {1:d},\n'
                                          'events rejected by cut = {2:d}'
                                    .format(len(array_end_times), number_time_pass_reco,
                                            number_time_rejected_reco))
plt.vlines(time_cut_min, ymin=0.0, ymax=(max(values8)+max(values8)*0.1), colors='k', linestyles='solid',
           label='$t_{min}$'+' = {0:.1f} ns'.format(time_cut_min))
plt.vlines(time_cut_max, ymin=0.0, ymax=(max(values8)+max(values8)*0.1), colors='k', linestyles='dashed',
           label='$t_{max}$'+' = {0:.1f} ms'.format(time_cut_max/1000000.0))
plt.xlim(xmin=50.0, xmax=(max_time_total-bin_width))
plt.xscale('log')
plt.xlabel("end times of delayed pulses in ns")
plt.ylabel("entries per bin (bin-width = {0:.0f} ns)".format(bin_width))
plt.title("End times of delayed pulses of {0} events".format(event_type))
plt.legend()
plt.grid()
plt.savefig(output_path_del + "histo_end_times_{0}_{1:.0f}ns_to_{2:.0f}ms_log_{3:d}.png"
            .format(event_type, time_cut_min, time_cut_max/1000000.0, analysis_part))
plt.close()

""" build 1D histogram of the end times of the delayed pulses: """
h10 = plt.figure(10, figsize=(11, 6))
values10, bins10, patches10 = plt.hist(np.asarray(array_end_times)/1000.0, bins_times_9, align='mid', histtype='step',
                                       linewidth='1.5',
                                       label='entries = {0:d},\nmean = {1:.0f} $\\mu$s'
                                       .format(len(array_end_times), np.mean(array_end_times)/1000.0))
plt.xlim(xmin=0.0, xmax=2000)
plt.xlabel("end times of delayed pulses in $\\mu$s")
plt.ylabel("entries per bin (bin-width = {0:.0f} $\\mu$s)".format(bin_width_9))
plt.title("End times of delayed pulses of {0} events".format(event_type))
plt.legend()
plt.grid()
plt.savefig(output_path_del + "histo_end_times_{0}_{1:.0f}ns_to_{2:.0f}ms_{3:d}.png"
            .format(event_type, time_cut_min, time_cut_max/1000000.0, analysis_part))
plt.close()

""" build 1D histogram of NeutronCaptureT from MC truth the delayed pulses with x axis in log scale: """
h11 = plt.figure(11, figsize=(11, 6))
values11, bins11, patches11 = plt.hist(array_NeutronCaptureT, bins_times, align='mid', histtype='step',
                                       linewidth='1.5',
                                       label='entries = {0:d}'.format(len(array_NeutronCaptureT)))
plt.vlines(time_cut_min, ymin=0.0, ymax=(max(values11)+max(values11)*0.1), colors='k', linestyles='solid',
           label='$t_{min}$'+' = {0:.1f} ns'.format(time_cut_min))
plt.vlines(time_cut_max, ymin=0.0, ymax=(max(values11)+max(values11)*0.1), colors='k', linestyles='dashed',
           label='$t_{max}$'+' = {0:.1f} ms'.format(time_cut_max/1000000.0))
plt.xlim(xmin=50.0, xmax=(max_time_total-bin_width))
plt.xscale('log')
plt.xlabel("Neutron capture time in ns")
plt.ylabel("entries per bin (bin-width = {0:.0f} ns)".format(bin_width))
plt.title("Neutron capture time of MC truth of {0} events".format(event_type))
plt.legend()
plt.grid()
plt.savefig(output_path_del + "histo_NeutronCaptureT_{0}_{1:.0f}ns_to_{2:.0f}ms_log_{3:d}.png"
            .format(event_type, time_cut_min, time_cut_max/1000000.0, analysis_part))
plt.close()

""" build 1D histogram of NeutronCaptureT from MC truth the delayed pulses: """
h12 = plt.figure(12, figsize=(11, 6))
values12, bins12, patches12 = plt.hist(np.asarray(array_NeutronCaptureT)/1000.0, bins_times_9,
                                       align='mid', histtype='step', linewidth='1.5',
                                       label='entries = {0:d},\nmean = {1:.0f} $\\mu$s'
                                       .format(len(array_NeutronCaptureT), np.mean(array_NeutronCaptureT)/1000.0))
plt.xlim(xmin=0.0, xmax=2000)
plt.xlabel("Neutron capture time in $\\mu$s")
plt.ylabel("entries per bin (bin-width = {0:.0f} $\\mu$s)".format(bin_width_9))
plt.title("Neutron capture time of MC truth of {0} events".format(event_type))
plt.legend()
plt.grid()
plt.savefig(output_path_del + "histo_NeutronCaptureT_{0}_{1:.0f}ns_to_{2:.0f}ms_{3:d}.png"
            .format(event_type, time_cut_min, time_cut_max/1000000.0, analysis_part))
plt.close()

""" save different number of events of the time cut to txt file: """
np.savetxt(output_path_del + "numbers_time_cut_{0}_{1:.0f}ns_to_{2:.0f}ms_{3:d}.txt"
           .format(event_type, time_cut_min, time_cut_max/1000000.0, analysis_part),
           np.array([number_total_events, number_without_particles, number_total_events - number_without_particles,
                     number_time_pass_reco, number_time_rejected_reco,
                     number_time_pass_MC, number_time_rejected_MC,
                     number_time_toomuch, number_time_tooless]), fmt='%i',
           header='number of events from analyze_prompt_delayed_cut_v2.py ({0}):'
                  '\nnumber of events of the time cut.'
                  '\nanalyzed files: user_{1}_{4:d}.root to user_{1}_{5:d}.root'
                  '\n{1} events are analyzed for t_min = {2:.2f} ns to t_max = {3:.2f} ns.'
                  '\nValues below:'
                  '\nnumber of total events,'
                  '\nnumber of events without initial particles,'
                  '\nnumber of events before cut (these events were analyzed),'
                  '\nnumber of events that pass time cut on begin/end time,'
                  '\nnumber of events that are rejected by time cut on begin/end time,'
                  '\nnumber of events that pass time cut on NeutronCaptureT,'
                  '\nnumber of events that are rejected by time cut on NeutronCaptureT,'
                  '\nnumber of events counted too much (begin/end time pass, but NeutronCaptureT is rejected),'
                  '\nnumber of events counted too less (begin/end time is rejected, but NeutronCaptureT pass):'
                  .format(now, event_type, time_cut_min, time_cut_max, start_number, stop_number))

""" build histogram from array_multiplicity: """
h3 = plt.figure(3, figsize=(11, 6))
bin_width = 1.0
bins_multiplicity = np.arange(0, 10+bin_width, bin_width)
values3, bins3, patches3 = plt.hist(array_multiplicity, bins_multiplicity, align='left', histtype='step',
                                    linewidth='1.5',
                                    label='entries = {0:d},\n'
                                          'events passing cut = {1:d},\n'
                                          'events rejected by cut = {2:d}'
                                    .format(number_time_pass_reco, number_multiplicity_pass_reco,
                                            number_multiplicity_rejected_reco))
plt.vlines(multiplicity, ymin=0.0, ymax=(max(values3)+max(values3)*0.1), colors='k', linestyles='solid',
           label='multiplicity cut ({0:d} signal in delayed time window)'.format(multiplicity))
plt.xlim(xmin=0.0, xmax=10)
plt.xticks(bins_multiplicity)
plt.xlabel("number of signal in delayed time window")
plt.ylabel("number of events per bin")
plt.title("Number of signals in delayed time window for {0} events".format(event_type))
plt.legend()
plt.grid()
plt.savefig(output_path_del + "histo_multiplicity_{0}_mult{1:d}_{2:d}.png"
            .format(event_type, multiplicity, analysis_part))
plt.close()

""" save different number of events of the multiplicity cut to txt file: """
np.savetxt(output_path_del + "numbers_multiplicity_cut_{0}_mult{1:d}_{2:d}.txt"
           .format(event_type, multiplicity, analysis_part),
           np.array([number_total_events, number_without_particles, number_time_pass_reco,
                     number_multiplicity_pass_reco, number_multiplicity_rejected_reco,
                     number_multiplicity_pass_MC, number_multiplicity_rejected_MC,
                     number_multiplicity_toomuch, number_multiplicity_tooless]), fmt='%i',
           header='number of events from analyze_prompt_delayed_cut_v2.py ({0}):'
                  '\nnumber of events of the multiplicity cut.'
                  '\nanalyzed files: user_{1}_{3:d}.root to user_{1}_{4:d}.root'
                  '\n{1} events are analyzed for multiplicity = {2:d} (exact {2:d} signal in delayed time window).'
                  '\nValues below:'
                  '\nnumber of total events,'
                  '\nnumber of events without initial particles,'
                  '\nnumber of events that pass the time cut+,'
                  '\nnumber of events that pass multiplicity cut based on begin/end time,'
                  '\nnumber of events that are rejected by multiplicity cut based on begin/end time,'
                  '\nnumber of events that pass multiplicity cut based on NeutronCaptureT,'
                  '\nnumber of events that are rejected by multiplicity cut based on NeutronCaptureT,'
                  '\nnumber of events counted too much (multiplicity based on begin/end time pass,'
                  ' but based on NeutronCaptureT is rejected),'
                  '\nnumber of events counted too less (multiplicity based on begin/end time is rejected, '
                  'but based on NeutronCaptureT pass):'
                  .format(now, event_type, multiplicity, start_number, stop_number))

""" build histogram from array_Evis_delayed_energy: """
h13 = plt.figure(13, figsize=(11, 6))
bin_width = 0.05
bins_Evis_del = np.arange(0, 13.0+bin_width, bin_width)
values14, bins14, patches14 = plt.hist(array_Evis_delayed_energy, bins_Evis_del, align='mid', histtype='step',
                                       linewidth='1.5',
                                       label='entries = {0:d},\n'
                                             'events passing cut = {1:d},\n'
                                             'events rejected by cut = {2:d}'
                                       .format(number_multiplicity_pass_reco, number_delayed_energy_pass_reco,
                                               number_delayed_energy_rejected_reco))
plt.vlines(min_PE_delayed_MeV, ymin=0.0, ymax=(max(values14)+max(values14)*0.1), colors='k', linestyles='solid',
           label='lower edge of delayed energy window = {0:.1f} MeV'.format(min_PE_delayed_MeV))
plt.vlines(max_PE_delayed_MeV, ymin=0.0, ymax=(max(values14)+max(values14)*0.1), colors='k', linestyles='dashed',
           label='upper edge of delayed energy window = {0:.1f} MeV'.format(max_PE_delayed_MeV))
plt.xlim(xmin=0.0, xmax=13.0)
plt.xlabel("visible energy of delayed signal in MeV")
plt.ylabel("number of events per bin (bin-width = {0:.2f} MeV)".format(bin_width))
plt.title("Visible energy of delayed signals of {0} events".format(event_type))
plt.legend()
plt.grid()
plt.savefig(output_path_del + "histo_delayed_energy_{0}_{1:.0f}keV_to_{2:.0f}keV_{3}.png"
            .format(event_type, min_PE_delayed_MeV*1000, max_PE_delayed_MeV*1000, analysis_part))
plt.close()

""" build histogram from array_Evis_delayed_energy with Gaussian Fit: """
h4 = plt.figure(4, figsize=(11, 6))
values4, bins4, patches4 = plt.hist(array_Evis_delayed_energy, bins_Evis_del, align='mid', histtype='step',
                                    linewidth='1.5',
                                    label='entries = {0:d},\n'
                                          'events passing cut = {1:d},\n'
                                          'events rejected by cut = {2:d}'
                                    .format(number_multiplicity_pass_reco, number_delayed_energy_pass_reco,
                                            number_delayed_energy_rejected_reco))
# do gaussian fit on the energy:
popt_4, pcov_4 = curve_fit(gaus, bins4[:-1], values4)
energy_range = np.arange(0, 13.0+0.01, 0.01)
plt.plot(energy_range, gaus(energy_range, *popt_4), "r-", label="Gaussian fit: $\\mu$ = {0:.2f} MeV, "
                                                                "$\\sigma$ = {1:.2f} MeV"
         .format(popt_4[1], popt_4[2]))
plt.vlines(min_PE_delayed_MeV, ymin=0.0, ymax=(max(values4)+max(values4)*0.1), colors='k', linestyles='solid',
           label='lower edge of delayed energy window = {0:.1f} MeV'.format(min_PE_delayed_MeV))
plt.vlines(max_PE_delayed_MeV, ymin=0.0, ymax=(max(values4)+max(values4)*0.1), colors='k', linestyles='dashed',
           label='upper edge of delayed energy window = {0:.1f} MeV'.format(max_PE_delayed_MeV))
plt.xlim(xmin=0.0, xmax=13.0)
plt.xlabel("visible energy of delayed signal in MeV")
plt.ylabel("number of events per bin (bin-width = {0:.2f} MeV)".format(bin_width))
plt.title("Visible energy of delayed signals of {0} events".format(event_type))
plt.legend()
plt.grid()
plt.savefig(output_path_del + "histo_delayed_energy_{0}_{1:.0f}keV_to_{2:.0f}keV_FIT_{3}.png"
            .format(event_type, min_PE_delayed_MeV*1000, max_PE_delayed_MeV*1000, analysis_part))
plt.close()

""" save different number of events of the delayed energy cut to txt file: """
np.savetxt(output_path_del + "numbers_delayed_energy_cut_{0}_{1:.0f}keV_to_{2:.0f}keV_{3:d}.txt"
           .format(event_type, min_PE_delayed_MeV*1000, max_PE_delayed_MeV*1000, analysis_part),
           np.array([number_total_events, number_without_particles, number_multiplicity_pass_reco,
                     number_delayed_energy_pass_reco, number_delayed_energy_rejected_reco,
                     efficiency_conversion]), fmt='%i',
           header='number of events from analyze_prompt_delayed_cut_v2.py ({0}):'
                  '\nnumber of events of the delayed energy cut.'
                  '\nanalyzed files: user_{1}_{4:d}.root to user_{1}_{5:d}.root'
                  '\n{1} events are analyzed for Evis_min = {2:.2f} MeV to Evis_max = {3:.2f} MeV.'
                  '\nValues below:'
                  '\nnumber of total events,'
                  '\nnumber of events without initial particles,'
                  '\nnumber of events that pass multiplicity cut (these events were analyzed),'
                  '\nnumber of events that pass delayed energy cut on Evis_converted,'
                  '\nnumber of events that are rejected by delayed energy cut on Evis_converted,'
                  '\nefficiency of conversion from nPE to MeV (N_before_conversion / N_after_conversion):'
                  .format(now, event_type, min_PE_delayed_MeV, max_PE_delayed_MeV, start_number, stop_number))

""" build histogram from array_distance_reco: """
h5 = plt.figure(5, figsize=(11, 6))
bin_width = 100
bins_dist = np.arange(0, 15000+bin_width, bin_width)
values5, bins5, patches5 = plt.hist(array_distance_reco, bins_dist, align='mid', histtype='step',
                                    linewidth='1.5',
                                    label='entries = {0:d},\n'
                                          'events passing cut = {1:d},\n'
                                          'events rejected by cut = {2:d}'
                                    .format(number_delayed_energy_pass_reco, number_distance_cut_pass_reco,
                                            number_distance_cut_rejected_reco))
plt.vlines(distance_cut, ymin=0.0, ymax=(max(values5)+max(values5)*0.1), colors='k', linestyles='solid',
           label='distance cut value = {0:.1f} mm'.format(distance_cut))
plt.xlim(xmin=0.0, xmax=15000)
plt.xlabel("distance between prompt and delayed signal in mm")
plt.ylabel("number of events per bin (bin-width = {0:d} mm)".format(bin_width))
plt.title("Distance between reconstructed positions of prompt and delayed signal of {0} events".format(event_type))
plt.legend()
plt.grid()
plt.savefig(output_path_del + "histo_distance_{0}_{1:.0f}mm_{2}.png"
            .format(event_type, distance_cut, analysis_part))
plt.close()

""" save array_filenumber_distance_cut and array_evtID_distance_cut to txt file: """
np.savetxt(output_path_del + 'filenumber_evtID_distance_cut_{0}_{1:.0f}ns_to_{2:.0f}ms_mult{3:d}_{4:.0f}keV_{5:.0f}keV_'
                             'dist{6:.0f}mm_{7:d}.txt'
           .format(event_type, time_cut_min, time_cut_max/1000000.0, multiplicity, min_PE_delayed_MeV*1000,
                   max_PE_delayed_MeV*1000, distance_cut, analysis_part),
           np.c_[array_filenumber_delayed_cut, array_evtID_delayed_cut], fmt='%i',
           header='filenumber | evtID of events that pass distance cut (and time, multiplicity and delayed energy cut)')

""" save different number of events of the distance cut to txt file: """
np.savetxt(output_path_del + "numbers_distance_cut_{0}_{1:.0f}mm_{2:d}.txt"
           .format(event_type, distance_cut, analysis_part),
           np.array([number_total_events, number_without_particles, number_delayed_energy_pass_reco,
                     number_distance_cut_pass_reco, number_distance_cut_rejected_reco,
                     number_distance_cut_pass_MC, number_distance_cut_rejected_MC,
                     number_distance_cut_toomuch, number_distance_cut_tooless]), fmt='%i',
           header='number of events from analyze_prompt_delayed_cut_v2.py ({0}):'
                  '\nnumber of events of the distance cut.'
                  '\nanalyzed files: user_{1}_{3:d}.root to user_{1}_{4:d}.root'
                  '\n{1} events are analyzed for distance cut = {2:.1f} mm.'
                  '\nValues below:'
                  '\nnumber of total events,'
                  '\nnumber of events without initial particles,'
                  '\nnumber of events that pass delayed energy cut (these events were analyzed),'
                  '\nnumber of events that pass distance cut on distance_reco,'
                  '\nnumber of events that are rejected by distance cut on distance_reco,'
                  '\nnumber of events that pass distance cut on distance_MC,'
                  '\nnumber of events that are rejected by distance cut on distance_MC,'
                  '\nnumber of events counted too much (distance_reco pass, but distance_MC is rejected),'
                  '\nnumber of events counted too less (distance_reco is rejected, but distance_MC pass):'
                  .format(now, event_type, distance_cut, start_number, stop_number))

""" build histogram from array_volume_delayed_position: """
h6 = plt.figure(6, figsize=(11, 6))
bin_width = 500
bins_radius = np.arange(0, 18000+bin_width, bin_width)
values6, bins6, patches6 = plt.hist(array_volume_delayed_position, bins_radius, align='mid', histtype='step',
                                    linewidth='1.5',
                                    label='entries = {0:d},\n'
                                          'events passing cut = {1:d},\n'
                                          'events rejected by cut = {2:d}'
                                    .format(number_distance_cut_pass_reco, number_volume_pass_delayed_reco,
                                            number_volume_rejected_delayed_reco))
plt.vlines(R_cut_delayed_mm, ymin=0.0, ymax=(max(values6)+max(values6)*0.1), colors='k', linestyles='solid',
           label='fiducial volume cut on delayed signal $R_{FV}$ = '+'{0:.0f} mm'.format(R_cut_delayed_mm))
plt.xlim(xmin=0.0, xmax=18000)
plt.xlabel("distance to detector center in mm")
plt.ylabel("number of events per bin (bin-width = {0:d} mm)".format(bin_width))
plt.title("Reconstructed position of delayed signal of {0} events".format(event_type))
plt.legend(loc='upper left')
plt.grid()
plt.savefig(output_path_del + "histo_delayed_volume_cut_{0}_{1:.0f}mm_{2}.png"
            .format(event_type, R_cut_delayed_mm, analysis_part))
plt.close()

""" save different number of events of the delayed volume cut to txt file: """
np.savetxt(output_path_del + "numbers_delayed_volume_cut_{0}_{1:.0f}mm_{2:d}.txt"
           .format(event_type, R_cut_delayed_mm, analysis_part),
           np.array([number_total_events, number_without_particles, number_distance_cut_pass_reco,
                     number_volume_pass_delayed_reco, number_volume_rejected_delayed_reco,
                     number_volume_pass_delayed_MC, number_volume_rejected_delayed_MC,
                     number_volume_leak_in, number_volume_leak_out]), fmt='%i',
           header='number of events from analyze_prompt_delayed_cut_v2.py ({0}):'
                  '\nnumber of events of the volume cut on the delayed signal.'
                  '\nanalyzed files: user_{1}_{3:d}.root to user_{1}_{4:d}.root'
                  '\n{1} events are analyzed for fiducial volume cut = {2:.1f} mm.'
                  '\nValues below:'
                  '\nnumber of total events,'
                  '\nnumber of events without initial particles,'
                  '\nnumber of events that pass distance cut (these events were analyzed),'
                  '\nnumber of events that pass delayed volume cut on r_reco_ncapture,'
                  '\nnumber of events that are rejected by delayed volume cut on r_reco_ncapture,'
                  '\nnumber of events that pass delayed volume cut on r_MC_ncapture,'
                  '\nnumber of events that are rejected by delayed volume cut on r_MC_ncapture,'
                  '\nnumber of leak-in events (counted too much) (initial position outside, but recon. position '
                  'inside fiducial volume),'
                  '\nnumber of leak-out events (counted too less) (initial position inside, but recon. position '
                  'outside fiducial volume):'
                  .format(now, event_type, R_cut_delayed_mm, start_number, stop_number))

""" save array_filenumber_delayed_cut and array_evtID_delayed_cut to txt file: """
np.savetxt(output_path_del + 'filenumber_evtID_delayed_cut_{0}_{1:.0f}ns_to_{2:.0f}ms_mult{3:d}_{4:.0f}keV_{5:.0f}keV_'
                             'dist{6:.0f}mm_R{7:.0f}mm_{8:d}.txt'
           .format(event_type, time_cut_min, time_cut_max/1000000.0, multiplicity, min_PE_delayed_MeV*1000,
                   max_PE_delayed_MeV*1000, distance_cut, R_cut_delayed_mm, analysis_part),
           np.c_[array_filenumber_delayed_cut, array_evtID_delayed_cut], fmt='%i',
           header='filenumber | evtID of events that pass delayed cut')

""" save array_filenumber_delayed_cut_MCtruth and array_evtID_delayed_cut_MCtruth to txt file: """
np.savetxt(output_path_del + 'filenumber_evtID_delayed_cut_MCtruth_{0}_{1:.0f}ns_to_{2:.0f}ms_mult{3:d}_{4:.0f}keV_'
                             '{5:.0f}keV_dist{6:.0f}mm_R{7:.0f}mm_{8:d}.txt'
           .format(event_type, time_cut_min, time_cut_max/1000000.0, multiplicity, min_PE_delayed_MeV*1000,
                   max_PE_delayed_MeV*1000, distance_cut, R_cut_delayed_mm, analysis_part),
           np.c_[array_filenumber_delayed_cut_MCtruth, array_evtID_delayed_cut_MCtruth], fmt='%i',
           header='filenumber | evtID of events that pass delayed cut')





