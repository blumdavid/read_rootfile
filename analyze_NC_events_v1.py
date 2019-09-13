""" Script to analyze the NC events:

    Version 1:  the efficiencies are calculated depending on each other and no distributions of the cut parameters are
                saved.
                -> use Version 2: here the efficiencies are calculated independent of each other!

    Apply cut on the NC events to select IBD-like signals only.

    And get the energy and hittime distribution of the prompt signal of atmospheric NC neutrino background that
    are simulated with JUNO detector simulation.

    1.  read all NC events (user_atmoNC_0.root to user_atmoNC_999.root, each file contains 100 events.
        Therefore, 100000 events are analyzed.)

    2.  Start analyzing the events:
        2.1.    first analysis: do volume cut on the reconstructed position of the event
                (you only need InitX,Y,Z and sum of Qedep to get the reconstructed position).
                Only analyze event that pass volume cut further.

    3.  Calculate hittime distribution (with time-of-flight correction and PMT time resolution) for each event:
        Procedure to get the hittime distribution with vertex reconstruction and time smearing of PMTs (same procedure
        like in script hittime_distribution_positron.py):

        2.1.    calculate time of flight:
            2.1.1   for every photon, that hits a PMT (20inch and 3inch), take the PMT position (via PMT ID from file
                    PMT_position.root) and calculate the time-of-flight distance with the reconstructed position from
                    file evtID_preselected_{}.txt
            2.1.2   with the time-of-flight distance, calculate the time-of-flight of this photon from production to
                    PMT by considering an effective speed of light in the LS.

        2.2.    consider TTS of PMTs:
            2.2.1   for every photon, that hits a PMT (20inch and 3inch), take the time resolution (sigma) of the PMT
                    (via PMT ID either from file PmtData.root for the 20inch PMTs or set TTS = 5 ns for 3inch PMTs.)
            2.2.2   the TTS of the PMT is FWHM. Therefore calculate sigma from TTS (FWHM = 2*sqrt(2*ln(2)) * sigma).
            2.2.3   smear hittime of detsim with gaussian of sigma (time resolution) around the value of detsim hittime
                    to get the smeared hittime

        2.3.    for every photon, calculate the 'real' hittime (= smeared hittime - time_of_flight) and store it in
                array

        2.4.    Do points 2.1 to 2.3 for every photon. Then you get the correct hittime of this event. Build histogram
                with correct hittimes.

    3.  Take the prompt signal of the corrected hittime histogram and do a cut on the prompt signal:
        use function conversion_npe_to_evis() and convert the number of pe of the prompt signal to visible energy in
        MeV.
        Then apply energy resolution on the visible energy and do a cut on the visible energy (10 MeV to 100 MeV).
        Only analyze events further that pass prompt energy cut.

    4.  Analyze delayed signal:
        Take the delayed signal of the corrected hittime histogram and do several cuts on the delayed signal (time cut,
        delayed energy cut, multiplicity cut and distance cut):

        4.1.    analyze delayed signal with function analyze_delayed_signal(). With this you get:
            4.1.1.  if there is a delayed signal in the event
            4.1.2.  the number of pe in this delayed signal
            4.1.3.  the time, where signal begins and ends

            for all delayed signals in one event (there could also be several delayed signals in one event)

        4.2.    do time cut: begin and end time must be within delayed time window.
                Only analyze events that pass time cut further.
        4.3.    do delayed energy cut: take number of pe of each delayed signal, consider energy smearing to get the
                reconstructed energy and check if reconstructed energy is within delayed energy window.
                Only analyze events that pass delayed energy cut further.
        4.4.    do neutron multiplicity cut: check if there is only 1 delayed signal in the event that pass time and
                delayed energy cut.
                Only analyze events that pass multiplicity cut further.
        4.5.    do distance cut: distance between reconstructed position of prompt signal and reconstructed position of
                delayed signal must be smaller than cut-distance (normally 1.5 m).

        For the delayed cut, also the efficiency of the efficiency of the cuts must be determine (how many events are
        rejected/ pass falsely).

    5.  For events that pass all cuts:
        5.1.    Save the hittime histogram of the prompt signal to txt file and png file for further analysis with
                script pulse_shape_analysis_v1.py.

        5.2.    Save evtID and visible energy of prompt signal of each event in txt file.

"""
import datetime
import ROOT
import sys
import NC_background_functions
import numpy as np
from matplotlib import pyplot as plt

# get the date and time, when the script was run:
date = datetime.datetime.now()
now = date.strftime("%Y-%m-%d %H:%M")

""" file information: """
# first file to be read:
start_number = 900
# last file to be read:
stop_number = 999
# number of entries in the input files:
Number_entries_input = 100
# set the path of the input root files:
input_path = "/local/scratch1/pipc51/astro/blum/detsim_output_data/"
# set the path of the output, where the txt file with visible energy of events that pass cuts is saved:
output_path = "/home/astro/blum/juno/atmoNC/data_NC/output_detsim/"

""" define parameters depending on the cuts: """
# fiducial volume cut in mm:
R_cut_mm = 16000
# minimum of prompt energy in MeV:
prompt_energy_min = 10.0
# maximum of prompt energy in MeV:
prompt_energy_max = 100.0
# set time window of whole signal in ns:
min_time_total = -50
max_time_total = 1001000
# time window of delayed signal in ns:
time_cut_min = 500
time_cut_max = 1000000
# Set bin-width of hittime histogram in ns:
binwidth = 5.0
# min and max number of PE for delayed energy cut:
min_PE_delayed = 2500
max_PE_delayed = 3400
# Set threshold of number of PE per bin for possible delayed signal (bin-width = 5 ns):
threshold1_del = 50
# set threshold2 of number of PEs per bin (signal peak is summed as long as nPE is above threshold2):
threshold2_del = 0
# distance cut between prompt and delayed signal in mm:
distance_cut = 1500

# number of analyzed events (total number of events):
number_total_events = (stop_number + 1 - start_number) * Number_entries_input

""" preallocate number of events for different scenarios: """
# number of events that pass all cuts:
number_IBDlike_events = 0
# number of events that are rejected by cuts:
number_rejected = 0
# number of events without initial particles:
number_without_particles = 0
""" volume """
# number of events that pass volume cut on prompt position:
number_volume_pass_prompt = 0
# number of events that are rejected by volume cut on prompt position:
number_volume_rejected_prompt = 0
# number of events that are rejected by volume cut on neutron capture position:
number_volume_rejected_delayed = 0
""" prompt energy """
# number of events that pass prompt energy cut (and volume cut):
number_prompt_energy_pass = 0
# number of events that are rejected by prompt energy cut:
number_prompt_energy_rejected = 0
# number of events that are rejected by prompt energy cut (E < prompt_energy_min):
number_prompt_energy_rejected_min = 0
# number of events that are rejected by prompt energy cut (E > prompt_energy_max):
number_prompt_energy_rejected_max = 0
# number of events, where smeared energy pass the cut, but real energy would be rejected:
number_prompt_energy_toomuch = 0
# number of events, where smeared energy is rejected by cut, but real energy would pass the cut:
number_prompt_energy_tooless = 0
""" time """
# number of events that pass time cut (and volume and prompt energy cut):
number_time_pass = 0
# number of events that are rejected by time cut:
number_time_rejected = 0
# number of events that pass time cut, but would be rejected by nCaptureT (events counted too much):
number_time_toomuch = 0
# number of events that are rejected by time cut, but nCaptureT would pass (events counted too less):
number_time_tooless = 0
""" delayed energy """
# number of events that pass delayed energy cut (and volume, prompt energy and time cut):
number_delayed_energy_pass = 0
# number of events that are rejected by delayed energy cut:
number_delayed_energy_rejected = 0
# number of events, where smeared energy pass the cut, but real energy would be rejected:
number_delayed_energy_toomuch = 0
# number of events, where smeared energy is rejected by cut, but real energy would pass the cut:
number_delayed_energy_tooless = 0
""" multiplicity """
# number of events that pass neutron multiplicity cut (and volume, prompt energy, time and delayed energy cut:
number_n_mult_pass = 0
# number of events that do not pass neutron multiplicity cut:
number_n_mult_rejected = 0
# number of events, where smeared energy pass multiplicity cut, but real energy would be rejected by multiplicity cut:
number_n_mult_toomuch = 0
# number of events, where smeared energy is rejected by multiplicity cut, but real energy would pass multiplicity cut:
number_n_mult_tooless = 0
""" distance """
# number of events that pass distance cut (and volume, prompt energy, time, delayed energy and multiplicity cut):
number_distance_cut_pass = 0
# number of events that are rejected by distance cut:
number_distance_cut_rejected = 0
# number of events, where reconstructed distance pass the cut, but "real" (MC truth) distance is rejected:
number_distance_cut_toomuch = 0
# number of events, where reconstructed distance is rejected by cut, but "real" distance pass:
number_distance_cut_tooless = 0

""" load position of the PMTs and corresponding PMT ID from file PMT_position.root: """
file_PMT_position = "/home/astro/blum/juno/atmoNC/PMT_information/PMT_position.root"
# array with PMT ID and corresponding x, y, z position in mm:
pmtID_pos_file, x_pos_pmt, y_pos_pmt, z_pos_pmt = NC_background_functions.get_pmt_position(file_PMT_position)

""" load 'time resolution' in ns of the 20 inch PMTs and corresponding PMT ID from file PmtData.root: """
file_PMT_time = "/home/astro/blum/juno/atmoNC/PMT_information/PmtData.root"
# array with PMT ID and corresponding sigma in ns:
pmtID_time_file, sigma_time_20inch = NC_background_functions.get_20inchpmt_tts(file_PMT_time)
# set TTS (FWHM) of the 3inch PMTs in ns:
tts_3inch = 5.0
# calculate time resolution (sigma) for the 3inch PMTs in ns:
sigma_time_3inch = tts_3inch / (2 * np.sqrt(2 * np.log(2)))
# set effective speed of light in the liquid scintillator in mm/ns (see page 7 of c_effective_JUNO-doc-3144-v2.pdf in
# folder /home/astro/blum/PhD/paper/Pulse_Shape_Discrimination/). Effective refraction index in LS n_eff = 1.54.
# c/n_eff = 299792458 m / 1.54 s ~ 194670427 m/s = 194670427 * 10**(-6) mm/ns ~ 194.67 mm/ns:
c_effective = 194.67

# loop over the files that are read:
for filenumber in range(start_number, stop_number+1):

    # file name of the input file:
    input_name = input_path + "user_atmoNC_{0:d}.root".format(filenumber)
    print("------------------------------------------------------------------------------------")
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

    """ preallocate arrays to be saved for each file: """
    # evtID of the events, that pass all cuts:
    eventID_pass_array = np.array([])
    # visible energy in MeV of events, that pass all cuts:
    Evis_pass_array = np.array([])

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
            number_rejected += 1
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

        # get current event in TTree:
        rtree_prmtrkdep.GetEntry(event)
        # get number of initial particles in prmtrkdep:
        nInitParticles_prmtrkdep = int(rtree_geninfo.GetBranch('nInitParticles').GetLeaf('nInitParticles').GetValue())

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
        r_reconstructed = np.sqrt(x_reconstructed**2 + y_reconstructed**2 + z_reconstructed**2)

        """ volume cut on prompt signal """
        # apply volume cut:
        if r_reconstructed >= R_cut_mm:
            # event is rejected by volume cut:
            number_volume_rejected_prompt += 1
            number_rejected += 1
            # go to next event:
            continue
        else:
            # event passes volume cut on prompt signal:
            number_volume_pass_prompt += 1

        """ calculate the real hittime distribution (time of flight correction with reconstructed position and time 
        smearing with TTS for each hit): """
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

            # get position of the PMT with specific pmtID (pmtID is ascending number from 0 to 17738 (17739 large PMTs)
            # and from 300000 to 336571 (36572 small PMTs)).
            # For large PMTs -> For 20inch PMTs, the pmtID is equal to index of x,y,z_pos_pmt array.
            # For small PMTs -> For 3inch PMTs, the pmtID - (300000 - 17739) is equal to index of x,y,z_pos_pmt array.
            # check if PMT is 20 inch or 3inch (pmtID < 20000 means 20inch PMT):
            if pmtID < 20000:
                # 20inch PMT:
                # get PMT position in mm from arrays:
                x_pmt = x_pos_pmt[pmtID]
                y_pmt = y_pos_pmt[pmtID]
                z_pmt = z_pos_pmt[pmtID]

            elif 20000 < pmtID < 40000:
                # there are some PMTs with ID around 30000 (user_atmoNC_7.root, event=32: 30637, 30276, 30573, 30561,
                # 30377) -> PMTs with ID above 30000 are Water Pool PMTs!!
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

            # calculate time of flight in ns:
            time_of_flight = distance_tof / c_effective

            """ time resolution of PMT: """
            # get time resolution of PMT with specific pmtID (pmtID is ascending number from 0 to 17738 (17739 large
            # PMTs)) -> For 20inch PMTs, the pmtID is equal to index of sigma_time_20inch array.
            # check if PMT is 20 inch or 3inch (pmtID < 20000 means 20inch PMT):
            if pmtID < 20000:
                # 20inch PMT:
                # get time resolution (sigma) of PMT in ns from array:
                sigma_pmt = sigma_time_20inch[pmtID]

            elif 20000 < pmtID < 40000:
                # there are some PMTs with ID around 30000 (user_atmoNC_7.root, event=32: 30637, 30276, 30573, 30561,
                # 30377) -> PMTs with ID above 30000 are Water Pool PMTs!!
                # go to next photon:
                continue

            else:
                # 3inch PMT:
                sigma_pmt = sigma_time_3inch

            # consider time resolution of PMT by generating normal distributed random number with mu = hit_time and
            # sigma = sigma_pmt (only the hit_time at the PMT must be smeared, not the time-of-flight):
            hittime_tts = np.random.normal(hit_time, sigma_pmt)

            # calculate the 'real' hittime of the photon in ns:
            hittime_real = hittime_tts - time_of_flight

            # append hittime to array:
            hittime_array.append(hittime_real)

        # hittime_array contains now the corrected hittimes of all PMTs in ns!

        """ analyze prompt signal: """
        # build histogram, where hittimes are saved:
        # set bin-edges of hittime histogram in ns in whole time window:
        bins_hittime = np.arange(min_time_total, max_time_total + 2 * binwidth, binwidth)
        # build hittime histogram:
        npe_per_hittime, bin_edges_hittime = np.histogram(hittime_array, bins_hittime)

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
            for index in range(index_time_cut_min, index_time_cut_min + 200):
                if npe_per_hittime[index] == npe_per_hittime[index + 1] == npe_per_hittime[index + 2] == 0:
                    index_max_hittime_prompt = index
                    break

        # calculate nPE as function of hittime only for prompt time window (last index should be included):
        npe_per_hittime_prompt = npe_per_hittime[index_min_hittime_prompt:index_max_hittime_prompt + 1]
        # bin edges of hittime histogram only for prompt time window:
        bins_hittime_prompt = bin_edges_hittime[index_min_hittime_prompt:index_max_hittime_prompt + 1]

        # get the minimum and maximum time of the prompt signal time window in ns:
        min_time_prompt = bins_hittime_prompt[0]
        max_time_prompt = bins_hittime_prompt[-1]

        # sum up the values of npe_per_hittime_prompt to get the total number of pe of the prompt signal:
        number_pe_prompt = np.sum(npe_per_hittime_prompt)

        # convert the total number of pe to quenched deposited energy in MeV:
        Qedep_real = NC_background_functions.conversion_npe_to_evis(number_pe_prompt)

        if Qedep_real == 0:
            # no prompt signal:
            # set Qedep_smeared = 0
            Qedep_smeared = 0.0
        else:
            # smear Qedep_real with energy resolution:
            # get the value of sigma of energy resolution for value of Qedep_real:
            sigma_energy_resolution = NC_background_functions.energy_resolution(Qedep_real)
            # generate normal distributed random number with mean = Qedep_real and sigma = sigma_energy_resolution:
            Qedep_smeared = np.random.normal(Qedep_real, sigma_energy_resolution)

        """ prompt energy cut: """
        # check efficiency because of energy resolution:
        if prompt_energy_min <= Qedep_smeared <= prompt_energy_max and (Qedep_real < prompt_energy_min or
                                                                        Qedep_real > prompt_energy_max):
            # smeared energy passes cut, but "real" energy would be rejected (event is counted too much):
            number_prompt_energy_toomuch += 1
        elif prompt_energy_min < Qedep_real < prompt_energy_max and (Qedep_smeared < prompt_energy_min or
                                                                     Qedep_smeared > prompt_energy_max):
            # smeared energy is rejected by cut, but real energy would pass cut (event is counted too less):
            number_prompt_energy_tooless += 1

        # do prompt energy cut
        if prompt_energy_min <= Qedep_smeared <= prompt_energy_max:
            # event passes prompt energy cut:
            number_prompt_energy_pass += 1
        else:
            # event is rejected by prompt energy cut:
            number_prompt_energy_rejected += 1
            number_rejected += 1

            if Qedep_smeared < prompt_energy_min:
                number_prompt_energy_rejected_min += 1
            else:
                number_prompt_energy_rejected_max += 1

            # event rejected -> go to next event:
            continue

        """ analyze delayed signal: """
        # get index, where delayed time window starts:
        index_min_hittime_delayed = index_max_hittime_prompt+2
        # get time corresponding to index_min_hittime_delayed in ns:
        min_time_delayed = bins_hittime[index_min_hittime_delayed]

        # get npe_per_hittime and bins_hittime for the delayed time window (from index_max_hittime_prompt+1 to
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
            # is_delayed_signal (=0 if no signal, =1 if there is delayed signal with
            # min_PE_delayed < nPE < max_PE_delayed), index_first_del (index after delayed
            # signal and start of next analysis), num_pe_delayed (nPE of delayed signal), begin_pulse (time, where
            # signal begins in ns), end_pulse (time, where signal ends in ns):
            is_delayed_signal, index_first_del, num_pe_delayed, begin_pulse, end_pulse = \
                NC_background_functions.analyze_delayed_signal(npe_per_hittime_delayed, bins_hittime_delayed,
                                                               index_first_del, threshold1_del, threshold2_del,
                                                               min_PE_delayed, max_PE_delayed, event_id)

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

        """ time cut: """
        # preallocate array, where indices of number_pe_delayed_array, begin_pulse_array, end_pulse_array, that pass
        # time cut are stored (important for further cuts):
        array_index_pass_time_cut = np.array([])

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
                    break

            if flag_pass_NeutronCaptureT == 0:
                flag_rejected_NeutronCaptureT = 1

        # Are there signals with begin_pulse and end_pulse inside delayed time window (min_time_delayed to
        # time_cut_max)?
        if len(number_pe_delayed_array) == 0:
            # no signal in delayed time window (either there is no delayed signal at all or the delayed signal is before
            # min_time_delayed and therefore inside prompt time window (overlap with prompt signal)).
            # event is rejected by time cut:
            flag_rejected_begin_end_pulse = 1

            # print("---------------ERROR: no delayed signal in event {0:d}, file = {1:d}".format(event_id, filenumber))
            h1 = plt.figure(1)
            plt.step(bins_hittime[:-1], npe_per_hittime)
            plt.xlabel("hit-time in ns")
            plt.ylabel("number of p.e. per bin (bin-width = {0:0.2f} ns)".format(binwidth))
            plt.title("Hit-time distribution of event {0:d}\nNO delayed signal".format(event_id))
            plt.xlim(xmin=min_time_total, xmax=min_time_delayed+1000)
            plt.legend()
            plt.grid()
            plt.savefig(output_path + "no_delayed_signal/file{1:d}_evt{0:d}_no_delayed_signal.png"
                        .format(event_id, filenumber))
            plt.close()

        else:
            # there is at least one signal in delayed time window. check, if begin_pulse and end_pulse of one signal is
            # inside time window (min_time_delayed to time_cut_max).
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
            # event is rejected by begin and end time, but would pass NeutronCaptureT time cut -> event counted too less
            number_time_tooless += 1

        # check time cut with begin time and end time:
        if flag_pass_begin_end_pulse == 1:
            # event pass time cut:
            number_time_pass += 1
        elif flag_rejected_begin_end_pulse == 1:
            # event is rejected by time cut:
            number_time_rejected += 1
            number_rejected += 1
            # go to next event:
            continue

        """ delayed energy cut: """
        # set flags to determine delayed energy cut efficiency:
        flag_pass_smeared = 0
        flag_rejected_smeared = 0
        flag_pass_real = 0
        flag_rejected_real = 0

        # loop over array_index_pass_time_cut:
        for index in range(len(array_index_pass_time_cut)):

            # get nPE of delayed signal corresponding to array_index_pass_time_cut[index] ("real" energy in PE):
            nPE_del_real = number_pe_delayed_array[int(array_index_pass_time_cut[index])]

            # get sigma of normal distribution because of energy resolution:
            sigma_energy_resolution_pe = NC_background_functions.energy_resolution_pe(nPE_del_real)
            # smear nPE_del_real with sigma_energy_resolution_pe to get reconstructed energy:
            nPE_del_smeared = np.random.normal(nPE_del_real, sigma_energy_resolution_pe)

            if min_PE_delayed <= nPE_del_real <= max_PE_delayed:
                # signal with "real" energy would pass delayed energy cut:
                flag_pass_real += 1
            else:
                # signal with "real" energy would be rejected by delayed energy cut:
                flag_rejected_real += 1

            if min_PE_delayed <= nPE_del_smeared <= max_PE_delayed:
                # signal with smeared energy pass delayed energy cut
                flag_pass_smeared += 1

            else:
                # signal with smeared energy is rejected by delayed energy cut:
                flag_rejected_smeared += 1

        # check efficiency:
        if flag_pass_smeared > 0 and flag_pass_real == 0:
            # event with smeared energy passes energy cut, but real energy would be rejected -> event counted too much:
            number_delayed_energy_toomuch += 1
        elif flag_pass_smeared == 0 and flag_pass_real > 0:
            # event with smeared energy is rejected by energy cut, but real energy would pass -> event counted too less
            number_delayed_energy_tooless += 1

        # do delayed energy cut:
        if flag_pass_smeared > 0:
            # at least on signal that pass delayed energy cut:
            number_delayed_energy_pass += 1
        else:
            # no signal pass delayed energy cut:
            number_delayed_energy_rejected += 1
            number_rejected += 1
            # go to next event:
            continue

        """ neutron multiplicity cut: """
        # check efficiency:
        # Are there events with only 1 signal with smeared energy, but more than 1 signal with real energy?
        # Are there events with more than 1 signal with smeared energy, but only 1 signal with real energy?
        if flag_pass_smeared == 1 and flag_pass_real > 1:
            # event with smeared energy pass multiplicity cut, but real energy would be rejected -> events counted too
            # much:
            number_n_mult_toomuch += 1
        elif flag_pass_smeared > 1 and flag_pass_real == 1:
            # event with smeared energy is rejected by multiplicity cut, but real energy would pass -> events counted
            # too less:
            number_n_mult_tooless += 1

        # do neutron multiplicity cut:
        if flag_pass_smeared == 1:
            # only 1 signal in delayed time window:
            number_n_mult_pass += 1
        else:
            # more than 1 signal in delayed time window:
            number_n_mult_rejected += 1
            number_rejected += 1
            # go to next event
            continue

        """ distance cut: """
        # preallocate number of signals with NeutronCaptureT in delayed time window:
        number_NeutronCaptureT = 0
        index_NeutronCaptureT = 0
        # loop over entries in NeutronCaptureT_array:
        for index in range(len(NeutronCaptureT_array)):
            # check if NeutronCaptureT inside delayed time window:
            if min_time_delayed <= NeutronCaptureT_array[index] <= time_cut_max:
                number_NeutronCaptureT += 1
                index_NeutronCaptureT = index

        # check NeutronN:
        if NeutronN == 0 and number_NeutronCaptureT == 0:
            print("+-+-+- ERROR: NeutronN = 0 when applying distance cut (event = {0:d}, file = {1:d})"
                  .format(event, filenumber))
        elif NeutronN > 1 and number_NeutronCaptureT > 1:
            print("----++++ ERROR: more than 1 nCapture inside delayed time window when applying distance cut "
                  "((event = {0:d}, file = {1:d})"
                  .format(event, filenumber))
        else:
            # 3 possibilities left:
            # 1.:   NeutronN = 1 and number_NeutronCaptureT > 1 -> not possible because 1 neutron can not be captured
            #       2 times
            # 2.:   NeutronN > 1 and number_NeutronCaptureT = 1 -> no problem because only one nCapture specified
            #       by index_NeutronCaptureT
            # 3.:   NeutronN = 1 and number_NeutronCaptureT = 1 -> no problem because only one nCapture specified
            #       by index_NeutronCaptureT

            # set Qedep of delayed signal to 2.2 MeV (neutron capture on Hydrogen):
            Qedep_capture = 2.2

            # get the start position of neutron capture in mm:
            x_ncapture = float(rtree_ncapture.GetBranch("NCStartX").GetLeaf("NCStartX").GetValue(index_NeutronCaptureT))
            y_ncapture = float(rtree_ncapture.GetBranch("NCStartY").GetLeaf("NCStartY").GetValue(index_NeutronCaptureT))
            z_ncapture = float(rtree_ncapture.GetBranch("NCStartZ").GetLeaf("NCStartZ").GetValue(index_NeutronCaptureT))

            # do vertex reconstruction of neutron capture position with function position_smearing():
            x_reco_ncapture = NC_background_functions.position_smearing(x_ncapture, Qedep_capture)
            y_reco_ncapture = NC_background_functions.position_smearing(y_ncapture, Qedep_capture)
            z_reco_ncapture = NC_background_functions.position_smearing(z_ncapture, Qedep_capture)

            # calculate distance of reconstructed neutron capture position to detector center in mm:
            r_reconstructed_ncapture = np.sqrt(x_reco_ncapture**2 + y_reco_ncapture**2 + z_reco_ncapture**2)

            # calculate distance between real initial position and real nCapture position in mm:
            distance_real = np.sqrt((init_x - x_ncapture)**2 + (init_y - y_ncapture)**2 + (init_z - z_reconstructed)**2)

            # calculate distance between reconstructed initial position and reconstructed nCapture position in mm:
            distance_reco = np.sqrt((x_reconstructed - x_reco_ncapture)**2 + (y_reconstructed - y_reco_ncapture)**2 +
                                    (z_reconstructed - z_reco_ncapture)**2)

            # check efficiency of distance cut:
            if distance_reco < distance_cut and distance_real >= distance_cut:
                # reco distance pass the cut, but real distance would be rejected -> event is counted too much:
                number_distance_cut_toomuch += 1
            elif distance_reco >= distance_cut and distance_real < distance_cut:
                # reco distance is rejected by cut, but real distance would pass -> event is counted too less:
                number_distance_cut_tooless += 1

            # do distance cut:
            if distance_reco < distance_cut:
                # event pass distance cut:
                number_distance_cut_pass += 1
            else:
                # event is rejected by distance cut:
                number_distance_cut_rejected += 1
                number_rejected += 1
                continue

            """ volume cut on delayed signal """
            # You only reach this line if the event passes all cuts (volume, prompt energy, time, delayed energy,
            # multiplicity and distance cut).
            # As an additional volume cut, the reconstructed position of neutron capture must also be considered,
            # because prompt AND delayed signal must be within fiducial volume:
            if r_reconstructed_ncapture >= R_cut_mm:
                # delayed signal is outside fiducial volume -> event is rejected:
                number_volume_rejected_delayed += 1
                number_rejected += 1
                # go to next event:
                continue

        # increment number_IBDlike_events and append event_id and Qedep_smeared of event that passes all cuts to arrays
        # -> represent the IBD-like events:
        number_IBDlike_events += 1
        eventID_pass_array = np.append(eventID_pass_array, event_id)
        Evis_pass_array = np.append(Evis_pass_array, Qedep_smeared)

        # save hittime distribution of the IBD-like event to png and txt file:
        h2 = plt.figure(2)
        plt.step(bins_hittime_prompt, npe_per_hittime_prompt, label="number of pe = {0:d}\n"
                                                                    "visible energy = {1:.3f} MeV"
                 .format(number_pe_prompt, Qedep_smeared))
        plt.xlabel("hit-time in ns")
        plt.ylabel("number of p.e. per bin (bin-width = {0:0.2f} ns)".format(binwidth))
        plt.title("Hit-time distribution of prompt time window of event {0:d}".format(event_id))
        plt.xlim(xmin=min_time_prompt, xmax=max_time_prompt)
        plt.legend()
        plt.grid()
        plt.savefig(output_path + "file{1:d}_evt{0:d}_prompt_signal.png".format(event_id, filenumber))
        plt.close()

        # save npe_per_hittime_prompt to txt file:
        # build list, where 0th entry is start-hittime in ns, 1st entry is last-hittime in ns, 2nd entry is binwidth in
        # ns and the following entries are nPE of each hittime-bin of prompt signal:
        npe_per_hittime_prompt_save = [min_time_prompt, max_time_prompt, binwidth]
        npe_per_hittime_prompt_save.extend(npe_per_hittime_prompt)
        np.savetxt(output_path + "file{0:d}_evt{1:d}_prompt_signal.txt".format(filenumber, event_id),
                   npe_per_hittime_prompt_save, fmt='%1.2f',
                   header="Number of pe as function of the hittime of the prompt signal (time-of-flight correction "
                          "and TTS smearing) of file user_atmoNC_{0:d}.root,"
                          "\nevent {1:d}, {2}:"
                          "\ntime window of hittime: from {3:.3f} ns to {4:.3f} ns with bin-width = {5:0.3f} ns,"
                          "\nEnergy cut on prompt signal is applied: {6:0.1f} MeV <= E_vis <= {7:0.1f} MeV,"
                          "\nConversion function E_vis = 0.0007483 * nPE:"
                   .format(filenumber, event_id, now, min_time_prompt, max_time_prompt, binwidth, prompt_energy_min,
                           prompt_energy_max))

    """ save event ID of IBD-like events to txt file together with information about cut parameters: """
    np.savetxt(output_path + "evtID_IBDlike_{0:d}.txt".format(filenumber),
               eventID_pass_array, fmt='%i',
               header="event ID's of IBD-like NC events that pass all cuts (volume cut, prompt energy cut, time cut,\n"
                      "delayed energy cut, neutron multiplicity cut and distance cut) with script "
                      "analyze_NC_events_v1.py ({0}):\n"
                      "{1:d} events analyzed from file: user_atmoNC_{2:d}.root {3:d}, IBD-like events after all cuts.\n"
                      "Cut Parameters:\n"
                      "volume cut: {4:d} mm (reconstructed initial position and reconstructed nCapture position "
                      "(smeared by vertex resolution sigma = 120mm/sqrt(E[MeV])),\n"
                      "prompt energy cut: {5:0.1f} MeV <= Qedep_smeared <= {6:0.1f} MeV (Qedep_smeared converted from "
                      "nPE of hittime distribution and smeared with energy resolution),\n"
                      "time cut (begin and end of delayed signal inside delayed time window from around {7:0.1f} ns to "
                      "{8:0.1f} ns),\n"
                      "delayed energy cut: {9:0.1f} PE <= nPE_del_smeared <= {10:0.1f} PE (nPE_del_smeared from hittime"
                      " distribution in delayed time window and smeared with energy resolution),\n"
                      "neutron multiplicity cut: only 1 delayed signal with correct delayed energy in delayed time "
                      "window,\n"
                      "distance cut (distance between reconstructed initial position and reconstructed neutron "
                      "capture position): {11:d} mm:"
               .format(now, Number_entries_input, filenumber, number_IBDlike_events, R_cut_mm, prompt_energy_min,
                       prompt_energy_max, time_cut_min, time_cut_max, min_PE_delayed, max_PE_delayed,
                       distance_cut))

    """ save visible energy of prompt signal of IBD-like events to txt file: """
    np.savetxt(output_path + "Evis_file{0:d}.txt".format(filenumber), Evis_pass_array, fmt="%.3f",
               header="Visible energy (Qedep_smeared converted from nPE of hittime distribution and smeared with "
                      "energy resolution)\n"
                      "of prompt signal for IBD-like NC events from file user_atmoNC_{0:d}.root ({1}).\n"
                      "Qedep_smeared is analyzed with script analyze_NC_events_v1.py.\n"
                      "Time window of prompt signal is defined from {2:0.2f} ns to around {3:0.2f} ns,\n"
                      "Energy cut on prompt signal is applied: {4:0.1f} MeV <= E_vis <= {5:0.1f} MeV,\n"
                      "Conversion function E_vis = 0.0007483 * nPE:"
               .format(filenumber, now, min_time_prompt, time_cut_min, prompt_energy_min, prompt_energy_max))

    """ save information about the number of events in different cases to txt file for each analyzed file: """
    np.savetxt(output_path + "information_number_of_events/numbers_file{0:d}.txt".format(filenumber),
               np.array([number_total_events, number_IBDlike_events, number_rejected, number_without_particles,
                         number_volume_pass_prompt, number_volume_rejected_prompt, number_volume_rejected_delayed,
                         number_prompt_energy_pass, number_prompt_energy_rejected, number_prompt_energy_rejected_min,
                         number_prompt_energy_rejected_max, number_prompt_energy_toomuch, number_prompt_energy_tooless,
                         number_time_pass, number_time_rejected, number_time_toomuch, number_time_tooless,
                         number_delayed_energy_pass, number_delayed_energy_rejected, number_delayed_energy_toomuch,
                         number_delayed_energy_tooless, number_n_mult_pass, number_n_mult_rejected,
                         number_n_mult_toomuch, number_n_mult_tooless, number_distance_cut_pass,
                         number_distance_cut_rejected, number_distance_cut_toomuch, number_distance_cut_tooless]),
               fmt="%i",
               header="Information about number of events for different cases of user_atmoNC_{0:d}.root {1}:"
                      "\nnumber of analyzed events, number of events that pass all cuts (IBD-like events),"
                      "\nnumber of events that are rejected by cuts, number of events without initial particles,"
                      "\nnumber of events that pass volume cut on prompt position, "
                      "number of events that are rejected by volume cut on prompt position,"
                      "\nnumber of events that are rejected by volume cut on delayed position, "
                      "number of events that pass prompt energy cut (and volume cut),"
                      "\nnumber of events that are rejected by prompt energy cut, "
                      "number of events that are rejected by prompt energy cut (E < prompt_energy_min),"
                      "\nnumber of events that are rejected by prompt energy cut (E > prompt_energy_max), "
                      "number of events, where smeared energy pass the cut, but real energy would be rejected "
                      "(counted too much),"
                      "\nnumber of events, where smeared energy is rejected by cut, but real energy would pass the cut "
                      "(counted too less), number of events that pass time cut (and volume and prompt energy cut),"
                      "\nnumber of events that are rejected by time cut, number of events that pass time cut, but "
                      "would be rejected by nCaptureT (events counted too much),"
                      "\nnumber of events that are rejected by time cut, but nCaptureT would pass (events counted too "
                      "less), number of events that pass delayed energy cut (and volume, prompt energy and time cut),"
                      "\nnumber of events that are rejected by delayed energy cut, number of events, where smeared "
                      "energy pass the cut, but real energy would be rejected (counted too much),"
                      "\nnumber of events, where smeared energy is rejected by cut, but real energy would pass the cut "
                      "(counted too less), number of events that pass neutron multiplicity cut (and volume, "
                      "prompt energy, time and delayed energy cut),"
                      "\nnumber of events that do not pass neutron multiplicity cut, number of events, where smeared "
                      "energy pass multiplicity cut, but real energy would be rejected by multiplicity cut (counted "
                      "too much),"
                      "\nnumber of events, where smeared energy is rejected by multiplicity cut, but real energy would "
                      "pass multiplicity cut (counted too less), number of events that pass distance cut (and volume, "
                      "prompt energy, time, delayed energy and multiplicity cut),"
                      "\nnumber of events that are rejected by distance cut, number of events, where reconstructed "
                      "distance pass the cut, but real (MC truth) distance is rejected (counted too much),"
                      "\nnumber of events, where reconstructed distance is rejected by cut, but real distance pass "
                      "(counted too less):"
               .format(filenumber, now))

# print information about the different cut efficiencies:
print("\n-----------------------------------------------------------")
print("results from user_atmoNC_{0:d}.root to user_atmoNC_{1:d}.root\n".format(start_number, stop_number))
print("\n------------")
print("number of analyzed events = {0:d}".format(number_total_events))
print("number of events that pass all cuts = {0:d}".format(number_IBDlike_events))
print("number of events that are rejected by cuts = {0:d}".format(number_rejected))
print("\n------------")
print("number of events without initial particles = {0:d}".format(number_without_particles))
print("\n------------")
print("volume cut:")
print("number of events that pass volume cut on prompt position = {0:d}".format(number_volume_pass_prompt))
print("number of events that are rejected by volume cut on prompt position = {0:d}"
      .format(number_volume_rejected_prompt))
print("number of events that are rejected by volume cut on delayed position = {0:d}"
      .format(number_volume_rejected_delayed))
print("\n------------")
print("prompt energy cut:")
print("number of events that pass prompt energy cut (and volume cut) = {0:d}".format(number_prompt_energy_pass))
print("number of events that are rejected by prompt energy cut = {0:d}".format(number_prompt_energy_rejected))
print("number of events that are rejected by prompt energy cut (E < prompt_energy_min) = {0:d}"
      .format(number_prompt_energy_rejected_min))
print("number of events that are rejected by prompt energy cut (E > prompt_energy_max) = {0:d}"
      .format(number_prompt_energy_rejected_max))
print("number of events, where smeared energy pass the cut, but real energy would be rejected (counted too much) "
      "= {0:d}".format(number_prompt_energy_toomuch))
print("number of events, where smeared energy is rejected by cut, but real energy would pass the cut (counted too less)"
      " = {0:d}".format(number_prompt_energy_tooless))
print("\n------------")
print("time cut:")
print("number of events that pass time cut (and volume and prompt energy cut) = {0:d}".format(number_time_pass))
print("number of events that are rejected by time cut = {0:d}".format(number_time_rejected))
print("number of events that pass time cut, but would be rejected by nCaptureT (events counted too much) = {0:d}"
      .format(number_time_toomuch))
print("number of events that are rejected by time cut, but nCaptureT would pass (events counted too less) = {0:d}"
      .format(number_time_tooless))
print("\n------------")
print("delayed energy cut:")
print("number of events that pass delayed energy cut (and volume, prompt energy and time cut) = {0:d}"
      .format(number_delayed_energy_pass))
print("number of events that are rejected by delayed energy cut = {0:d}".format(number_delayed_energy_rejected))
print("number of events, where smeared energy pass the cut, but real energy would be rejected (counted too much) = "
      "{0:d}".format(number_delayed_energy_toomuch))
print("number of events, where smeared energy is rejected by cut, but real energy would pass the cut (counted too less)"
      " = {0:d}".format(number_delayed_energy_tooless))
print("\n------------")
print("number of events that pass neutron multiplicity cut (and volume, prompt energy, time and delayed energy cut) "
      "= {0:d}".format(number_n_mult_pass))
print("number of events that do not pass neutron multiplicity cut = {0:d}".format(number_n_mult_rejected))
print("number of events, where smeared energy pass multiplicity cut, but real energy would be rejected by multiplicity "
      "cut (counted too much) = {0:d}".format(number_n_mult_toomuch))
print("number of events, where smeared energy is rejected by multiplicity cut, but real energy would pass multiplicity "
      "cut (counted too less) = {0:d}".format(number_n_mult_tooless))
print("\n------------")
print("number of events that pass distance cut (and volume, prompt energy, time, delayed energy and multiplicity cut) "
      "= {0:d}".format(number_distance_cut_pass))
print("number of events that are rejected by distance cut = {0:d}".format(number_distance_cut_rejected))
print("number of events, where reconstructed distance pass the cut, but real (MC truth) distance is rejected (counted "
      "too much) = {0:d}".format(number_distance_cut_toomuch))
print("number of events, where reconstructed distance is rejected by cut, but real distance pass (counted too less) "
      "= {0:d}".format(number_distance_cut_tooless))


