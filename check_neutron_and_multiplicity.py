""" Script to check neutrons simulated with JUNO detsim (neutron initial momentum uniformly distributed from 0.001 MeV
    to 30 MeV within radius of R<16m)

    To check the delayed cut (delayed energy cut, neutron multiplicity cut, time cut and distance cut):

    Four important things can be checked with this script:

    1.  Efficiency of the time cut between prompt and delayed signal. In this case there is no prompt signal, but,
        because it is the prompt signal, it is assumed to be at very early hittimes (around 0 ns).

    2.  Efficiency of the delayed energy cut:
        There are some problems with the first investigation of the delayed energy cut (based on 1.9 MeV, 2.2 MeV and
        2.5 MeV gammas), because many delayed signals of NC events have to small nPE compared to the result of the
        gamma simulation (see script OLD_check_delayed_energy.py and folder /output_gamma_2_2_MeV)

    3.  Efficiency of the neutron multiplicity cut (momentum of 0.001 MeV corresponds to minimum kinetic energy of a
        neutron after inverse beta decay with a neutrino with initial energy of 10 MeV; momentum of 28 MeV corresponds
        to the maximum kinetic energy a neutron can get after inverse beta decay with a neutrino with initial energy
        of 115 MeV)
        -> Can such neutrons produce further neutrons via scattering, which are also captured on H?
        -> Is it possible that such a neutron is not captured by H?
        These two questions are important for the efficiency of the neutron multiplicity cut (only 1 neutron capture
        in time window 500 ns to 1 ms)

    4.  Efficiency of the distance cut between initial position (position, where prompt signal deposits its energy) and
        position of the neutron capture. Both positions (initial position and start position of nCapture) must be
        smeared

    So, with this script, you can calculate the efficiency of the delayed cut for real IBD events.

    The efficiency of the preselection, which is done in preselection_detsim_user.py, and the efficiency of the
    delayed cut for IBD-like events (NC events) must be determined in another way.

"""
import datetime
import ROOT
import sys
import NC_background_functions
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

# get the date and time, when the script was run:
date = datetime.datetime.now()
now = date.strftime("%Y-%m-%d %H:%M")

# set the path of the input files:
input_path = "/home/astro/blum/juno/atmoNC/data_NC/output_neutron_multiplicity/"

# set path, where results should be saved:
output_path = input_path + "results/"

# start file and stop file:
file_start = 0
file_end = 500
# file_end = 13
# number of events per file (for simulation with optical processes):
number_evts_per_file = 100

""" define time window and radius cut: """
# set time window of whole signal in ns (start at 200 ns to avoid the analysis of the small prompt signal):
min_time = 200
max_time = 1001000
# set maximum of time window of delayed signal (neutron capture window) in ns:
max_time_ncap = 1000000
# set time in ns, where the neutron capture signal should start:
time_limit = 500
# Set bin-width of hittime histogram in ns:
binwidth = 5.0
# set the radius for the volume cut in mm:
r_cut = 16000
# set the distance for the distance cut in mm:
distance_cut = 1500

""" thresholds and cuts for ncapture signal: """
# Set threshold of number of PE per bin for possible delayed signal (bin-width = 5 ns):
threshold1_del = 50
# set threshold2 of number of PEs per bin (signal peak is summed as long as nPE is above threshold2):
threshold2_del = 0
# min and max number of PE for delayed energy cut (delayed energy cut: 1.9 MeV / 0.0007483 = 2573 PE ~ 2500 PE,
# 2.5 MeV / 0.0007384 = 3385 PE ~ 3400 PE):
min_PE_delayed = 2400
max_PE_delayed = 3400

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

""" preallocate variables: """
# preallocate the total number of events simulated:
number_of_events_total = 0
# preallocate the number of events inside fiducial volume (defined by r_cut) (equal to number of initial neutrons):
number_of_events = 0
# preallocate the number of events with at least one signal pulse in the delayed time window (from around 500 ns to
# max_time):
number_timecut_pass = 0
# preallocate the number of events with NO signal pulse in the delayed time window (from around 500 ns to max_time):
number_timecut_rejected = 0
# preallocate the number of events, where delayed pulse begins before time_limit:
number_timecut_rejected_min = 0
# preallocate the number of events, where delayed pulse ends after max_time_ncap:
number_timecut_rejected_max = 0
# preallocate number of events with at least on signal pulse that pass the delayed energy cut in the delayed time window
# (min_PE_delayed < E_d < max_PE_delayed):
number_delayed_energy_pass = 0
# preallocate number of events with NO signal pulse that pass the delayed energy cut in the delayed time window:
number_delayed_energy_rejected = 0
# preallocate number of events, where smeared energy pass the cut, but real energy would be rejected:
number_delayed_energy_toomuch = 0
# preallocate number of events, where smeared energy is rejected by cut, but real energy would pass the cut:
number_delayed_energy_tooless = 0
# preallocate number of events with only 1 signal pulse in delayed time window, that pass delayed energy cut:
number_n_mult_pass = 0
# preallocate number of events that do not pass neutron multiplicity cut:
number_n_mult_rejected = 0
# preallocate number of events with distance between reco. initial position to reco. nCapture position below
# distance_cut, that pass all cuts above:
number_distance_cut_pass = 0
# preallocate number of events, that pass all cuts above, except of the distance cut:
number_distance_cut_rejected = 0
# preallocate number of events with distance between real initial position to real nCapture position below distance cut,
# that pass all cuts above:
number_distance_cut_pass_MCtruth = 0
# preallocate number of events, that pass all cuts above, except the distance cut of the real MC truth position:
number_distance_cut_rejected_MCtruth = 0

number_pe0_qedep2 = 0
number_pe3000_qedep2 = 0
number_pe0_qedep5 = 0
number_pe6500_qedep5 = 0
number_pe15000_qedep11 = 0
number_no_ncapture = 0

# preallocate list, where number of pe (directly from nPhotons) is saved:
array_npe = []
# preallocate list, where number of pe from analyzing the corrected hittime distribution is saved:
array_npe_from_hittime = []
# preallocate list, where edep in MeV is saved:
array_edep = []
# preallocate list, where Qedep in MeV is saved:
array_Qedep = []

""" Analyze the file user_neutron_multiplicity_{}.root: """
# loop over files:
for filenumber in range(file_start, file_end+1, 1):
    # load file:
    rfile = ROOT.TFile(input_path + "user_neutron_multiplicity_{0:d}.root".format(filenumber))
    # print("... read {0}...".format(rfile))

    # get evt tree from TFile:
    rtree_evt = rfile.Get("evt")
    # get geninfo tree from TFile:
    rtree_geninfo = rfile.Get("geninfo")
    # get prmtrkdep tree from TFile:
    rtree_prmtrkdep = rfile.Get("prmtrkdep")
    # get nCapture tree from TFile:
    rtree_ncapture = rfile.Get("nCapture")

    # get the number of events in the 'geninfo' Tree:
    num_evts = rtree_geninfo.GetEntries()
    if num_evts != number_evts_per_file:
        sys.exit("ERROR: number of events differ in file {0}".format(rfile))

    # add num_evts to number_of_events_total:
    number_of_events_total += num_evts

    # loop over events:
    for event in range(num_evts):
        # set flag, that event passes time cut (time_limit <= nCaptureT <= max_time):
        flag_pass_timecut = False

        # get current event:
        rtree_prmtrkdep.GetEntry(event)
        # get deposit energy of initial neutron in MeV (neutrons deposit most of the energy while being captured):
        Qedep_capture = float(rtree_prmtrkdep.GetBranch('Qedep').GetLeaf('Qedep').GetValue())

        # get current event:
        rtree_geninfo.GetEntry(event)
        # get initial x,y,z position in mm:
        x_init = float(rtree_geninfo.GetBranch('InitX').GetLeaf('InitX').GetValue())
        y_init = float(rtree_geninfo.GetBranch('InitY').GetLeaf('InitY').GetValue())
        z_init = float(rtree_geninfo.GetBranch('InitZ').GetLeaf('InitZ').GetValue())

        # do vertex reconstruction with function position_smearing() for distance cut:
        # Smear x,y and z position of the initial position (returns reconstructed position in mm)
        # (for Qedep use random number from uniform distribution between 10 MeV and 100 MeV. This represents the prompt
        # energy of a positron like in a real IBD event, since in user_neutron_multiplicity.root only the neutron is
        # simulated):
        Qedep_init = np.random.uniform(10, 100)
        x_reco_init = NC_background_functions.position_smearing(x_init, Qedep_init)
        y_reco_init = NC_background_functions.position_smearing(y_init, Qedep_init)
        z_reco_init = NC_background_functions.position_smearing(z_init, Qedep_init)

        # get nCapture tree:
        rtree_ncapture.GetEntry(event)
        # get number of neutron captures:
        NeutronN = int(rtree_ncapture.GetBranch('NeutronN').GetLeaf('NeutronN').GetValue())

        # set variables to zeros:
        nCaptureT = 0.0
        x_reco_ncapture = -17000
        y_reco_ncapture = 0
        z_reco_ncapture = 0

        # check NeutronN:
        if NeutronN < 1:
            # no neutron capture in event:
            number_no_ncapture += 1
            print("-------------no neutron capture in event {0:d} in file {1}".format(event, rfile))
        elif NeutronN > 1:
            # more than one neutron capture in event:
            print("+++++++++++++more than 1 neutron capture in event {0:d} in file {1}".format(event, rfile))
        else:
            # NeutronN == 1 -> 1 neutron capture in event.

            # check if captured neutron was the initial neutron:
            NeutronTrkid = int(rtree_ncapture.GetBranch('NeutronTrkid').GetLeaf('NeutronTrkid').GetValue(0))
            if NeutronTrkid != 1:
                print("captured neutron is not initial neutron (event = {0:d}, file {1}".format(event, rfile))

            # check neutron capture time in ns:
            nCaptureT = float(rtree_ncapture.GetBranch("NeutronCaptureT").GetLeaf("NeutronCaptureT").GetValue(0))

            # get the start position of neutron capture in mm:
            x_ncapture = float(rtree_ncapture.GetBranch("NCStartX").GetLeaf("NCStartX").GetValue(0))
            y_ncapture = float(rtree_ncapture.GetBranch("NCStartY").GetLeaf("NCStartY").GetValue(0))
            z_ncapture = float(rtree_ncapture.GetBranch("NCStartZ").GetLeaf("NCStartZ").GetValue(0))

            # do vertex reconstruction of neutron capture position with function position_smearing():
            x_reco_ncapture = NC_background_functions.position_smearing(x_ncapture, Qedep_capture)
            y_reco_ncapture = NC_background_functions.position_smearing(y_ncapture, Qedep_capture)
            z_reco_ncapture = NC_background_functions.position_smearing(z_ncapture, Qedep_capture)

        # calculate distance of neutron capture to detector center in mm:
        r_reco_ncapture = np.sqrt(x_reco_ncapture**2 + y_reco_ncapture**2 + z_reco_ncapture**2)

        # check volume cut:
        if r_reco_ncapture >= r_cut:
            # event outside fiducial volume:
            continue
        else:
            # event inside fiducial volume:
            number_of_events += 1

        """ calculate the real hittime distribution (time of flight correction with reconstructed position of neutron 
        capture and time smearing with TTS for each hit): """
        # get current event:
        rtree_evt.GetEntry(event)
        # get number of photons:
        n_photons = int(rtree_evt.GetBranch("nPhotons").GetLeaf("nPhotons").GetValue())

        # preallocate empty array to build default hittime-histogram:
        hittime_array = []

        for index in range(n_photons):
            # get number of PE per photon:
            nPE = int(rtree_evt.GetBranch("nPE").GetLeaf("nPE").GetValue(index))
            if nPE != 1:
                # more than 1 PE per photon
                sys.exit("ERROR: more than 1 PE per photon (event {0:d}, file {1})".format(event, rfile))

            # get the pmtID of the hit PMT:
            pmtID = int(rtree_evt.GetBranch('pmtID').GetLeaf('pmtID').GetValue(index))

            """ time of flight correction: """
            # get hittime of PMT from tree in ns:
            hittime = float(rtree_evt.GetBranch('hitTime').GetLeaf('hitTime').GetValue(index))

            # get position of the PMT with specific pmtID (pmtID is ascending number from 0 to 17738 (17739 large PMTs)
            # and from 300000 to 336571 (36572 small PMTs)).
            # For large PMTs -> For 20inch PMTs, the pmtID is equal to index of x,y,z_pos_pmt array.
            # For small PMTs -> For 3inch PMTs, the pmtID - (300000 - 17739) is equal to index of x,y,z_pos_pmt array.
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

            # calculate distance between reconstructed position neutron capture and position of PMT (in mm):
            distance_tof = np.sqrt((x_reco_ncapture - x_pmt)**2 + (y_reco_ncapture - y_pmt)**2 +
                                   (z_reco_ncapture - z_pmt)**2)

            # calculate time of flight in ns:
            time_of_flight = distance_tof / c_effective

            """ time resolution of PMT: """
            # get time resolution of PMT with specific pmtID (pmtID is ascending number from 0 to 17738 (17739 large
            # PMTs)) -> For 20inch PMTs, the pmtID is equal to index of sigma_time_20inch array.
            # check if PMT is 20 inch or 3inch (pmtID < 50000 means 20inch PMT):
            if pmtID < 50000:
                # 20inch PMT:
                # get time resolution (sigma) of PMT in ns from array:
                sigma_pmt = sigma_time_20inch[pmtID]

            else:
                # 3inch PMT:
                sigma_pmt = sigma_time_3inch

            # consider time resolution of PMT by generating normal distributed random number with mu = hittime and
            # sigma = sigma_pmt (only the hittime at the PMT must be smeared, not the time-of-flight):
            hittime_tts = np.random.normal(hittime, sigma_pmt)

            """ calculate the 'real' hittime of the photon in ns: """
            hittime_real = hittime_tts - time_of_flight

            # append real hittime to array:
            hittime_array.append(hittime_real)

        # analyze hittime distribution (the same way like in prompt_signal_preselected_evts.py):
        # build histogram, where hittimes are saved:
        # set bin-edges of hittime histogram in ns (whole time window from min_time to max_time):
        bins_hittime = np.arange(min_time, max_time + 2 * binwidth, binwidth)
        # build hittime histogram:
        npe_per_hittime_all, bin_edges_hittime_all = np.histogram(hittime_array, bins_hittime)

        # analyze the whole time window (not only neutron capture time window) for the time cut:
        num_pulses, index_test1, num_of_ncaptures, begin_time_pulse, end_time_pulse = \
            NC_background_functions.analyze_delayed_signal(npe_per_hittime_all, bin_edges_hittime_all, 0,
                                                           threshold1_del, threshold2_del, min_PE_delayed,
                                                           max_PE_delayed, event)

        # do time cut:
        if begin_time_pulse <= time_limit:
            # pulse begins before time_limit:
            number_timecut_rejected += 1
            number_timecut_rejected_min += 1
        elif end_time_pulse >= max_time_ncap:
            # pulse ends after time_limit:
            number_timecut_rejected += 1
            number_timecut_rejected_max += 1
        else:
            # event passes time cut:
            flag_pass_timecut = True
            number_timecut_pass += 1

        # compare time cut done with begin_time_pulse and end_time_pulse with the time cut done with nCapture T:
        if begin_time_pulse != 0 and begin_time_pulse <= time_limit and time_limit < nCaptureT < max_time_ncap:
            print("----- event is rejected by time cut (begin_pulse = {0:.0f} ns), but would pass nCaptureT "
                  "cut (nCaptureT = {1:.0f} ns)".format(begin_time_pulse, nCaptureT))
        elif begin_time_pulse != 0 and time_limit < begin_time_pulse < max_time_ncap and nCaptureT <= time_limit:
            print("+++++ event pass time cut (begin_pulse = {0:.0f} ns), but would be rejected by nCaptureT cut "
                  "(nCaptureT = {1:.0f} ns)".format(begin_time_pulse, nCaptureT))

        if end_time_pulse != 2000000 and end_time_pulse >= max_time_ncap and time_limit < nCaptureT < max_time_ncap:
            print("------------ event is rejected by time cut (end_pulse = {0:.0f} ns), but would pass nCaptureT cut "
                  "(nCaptureT = {1:.0f} ns)".format(end_time_pulse, nCaptureT))
        elif end_time_pulse != 2000000 and time_limit < end_time_pulse < max_time_ncap and nCaptureT >= max_time_ncap:
            print("++++++++++++ event pass time cut (end_pulse = {0:.0f} ns), but would be rejected by nCaptureT cut "
                  "(nCaptureT = {1:.0f} ns)".format(end_time_pulse, nCaptureT))

        #################################################
        # plt.plot(bin_edges_hittime_all[:-1], npe_per_hittime_all)
        # plt.show()

        # get index of bins_hittime corresponding to time_limit_prompt:
        index_time_limit = int((time_limit - min_time) / binwidth)
        # get index of bins_hittime corresponding to max_time_ncap:
        index_max_time_ncap = int((max_time_ncap - min_time) / binwidth)

        # take only hittime histogram from index_time_limit to index_max_time_ncap:
        bin_edges_hittime = bin_edges_hittime_all[index_time_limit:(index_max_time_ncap)]
        npe_per_hittime = npe_per_hittime_all[index_time_limit:(index_max_time_ncap)]
        index_test = 0

        ############################################
        # plt.plot(bin_edges_hittime, npe_per_hittime)
        # plt.show()

        number_nCapture_pass_e_cut = 0
        number_pe_ncapture = 0
        number_nCapture_pass_e_cut_smeared = 0
        number_pe_ncapture_smeared = 0

        # analyze neutron capture signal (num_n_captures: number of signal with correct energy in event (0 or 1);
        # index_test: index, where analysis of hittime starts;
        # number_pe_ncapture: number of pe in neutron capture signal peak):
        while index_test < len(npe_per_hittime):
            # analyze delayed signal until you reach the end of the time window:
            num_n_captures, index_test, num_pe_ncapture, begin_pulse, end_pulse = \
                NC_background_functions.analyze_delayed_signal(npe_per_hittime, bin_edges_hittime, index_test,
                                                               threshold1_del,
                                                               threshold2_del, min_PE_delayed, max_PE_delayed, event)

            number_nCapture_pass_e_cut += num_n_captures
            number_pe_ncapture += num_pe_ncapture

        if number_pe_ncapture != 0:
            # apply energy resolution on the number of pe from ncapture:
            # get sigma for this certain number of pe:
            sigma_nPE = NC_background_functions.energy_resolution_pe(number_pe_ncapture)
            # generate normal distributed random number with mean = number_pe_ncapture and sigma = sigma_nPE:
            number_pe_ncapture_smeared = np.random.normal(number_pe_ncapture, sigma_nPE)
            # check if smeared energy would pass the delayed energy cut:
            if min_PE_delayed < number_pe_ncapture_smeared < max_PE_delayed:
                number_nCapture_pass_e_cut_smeared = 1
            else:
                number_nCapture_pass_e_cut_smeared = 0

        # do event pass time cut:
        if flag_pass_timecut:

            # check, if there is at least 1 n capture that pass the delayed energy cut in the delayed time window:
            if number_nCapture_pass_e_cut_smeared > 0:
                number_delayed_energy_pass += 1

                # check, if there is only 1 signal pulse in delayed time window, that pass delayed energy cut:
                if number_nCapture_pass_e_cut_smeared == 1:
                    # only 1 neutron capture in time window with correct energy:
                    number_n_mult_pass += 1

                    # check, if distance between reco. initial position to reco. nCapture position below
                    # distance_cut of event, that passes all cuts above (n-mult, time, energy):
                    # calculate distance between reco. initial position and reco. nCapture position:
                    distance = np.sqrt((x_reco_init - x_reco_ncapture) ** 2 + (y_reco_init - y_reco_ncapture) ** 2 +
                                       (z_reco_init - z_reco_ncapture) ** 2)

                    if distance < distance_cut:
                        number_distance_cut_pass += 1
                    else:
                        number_distance_cut_rejected += 1

                    # check, if distance between real initial position to real nCapture position (from MC truth) is
                    # below distance_cut of event:
                    # calculate distance between real initial and real nCapture position:
                    distance_MCtruth = np.sqrt((x_init - x_ncapture)**2 + (y_init - y_ncapture)**2 +
                                               (z_init - z_ncapture)**2)
                    if distance_MCtruth < distance_cut:
                        number_distance_cut_pass_MCtruth += 1
                    else:
                        number_distance_cut_rejected_MCtruth += 1

                else:
                    number_n_mult_rejected += 1

            else:
                number_delayed_energy_rejected += 1

            if number_nCapture_pass_e_cut_smeared == 0 and number_nCapture_pass_e_cut == 1:
                # smeared energy is rejected, but real energy would pass:
                number_delayed_energy_tooless += 1
            elif number_nCapture_pass_e_cut_smeared == 1 and number_nCapture_pass_e_cut == 0:
                # smeared energy pass, but real energy would be rejected:
                number_delayed_energy_toomuch += 1

        if 2 < Qedep_capture < 3 and number_pe_ncapture_smeared < 2:
            number_pe0_qedep2 += 1
        elif 2 < Qedep_capture < 3 and 2000 < number_pe_ncapture_smeared < 4000:
            number_pe3000_qedep2 += 1
        elif 4 < Qedep_capture < 6 and number_pe_ncapture_smeared < 2:
            number_pe0_qedep5 += 1
        elif 4 < Qedep_capture < 6 and 5000 < number_pe_ncapture_smeared < 8000:
            number_pe6500_qedep5 += 1
        elif 8 < Qedep_capture < 12 and 10000 < number_pe_ncapture_smeared < 20000:
            number_pe15000_qedep11 += 1
        else:
            print("rfile {0}, event = {1:d}".format(rfile, event))

        # append number_pe_ncapture to array:
        array_npe_from_hittime.append(number_pe_ncapture_smeared)
        # append values to list:
        array_npe.append(n_photons)
        array_Qedep.append(Qedep_capture)

print("\ntotal number of events = {0:d}".format(number_of_events_total))
print("\nnumber of events within volume (r < {0:.0f} mm) = {1:d}".format(r_cut, number_of_events))
print("\nnumber of events that pass the time cut = {0:d}".format(number_timecut_pass))
print("number of events that are rejected by time cut = {0:d}".format(number_timecut_rejected))
print("\nnumber of events that pass time and delayed energy cut (min={1:.0f}, max={2:.0f}) = {0:d}"
      .format(number_delayed_energy_pass, min_PE_delayed, max_PE_delayed))
print("number of events that are rejected by delayed energy cut (but pass time cut) = {0:d}"
      .format(number_delayed_energy_rejected))
print("number of events, that falsely pass the delayed energy cut (counted too much) (smeared energy pass, but real "
      "energy is rejected) = {0:d}".format(number_delayed_energy_toomuch))
print("number of events, that are falsely rejected by delayed energy cut (counted too less) (smeared energy rejected, "
      "but real energy pass) = {0:d}".format(number_delayed_energy_tooless))
print("\nnumber of events that pass the neutron multiplicity cut (and also time and energy cut) = {0:d}"
      .format(number_n_mult_pass))
print("number of events that are rejected by neutron multiplicity cut (but pass time and energy cut) = {0:d}"
      .format(number_n_mult_rejected))
print("\nnumber of events that pass distance cut (and also multiplicity, time and energy cut) = {0:d}"
      .format(number_distance_cut_pass))
print("number of events that are rejected by distance cut (but pass multiplicity, time and energy cut) = {0:d}"
      .format(number_distance_cut_rejected))
print("\nnumber of events that pass distance cut with MC truth position (and also multiplicity, time and energy cut) = "
      "{0:d}".format(number_distance_cut_pass_MCtruth))
print("number of events that are rejected by distance cut with MC truth position (but pass multiplicity, time and "
      "energy cut) = {0:d}".format(number_distance_cut_rejected_MCtruth))

print("\nnumber_pe0_qedep2 = {0:d}".format(number_pe0_qedep2))
print("number_pe3000_qedep2 = {0:d}".format(number_pe3000_qedep2))
print("number_pe0_qedep5 = {0:d}".format(number_pe0_qedep5))
print("number_pe6500_qedep5 = {0:d}".format(number_pe6500_qedep5))
print("number_pe15000_qedep11 = {0:d}".format(number_pe15000_qedep11))
print("\nnumber of events with 0 neutron captures = {0:d}".format(number_no_ncapture))

h1 = plt.figure(1, figsize=(15, 8))
plt.plot(array_npe_from_hittime, array_Qedep, "xb", label="entries = {0:.0f}".format(len(array_Qedep)))
plt.vlines(min_PE_delayed, ymin=0.0, ymax=11.0, colors="r", linestyles="dashed",
           label="min. nPE = {0:.0f} PE".format(min_PE_delayed))
plt.vlines(max_PE_delayed, ymin=0.0, ymax=11.0, colors="r", linestyles="dotted",
           label="max. nPE = {0:.0f} PE".format(max_PE_delayed))
plt.ylim(ymin=0.0, ymax=11.0)
plt.xlabel("number of p.e. per event (calculated from hittime distribution)")
plt.ylabel("visible energy (quenched deposited energy) in MeV")
plt.title("Correlation of number of p.e. to energy for captured neutrons in JUNO detector\n"
          "(within time window between {0:.0f} ns and {1:.0f} ms)".format(time_limit, max_time/1000000))
plt.grid()
plt.legend()

h2 = plt.figure(2, figsize=(15, 8))
plt.plot(array_npe_from_hittime, array_Qedep, "xb", label="entries = {0:.0f}".format(number_pe3000_qedep2))
plt.vlines(min_PE_delayed, ymin=0.0, ymax=11.0, colors="r", linestyles="dashed",
           label="min. nPE = {0:.0f} PE".format(min_PE_delayed))
plt.vlines(max_PE_delayed, ymin=0.0, ymax=11.0, colors="r", linestyles="dotted",
           label="max. nPE = {0:.0f} PE".format(max_PE_delayed))
plt.xlim(xmin=2100, xmax=3600)
plt.ylim(ymin=1.5, ymax=3.0)
plt.xlabel("number of p.e. per event (calculated from hittime distribution)")
plt.ylabel("visible energy (quenched deposited energy) in MeV")
plt.title("Correlation of number of p.e. to energy for captured neutrons on H\n"
          "(within time window between {0:.0f} ns and {1:.0f} ms)".format(time_limit, max_time/1000000))
plt.grid()
plt.legend()

plt.show()




