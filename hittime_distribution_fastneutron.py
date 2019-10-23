""" Script to read the root files from neutron simulation (representing fast neutron background) and saved their
    hittime distributions to txt file.

    These hittime distributions can then be analyzed further with pulse_shape_analysis_v1.py as reference to the prompt
    signal of IBD-like NC events (to compare hittime distributions of neutron and NC events).

    As simulated neutrons, the files user_neutron_10_MeV_{}.root, user_neutron_100_MeV_{}.root,
    user_neutron_300_MeV_{}.root and user_neutron_500_MeV_{}.root (from 'conversion_nPE_MeV') are used.

    Files user_neutron_10_MeV_0.root to user_neutron_10_MeV_99.root can not be used, because they deposit too less
    energy in the detector and therefore the hittime distribution has only a few ten entries (Qedep only from 0.2 MeV
    to 5 MeV).
    For files user_neutron_100_MeV_0.root to user_neutron_100_MeV_99.root, Qedep is from 1.2 MeV to 7.4 MeV.
    For files user_neutron_300_MeV_0.root to user_neutron_300_MeV_99.root, Qedep is from 1.5 MeV to 40 MeV.
    For files user_neutron_500_MeV_0.root to user_neutron_500_MeV_99.root, Qedep is from 0.0 MeV to 113 MeV.

    This is the same script like hittime_distribution_positron.py.

    Procedure to get the hittime distribution with vertex reconstruction and time smearing of PMTs:

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

"""
import datetime
import ROOT
import sys
import NC_background_functions
import numpy as np
from matplotlib import pyplot as plt

# get the date and time, when the script was run:
date = datetime.datetime.now()
NOW = date.strftime("%Y-%m-%d %H:%M")

""" define time window and bin width: """
# set time window of whole signal in ns:
min_time = -50
max_time = 1000000
# set time in ns, where the prompt signal should be 0:
time_limit_prompt = 500
# Set bin-width of hittime histogram in ns:
binwidth = 5.0

""" set parameter for volume cut (must be the same like for NC events): """
# cut-radius in mm:
radius_cut = 16000

# path, where root files of positron simulation are saved:
input_path_neutron = "/local/scratch1/pipc51/astro/blum/conversion_nPE_MeV/neutron_output/"
# path, where hittime distributions (png and txt) are saved:
output_path_neutron = "/home/astro/blum/PhD/work/MeVDM_JUNO/fast_neutrons/hittimes/"

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
# preallocate number of events that are analyzed (pass volume cut):
number_analyzed = 0

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

# loop over root files with neutron simulation:
for index in range(first_file_neutron, last_file_neutron + 1, 1):
    # load user_positron_{}.root file:
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

    print("... read {0}...".format(rfile))

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

        # print("\nanalyze event {0:d}".format(evt_id))

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
        if 10.0 < qedep_prmtrkdep < 100.0:
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

        # check if event passes the volume cut:
        if r_reconstructed >= radius_cut:
            # event is rejected by volume cut.
            print("file {0:d}, event = {1:d}: r_init = {2:0.2f} mm".format(index, event, r_reconstructed))
            # go to next event
            continue
        else:
            # event passes volume cut. increment number_analyzed:
            number_analyzed += 1

        """ calculate the real hittime distribution (time of flight correction with reconstructed position and time 
        smearing with TTS for each hit): """
        # get event of 'evt'-tree:
        rtree_evt.GetEntry(event)
        # get evtID of the tree and compare with event:
        evt_id = int(rtree_evt.GetBranch('evtID').GetLeaf('evtID').GetValue())
        if evt_id != event:
            sys.exit("ERROR: evtID of tree ({0:d}) != {1:d}".format(evt_id, event))

        # get number of photons of this event:
        n_photons = int(rtree_evt.GetBranch('nPhotons').GetLeaf('nPhotons').GetValue())

        # preallocate list, where corrected (real) hittimes are saved:
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

            # calculate distance between reconstructed position of event and position of PMT (in mm):
            distance_tof = np.sqrt((x_reconstructed - x_pmt)**2 + (y_reconstructed - y_pmt)**2 +
                                   (z_reconstructed - z_pmt)**2)

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
            if hittime_real < min_time:
                print("------")
                print(hittime_real)
                print(pmtID)
                print(sigma_pmt)

            # append real hittime to array:
            hittime_array.append(hittime_real)

        """ analyze prompt signal: """
        # build histogram, where hittimes are saved:
        # set bin-edges of hittime histogram in ns:
        bins_hittime = np.arange(min_time, max_time + 2 * binwidth, binwidth)
        # build hittime histogram:
        npe_per_hittime, bin_edges_hittime = np.histogram(hittime_array, bins_hittime)

        # get index of bins_hittime corresponding to min_time (should be index = 0):
        index_min_hittime_prompt = 0

        # Where does prompt signal end?
        # get index of bins_hittime corresponding to t_limit:
        index_time_limit_prompt = int((time_limit_prompt + np.abs(min_time)) / binwidth)
        # check if npe_per_hittime (and following two bins) are 0 for this index:
        if (npe_per_hittime[index_time_limit_prompt] == npe_per_hittime[index_time_limit_prompt+1] ==
                npe_per_hittime[index_time_limit_prompt+2] == 0):
            # prompt signal already 0:
            index_max_hittime_prompt = index_time_limit_prompt
        else:
            # prompt signal not yet 0.
            # loop over npe_per_hittime from index_time_limit_prompt until npe_per_hittime (and following two bins)
            # are 0:
            for index2 in range(index_time_limit_prompt, index_time_limit_prompt + 200):
                if npe_per_hittime[index2] == npe_per_hittime[index2+1] == npe_per_hittime[index2+2] == 0:
                    index_max_hittime_prompt = index2
                    break

        # calculate nPE as function of hittime only for prompt time window (from min_hittime_prompt to
        # max_hittime_prompt+1 (last index should be included)):
        npe_per_hittime_prompt = npe_per_hittime[index_min_hittime_prompt:index_max_hittime_prompt + 1]
        # bin edges of hittime histogram only for prompt time window:
        bins_hittime_prompt = bin_edges_hittime[index_min_hittime_prompt:index_max_hittime_prompt + 1]

        # get the minimum and maximum time of the prompt signal time window in ns:
        min_time_prompt = bins_hittime_prompt[0]
        max_time_prompt = bins_hittime_prompt[-1]

        # sum up the values of npe_per_hittime_prompt to get the total number of pe of the prompt signal:
        number_pe_prompt = np.sum(npe_per_hittime_prompt)

        """ save hittime distribution in png file """
        h1 = plt.figure(1)
        plt.step(bins_hittime_prompt, npe_per_hittime_prompt, label="number of pe = {0:d}".format(number_pe_prompt))
        plt.xlabel("hit-time in ns")
        plt.ylabel("number of p.e. per bin (bin-width = {0:0.2f} ns)".format(binwidth))
        plt.title("Hit-time distribution of prompt time window of event {0:d}".format(evt_id))
        plt.xlim(xmin=min_time_prompt, xmax=max_time_prompt)
        plt.legend()
        plt.grid()
        plt.savefig(output_path_neutron + "file{1:d}_evt{0:d}_neutron.png".format(evt_id, index))
        plt.close()
        # plt.show()

        # save npe_per_hittime_prompt to txt file:
        # build list, where 0th entry is start-hittime in ns, 1st entry is last-hittime in ns, 2nd entry is binwidth
        # in ns and the following entries are nPE of each hittime-bin of prompt signal:
        npe_per_hittime_prompt_save = [min_time_prompt, max_time_prompt, binwidth]
        npe_per_hittime_prompt_save.extend(npe_per_hittime_prompt)
        np.savetxt(output_path_neutron + "file{0:d}_evt{1:d}_neutron.txt".format(index, evt_id),
                   npe_per_hittime_prompt_save, fmt='%1.2f',
                   header="Number of pe as function of the corrected hittime (time-of-flight correction and TTS "
                          "smearing) of the prompt neutron signal of file "
                          "user_positron_{0:d}.root,"
                          "\nevent = {1:d}, (analyzed with hittime_distribution_neutron.py, {2}):"
                          "\ntime window of hittime: from {3:.3f} ns to {4:.3f} ns with bin-width = {5:0.3f} ns:"
                   .format(index, evt_id, NOW, min_time_prompt, max_time_prompt, binwidth))

print("total number of events = {0:d}".format(number_evts_total))
print("number of analyzed events = {0:d}".format(number_analyzed))
