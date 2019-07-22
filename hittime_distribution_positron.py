""" Script to read the root files from positron simulation and saved their hittime distributions to txt file.

    These hittime distributions can then be analyzed further with pulse_shape_analysis.py as reference to the prompt
    signal of IBD-like NC events (to compare hittime distributions of positrons and NC events).

    Procedure to get the hittime distribution with vertex reconstruction and time smearing of PMTs:

    1.  apply same cuts like on prompt signals of NC events:
        1.1 energy cut on prompt signal: only positrons with energy from 10 MeV to 100 MeV (uniformly distributed)
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


def get_hittime_from_rootfile_fixenergy(input_path, output_path, first_file, last_file, num_evts, kin_energy, min_t,
                                        max_t, t_limit, bin_width, radius, now):
    """
    function to read events of the root files of positron simulation and save hittime distribution to png and txt file.

    function is used for positron with fixed kinetic energy.

    IMPORTANT: time of flight correction is NOT correct!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    :param input_path: path, where root-files are saved
    :param output_path: path, where png and txt files of hittime distribution are saved
    :param first_file: number of the first file to read
    :param last_file: number of the last file to read
    :param num_evts: number of events per file
    :param kin_energy: kinetic energy of positron in MeV
    :param min_t: minimum of time window of whole signal in ns
    :param max_t: maximum of time window of whole signal in ns
    :param t_limit: time in ns, where prompt signal should be 0
    :param bin_width: bin-width of hittime-distribution in ns
    :param radius: cut radius in mm
    :param now: actual time
    :return:
    """
    # preallocate array, where total nPE of prompt signal per event of all files is stored:
    num_pe_total = np.array([])
    # number of events that are analyzed (pass volume cut):
    number_analyzed = 0

    # loop over root files with positron simulation:
    for index in range(first_file, last_file + 1, 1):
        # load user_positron_{}.root file:
        rfile = ROOT.TFile(input_path + "user_positron_{0:d}_MeV_{1:d}.root".format(kin_energy, index))
        print("... read {0}...".format(rfile))

        # get the "evt"-TTree from the TFile:
        rtree_evt = rfile.Get("evt")
        # get geninfo tree from TFile:
        rtree_geninfo = rfile.Get("geninfo")
        # get the number of events in the 'evt' Tree:
        num_events_evt = rtree_evt.GetEntries()
        # check number of events:
        if num_events_evt != num_evts:
            sys.exit("ERROR: number of events in root file ({0:d}) != {1:d}"
                     .format(num_events_evt, num_evts))

        # loop over the events:
        for event in range(num_evts):
            # check volume cut:
            rtree_geninfo.GetEntry(event)
            # get number of initial particles:
            n_init_particles = int(rtree_geninfo.GetBranch('nInitParticles').GetLeaf('nInitParticles').GetValue())
            if n_init_particles != 1:
                sys.exit("ERROR: more than 1 initial particles in event {0:d}".format(event))
            # get initial x, y, z position:
            x_init = float(rtree_geninfo.GetBranch('InitX').GetLeaf('InitX').GetValue())
            y_init = float(rtree_geninfo.GetBranch('InitY').GetLeaf('InitY').GetValue())
            z_init = float(rtree_geninfo.GetBranch('InitZ').GetLeaf('InitZ').GetValue())
            # calculate distance to center:
            r_init = np.sqrt(x_init**2 + y_init**2 + z_init**2)

            if r_init >= radius:
                print("file {0:d}, event = {1:d}: r_init = {2:0.2f} mm".format(index, event, r_init))
                continue

            # get event of 'evt'-tree:
            rtree_evt.GetEntry(event)
            # get evtID of the tree and compare with event:
            evt_id = int(rtree_evt.GetBranch('evtID').GetLeaf('evtID').GetValue())
            if evt_id != event:
                sys.exit("ERROR: evtID of tree ({0:d}) != {1:d}".format(evt_id, event))

            print("\nanalyze event {0:d}".format(evt_id))

            # increment number_analyzed:
            number_analyzed += 1

            # get number of photons of this event:
            n_photons = int(rtree_evt.GetBranch('nPhotons').GetLeaf('nPhotons').GetValue())

            # preallocate empty array to build default hittime-histogram:
            hittime_array = []

            # loop over every photon in the event:
            for index1 in range(n_photons):

                # get PMT ID, where photon is absorbed:
                pmt_id = int(rtree_evt.GetBranch('pmtID').GetLeaf('pmtID').GetValue(index1))

                # only 20 inch PMTs (PMT ID of 20 inch PMTs are below 21000, PMT ID of 3 inch PMTs start at 290000):
                if pmt_id < 25000:
                    # get nPE for this photon:
                    n_pe = int(rtree_evt.GetBranch('nPE').GetLeaf('nPE').GetValue(index1))
                    # check, if photon produces only 1 PE:
                    if n_pe != 1:
                        print("{1:d} PE for 1 photon in event {0:d} in file user_positron_{3:d}_MeV_{2:d}.root"
                              .format(evt_id, n_pe, index, kin_energy))

                    # get hittime of this photon:
                    hit_time = float(rtree_evt.GetBranch('hitTime').GetLeaf('hitTime').GetValue(index1))

                    # append hittime to array:
                    hittime_array.append(hit_time)

                else:
                    continue

            """ analyze prompt signal: """
            # build histogram, where hittimes are saved:
            # set bin-edges of hittime histogram in ns:
            bins_hittime = np.arange(min_t, max_t + 2 * bin_width, bin_width)
            # build hittime histogram:
            npe_per_hittime, bin_edges_hittime = np.histogram(hittime_array, bins_hittime)

            # get index of bins_hittime corresponding to min_time (should be index = 0):
            index_min_hittime_prompt = int(min_t / bin_width)

            # Where does prompt signal end?
            # get index of bins_hittime corresponding to t_limit:
            index_time_limit_prompt = int(t_limit / bin_width)
            # check if npe_per_hittime is 0 for this index:
            if npe_per_hittime[index_time_limit_prompt] == 0:
                # prompt signal already 0:
                index_max_hittime_prompt = index_time_limit_prompt
            else:
                # prompt signal not yet 0.
                # loop over npe_per_hittime from index_time_limit_prompt until npe_per_hittime is 0:
                for index2 in range(index_time_limit_prompt, index_time_limit_prompt+200):
                    if npe_per_hittime[index2] == 0:
                        index_max_hittime_prompt = index2
                        break

            # calculate nPE as function of hittime only for prompt time window (from min_hittime_prompt to
            # max_hittime_prompt+1 (last index should be included)):
            npe_per_hittime_prompt = npe_per_hittime[index_min_hittime_prompt:index_max_hittime_prompt+1]
            # bin edges of hittime histogram only for prompt time window:
            bins_hittime_prompt = bin_edges_hittime[index_min_hittime_prompt:index_max_hittime_prompt+1]

            # get the minimum and maximum time of the prompt signal time window in ns:
            min_time_prompt = bins_hittime_prompt[0]
            max_time_prompt = bins_hittime_prompt[-1]

            # sum up the values of npe_per_hittime_prompt to get the total number of pe of the prompt signal:
            number_pe_prompt = np.sum(npe_per_hittime_prompt)

            # append number of pe to array:
            num_pe_total = np.append(num_pe_total, number_pe_prompt)

            h1 = plt.figure(1)
            plt.step(bins_hittime_prompt, npe_per_hittime_prompt, label="number of pe = {0:d}".format(number_pe_prompt))
            plt.xlabel("hit-time in ns")
            plt.ylabel("number of p.e. per bin (bin-width = {0:0.2f} ns)".format(bin_width))
            plt.title("Hit-time distribution of prompt time window of event {0:d}".format(evt_id))
            plt.xlim(xmin=min_time_prompt, xmax=max_time_prompt)
            plt.legend()
            plt.grid()
            plt.savefig(output_path + "file{1:d}_evt{0:d}_positron_{2:d}_MeV.png".format(evt_id, index, kin_energy))
            plt.close()
            # plt.show()

            # save npe_per_hittime_prompt to txt file:
            # build list, where 0th entry is start-hittime in ns, 1st entry is last-hittime in ns, 2nd entry is binwidth
            # in ns and the following entries are nPE of each hittime-bin of prompt signal:
            npe_per_hittime_prompt_save = [min_time_prompt, max_time_prompt, bin_width]
            npe_per_hittime_prompt_save.extend(npe_per_hittime_prompt)
            np.savetxt(output_path + "file{0:d}_evt{1:d}_positron_{2:d}_MeV.txt".format(index, evt_id, kin_energy),
                       npe_per_hittime_prompt_save, fmt='%1.2f',
                       header="Number of pe as function of the hittime of the prompt positron signal of file "
                              "user_positron_{6:d}_MeV_{0:d}.root,"
                              "\nevent = {1:d}, (analyzed with hittime_distribution_positron.py, {2}):"
                              "\ntime window of hittime: from {3:.3f} ns to {4:.3f} ns with bin-width = {5:0.3f} ns:"
                       .format(index, evt_id, now, min_time_prompt, max_time_prompt, bin_width, kin_energy))

    return num_pe_total, number_analyzed


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
input_path_positron = "/local/scratch1/pipc51/astro/blum/positron_output/"
# path, where hittime distributions (png and txt) are saved:
output_path_positron = "/home/astro/blum/juno/atmoNC/data_NC/output_PSD/positron_hittime/"

""" analyze positron hittime distribution for kinetic energy from 10 MeV to 100 MeV (uniformly distributed): """
# first file of positron simulation:
first_file_positron = 0
# last file of positron simulation:
last_file_positron = 99
# number of events per file:
number_evts_per_file = 100
# total number of positron events:
number_evts_total = (last_file_positron - first_file_positron + 1) * number_evts_per_file
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

# loop over root files with positron simulation:
for index in range(first_file_positron, last_file_positron + 1, 1):
    # load user_positron_{}.root file:
    rfile = ROOT.TFile(input_path_positron + "user_positron_{0:d}.root".format(index))
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
        plt.savefig(output_path_positron + "file{1:d}_evt{0:d}_positron.png".format(evt_id, index))
        plt.close()
        # plt.show()

        # save npe_per_hittime_prompt to txt file:
        # build list, where 0th entry is start-hittime in ns, 1st entry is last-hittime in ns, 2nd entry is binwidth
        # in ns and the following entries are nPE of each hittime-bin of prompt signal:
        npe_per_hittime_prompt_save = [min_time_prompt, max_time_prompt, binwidth]
        npe_per_hittime_prompt_save.extend(npe_per_hittime_prompt)
        np.savetxt(output_path_positron + "file{0:d}_evt{1:d}_positron.txt".format(index, evt_id),
                   npe_per_hittime_prompt_save, fmt='%1.2f',
                   header="Number of pe as function of the corrected hittime (time-of-flight correction and TTS "
                          "smearing) of the prompt positron signal of file "
                          "user_positron_{0:d}.root,"
                          "\nevent = {1:d}, (analyzed with hittime_distribution_positron.py, {2}):"
                          "\ntime window of hittime: from {3:.3f} ns to {4:.3f} ns with bin-width = {5:0.3f} ns:"
                   .format(index, evt_id, NOW, min_time_prompt, max_time_prompt, binwidth))

print("total number of events = {0:d}".format(number_evts_total))
print("number of analyzed events = {0:d}".format(number_analyzed))

""" analyze 10 MeV positron hittime distributions """
# # event per root file:
# number_evts_positron = 10
#
# # first file of positron simulation for 10 MeV:
# first_file_10 = 0
# # last file of positron simulation for 10 MeV:
# last_file_10 = 99
# # kinetic energy of positrons in MeV:
# energy_positron_10 = 10
# # total number of positron events with 10 MeV:
# number_evts_total_10 = (last_file_10 - first_file_10 + 1) * number_evts_positron
#
# # get hittime distribution (return values: number of pe of each event, number of events that are analyzed):
# number_pe_10, number_analyzed_10 = get_hittime_from_rootfile_fixenergy(input_path_positron, output_path_positron,
#                                                                        first_file_10, last_file_10,
#                                                                        number_evts_positron, energy_positron_10,
#                                                                        min_time, max_time, time_limit_prompt,
#                                                                        binwidth, radius_cut, NOW)
#
# # check the distribution of total nPE for all events with positrons of 10 MeV kinetic energy:
# # calculate mean of number of PE:
# mean_10_MeV = np.mean(number_pe_10)
# # display number_pe_10 in histogram:
# h1 = plt.figure(1, figsize=(15, 8))
# n_PE_10MeV, bins_10MeV, patches1 = plt.hist(number_pe_10, align='mid', bins=100,
#                                             label="{0:d} positrons with kinetic energy = {1:d} MeV"
#                                             .format(number_evts_total_10, energy_positron_10))
# plt.vlines(mean_10_MeV, ymin=0, ymax=max(n_PE_10MeV), label="mean = {0:0.2f} nPE".format(mean_10_MeV))
# plt.xlabel("number of PE (per positron)", fontsize=13)
# plt.ylabel("entries per bin", fontsize=13)
# plt.title("Number of PE of 10 MeV positrons", fontsize=18)
# plt.legend()
# plt.grid()
# plt.savefig(output_path_positron + "hist_nPE_10_MeV.png")
# plt.close()

""" analyze 100 MeV positron hittime distributions """
# # event per root file:
# number_evts_positron = 10
#
# # first file of positron simulation for 100 MeV:
# first_file_100 = 0
# # last file of positron simulation for 100 MeV:
# last_file_100 = 99
# # kinetic energy of positrons in MeV:
# energy_positron_100 = 100
# # total number of positron events with 10 MeV:
# number_evts_total_100 = (last_file_100 - first_file_100 + 1) * number_evts_positron
#
# number_pe_100, number_analyzed_100 = get_hittime_from_rootfile_fixenergy(input_path_positron, output_path_positron,
#                                                                          first_file_100, last_file_100,
#                                                                          number_evts_positron, energy_positron_100,
#                                                                          min_time, max_time, time_limit_prompt,
#                                                                          binwidth, radius_cut, NOW)
#
# # check the distribution of total nPE for all events with positrons of 100 MeV kinetic energy:
# # calculate mean of number of PE:
# mean_100_MeV = np.mean(number_pe_100)
# # display number_pe_100 in histogram:
# h2 = plt.figure(2, figsize=(15, 8))
# n_PE_100MeV, bins_100MeV, patches2 = plt.hist(number_pe_100, align='mid', bins=100,
#                                               label="{0:d} positrons with kinetic energy = {1:d} MeV"
#                                               .format(number_evts_total_100, energy_positron_100))
# plt.vlines(mean_100_MeV, ymin=0, ymax=max(n_PE_100MeV), label="mean = {0:0.2f} nPE".format(mean_100_MeV))
# plt.xlabel("number of PE (per positron)", fontsize=13)
# plt.ylabel("entries per bin", fontsize=13)
# plt.title("Number of PE of 100 MeV positrons", fontsize=18)
# plt.legend()
# plt.grid()
# plt.savefig(output_path_positron + "hist_nPE_100_MeV.png")
# plt.close()

