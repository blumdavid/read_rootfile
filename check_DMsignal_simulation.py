""" script to check the DM signal simulated with JUNO detsim simulation to check the calculation of the DM signal
    done in script gen_spectrum_v3.py in functions darkmatter_signal_v3() and ibd_kinematics():


"""
import ROOT
import datetime
from NC_background_functions import conversion_npe_to_evis, get_pmt_position, get_20inchpmt_tts, energy_resolution
import numpy as np
from matplotlib import pyplot as plt

# get the date and time, when the script was run:
date = datetime.datetime.now()
now = date.strftime("%Y-%m-%d %H:%M")

# set flag, if pulse shape should be read or not:
flag_read_pulseshape_from_file = True

# set DM mass, which was simulated, in MeV
DM_mass = 60.0
# set minimum of radius in mm:
radius_cut_min = 0.0
# set maximum radius of fiducial volume cut in mm:
radius_cut = 17000.0

# set the path, where user-root files of the simulation are saved:
# path_input = "/local/scratch1/pipc51/astro/blum/DMsignal60MeV/test_CDcenter/"
path_input = "/local/scratch1/pipc51/astro/blum/DMsignal60MeV/"

# set path, where output should be saved:
path_output = "/home/astro/blum/PhD/work/MeVDM_JUNO/gen_spectrum_v2/DM_signal_detsim/"

# set name of output file (txt and png file have same name):
# filename_output = "DMsignal_{0:.0f}MeV_CDcenter".format(DM_mass)
filename_output = "DMsignal_{0:.0f}MeV_R{1:.0f}mTo{2:.0f}m".format(DM_mass, radius_cut_min/1000.0, radius_cut/1000.0)

# set path, where pulse shapes should be saved:
# path_pulse_shape = path_output + "pulse_shape_CDcenter/"
path_pulse_shape = path_output + "pulse_shape_R17m/"

# define array for visible energy in MeV:
min_E = 10.0
max_E = 100.0
interval_E = 0.5
energy = np.arange(min_E, max_E + interval_E, interval_E)

# number of first file to read:
first_file = 0
# number of last file to read:
last_file = 99
# number of events per file:
number_per_file = 100

# define start and end time of the prompt signal in ns:
start_time_prompt = -50.0
end_time_prompt = 500.0

# preallocate array, where number of pe of prompt signal are stored:
array_PE_prompt = []
# preallocate array, where total number of pe is stored:
array_PE_total = []
# preallocate array, where visible energy of prompt signal are stored:
array_Evis = []
# preallocate array, where distance to center of initial position is stored in mm:
array_Rinit = []
# preallocate array, where deposited energy of event in MeV is stored:
array_edep = []
# preallocate array, where quenched deposited energy in MeV is stored:
array_Qedep = []
# preallocate array, where smeared quenched deposited energy in MeV is stored:
array_Qedep_smeared = []
# preallocate the number of analyzed events:
events_analyzed = 0

""" load position of the PMTs and corresponding PMT ID from file PMT_position.root: """
file_PMT_position = "/home/astro/blum/juno/atmoNC/PMT_information/PMT_position.root"
# array with PMT ID and corresponding x, y, z position in mm:
pmtID_pos_file, x_pos_pmt, y_pos_pmt, z_pos_pmt = get_pmt_position(file_PMT_position)

""" load 'time resolution' in ns of the 20 inch PMTs and corresponding PMT ID from file PmtData.root: """
file_PMT_time = "/home/astro/blum/juno/atmoNC/PMT_information/PmtData.root"
# array with PMT ID and corresponding sigma in ns:
pmtID_time_file, sigma_time_20inch = get_20inchpmt_tts(file_PMT_time)
# set TTS (FWHM) of the 3inch PMTs in ns:
tts_3inch = 5.0
# calculate time resolution (sigma) for the 3inch PMTs in ns:
sigma_time_3inch = tts_3inch / (2 * np.sqrt(2 * np.log(2)))
# set effective speed of light in the liquid scintillator in mm/ns (see page 7 of c_effective_JUNO-doc-3144-v2.pdf in
# folder /home/astro/blum/PhD/paper/Pulse_Shape_Discrimination/). Effective refraction index in LS n_eff = 1.54.
# c/n_eff = 299792458 m / 1.54 s ~ 194670427 m/s = 194670427 * 10**(-6) mm/ns ~ 194.67 mm/ns:
c_effective = 194.67

# loop over user-root:
for filenumber in range(first_file, last_file+1, 1):

    print("... analyzed file {0:d}...".format(filenumber))

    # input user-root file:
    input_file = path_input + "user_DMsignal{0:.0f}MeV_{1:d}.root".format(DM_mass, filenumber)
    # load the ROOT file:
    rfile = ROOT.TFile(input_file)
    # get the "evt"-TTree from the TFile:
    rtree_evt = rfile.Get("evt")
    # get 'geninfo' tree from TFile:
    rtree_geninfo = rfile.Get("geninfo")
    # get 'prmtrkdep' tree from TFile:
    rtree_prmtrkdep = rfile.Get("prmtrkdep")

    # loop ever each event:
    for event in range(number_per_file):

        # get the current event in geninfo tree:
        rtree_geninfo.GetEntry(event)

        # get InitX, InitY and InitZ from geninfo tree in mm:
        InitX = float(rtree_geninfo.GetBranch('InitX').GetLeaf('InitX').GetValue(0))
        InitY = float(rtree_geninfo.GetBranch('InitY').GetLeaf('InitY').GetValue(0))
        InitZ = float(rtree_geninfo.GetBranch('InitZ').GetLeaf('InitZ').GetValue(0))

        # calculate distance to detector center in mm
        r_init = np.sqrt(InitX**2 + InitY**2 + InitZ**2)

        # do fiducial volume cut:
        if r_init > radius_cut or r_init < radius_cut_min:
            continue

        # increment events_analyzed:
        events_analyzed += 1

        # get event in prmtrkdep tree:
        rtree_prmtrkdep.GetEntry(event)
        # preallocate Qedep of initial particles:
        Qedep = 0.0
        # get number of initial particles:
        nInitParticles = int(rtree_prmtrkdep.GetBranch("nInitParticles").GetLeaf("nInitParticles").GetValue())
        # loop over initial particles
        for index in range(nInitParticles):
            # get Qedep of initial particle in MeV:
            Qedep_value = float(rtree_prmtrkdep.GetBranch("Qedep").GetLeaf("Qedep").GetValue(index))
            # add it to Qedep:
            Qedep += Qedep_value

        # append Qedep to array:
        array_Qedep.append(Qedep)

        # smear Qedep with energy resolution:
        if Qedep == 0.0:
            array_Qedep_smeared.append(0.0)
        else:
            # get sigma:
            sigma = energy_resolution(Qedep)
            # smear Qedep with sigma:
            array_Qedep_smeared.append(np.random.normal(Qedep, sigma))

        # get event in evt tree:
        rtree_evt.GetEntry(event)

        # get deposited energy of event in MeV:
        edep = float(rtree_evt.GetBranch("edep").GetLeaf("edep").GetValue())
        # append edep to array:
        array_edep.append(edep)

        # get total number of PE of event:
        total_PE = int(rtree_evt.GetBranch("totalPE").GetLeaf("totalPE").GetValue())
        # append total_PE to array:
        array_PE_total.append(total_PE)

        if not flag_read_pulseshape_from_file:

            # preallocate array, where hittimes are stored:
            hittime_array = []

            # get total number of photons:
            nPhotons = int(rtree_evt.GetBranch("nPhotons").GetLeaf("nPhotons").GetValue())

            # loop over nPhotons:
            for index in range(nPhotons):

                # get hittime of photon in ns:
                hittime = float(rtree_evt.GetBranch("hitTime").GetLeaf("hitTime").GetValue(index))
                nPE = int(rtree_evt.GetBranch("nPE").GetLeaf("nPE").GetValue(index))

                # get the pmtID of the hit PMT:
                pmtID = int(rtree_evt.GetBranch('pmtID').GetLeaf('pmtID').GetValue(index))

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
                distance_tof = np.sqrt((InitX - x_pmt)**2 + (InitY - y_pmt)**2 + (InitZ - z_pmt)**2)

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
                    # there are some PMTs with ID around 30000 (user_atmoNC_7.root, event=32: 30637, 30276, 30573,30561,
                    # 30377) -> PMTs with ID above 30000 are Water Pool PMTs!!
                    # go to next photon:
                    continue

                else:
                    # 3inch PMT:
                    sigma_pmt = sigma_time_3inch

                # consider time resolution of PMT by generating normal distributed random number with mu = hittime and
                # sigma = sigma_pmt (only the hittime at the PMT must be smeared, not the time-of-flight):
                hittime_tts = np.random.normal(hittime, sigma_pmt)

                # calculate the 'real' hittime of the photon in ns:
                hittime_real = hittime_tts - time_of_flight

                # check if hittime is within start_time_prompt to end_time_prompt*10 (-50 ns to 5000 ns):
                if start_time_prompt < hittime_real < end_time_prompt*10:
                    # append hittime_real to hittime array:
                    hittime_array.append(hittime_real)
                else:
                    continue

            # only take prompt time window of hittime_array:
            # bin-width of hittime histogram:
            bin_width = 5.0
            # build histogram of hittime_array:
            # set bin-edges of hittime histogram in ns in whole time window:
            bins_hittime = np.arange(start_time_prompt, end_time_prompt*10 + bin_width, bin_width)
            # build hittime histogram:
            npe_per_hittime, bin_edges_hittime = np.histogram(hittime_array, bins_hittime)

            # get index of bins_hittime, where prompt time window ends:
            index_time_cut = int((end_time_prompt + np.abs(start_time_prompt)) / bin_width)

            # check if npe_per_hittime (and the following two bins) are 0 for this index:
            if (npe_per_hittime[index_time_cut] == npe_per_hittime[index_time_cut + 1]
                    == npe_per_hittime[index_time_cut + 2] == 0):
                # prompt signal already 0:
                index_max_hittime_prompt = index_time_cut
            else:
                # prompt signal not yet 0.
                # loop over npe_per_hittime from index_time_cut_min until npe_per_hittime (and the following two bins)
                # are 0:
                for index in range(index_time_cut, index_time_cut + 500):
                    if npe_per_hittime[index] == npe_per_hittime[index + 1] == npe_per_hittime[index + 2] == 0:
                        index_max_hittime_prompt = index
                        break

            # hittime histogram for prompt time window:
            bins_hittime = bins_hittime[0:index_max_hittime_prompt + 1]
            npe_per_hittime = npe_per_hittime[0:index_max_hittime_prompt + 1]

        else:
            # read pulse shape from file:
            file_name = path_pulse_shape + "file{0:d}_evt{1:d}_prompt_signal.txt".format(filenumber, event)
            # load file from txt file:
            npe_from_file = np.loadtxt(file_name)

            # bin-width from file:
            bin_width = npe_from_file[2]
            # set bins of hittime histogram in ns:
            bins_hittime = np.arange(npe_from_file[0], npe_from_file[1]+bin_width, bin_width)
            # get hittime histogram from file:
            npe_per_hittime = npe_from_file[3:]

        # calculate sum of npe_per_hittime to get number of pe in prompt time window:
        nPE_event = np.sum(npe_per_hittime)
        # append nPE_event to array:
        array_PE_prompt.append(nPE_event)

        # convert nPE to MeV:
        Evis = conversion_npe_to_evis(nPE_event)

        # append E_vis to array:
        array_Evis.append(Evis)
        # append r_init to array:
        array_Rinit.append(r_init)

        if not flag_read_pulseshape_from_file:
            # save corrected pulse shape of prompt signal to txt file:
            npe_per_hittime_save_prompt = [bins_hittime[0], bins_hittime[index_max_hittime_prompt], bin_width]
            npe_per_hittime_save_prompt.extend(npe_per_hittime)
            np.savetxt(path_pulse_shape + "file{0:d}_evt{1:d}_prompt_signal.txt".format(filenumber, event),
                       npe_per_hittime_save_prompt, fmt='%1.2f',
                       header="Pulse shape of prompt signal: Number of pe as function of the time "
                              "(time-of-flight correction and TTS smearing) of file user_DMsignal_{0:d}.root,"
                              "\nevent {1:d}, {2}:"
                              "\ntime window of pulse shape: from {3:.3f} ns ns with bin-width = {4:0.3f} "
                              "ns,"
                       .format(filenumber, event, now, start_time_prompt, bin_width))

            if filenumber < 1:
                # save figure of pulse shape to png file:
                h1 = plt.figure(1, figsize=(15, 8))
                plt.step(bins_hittime, npe_per_hittime, label="number of pe = {0:.0f}\nE_vis = {1:.2f} MeV"
                                                              "\nR = {2:.1f} mm".format(nPE_event, Evis, r_init))
                plt.xlabel("time in ns")
                plt.ylabel("nPE per bin (bin-width = {0:.2f} ns)".format(bin_width))
                plt.title("Pulse shape of prompt signal in JUNO detector")
                plt.legend()
                plt.grid()
                plt.savefig(path_pulse_shape + "file{0:d}_evt{1:d}_prompt_signal.png".format(filenumber, event))
                plt.close()

# display DM signal in histogram:
h1 = plt.figure(1, figsize=(15, 8))
plt.hist(array_Evis, bins=energy, histtype="step", align='mid', color="b", linewidth=1.5,
         label="entries = {0:d}".format(events_analyzed))
plt.xlabel("visible energy of positron in MeV")
plt.ylabel("events")
plt.title("DM mass = {0:.0f} MeV".format(DM_mass))
plt.legend()
plt.grid()
plt.savefig(path_output + filename_output + ".png")
plt.close()

# display array_Evis and array_Rinit in 2D histogram:
h2 = plt.figure(2, figsize=(15, 8))
plt.hist2d(array_Evis, array_Rinit, bins=[energy, np.arange(0, 17700+400, 1000)], cmap="Reds")
plt.xlabel("visible energy of positron in MeV")
plt.ylabel("distance to detector center in mm")
plt.title("DM mass = {0:.0f} MeV".format(DM_mass))
plt.colorbar()
plt.grid()
plt.savefig(path_output + filename_output + "_2D.png")
plt.close()

# display edep of event in histogram:
h3 = plt.figure(3, figsize=(15, 8))
plt.hist(array_edep, bins=energy, histtype="step", align='mid', color="b", linewidth=1.5,
         label="entries = {0:d}".format(events_analyzed))
plt.xlabel("total deposited energy in MeV")
plt.ylabel("events")
plt.title("DM mass = {0:.0f} MeV".format(DM_mass))
plt.legend()
plt.grid()
plt.savefig(path_output + filename_output + "_edep.png")
plt.close()

# display total PE of event in histogram:
h4 = plt.figure(4, figsize=(15, 8))
bins_PE = np.arange(min(array_PE_total), max(array_PE_total), 500)
plt.hist(array_PE_total, bins=bins_PE, histtype="step", align='mid', color="r", linewidth=1.5,
         label="total number of p.e.\n(entries = {0:d})".format(events_analyzed))
plt.hist(array_PE_prompt, bins=bins_PE, histtype="step", align='mid', color="b", linewidth=1.5,
         label="number of p.e. of prompt signal\n(entries = {0:d})".format(events_analyzed))
plt.xlabel("number of PE")
plt.ylabel("events")
plt.title("DM mass = {0:.0f} MeV".format(DM_mass))
plt.legend()
plt.grid()
plt.savefig(path_output + filename_output + "_nPE.png")
plt.close()

# display edep of event in histogram:
h5 = plt.figure(5, figsize=(15, 8))
Qedep_hist, bins_hist, patches = plt.hist(array_Qedep, bins=energy, histtype="step", align='mid', color="b",
                                          linewidth=1.5, label="entries = {0:d}".format(events_analyzed))
plt.xlabel("total quenched deposited energy in MeV")
plt.ylabel("events")
plt.title("DM mass = {0:.0f} MeV".format(DM_mass))
plt.legend()
plt.grid()
plt.savefig(path_output + filename_output + "_Qedep.png")
plt.close()

# display edep of event in histogram:
h6 = plt.figure(6, figsize=(15, 8))
Qedep_smeared_hist, bins_hist1, patches1 = plt.hist(array_Qedep_smeared, bins=energy, histtype="step", align='mid',
                                                    color="b", linewidth=1.5, label="entries = {0:d}"
                                                    .format(events_analyzed))
plt.xlabel("total quenched deposited energy smeared with energy resolution in MeV")
plt.ylabel("events")
plt.title("DM mass = {0:.0f} MeV".format(DM_mass))
plt.legend()
plt.grid()
plt.savefig(path_output + filename_output + "_Qedep_smeared.png")
plt.close()

# plt.show()

# save histogram to txt file:
np.savetxt(path_output + filename_output + "_Qedep.txt", Qedep_hist, fmt='%i',
           header="Energy spectrum (number events per bin vs quenched deposited energy) for DM signal of mass = "
                  "{0:.1f} MeV\n"
                  "simulated with JUNO detsim (analyzed with check_DMsignal_simulation.py):\n"
                  "Input root files: user_DMsignal{0:.0f}MeV_{1:d}.root to user_DMsignal{0:.0f}MeV_{2:d}.root,\n"
                  "Number of analyzed events = {3:d}, only events within R={7:.1f}mm are analyzed,\n"
                  "E_visible_min = {4:.1f} MeV, E_visible_max = {5:.1f} MeV, interval_E_vis = {6:.2f} MeV,\n"
                  "Qedep of prmtrkdep tree of user-root files are taken (in MeV):"
           .format(DM_mass, first_file, last_file, events_analyzed, min_E, max_E, interval_E, radius_cut))

# save histogram to txt file:
np.savetxt(path_output + filename_output + "_Qedep_smeared.txt", Qedep_smeared_hist, fmt='%i',
           header="Energy spectrum (number events per bin vs quenched deposited energy) for DM signal of mass = "
                  "{0:.1f} MeV\n"
                  "simulated with JUNO detsim (analyzed with check_DMsignal_simulation.py):\n"
                  "Input root files: user_DMsignal{0:.0f}MeV_{1:d}.root to user_DMsignal{0:.0f}MeV_{2:d}.root,\n"
                  "Number of analyzed events = {3:d}, only events within R={7:.1f}mm are analyzed,\n"
                  "E_visible_min = {4:.1f} MeV, E_visible_max = {5:.1f} MeV, interval_E_vis = {6:.2f} MeV,\n"
                  "Qedep of prmtrkdep tree of user-root files are taken and smeared with energy resolution (in MeV):"
           .format(DM_mass, first_file, last_file, events_analyzed, min_E, max_E, interval_E, radius_cut))
















