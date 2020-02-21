""" Script to check conversion from nPE to MeV of positron, which were simulated with tut_detsim.py of JUNO offline
    version J18v1r1-pre1.

    This is a cross-check to script check_conversion_npe_mev.py, where the conversion from nPE to MeV of neutron and
    protons is calculated.

    The conversion factor of positrons should be equal to the conversion factor of neutrons and protons.

    More information: info_conversion_proton_neutron.odt (/home/astro/blum/juno/atmoNC/data_NC/conversion_nPE_MeV/)
"""
import datetime
import ROOT
import sys
from NC_background_functions import energy_resolution
import numpy as np
from matplotlib import pyplot as plt
from decimal import Decimal
from matplotlib.colors import LogNorm

# get the date and time, when the script was run:
date = datetime.datetime.now()
now = date.strftime("%Y-%m-%d %H:%M")

# set the path of the input files (in this folder the already analyzed pulse shapes of the prompt signal of the IBD
# events are saved):
input_path_PE = "/home/astro/blum/juno/IBD_events/hittimes/"
# path, where the root files of the IBD events are saved to get the visible energy of the event:
input_path_Qedep = "/local/scratch1/pipc51/astro/blum/IBD_hepevt/"

# set path, where results should be saved:
output_path = "/home/astro/blum/juno/atmoNC/data_NC/conversion_nPE_MeV/"

# set the number of the first file and number of the last file that should be read:
start_number = 0
stop_number = 199
# number of entries in the input files:
Number_entries_input = 100
# total number of events:
number_positron = (stop_number - start_number + 1) * Number_entries_input

# set the radius for the volume cut in mm:
r_cut = 17700

# preallocate array, where number of PE of the prompt signal (positron) is stored:
array_number_pe = np.array([])
# preallocate array, where visible energy of the positron in MeV is stored:
array_Qedep = np.array([])

# preallocate arrays for nPE and Qedep for different radii:
array_Npe_0_8m = np.array([])
array_Npe_8_12m = np.array([])
array_Npe_12_14m = np.array([])
array_Npe_14_16m = np.array([])
array_Npe_16_17m = np.array([])
array_Npe_17_17_7m = np.array([])
array_Qedep_0_8m = np.array([])
array_Qedep_8_12m = np.array([])
array_Qedep_12_14m = np.array([])
array_Qedep_14_16m = np.array([])
array_Qedep_16_17m = np.array([])
array_Qedep_17_17_7m = np.array([])

# loop over all files:
for filenumber in range(start_number, stop_number+1, 1):

    # load ROOT file:
    rfile = ROOT.TFile(input_path_Qedep + "user_IBD_hepevt_{0:d}.root".format(filenumber))

    # get the "geninfo"-TTree from the TFile:
    rtree_geninfo = rfile.Get("geninfo")
    # get the number of events in the geninfo Tree:
    number_events_geninfo = rtree_geninfo.GetEntries()

    # get the 'prmtrkdep' tree:
    rtree_prmtrkdep = rfile.Get('prmtrkdep')
    # get number of events in the tree:
    number_events_prmtrkdep = rtree_prmtrkdep.GetEntries()

    # check if number of events are equal in both trees:
    if number_events_geninfo == number_events_prmtrkdep:
        number_events = number_events_geninfo
    else:
        sys.exit("ERROR: number of events in t Trees are NOT equal!!")

    # check if number_events is equal to number_entries_input (if not, the detector simulation was incorrect!!):
    if number_events != Number_entries_input:
        sys.exit("ERROR: number of events are not equal to {0:d} -> Detector Simulation not correct!"
                 .format(Number_entries_input))

    # loop over every event, i.e. every entry, in the TTree:
    for event in range(number_events):

        """ read 'geninfo' tree to check initial energy and to apply volume cut: """
        # get the current event in the TTree:
        rtree_geninfo.GetEntry(event)

        # get event ID of geninfo-tree:
        evt_id_geninfo = int(rtree_geninfo.GetBranch('evtID').GetLeaf('evtID').GetValue())

        # get position of 0th initial particle in x, y, z in mm (positions of the initial particles are equal):
        x_init = float(rtree_geninfo.GetBranch('InitX').GetLeaf('InitX').GetValue(0))
        y_init = float(rtree_geninfo.GetBranch('InitY').GetLeaf('InitY').GetValue(0))
        z_init = float(rtree_geninfo.GetBranch('InitZ').GetLeaf('InitZ').GetValue(0))

        # calculate the radius in mm:
        r_init = np.sqrt(x_init**2 + y_init**2 + z_init**2)

        if r_init > r_cut:
            # apply volume cut:
            continue

        """ read 'prmtrkdep' tree to check quenched deposit energy: """
        # get the current event in the TTree:
        rtree_prmtrkdep.GetEntry(event)

        # get number of particles:
        n_part = int(rtree_prmtrkdep.GetBranch('nInitParticles').GetLeaf('nInitParticles').GetValue())

        # loop over the initial particles to get Qedep of the positron:
        for index in range(n_part):

            # get PDGID of the particle:
            PDGID = int(rtree_prmtrkdep.GetBranch('PDGID').GetLeaf('PDGID').GetValue(index))

            # check if particle is positron:
            if PDGID == -11:

                # get quenched energy of positron in MeV:
                qedep_value = float(rtree_prmtrkdep.GetBranch('Qedep').GetLeaf('Qedep').GetValue(index))

                # consider energy resolution of detector:
                if qedep_value > 0:
                    # get the value of sigma of energy resolution for value of qedep_value:
                    sigma_energy = energy_resolution(qedep_value)
                    # generate normal distributed random number with mean = qedep_value and sigma = sigma_energy:
                    qedep_value = np.random.normal(qedep_value, sigma_energy)

                array_Qedep = np.append(array_Qedep, qedep_value)

                if r_init < 8000:
                    array_Qedep_0_8m = np.append(array_Qedep_0_8m, qedep_value)
                elif 8000 <= r_init < 12000:
                    array_Qedep_8_12m = np.append(array_Qedep_8_12m, qedep_value)
                elif 12000 <= r_init < 14000:
                    array_Qedep_12_14m = np.append(array_Qedep_12_14m, qedep_value)
                elif 14000 <= r_init < 16000:
                    array_Qedep_14_16m = np.append(array_Qedep_14_16m, qedep_value)
                elif 16000 <= r_init < 17000:
                    array_Qedep_16_17m = np.append(array_Qedep_16_17m, qedep_value)
                else:
                    array_Qedep_17_17_7m = np.append(array_Qedep_17_17_7m, qedep_value)

        """ get number of PE of prompt signal from the hittimes txt files: """
        prompt_signal = np.loadtxt(input_path_PE + "file{0:d}_evt{1:d}_prompt_signal.txt".format(filenumber, event))

        # sum up prompt signal from index 6 to end (index 0 to 2 defines reconstructed position and index 3 to 5
        # defines the time window of the prompt signal):
        nPE = np.sum(prompt_signal[6:])

        # append nPE of the prompt signal to array:
        array_number_pe = np.append(array_number_pe, nPE)

        if r_init < 8000:
            array_Npe_0_8m = np.append(array_Npe_0_8m, nPE)
        elif 8000 <= r_init < 12000:
            array_Npe_8_12m = np.append(array_Npe_8_12m, nPE)
        elif 12000 <= r_init < 14000:
            array_Npe_12_14m = np.append(array_Npe_12_14m, nPE)
        elif 14000 <= r_init < 16000:
            array_Npe_14_16m = np.append(array_Npe_14_16m, nPE)
        elif 16000 <= r_init < 17000:
            array_Npe_16_17m = np.append(array_Npe_16_17m, nPE)
        else:
            array_Npe_17_17_7m = np.append(array_Npe_17_17_7m, nPE)

# array_number_pe contains nPE and array_Qedep contains Qedep of all events inside r_cut!

""" do linear fit """
# do linear fit with np.linalg.lstsq:
# The model is y = a * x; x = array_number_pe, y = array_qedep
# x needs to be a column vector instead of a 1D vector for this, however.
array_number_pe_columnvector = array_number_pe[:, np.newaxis]
# first value of output is slope of linear fit (fir_result is array with one entry):
fit_result = np.linalg.lstsq(array_number_pe_columnvector, array_Qedep, rcond=None)[0]
# take first entry of fit_result:
fit_result = fit_result[0]
# set x axis for linear fit:
fit_x_axis = np.arange(0, max(array_number_pe), 100)
# set y axis for linear fit:
fit_y_axis = fit_result * fit_x_axis

""" do linear fit for different radii: """
# x needs to be a column vector instead of a 1D vector for this, however.
array_Npe_c_0_8m = array_Npe_0_8m[:, np.newaxis]
# first value of output is slope of linear fit (fir_result is array with one entry):
fit_result_0_8m = np.linalg.lstsq(array_Npe_c_0_8m, array_Qedep_0_8m, rcond=None)[0]
# take first entry of fit_result:
fit_result_0_8m = fit_result_0_8m[0]
# set y axis for linear fit:
fit_y_axis_0_8m = fit_result_0_8m * fit_x_axis

# x needs to be a column vector instead of a 1D vector for this, however.
array_Npe_c_8_12m = array_Npe_8_12m[:, np.newaxis]
# first value of output is slope of linear fit (fir_result is array with one entry):
fit_result_8_12m = np.linalg.lstsq(array_Npe_c_8_12m, array_Qedep_8_12m, rcond=None)[0]
# take first entry of fit_result:
fit_result_8_12m = fit_result_8_12m[0]
# set y axis for linear fit:
fit_y_axis_8_12m = fit_result_8_12m * fit_x_axis

# x needs to be a column vector instead of a 1D vector for this, however.
array_Npe_c_12_14m = array_Npe_12_14m[:, np.newaxis]
# first value of output is slope of linear fit (fir_result is array with one entry):
fit_result_12_14m = np.linalg.lstsq(array_Npe_c_12_14m, array_Qedep_12_14m, rcond=None)[0]
# take first entry of fit_result:
fit_result_12_14m = fit_result_12_14m[0]
# set y axis for linear fit:
fit_y_axis_12_14m = fit_result_12_14m * fit_x_axis

# x needs to be a column vector instead of a 1D vector for this, however.
array_Npe_c_14_16m = array_Npe_14_16m[:, np.newaxis]
# first value of output is slope of linear fit (fir_result is array with one entry):
fit_result_14_16m = np.linalg.lstsq(array_Npe_c_14_16m, array_Qedep_14_16m, rcond=None)[0]
# take first entry of fit_result:
fit_result_14_16m = fit_result_14_16m[0]
# set y axis for linear fit:
fit_y_axis_14_16m = fit_result_14_16m * fit_x_axis

# x needs to be a column vector instead of a 1D vector for this, however.
array_Npe_c_16_17m = array_Npe_16_17m[:, np.newaxis]
# first value of output is slope of linear fit (fir_result is array with one entry):
fit_result_16_17m = np.linalg.lstsq(array_Npe_c_16_17m, array_Qedep_16_17m, rcond=None)[0]
# take first entry of fit_result:
fit_result_16_17m = fit_result_16_17m[0]
# set y axis for linear fit:
fit_y_axis_16_17m = fit_result_16_17m * fit_x_axis

# x needs to be a column vector instead of a 1D vector for this, however.
array_Npe_c_17_17_7m = array_Npe_17_17_7m[:, np.newaxis]
# first value of output is slope of linear fit (fir_result is array with one entry):
fit_result_17_17_7m = np.linalg.lstsq(array_Npe_c_17_17_7m, array_Qedep_17_17_7m, rcond=None)[0]
# take first entry of fit_result:
fit_result_17_17_7m = fit_result_17_17_7m[0]
# set y axis for linear fit:
fit_y_axis_17_17_7m = fit_result_17_17_7m * fit_x_axis

# """ plot Qedep as function of nPE: """
# h1 = plt.figure(1, figsize=(15, 8))
# plt.plot(array_number_pe, array_Qedep, "rx", label="positron ({0:d} entries)".format(len(array_number_pe)))
# plt.xlabel("number of p.e.")
# plt.ylabel("visible energy in JUNO detector (in MeV)")
# plt.title("Visible energy vs. number of p.e.")
# plt.legend()
# plt.grid()
# plt.savefig(output_path + "qedep_vs_nPE_positron.png")
#
# """ plot Qedep as function of nPE with fit: """
# h2 = plt.figure(2, figsize=(15, 8))
# plt.plot(array_number_pe, array_Qedep, "rx",
#          label="{0:d} entries".format(len(array_number_pe)))
# plt.plot(fit_x_axis, fit_y_axis, "b", label="linear fit: f(x) = {0:.3E} * x"
#          .format(fit_result))
# plt.xlabel("number of p.e.")
# plt.ylabel("visible energy in JUNO detector (in MeV)")
# plt.title("Visible energy vs. number of p.e.\n(with linear fit)")
# plt.legend()
# plt.grid()
# plt.savefig(output_path + "fit_qedep_vs_nPE_positron.png")
#
# """ display Qedep as function of nPE in 2D histogram: """
# h3 = plt.figure(3, figsize=(15, 8))
bins_edges_nPE = np.arange(0, max(array_number_pe), 2000)
bins_edges_Qedep = np.arange(0, max(array_Qedep)+2, 2)
# plt.hist2d(array_number_pe, array_Qedep, [bins_edges_nPE, bins_edges_Qedep], norm=LogNorm(),
#            cmap="rainbow")
# plt.xlabel("number of p.e.")
# plt.ylabel("visible energy in JUNO detector (in MeV)")
# plt.title("Visible energy vs. number of p.e.")
# plt.colorbar()
# plt.legend()
# plt.grid()
# plt.savefig(output_path + "hist2d_Qedep_vs_nPE_positron.png")
#
# """ display Qedep as function of nPE in 2D histogram with fit: """
# h4 = plt.figure(4, figsize=(15, 8))
# plt.hist2d(array_number_pe, array_Qedep, [bins_edges_nPE, bins_edges_Qedep], norm=LogNorm(),
#            cmap="rainbow")
# plt.plot(fit_x_axis, fit_y_axis, "k", label="{1:d} entries\nlinear fit: f(x) = {0:.3E} * x"
#          .format(fit_result, len(array_number_pe)))
# plt.xlabel("number of p.e.")
# plt.ylabel("visible energy in JUNO detector (in MeV)")
# plt.title("Visible energy vs. number of p.e.\nwith linear fit")
# plt.colorbar()
# plt.legend()
# plt.grid()
# plt.savefig(output_path + "hist2d_Qedep_vs_nPE_positron_fit.png")

""" display Qedep as function of nPE in 2D histogram with fit: """
h5 = plt.figure(5, figsize=(15, 8))
plt.hist2d(array_Npe_0_8m, array_Qedep_0_8m, [bins_edges_nPE, bins_edges_Qedep], norm=LogNorm(),
           cmap="rainbow")
plt.plot(fit_x_axis, fit_y_axis_0_8m, "k", label="{1:d} entries\nlinear fit: f(x) = {0:.3E} * x"
         .format(fit_result_0_8m, len(array_Npe_0_8m)))
plt.xlabel("number of p.e.")
plt.ylabel("visible energy in JUNO detector (in MeV)")
plt.title("Visible energy vs. number of p.e.\nwith linear fit for R < 8 m")
plt.colorbar()
plt.legend()
plt.grid()
plt.savefig(output_path + "hist2d_Qedep_vs_nPE_positron_fit_0_8m.png")

""" display Qedep as function of nPE in 2D histogram with fit: """
h6 = plt.figure(6, figsize=(15, 8))
plt.hist2d(array_Npe_8_12m, array_Qedep_8_12m, [bins_edges_nPE, bins_edges_Qedep], norm=LogNorm(),
           cmap="rainbow")
plt.plot(fit_x_axis, fit_y_axis_8_12m, "k", label="{1:d} entries\nlinear fit: f(x) = {0:.3E} * x"
         .format(fit_result_8_12m, len(array_Npe_8_12m)))
plt.xlabel("number of p.e.")
plt.ylabel("visible energy in JUNO detector (in MeV)")
plt.title("Visible energy vs. number of p.e.\nwith linear fit for 8 m < R < 12 m")
plt.colorbar()
plt.legend()
plt.grid()
plt.savefig(output_path + "hist2d_Qedep_vs_nPE_positron_fit_8_12m.png")

""" display Qedep as function of nPE in 2D histogram with fit: """
h7 = plt.figure(7, figsize=(15, 8))
plt.hist2d(array_Npe_12_14m, array_Qedep_12_14m, [bins_edges_nPE, bins_edges_Qedep], norm=LogNorm(),
           cmap="rainbow")
plt.plot(fit_x_axis, fit_y_axis_12_14m, "k", label="{1:d} entries\nlinear fit: f(x) = {0:.3E} * x"
         .format(fit_result_12_14m, len(array_Npe_12_14m)))
plt.xlabel("number of p.e.")
plt.ylabel("visible energy in JUNO detector (in MeV)")
plt.title("Visible energy vs. number of p.e.\nwith linear fit for 12 m < R < 14 m")
plt.colorbar()
plt.legend()
plt.grid()
plt.savefig(output_path + "hist2d_Qedep_vs_nPE_positron_fit_12_14m.png")

""" display Qedep as function of nPE in 2D histogram with fit: """
h8 = plt.figure(8, figsize=(15, 8))
plt.hist2d(array_Npe_14_16m, array_Qedep_14_16m, [bins_edges_nPE, bins_edges_Qedep], norm=LogNorm(),
           cmap="rainbow")
plt.plot(fit_x_axis, fit_y_axis_14_16m, "k", label="{1:d} entries\nlinear fit: f(x) = {0:.3E} * x"
         .format(fit_result_14_16m, len(array_Npe_14_16m)))
plt.xlabel("number of p.e.")
plt.ylabel("visible energy in JUNO detector (in MeV)")
plt.title("Visible energy vs. number of p.e.\nwith linear fit for 14 m < R < 16 m")
plt.colorbar()
plt.legend()
plt.grid()
plt.savefig(output_path + "hist2d_Qedep_vs_nPE_positron_fit_14_16m.png")

""" display Qedep as function of nPE in 2D histogram with fit: """
h9 = plt.figure(9, figsize=(15, 8))
plt.hist2d(array_Npe_16_17m, array_Qedep_16_17m, [bins_edges_nPE, bins_edges_Qedep], norm=LogNorm(),
           cmap="rainbow")
plt.plot(fit_x_axis, fit_y_axis_16_17m, "k", label="{1:d} entries\nlinear fit: f(x) = {0:.3E} * x"
         .format(fit_result_16_17m, len(array_Npe_16_17m)))
plt.xlabel("number of p.e.")
plt.ylabel("visible energy in JUNO detector (in MeV)")
plt.title("Visible energy vs. number of p.e.\nwith linear fit for 16 m < R < 17 m")
plt.colorbar()
plt.legend()
plt.grid()
plt.savefig(output_path + "hist2d_Qedep_vs_nPE_positron_fit_16_17m.png")

""" display Qedep as function of nPE in 2D histogram with fit: """
h10 = plt.figure(10, figsize=(15, 8))
plt.hist2d(array_Npe_17_17_7m, array_Qedep_17_17_7m, [bins_edges_nPE, bins_edges_Qedep], norm=LogNorm(),
           cmap="rainbow")
plt.plot(fit_x_axis, fit_y_axis_17_17_7m, "k", label="{1:d} entries\nlinear fit: f(x) = {0:.3E} * x"
         .format(fit_result_17_17_7m, len(array_Npe_17_17_7m)))
plt.xlabel("number of p.e.")
plt.ylabel("visible energy in JUNO detector (in MeV)")
plt.title("Visible energy vs. number of p.e.\nwith linear fit for 17 m < R < 17.7 m")
plt.colorbar()
plt.legend()
plt.grid()
plt.savefig(output_path + "hist2d_Qedep_vs_nPE_positron_fit_17_17_7m.png")

plt.close()








