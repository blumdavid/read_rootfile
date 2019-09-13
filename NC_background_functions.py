""" Script contains different functions which can be used to check the root-output file from DSNB-NC.exe:

    DSNB-NC.exe is the generator of the Neutral Current background from atmospheric neutrinos and antineutrinos
    (build by Jie Cheng).

    The interactions of the neutrino and antineutrinos with the liquid scintillator are simulated with GENIE (or NuWro).
    For the flux of the atmospheric neutrinos (nu_e, nu_e_bar, nu_mu, nu_mu_bar) the flux of Honda for JUNO site is used
    (same flux like in my simulation of atmospheric CC background, energy range from 100 MeV to 10 GeV).

    -> is the solar average used? or solar minimum? or solar maximum?

    For the interactions with the LS, the interactions on C12 (-> neutral current interactions) and on free protons
    (-> elastic scattering) are considered.

    The output of the GENIE simulation is saved in "genie_data.root".

    Then the deexcitation of the products from the NC interactions is simulated:
    The status of the residual nucleus is unknown. Therefor the deexcitation of the residual nucleus is assumed from
    a simple shell model
    -> 1/3 of the residual nuclei are excited
    -> 2/3 of the residual nuclei are NOT excited, but in the ground state
    (Jie is working on the implementation of a more complicated nuclear model, but it is not implemented in the DSNB-NC
    generator yet (10.10.2018))

    The deexitation of the residual nucleus is simulated with the TALYS software and saved in "*deexcitation.root" for
    all nuclei (Li7, Li8, Li9, C9, C10, C11, Be7, Be9, Be10, B8, B9, B10, B11)

    With these root-files, the deexcitation is calculated and all final particles with their PDG ID, momentum and
    mass printed in the terminal.
    These information together with the target ID, the channel ID, the deexcitation ID, the isotope ID and the energy
    of the incoming neutrino are saved in the root output file.

"""

import ROOT
from array import array
import sys
import numpy as np
from matplotlib import pyplot as plt
import xml.etree.ElementTree as ET


def pulse_shape(hittime, npe, tail_start, tail_end):
    """
    function to analyzed the hittime distribution of one event with the tail to total method (charge integration method)

    1. find the start of the distribution (so there are no values at the beginning with 0 nPE)
    2. normalize the hittime distribution to 1
    3. calculate the 'charge' of the whole distribution
    4. calculate the 'charge' of the tail of the distribution (define by tail_start and tail_end)
    5. calculate the ration between charge of tail and charge of total distribution


    :param hittime: hittime array, which defines the time window (array in ns)
    :param npe: number of pe per bin (corresponds to hittime array) (array)
    :param tail_start: defines the start value of the tail in ns
    :param tail_end: defines the stop value of the tail in ns
    :return:
    """
    """ find start of distribution: """
    # get maximum value of npe:
    # maximum_npe = np.max(npe)
    # calculate 10 % of maximum_pe:
    # start_condition = 0.1 * maximum_npe
    start_condition = 0

    # loop over npe until one entry is higher than start_condition:
    for index5 in range(len(npe)):
        if npe[index5] > start_condition:
            break

    # do 'time of flight' correction for hittime and npe:
    hittime = hittime[index5:]
    npe = npe[index5:]

    """ normalize hittime distribution to 1: """
    # calculate the integral of the whole hittime distribution:
    integral_npe = np.trapz(npe, hittime)
    # normalize npe distribution to 1:
    npe_norm = npe / integral_npe

    """ integral (charge) of total distribution: """
    # should be 1 because of normalization:
    integral_total = np.trapz(npe_norm, hittime)

    """ integral (charge) of the tail of the distribution: """
    # get the index of hittime, which correspond to tail_start:
    for index6 in range(len(hittime)):
        if hittime[index6] == tail_start:
            index_tail_start = index6
        elif hittime[index6] == tail_end:
            index_tail_end = index6
        else:
            continue

    # define the time window of the tail of the distribution:
    hittime_tail = hittime[index_tail_start:index_tail_end+1]
    # get the corresponding npe_norm array:
    npe_tail = npe_norm[index_tail_start:index_tail_end+1]

    # integrate the tail of the distribution:
    integral_tail = np.trapz(npe_tail, hittime_tail)

    """ tail to total ratio: """
    # calculate the ratio between integral od tail of distribution and of total distribution:
    tot_ratio = integral_tail / integral_total

    return tot_ratio, npe_norm


def tot_efficiency(tot_values_positron, tot_values_nc, eff_nc):
    """
    function to get the efficiency, how many positron events are cut away, depending on eff_nc
    :param tot_values_positron: list/array of all tail-to-total values of positron events
    :param tot_values_nc: list/array of all tail-to-total values of NC events
    :param eff_nc: efficiency, how many NC events should be cut away (in percent). For example: eff_nc=99% -> 99% of all
    NC events are cut away
    :return:
    """
    # calculate the percentile of tot_values_nc defined by (100-eff_nc) (value of tot, where (100-eff_nc) % of all
    # tot-values are smaller)
    tot_cut_value = np.percentile(tot_values_nc, 100 - eff_nc)

    # use tot_cut_value to calculate, how many positron events are cut away:
    num_positron_cut = 0
    # loop over tot_values_positron:
    for index in range(len(tot_values_positron)):
        if tot_values_positron[index] >= tot_cut_value:
            # positron events is cut away by PSD:
            num_positron_cut += 1
        else:
            continue

    # get total number of positron events:
    num_pos_total = len(tot_values_positron)

    # calculate efficiency, how many positron events are cut away by PSD based on eff_nc (in percent):
    eff_positron = float(num_positron_cut) / float(num_pos_total) * 100

    return eff_positron, tot_cut_value



def energy_resolution(e_vis):
    """ 'energy resolution' of the JUNO detector for energies given in MeV (same function like energy_resolution() in
        /home/astro/blum/PhD/work/MeVDM_JUNO/source/gen_spectrum_functions.py):

        INFO-me: in detail: energy_resolution returns the width sigma of the gaussian distribution.
        The real energy of the neutrino is smeared by a gaussian distribution characterized by sigma.

        :param e_vis: visible energy in MeV (float or np.array of float)

        :return sigma/width of the gaussian distribution in MeV (float or np.array of float)
        """
    # parameters to describe the energy resolution in percent (maximum values of table 13-4, page 196, PhysicsReport):
    # TODO-me: p0, p1 and p2 are determined by the maximum values of the parameters from table 13-4
    # p0: is the leading term dominated by the photon statistics (in percent):
    p0 = 2.8
    # p1 and p2 come from detector effects such as PMT dark noise, variation of the PMT QE and the
    # reconstructed vertex smearing (in percent):
    p1 = 0.26
    p2 = 0.9
    # energy resolution defined as sigma/E_visible in percent, 3-parameter function (page 195, PhysicsReport) (float):
    energy_res = np.sqrt((p0 / np.sqrt(e_vis))**2 + p1**2 + (p2 / e_vis)**2)
    # sigma or width of the gaussian distribution in MeV * percent (float):
    sigma_resolution = energy_res * e_vis
    # sigma in MeV (float):
    sigma_resolution = sigma_resolution / 100

    return sigma_resolution


def energy_resolution_pe(npe):
    """ 'energy resolution' of the JUNO detector for energies given in nPE (same function like energy_resolution() in
        /home/astro/blum/PhD/work/MeVDM_JUNO/source/gen_spectrum_functions.py):

        INFO-me: in detail: energy_resolution returns the width sigma of the gaussian distribution.
        The real energy of the neutrino is smeared by a gaussian distribution characterized by sigma.

        :param npe: number of photo-electrons (float or np.array of float)

        :return sigma/width of the gaussian distribution in nPE (float or np.array of float)
        """
    # parameters to describe the energy resolution in percent (maximum values of table 13-4, page 196, PhysicsReport):
    # TODO-me: p0, p1 and p2 are determined by the maximum values of the parameters from table 13-4
    # p0: is the leading term dominated by the photon statistics (in percent):
    p0 = 2.8
    # p1 and p2 come from detector effects such as PMT dark noise, variation of the PMT QE and the
    # reconstructed vertex smearing (in percent):
    p1 = 0.26
    p2 = 0.9

    # conversion factor a (E_vis = a * nPE) (units: MeV/PE):
    # INFO-me: must be same value as in function conversion_npe_to_evis()
    a = 0.0007483

    # convert p0, p1, p2 to parameters corresponding to nPE:
    p0_new = p0 / np.sqrt(a)
    p1_new = p1
    p2_new = p2 / a

    # sigma (width of the gaussian distribution) defined as sigma = resolution * nPE (page 195, PhysicsReport) in
    # units of nPE * percent (float):
    sigma_resolution = np.sqrt((p0_new / np.sqrt(npe)) ** 2 + p1_new ** 2 + (p2_new / npe) ** 2) * npe
    # sigma in nPE (float):
    sigma_resolution = sigma_resolution / 100

    return sigma_resolution


def analyze_delayed_signal_v2(npe_per_time, bins_time, first_index, threshold, threshold2, min_pe, evt):
    """
    function to analyze the time window, where a delayed signal could be. This time window is check for a possible
    delayed signal (signal has to be greater than 'threshold' and nPE of signal peak must be above min_pe)

    :param npe_per_time: number of pe as function of hittime for time window of delayed signal (array of integer)
    :param bins_time: array of bins, which contains the information about the hittime in ns (array of float)
    :param first_index: first index of npe_per_time, that should be analyzed (integer)
    :param threshold: threshold in nPE per bin, that a delayed signal must have (integer)
    :param threshold2: threshold in nPE per bin, that specifies the minimal nPE corresponding to a peak
    :param min_pe: lower threshold of PE for delayed signal
    :param evt: event ID of the event

    :return:
    """
    # preallocate number of delayed signal and index_after_peak1 and number of pe in the delayed signal:
    num_del_signal = 0
    index_after_peak1 = len(npe_per_time)
    num_pe_del = 0
    begin_peak = 0
    end_peak = 2000000

    # loop over values of the histogram bins of delayed window:
    for index4 in range(first_index, len(npe_per_time)):
        # check if number of PE in this bin is above the threshold:
        if npe_per_time[index4] > threshold:
            # possible delayed signal (signal in delayed window):

            # check if index4 = first_value -> first value of npe_per_time above threshold -> break from for loop and
            # return index_after_peak1 = index4 + 1:
            if index4 == first_index:
                index_after_peak1 = index4 + 1
                print("WARNING: npe_per_time[{1:d}] > threshold in evtID = {0:d}".format(evt, index4))
                break

            # calculate number of PEs in this signal peak:
            # preallocate sum of PE per bin in the signal peak:
            sum_pe_peak = 0

            # add nPE of npe_per_time[index4] to sum_pe_peak:
            sum_pe_peak = sum_pe_peak + npe_per_time[index4]

            # check previous bins:
            for num1 in range(1, 100, 1):
                # check if value in previous bins_time[index4 - num1] is above threshold2:
                if npe_per_time[index4 - num1] > threshold2:
                    # add nPE of this bin to sum_pe_peak:
                    sum_pe_peak = sum_pe_peak + npe_per_time[index4 - num1]

                else:
                    # hittime, when signal peak begins, in ns:
                    begin_peak = bins_time[index4 - num1]
                    break

            # check following bins:
            for num2 in range(1, 1000, 1):
                # check if (index4 + num2 + 2) is in the range of npe_per_time array:
                if (index4 + num2 + 2) >= len(npe_per_time):
                    print("Warning: iteration reaches last index of npe_per_time in evtID = {0:d}".format(evt))
                    break

                # check if value bins_time[index4 + num2] is <= threshold2 and if following bins_time[index4 + num2 + 1]
                # and bins_time[index4 + num2 + 2] are also <= threshold2:
                if (npe_per_time[index4 + num2] <= threshold2 and npe_per_time[index4 + num2 + 1] <= threshold2 and
                        npe_per_time[index4 + num2 + 2] <= threshold2):
                    # get hittime, when delayed signal peak ends, in ns:
                    end_peak = bins_time[index4 + num2]

                    # get index, where npe_hittime is below threshold2:
                    index_after_peak1 = index4 + num2
                    break
                else:
                    # npe_per_time[index4 + num2] or npe_per_time[index4 + num2 + 1] are above threshold2:
                    # add nPE of this bin to sum_pe_peak:
                    sum_pe_peak = sum_pe_peak + npe_per_time[index4 + num2]

            # print("number of pe in delayed signal = {0:d}".format(sum_pe_peak))

            # first peak is analyzed:
            # check, if number of PE in this signal peak is above threshold:
            if min_pe < sum_pe_peak:
                # PE of signal peak above threshold:
                # set delayed flag:
                num_del_signal = 1
                # set number of pe of delayed signal:
                num_pe_del = sum_pe_peak
            else:
                num_del_signal = 0

            # after analyzing the first peak, break out of for-loop (possible second delayed signal will be checked
            # below)
            break
        else:
            # go to next bin:
            continue

    return num_del_signal, index_after_peak1, num_pe_del, begin_peak, end_peak


def analyze_delayed_signal(npe_per_time, bins_time, first_index, threshold, threshold2, min_pe_delayed, max_pe_delayed,
                           evt):
    """
    function to analyze the time window, where a delayed signal could be. This time window is check for a possible
    delayed signal (signal has to be greater than 'threshold' and nPE of signal peak must be between 'min_pe_delayed'
    and 'max_pe_delayed')

    :param npe_per_time: number of pe as function of hittime for time window of delayed signal (array of integer)
    :param bins_time: array of bins, which contains the information about the hittime in ns (array of float)
    :param first_index: first index of npe_per_time, that should be analyzed (integer)
    :param threshold: threshold in nPE per bin, that a delayed signal must have (integer)
    :param threshold2: threshold in nPE per bin, that specifies the minimal nPE corresponding to a peak
    :param min_pe_delayed:  minimum number of PE for delayed energy cut of neutron capture on Hydrogen
                            (values from check_delayed_energy.py)
    :param max_pe_delayed:  maximum number of PE for delayed energy cut of neutron capture on Hydrogen
                            (values from check_delayed_energy.py)
    :param evt: event ID of the event

    :return:
    """
    # preallocate number of delayed signal and index_after_peak1 and number of pe in the delayed signal:
    num_del_signal = 0
    index_after_peak1 = len(npe_per_time)
    num_pe_del = 0
    begin_peak = 0
    end_peak = 2000000

    # loop over values of the histogram bins of delayed window:
    for index4 in range(first_index, len(npe_per_time)):
        # check if number of PE in this bin is above the threshold:
        if npe_per_time[index4] > threshold:
            # possible delayed signal (signal in delayed window):

            # check if index4 = first_value -> first value of npe_per_time above threshold -> break from for loop and
            # return index_after_peak1 = index4 + 1:
            if index4 == first_index:
                index_after_peak1 = index4 + 1
                print("WARNING: npe_per_time[{1:d}] > threshold in evtID = {0:d}".format(evt, index4))
                break

            # calculate number of PEs in this signal peak:
            # preallocate sum of PE per bin in the signal peak:
            sum_pe_peak = 0

            # add nPE of npe_per_time[index4] to sum_pe_peak:
            sum_pe_peak = sum_pe_peak + npe_per_time[index4]

            # check previous bins:
            for num1 in range(1, 100, 1):
                # check if value in previous bins_time[index4 - num1] is above threshold2:
                if npe_per_time[index4 - num1] > threshold2:
                    # add nPE of this bin to sum_pe_peak:
                    sum_pe_peak = sum_pe_peak + npe_per_time[index4 - num1]

                else:
                    # hittime, when signal peak begins, in ns:
                    begin_peak = bins_time[index4 - num1]
                    break

            # check following bins:
            for num2 in range(1, 1000, 1):
                # check if (index4 + num2 + 2) is in the range of npe_per_time array:
                if (index4 + num2 + 2) >= len(npe_per_time):
                    print("Warning: iteration reaches last index of npe_per_time in evtID = {0:d}".format(evt))
                    break

                # check if value bins_time[index4 + num2] is <= threshold2 and if following bins_time[index4 + num2 + 1]
                # and bins_time[index4 + num2 + 2] are also <= threshold2:
                if (npe_per_time[index4 + num2] <= threshold2 and npe_per_time[index4 + num2 + 1] <= threshold2 and
                        npe_per_time[index4 + num2 + 2] <= threshold2):
                    # get hittime, when delayed signal peak ends, in ns:
                    end_peak = bins_time[index4 + num2]

                    # get index, where npe_hittime is below threshold2:
                    index_after_peak1 = index4 + num2
                    break
                else:
                    # npe_per_time[index4 + num2] or npe_per_time[index4 + num2 + 1] are above threshold2:
                    # add nPE of this bin to sum_pe_peak:
                    sum_pe_peak = sum_pe_peak + npe_per_time[index4 + num2]

            # print("number of pe in delayed signal = {0:d}".format(sum_pe_peak))

            # set number of pe of delayed signal:
            num_pe_del = sum_pe_peak

            # first peak is analyzed:
            # check, if number of PE in this signal peak agree with delayed energy cut:
            if min_pe_delayed < sum_pe_peak < max_pe_delayed:
                # PE of signal peak agree with delayed energy cut
                # set delayed flag:
                num_del_signal = 1
            else:
                num_del_signal = 0
                # print("number of pe in delayed signal = {0:d}".format(sum_pe_peak))

            # after analyzing the first peak, break out of for-loop (possible second delayed signal will be checked
            # below)
            break
        else:
            # go to next bin:
            continue

    return num_del_signal, index_after_peak1, num_pe_del, begin_peak, end_peak


def position_smearing(pos_init, e_in_mev):
    """
    function to smear the initial position of an event with the vertex resolution.

    :param pos_init: initial x,y or z position in mm of the event
    :param e_in_mev: energy of the event in MeV (normally sum of quenched deposited energy of primary particles
    from prmtrkdep-tree is used)
    :return: reconstructed (smeared) x,y or z position of the event
    """
    # mean of the gaussian equal to initial position in mm:
    mu = pos_init
    # sigma of the gaussian given by the vertex resolution (sigma = 120 mm / sqrt(E[MeV])) (source YellowBook, p. 157):
    sigma = 120.0 / np.sqrt(e_in_mev)
    # generate random number from normal distribution defined by mu and sigma. This represents the reconstructed
    # position in mm:
    pos_recon = np.random.normal(mu, sigma)

    return pos_recon


def conversion_npe_to_evis(number_photo_electron):
    """
    Function to calculate the visible energy in MeV for a given number of p.e. for the prompt signal.
    This function is the result of linear fit from script check_conversion_npe_mev.py.

    :param number_photo_electron: number of photo-electrons of the prompt signal
    :return: quenched deposited energy (visible energy in MeV)
    """
    # TODO-me: parameters of the fit has to be checked!!!!!!!!!
    # first fit parameter (slope) in MeV/nPE:
    parameter_a = 0.0007483

    energy = parameter_a * number_photo_electron

    return energy


def get_pmt_position(filename):
    """
    function to read the file PMT_position.root to get the PMT position as function of the PMT ID.

    Position of all PMTs (20inch and 3inch) is saved.

    :param filename: path and file name of PMT_position.root
    :return:
    """
    # preallocate list for PMT ID:
    arr_pmt_id = []
    # preallocate list for PMT position in mm:
    arr_pmt_x = []
    arr_pmt_y = []
    arr_pmt_z = []

    # read root file:
    rfile = ROOT.TFile(filename)
    # get the pmtpos-tree from the TFile:
    rtree = rfile.Get("pmtpos")
    # get the number of events in the tree:
    number_pmts = rtree.GetEntries()

    # loop over PMTs:
    for index in range(number_pmts):
        rtree.GetEntry(index)

        # get PMT ID:
        pmt_id = int(rtree.GetBranch('pmtID').GetLeaf('pmtID').GetValue())
        arr_pmt_id.append(pmt_id)
        # get x-position of PMT in mm:
        pmt_x = float(rtree.GetBranch('x').GetLeaf('x').GetValue())
        arr_pmt_x.append(pmt_x)
        # get y-position of PMT in mm:
        pmt_y = float(rtree.GetBranch('y').GetLeaf('y').GetValue())
        arr_pmt_y.append(pmt_y)
        # get z-position of PMT in mm:
        pmt_z = float(rtree.GetBranch('z').GetLeaf('z').GetValue())
        arr_pmt_z.append(pmt_z)

    return arr_pmt_id, arr_pmt_x, arr_pmt_y, arr_pmt_z


def get_20inchpmt_tts(filename):
    """
    function to read the file PmtData.root (originally saved in folder
    /home/astro/blum/juno/JUNO-SOFT/data/Simulation/ElecSim/) to get the time spread (TTS)
    of the 20 inch PMTs as function of the PMT ID.

    TTS of the small 3inch PMTs has to be considered differently.

    :param filename: path and file name of PmtData.root
    :return:
    """
    # preallocate list for PMT ID:
    arr_pmt_id = []
    # preallocate list for PMT 'time resolution' in ns (must be sigma of normal distribution):
    arr_sigma_time = []

    # read root file:
    rfile = ROOT.TFile(filename)
    # get the pmtpos-tree from the TFile:
    rtree = rfile.Get("PmtData")

    # set the number of 20 inch PMTs in Central Detector by hand:
    number_pmts = 17739

    # INFO-me: in PmtData.root, 20000 large PMTs are saved. In detsim version J18v1r1-Pre1 (used for NC simulation),
    # INFO-me: only 17739 large PMTs are installed. -> Therefore number_pmts = 17739!

    # loop over 20inch PMTs:
    for index in range(number_pmts):
        # get event:
        rtree.GetEntry(index)

        # get PMT ID:
        pmt_id = int(rtree.GetBranch('pmtId').GetLeaf('m_pmtId').GetValue())
        arr_pmt_id.append(pmt_id)

        # get timeSpread of PMT in ns:
        timespread = float(rtree.GetBranch('timeSpread').GetLeaf('&m_timeSpread').GetValue())

        # The timespread (given by the TTS) defines the time resolution of the PMT, which is defined as FWHM.
        # (see Flyckt: "Photomultiplier Tubes, principles and applications", 2-11).
        # Therefore you have to calculate sigma of the normal distribution with the FWHM.
        # For normal distribution: FWHM = 2 * sqrt(2 * ln(2)) * sigma:
        # sigma in ns:
        sigma = timespread / (2 * np.sqrt(2 * np.log(2)))

        arr_sigma_time.append(sigma)

    return arr_pmt_id, arr_sigma_time


def check_neutron_cut(input_path, number_file, output_path, min_hittime, max_hittime, threshold, threshold2,
                      binwidth_hittime, min_pe_delayed, max_pe_delayed, number_entries_input, save_hittime):
    """
    function to read the user_atmoNC_.root file and check, if there are events with no neutron, but with a delayed
    signal.

    :param input_path: path to input root files from tut_detsim.py: user_atmoNC_{}.root
    :param number_file: number of the input root file: e.g. 2 -> user_atmoNC_2.root
    :param output_path: path, where plots of hittime should be saved
    :param min_hittime: minimum hittime in ns of a possible delayed signal
    :param max_hittime: maximum hittime in ns of a possible delayed signal
    :param threshold: threshold of number of PEs per bin for possible delayed signal
    :param threshold2: threshold2 of number of PEs per bin (signal peak is summed as long as nPE is above threshold2)
    :param binwidth_hittime: bin-width of the hittime histograms in ns
    :param min_pe_delayed: minimum number of PE for delayed energy cut (values from check_delayed_energy.py)
    :param max_pe_delayed: maximum number of PE for delayed energy cut (values from check_delayed_energy.py)
    :param number_entries_input: number of entries, that the input files should have (integer), normally = 100
    :param save_hittime: boolean variable, if True the hittime-histogram is saved
    :return:
    """

    # load the ROOT file:
    rfile = ROOT.TFile(input_path + "user_atmoNC_{0:d}.root".format(number_file))
    # get the "evt"-TTree from the TFile:
    rtree_evt = rfile.Get("evt")
    # get the "geninfo"-TTree from the TFile:
    rtree_geninfo = rfile.Get("geninfo")
    # get the "prmtrkdep"-TTree from the TFile:
    rtree_prmtrkdep = rfile.Get("prmtrkdep")

    # get the number of events in the 'evt' Tree:
    number_events_evt = rtree_evt.GetEntries()
    # get the number of events in the geninfo Tree:
    number_events_geninfo = rtree_geninfo.GetEntries()
    # get the number of events in the prmtrkdep Tree:
    number_events_prmtrkdep = rtree_prmtrkdep.GetEntries()
    if number_events_geninfo == number_events_prmtrkdep and number_events_geninfo == number_events_evt:
        number_events = number_events_geninfo
    else:
        sys.exit("ERROR: number of events in Trees are NOT equal!!")

    # check if number_events is equal to number_entries_input (if not, the detector simulation was incorrect!!):
    if number_events != number_entries_input:
        sys.exit("ERROR: number of events are not equal to {0:d} -> Detector Simulation not correct!"
                 .format(number_entries_input))

    # preallocate variables:
    # number of events with at least one neutron:
    number_neutron = 0
    # number of events with no neutron:
    number_no_neutron = 0
    # number of events with no neutron and no delayed signal:
    number_no_delayed = 0
    # number of events with no neutron, but with delayed signal:
    number_delayed = 0
    # number of possible delayed signal (delayed signal but not correct number of PE):
    number_possible_delayed = 0
    # number of possible second delayed signal (another possible delayed signal after first signal)
    # (only hittime is checked):
    number_possible_second_delayed = 0

    # loop over every event, i.e. every entry, in the TTree:
    # for event in range(47, 48, 1):
    for event in range(number_events):

        print("analyze event {0:d}".format(event))

        """ first read the 'geninfo' Tree: """
        # get the current event in the Tree:
        rtree_geninfo.GetEntry(event)

        # get the value of the event ID:
        evt_id_geninfo = int(rtree_geninfo.GetBranch('evtID').GetLeaf('evtID').GetValue())

        # get number of particles in the event:
        n_par_geninfo = int(rtree_geninfo.GetBranch('nInitParticles').GetLeaf('nInitParticles').GetValue())

        # check if there is at least 1 particle in the Tree:
        if n_par_geninfo == 0:
            number_no_neutron = number_no_neutron + 1
            number_no_delayed = number_no_delayed + 1
            # no particle in the event -> go to next event
            continue

        # preallocate arrays for PDG ID and ExitT:
        pdg_geninfo = np.array([])
        exit_time = np.array([])

        # set neutron flag:
        neutron_in_event = False

        # loop over number of particles in the event:
        for index1 in range(n_par_geninfo):
            # get PDG ID of the initial particles:
            pdgid_geninfo = int(rtree_geninfo.GetBranch('InitPDGID').GetLeaf('InitPDGID').GetValue(index1))

            # check if neutron:
            if pdgid_geninfo == 2112:
                # set neutron flag (boolean):
                neutron_in_event = True
                # close for loop over number of particles:
                break

            # append PDG ID:
            pdg_geninfo = np.append(pdg_geninfo, pdgid_geninfo)
            # get exit time in ns:
            exit_t = float(rtree_geninfo.GetBranch('ExitT').GetLeaf('ExitT').GetValue(index1))
            # append exit time:
            exit_time = np.append(exit_time, exit_t)

        if neutron_in_event:
            # increment number_neutron:
            number_neutron = number_neutron + 1
            # go to next event:
            continue
        else:
            # increment number_no_neutron:
            number_no_neutron = number_no_neutron + 1

        """ read prmtrkdep tree: """
        # get the current event in the Tree:
        rtree_prmtrkdep.GetEntry(event)

        # get evt ID of tree:
        evt_id_prm = int(rtree_prmtrkdep.GetBranch('evtID').GetLeaf('evtID').GetValue())

        # get number of particles in the event:
        n_par_prm = int(rtree_prmtrkdep.GetBranch('nInitParticles').GetLeaf('nInitParticles').GetValue())

        if n_par_prm == n_par_geninfo:
            n_par = n_par_prm
        else:
            sys.exit("Number of particles in 'geninfo' NOT equal to number of particles in 'prmtrkdep'")

        # preallocate arrays for edep and Qedep:
        edep_arr = np.array([])
        qedep_arr = np.array([])

        # loop over number of particles:
        for index2 in range(n_par):
            # get edep in MeV:
            edep = float(rtree_prmtrkdep.GetBranch('edep').GetLeaf('edep').GetValue(index2))
            # get Qedep in MeV:
            qedep = float(rtree_prmtrkdep.GetBranch('Qedep').GetLeaf('Qedep').GetValue(index2))

            # append edep and qedep to arrays:
            edep_arr = np.append(edep_arr, edep)
            qedep_arr = np.append(qedep_arr, qedep)

        """ read the "evt" Tree"""
        # get the current event in the TTree:
        rtree_evt.GetEntry(event)

        # get the value of the event ID:
        evt_id_evt = int(rtree_evt.GetBranch('evtID').GetLeaf('evtID').GetValue())

        if evt_id_evt == evt_id_geninfo and evt_id_evt == evt_id_prm:
            evt_id = evt_id_evt
        else:
            sys.exit("Event ID's of the three trees are NOT equal for 1 event!")

        # get number of photons of this event:
        n_photons = int(rtree_evt.GetBranch('nPhotons').GetLeaf('nPhotons').GetValue())
        print("number of photons = {0:d}".format(n_photons))

        # check, if n_photons = 0:
        if n_photons == 0:
            number_no_neutron = number_no_neutron + 1
            number_no_delayed = number_no_delayed + 1
            # no photons in event -> go to next event
            continue

        # preallocate variables:
        # number of pe in event:
        number_pe_event = 0

        # preallocate histogram, where hittimes are saved:
        # set bin-edges of hittime histogram in ns:
        bins_hittime = np.arange(0, max_hittime+2*binwidth_hittime, binwidth_hittime)
        # preallocate empty array to build default hittime-histogram:
        hittime_array = np.array([])
        # build default hittime histogram:
        npe_per_hittime, bin_edges_hittime = np.histogram(hittime_array, bins_hittime)

        # loop over every photon in the event:
        for index3 in range(n_photons):

            # get PMT ID, where photon is absorbed:
            pmt_id = int(rtree_evt.GetBranch('pmtID').GetLeaf('pmtID').GetValue(index3))

            # only 20 inch PMTs (PMT ID of 20 inch PMTs are below 21000, PMT ID of 3 inch PMTs start at 290000):
            if pmt_id < 25000:
                # get nPE for this photon:
                n_pe = int(rtree_evt.GetBranch('nPE').GetLeaf('nPE').GetValue(index3))
                # check, if photon produces only 1 PE:
                if n_pe != 1:
                    print("{1:d} PE for 1 photon in event {0:d} in file {2}".format(evt_id_evt, n_pe, rootfile_input))

                # add n_pe to number_pe_event:
                number_pe_event = number_pe_event + n_pe

                # get hittime of this photon:
                hit_time = float(rtree_evt.GetBranch('hitTime').GetLeaf('hitTime').GetValue(index3))

                # add hittime to default hittime histogram (build histogram with value 'hit_time' and bins_hittime;
                # take the values of the 'new' histogram and add them to default hittime histogram):
                npe_per_hittime += np.histogram(hit_time, bins_hittime)[0]

            else:
                continue

        """ check hittime histogram of the event: """
        # set possible delayed flag:
        possible_delayed_signal = False
        # set delayed flag:
        delayed_signal = False
        # set possible second delayed flag:
        possible_second_delayed = False

        # get index, where bin_edges_hittime = min_hittime, (min_hittime - 0) / binwidth_hittime,
        # e.g. (2000 - 0) / 5 = 400, bin_edges_hittime[400] = 2000.
        index_min_hittime = int(min_hittime / binwidth_hittime)

        # calculate nPE as function of hittime only for the delayed time window (from min_hittime to max_hittime):
        npe_per_hittime_window = npe_per_hittime[index_min_hittime:-1]
        # bin edges of hittime histogram only fro the delayed time window:
        bins_hittime_window = bin_edges_hittime[index_min_hittime:-1]

        # loop over values of the histogram bins of delayed window:
        for index4 in range(len(npe_per_hittime_window)):
            # check if number of PE in this bin is above the threshold:
            if npe_per_hittime_window[index4] > threshold:
                # possible delayed signal (signal in delayed window):

                # calculate number of PEs in this signal peak:
                # preallocate sum of PE per bin in the signal peak:
                sum_pe_peak = 0

                # add nPE of bins_hittime_window[index4] to sum_pe_peak:
                sum_pe_peak = sum_pe_peak + npe_per_hittime_window[index4]

                # check previous bins:
                for num in range(1, 100, 1):
                    # check if value in previous bins_hittime_window[index4 - num] is above threshold2:
                    if npe_per_hittime_window[index4 - num] > threshold2:
                        # add nPE of this bin to sum_pe_peak:
                        sum_pe_peak = sum_pe_peak + npe_per_hittime_window[index4 - num]

                    else:
                        # hittime, when signal peak begins, in ns:
                        begin_peak = bins_hittime_window[index4 - num]

                        break

                # check following bins:
                for num in range(1, 1000, 1):
                    # check if value in following bins_hittime_window[index4 + num] is above threshold2:
                    if npe_per_hittime_window[index4 + num] > threshold2:
                        # add nPE of this bin to sum_pe_peak:
                        sum_pe_peak = sum_pe_peak + npe_per_hittime_window[index4 + num]

                    else:
                        # hittime, when signal peak ends, in ns:
                        end_peak = bins_hittime_window[index4 + num]

                        # get index, where npe_hittime is below threshold2:
                        index_after_peak1 = index4 + num
                        break

                # first peak is analyzed:
                # check, if number of PE in this signal peak agree with delayed energy cut:
                if min_pe_delayed < sum_pe_peak < max_pe_delayed:
                    # PE of signal peak agree with delayed energy cut
                    # set delayed flag:
                    delayed_signal = True
                else:
                    # PE of signal peak do NOT agree with delayed energy cut:
                    # set possible delayed flag:
                    possible_delayed_signal = True

                # after analyzing the first peak, break out of for-loop (possible second delayed signal will be checked
                # below)
                break
            else:
                # go to next bin:
                continue

        # When there is a delayed or possible delayed signal in the first peak, check if there is also a second
        # possible peak in the event:
        if delayed_signal or possible_delayed_signal:
            # check if there is a possible second delayed signal

            # loop over values of histogram bins of delayed time window to get possible second peak:
            for index5 in range(index_after_peak1, len(npe_per_hittime_window), 1):
                # check if number of PE in this bin is above the threshold:
                if npe_per_hittime_window[index5] > threshold:
                    # possible second delayed signal:

                    # set possible second delayed signal flag:
                    possible_second_delayed = True

                    # possible second delayed signal after delayed or possible delayed signal
                    # (only time window is checked)
                    number_possible_second_delayed = number_possible_second_delayed + 1
                    break

                else:
                    # go to next bin:
                    continue

        # increment numbers:
        if delayed_signal:
            # 1 delayed signal (agree with time window and delayed energy cut)
            number_delayed = number_delayed + 1

            # plot hittime_event in histogram:
            if save_hittime:
                h1 = plt.figure(1, figsize=(15, 8))
                # plot histogram:
                plt.plot(bin_edges_hittime[0:-1], npe_per_hittime, drawstyle='steps',
                         label='Event Info:\nInitPDGID = {0}\nExitT = {1} ns\nedep = {2} MeV\nQedep = {3} MeV\n'
                               'Analysis:\nnPE in pulse = {4}\nstart hittime pulse = {5} ns\n'
                               'end hittime pulse = {6} ns'.format(pdg_geninfo, exit_time, edep_arr, qedep_arr,
                                                                   sum_pe_peak, begin_peak, end_peak))

                plt.xlabel("hit-time in ns", fontsize=13)
                # INFO-me: ylabel is only equal to number of PE, if nPE == 1 for all photons (1 PE each photon)
                plt.ylabel("number of PE per bin (bin-width = {0:=.1f} ns)".format(binwidth_hittime), fontsize=13)
                plt.title("hittime on 20inch PMTs for atmoNC events with delayed signal and without neutron\n"
                          "(agree with time and delayed energy cut)",
                          fontsize=18)
                plt.legend()
                plt.grid()
                plt.savefig(output_path + "delSig_hittime_atmoNC_{0:d}_evt{1:d}.png".format(number_file, evt_id))
                plt.close()

        elif possible_delayed_signal:
            # 1 possible delayed signal (agree with time window, but NOT with delayed energy cut)
            number_possible_delayed = number_possible_delayed + 1

            # plot hittime_event in histogram:
            if save_hittime:
                h1 = plt.figure(1, figsize=(15, 8))
                # plot histogram:
                plt.plot(bin_edges_hittime[0:-1], npe_per_hittime, drawstyle='steps',
                         label='Event Info:\nInitPDGID = {0}\nExitT = {1} ns\nedep = {2} MeV\nQedep = {3} MeV\n'
                               'Analysis:\nnPE in pulse = {4}\nstart hittime pulse = {5} ns\n'
                               'end hittime pulse = {6} ns'
                         .format(pdg_geninfo, exit_time, edep_arr, qedep_arr, sum_pe_peak, begin_peak, end_peak))

                plt.xlabel("hit-time in ns", fontsize=13)
                # INFO-me: ylabel is only equal to number of PE, if nPE == 1 for all photons (1 PE each photon)
                plt.ylabel("number of PE per bin (bin-width = {0:=.1f} ns)".format(binwidth_hittime), fontsize=13)
                plt.title("hittime on 20inch PMTs for atmoNC events with delayed signal and without neutron\n"
                          "(agree with time cut but not with delayed energy cut)",
                          fontsize=18)
                plt.legend()
                plt.grid()
                plt.savefig(output_path + "posDelSig_hittime_atmoNC_{0:d}_evt{1:d}.png".format(number_file, evt_id))
                plt.close()

        else:
            # no delayed and no possible delayed -> increment number_no_delayed:
            number_no_delayed = number_no_delayed + 1

    return number_events, number_neutron, number_no_neutron, number_no_delayed, number_delayed, \
           number_possible_delayed, number_possible_second_delayed


def conversion_npe_mev(rootfile_input, number_entries_input, radius_cut):
    """
    function to read user_proton_..._.root files from detsim simulation to convert neutron/proton energy of possible
    prompt signal from number of photo-electron (nPE) to MeV.

    :param rootfile_input: input root file from tut_detsim.py: e.g. user_proton_10_MeV_0.root (string)
    :param number_entries_input: number of entries, that the input files should have (int)
    :param radius_cut: radius, that defines the volume cut, in mm

    :return:
    """
    # load ROOT file:
    rfile = ROOT.TFile(rootfile_input)
    # get the "evt"-TTree from the TFile:
    rtree_evt = rfile.Get("evt")
    # get the number of events in the geninfo Tree:
    number_events_evt = rtree_evt.GetEntries()

    # get the "geninfo"-TTree from the TFile:
    rtree_geninfo = rfile.Get("geninfo")
    # get the number of events in the geninfo Tree:
    number_events_geninfo = rtree_geninfo.GetEntries()

    # get the 'prmtrkdep' tree:
    rtree_prmtrkdep = rfile.Get('prmtrkdep')
    # get number of events in the tree:
    number_events_prmtrkdep = rtree_prmtrkdep.GetEntries()

    # check if number of events are equal in both trees:
    if number_events_geninfo == number_events_evt and number_events_evt == number_events_prmtrkdep:
        number_events = number_events_geninfo
    else:
        sys.exit("ERROR: number of events in t Trees are NOT equal!!")

    # check if number_events is equal to number_entries_input (if not, the detector simulation was incorrect!!):
    if number_events != number_entries_input:
        sys.exit("ERROR: number of events are not equal to {0:d} -> Detector Simulation not correct!"
                 .format(number_entries_input))

    # preallocate array, where number of PE per event (e.g. per gamma) is saved (np.array of int):
    number_pe = np.array([])

    # preallocate array. where initial momentum per event is saved (in MeV):
    init_momentum = np.array([])
    # preallocate array, where deposit energy per event is saved (in MeV):
    edep = np.array([])
    # preallocate array, where quenched deposit energy (visible energy) is saved (in MeV):
    qedep = np.array([])

    # loop over every event, i.e. every entry, in the TTree:
    for event in range(number_events):

        """ read 'geninfo' tree to check initial energy and to apply volume cut: """
        # get the current event in the TTree:
        rtree_geninfo.GetEntry(event)

        # get event ID of geninfo-tree:
        evt_id_geninfo = int(rtree_geninfo.GetBranch('evtID').GetLeaf('evtID').GetValue())

        # get number of particles of this event:
        n_particles = int(rtree_geninfo.GetBranch('nInitParticles').GetLeaf('nInitParticles').GetValue())

        # check, that only 1 particle per event:
        if n_particles == 1:

            # get the initial position in x, y, z in mm:
            x_init = float(rtree_geninfo.GetBranch('InitX').GetLeaf('InitX').GetValue())
            y_init = float(rtree_geninfo.GetBranch('InitY').GetLeaf('InitY').GetValue())
            z_init = float(rtree_geninfo.GetBranch('InitZ').GetLeaf('InitZ').GetValue())

            # calculate the radius in mm:
            r_init = np.sqrt(x_init**2 + y_init**2 + z_init**2)

            if r_init > radius_cut:
                # apply volume cut:
                continue

            # get initial momenta in x, y, z in MeV:
            px_init = float(rtree_geninfo.GetBranch('InitPX').GetLeaf('InitPX').GetValue())
            py_init = float(rtree_geninfo.GetBranch('InitPY').GetLeaf('InitPY').GetValue())
            pz_init = float(rtree_geninfo.GetBranch('InitPZ').GetLeaf('InitPZ').GetValue())

            # total initial momentum in MeV:
            momentum = np.sqrt(px_init ** 2 + py_init ** 2 + pz_init ** 2)
            # append momentum to array:
            init_momentum = np.append(init_momentum, momentum)

        else:
            print("{0:d} particles in event {1:d} in file {2}".format(n_particles, evt_id_geninfo, rootfile_input))

        """ read 'evt' tree: """
        # get the current event in the TTree:
        rtree_evt.GetEntry(event)

        # get event ID of evt-tree:
        evt_id_evt = int(rtree_evt.GetBranch('evtID').GetLeaf('evtID').GetValue())

        # get number of photons of this event:
        n_photons = int(rtree_evt.GetBranch('nPhotons').GetLeaf('nPhotons').GetValue())

        """ preallocate variables: """
        # number of pe in event:
        number_pe_event = 0
        # hittime of photons in event in ns:
        hittime_event = []

        # loop over every photon in the event:
        for index in range(n_photons):

            # Read all PMTs (20inch AND 3inch):

            # get nPE for this photon:
            n_pe = int(rtree_evt.GetBranch('nPE').GetLeaf('nPE').GetValue(index))
            # check, if photon produces only 1 PE:
            if n_pe != 1:
                print("{1:d} pe for 1 photon in event {0:d} in file {2}".format(evt_id_evt, n_pe, rootfile_input))

            # add n_pe to number_pe_event:
            number_pe_event = number_pe_event + n_pe

            # get hittime of this photon:
            hit_time = float(rtree_evt.GetBranch('hitTime').GetLeaf('hitTime').GetValue(index))
            # append hit_time to hittime_event:
            hittime_event.append(hit_time)

        # append number of PE of this event to number_pe array:
        number_pe = np.append(number_pe, number_pe_event)

        """ read 'prmtrkdep' tree to check deposit energy and quenched deposit energy: """
        # get the current event in the TTree:
        rtree_prmtrkdep.GetEntry(event)

        # get number of particles:
        n_part = int(rtree_prmtrkdep.GetBranch('nInitParticles').GetLeaf('nInitParticles').GetValue())
        # check, that only 1 particle per event:
        if n_part == 1:

            # get deposit energy in MeV:
            edep_value = float(rtree_prmtrkdep.GetBranch('edep').GetLeaf('edep').GetValue())
            edep = np.append(edep, edep_value)

            # get quenched energy in MeV:
            qedep_value = float(rtree_prmtrkdep.GetBranch('Qedep').GetLeaf('Qedep').GetValue())
            # consider energy resolution of detector:
            if qedep_value > 0:
                # get the value of sigma of energy resolution for value of qedep_value:
                sigma_energy = energy_resolution(qedep_value)
                # generate normal distributed random number with mean = qedep_value and sigma = sigma_energy:
                qedep_value = np.random.normal(qedep_value, sigma_energy)
            qedep = np.append(qedep, qedep_value)

    return number_pe, hittime_event, init_momentum, edep, qedep


def read_gamma_delayed_signal(rootfile_input, number_entries_input, radius_cut, x_position_pmt, y_position_pmt,
                              z_position_pmt, sigma_t_20inch, sigma_t_3inch, c_eff, min_t, max_t, bin_width,
                              threshold1, threshold2):
    """
    function to read user_gamma_..._.root file from detsim simulation to convert gamma energy of possible delayed
    signal (from neutron capture) from MeV to number of PE.

    Do the exact same analysis of the signal like for the delayed signal in prompt_signal_preselected_evts.py.

    :param rootfile_input: input root file from tut_detsim.py: e.g. user_gamma_2_2_MeV_0.root (string)
    :param number_entries_input: number of entries, that the input files should have (integer), normally = 10
    :param radius_cut: radius, that defines the volume cut, in mm
    :param x_position_pmt: array of x position of the PMTs in mm
    :param y_position_pmt: array of y position of the PMTs in mm
    :param z_position_pmt: array of z position of the PMTs in mm
    :param sigma_t_20inch: array of time resolution (sigma) of the 20 inch PMTs in ns
    :param sigma_t_3inch: value of the time resolution (sigma) of the 3inch PMTs in ns
    :param c_eff: effective speed of light in liquid scintillator in mm/ns
    :param min_t: time, where the hittime histogram starts, in ns
    :param max_t: time, where the hittime histogram ends, in ns
    :param bin_width: binwidth of the hittime histogram in ns
    :param threshold1: threshold of number of PE per bin for signal (bin-width = 5 ns)
    :param threshold2: threshold2 of number of PEs per bin (signal peak is summed as long as nPE is above threshold2)

    :return:

    """
    # load ROOT file:
    rfile = ROOT.TFile(rootfile_input)
    # get the "evt"-TTree from the TFile:
    rtree_evt = rfile.Get("evt")
    # get the number of events in the geninfo Tree:
    number_events_evt = rtree_evt.GetEntries()

    # get the "geninfo"-TTree from the TFile:
    rtree_geninfo = rfile.Get("geninfo")
    # get the number of events in the geninfo Tree:
    number_events_geninfo = rtree_geninfo.GetEntries()

    # get "prmtrkdep" tree from root file:
    rtree_prmtrkdep = rfile.Get("prmtrkdep")
    # get number of events in prmtrkdep tree:
    number_events_prmtrkdep = rtree_prmtrkdep.GetEntries()

    # check if number of events are equal in both trees:
    if number_events_geninfo == number_events_evt and number_events_geninfo == number_events_prmtrkdep:
        number_events = number_events_geninfo
    else:
        sys.exit("ERROR: number of events in t Trees are NOT equal!!")

    # check if number_events is equal to number_entries_input (if not, the detector simulation was incorrect!!):
    if number_events != number_entries_input:
        sys.exit("ERROR: number of events are not equal to {0:d} -> Detector Simulation not correct!"
                 .format(number_entries_input))

    # preallocate array, where number of PE per event (e.g. per gamma) is saved (np.array of int):
    number_pe = np.array([])

    # preallocate number of analyzed events:
    number_analyzed = 0

    # loop over every event, i.e. every entry, in the TTree:
    for event in range(number_events):

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
        x_reconstructed = position_smearing(x_init, qedep_prmtrkdep)
        y_reconstructed = position_smearing(y_init, qedep_prmtrkdep)
        z_reconstructed = position_smearing(z_init, qedep_prmtrkdep)

        # calculate distance to detector center in mm:
        r_reconstructed = np.sqrt(x_reconstructed**2 + y_reconstructed**2 + z_reconstructed**2)

        # check if event passes the volume cut:
        if r_reconstructed >= radius_cut:
            # event is rejected by volume cut.
            print("file {0}, event = {1:d}: r_init = {2:0.2f} mm".format(rootfile_input, event, r_reconstructed))
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
            pmtid = int(rtree_evt.GetBranch('pmtID').GetLeaf('pmtID').GetValue(index1))

            """ time of flight correction: """
            # get hittime of PMT from tree in ns:
            hittime = float(rtree_evt.GetBranch('hitTime').GetLeaf('hitTime').GetValue(index1))

            # get position of the PMT with specific pmtID (pmtID is ascending number from 0 to 17738 (17739 large PMTs)
            # and from 300000 to 336571 (36572 small PMTs)).
            # For large PMTs -> For 20inch PMTs, the pmtID is equal to index of x,y,z_pos_pmt array.
            # For small PMTs -> For 3inch PMTs, the pmtID - (300000 - 17739) is equal to index of x,y,z_pos_pmt array.
            # check if PMT is 20 inch or 3inch (pmtID < 50000 means 20inch PMT):
            if pmtid < 50000:
                # 20inch PMT:
                # get PMT position in mm from arrays:
                x_pmt = x_position_pmt[pmtid]
                y_pmt = y_position_pmt[pmtid]
                z_pmt = z_position_pmt[pmtid]
            else:
                # 3inch PMT:
                # calculate index of pos_pmt array that correspond to pmtID of 3inch PMTs (for example:
                # first small PMT: 300000-282261 = 17739, last small PMT: 336571-282261 = 54310)
                index_3inch = pmtid - 282261
                # get PMT position in mm from arrays:
                x_pmt = x_position_pmt[index_3inch]
                y_pmt = y_position_pmt[index_3inch]
                z_pmt = z_position_pmt[index_3inch]

            # calculate distance between reconstructed position of event and position of PMT (in mm):
            distance_tof = np.sqrt((x_reconstructed - x_pmt)**2 + (y_reconstructed - y_pmt)**2 +
                                   (z_reconstructed - z_pmt)**2)

            # calculate time of flight in ns:
            time_of_flight = distance_tof / c_eff

            """ time resolution of PMT: """
            # get time resolution of PMT with specific pmtID (pmtID is ascending number from 0 to 17738 (17739 large
            # PMTs)) -> For 20inch PMTs, the pmtID is equal to index of sigma_time_20inch array.
            # check if PMT is 20 inch or 3inch (pmtID < 50000 means 20inch PMT):
            if pmtid < 50000:
                # 20inch PMT:
                # get time resolution (sigma) of PMT in ns from array:
                sigma_pmt = sigma_t_20inch[pmtid]

            else:
                # 3inch PMT:
                sigma_pmt = sigma_t_3inch

            # consider time resolution of PMT by generating normal distributed random number with mu = hittime and
            # sigma = sigma_pmt (only the hittime at the PMT must be smeared, not the time-of-flight):
            hittime_tts = np.random.normal(hittime, sigma_pmt)

            """ calculate the 'real' hittime of the photon in ns: """
            hittime_real = hittime_tts - time_of_flight
            if hittime_real < min_t:
                print("------")
                print(hittime_real)
                print(pmtid)
                print(sigma_pmt)

            # append real hittime to array:
            hittime_array.append(hittime_real)

        """ analyze signal: """
        # build histogram, where hittimes are saved:
        # set bin-edges of hittime histogram in ns:
        bins_hittime = np.arange(min_t, max_t + 2 * bin_width, bin_width)
        # build hittime histogram:
        npe_per_hittime, bin_edges_hittime = np.histogram(hittime_array, bins_hittime)

        # preallocate number of pe of this event:
        num_pe = 0

        # loop over values of the histogram bins:
        for index1 in range(len(npe_per_hittime)):
            # check if number of PE in this bin is above the threshold:
            if npe_per_hittime[index1] > threshold1:
                # calculate number of PEs in this signal peak:
                # preallocate sum of PE per bin in the signal peak:
                sum_pe_peak = 0

                # add nPE of npe_per_hittime[index1] to sum_pe_peak:
                sum_pe_peak = sum_pe_peak + npe_per_hittime[index1]

                # check previous bins:
                for num1 in range(1, 100, 1):
                    # check if value in previous bins_hittime[index1 - num1] is above threshold2:
                    if npe_per_hittime[index1 - num1] > threshold2:
                        # add nPE of this bin to sum_pe_peak:
                        sum_pe_peak = sum_pe_peak + npe_per_hittime[index1 - num1]

                    else:
                        break

                # check following bins:
                for num2 in range(1, 1000, 1):
                    # check if (index1 + num2 + 2) is in the range of npe_per_time array:
                    if (index1 + num2 + 2) >= len(npe_per_hittime):
                        print("Warning: iteration reaches last index of npe_per_time in evtID = {0:d}".format(evt))
                        break

                    # check if value bins_hittime[index1 + num2] is <= threshold2 and if following
                    # bins_hittime[index1 + num2 + 1] and bins_hittime[index1 + num2 + 2] are also <= threshold2:
                    if (npe_per_hittime[index1 + num2] <= threshold2 and
                            npe_per_hittime[index1 + num2 + 1] <= threshold2 and
                            npe_per_hittime[index1 + num2 + 2] <= threshold2):
                        break
                    else:
                        # npe_per_hittime[index1 + num2] or npe_per_hittime[index1 + num2 + 1] are above threshold2:
                        # add nPE of this bin to sum_pe_peak:
                        sum_pe_peak = sum_pe_peak + npe_per_hittime[index1 + num2]

                # print("number of pe in delayed signal = {0:d}".format(sum_pe_peak))

                # set number of pe of delayed signal:
                num_pe = sum_pe_peak
                break
            else:
                # go to next bin:
                continue

        # append the number of pe of one event to the array:
        number_pe = np.append(number_pe, num_pe)

    return number_pe, number_analyzed


def read_sample_detsim_user(rootfile_input, r_cut, e_prompt_min, e_prompt_max, e_delayed_min, e_delayed_max,
                            time_cut_min, time_cut_max, distance_cut, time_resolution, number_entries_input):
    """
    function to read the sample_detsim_user.root file and to get visible energy of the prompt signal of
    the IBD-like signals.

    IMPORTANT: the cuts are made with the parameters from 'geninfo' and 'prmtrkdep' tree of user-output!!!

    :param rootfile_input: input root file from tut_detsim.py: sample_detsim_user.root
    :param r_cut: specifies fiducial volume cut, radius is mm, normally r < 17 m = 17000 mm
    :param e_prompt_min: minimal prompt energy from energy cut in MeV, normally e_prompt_min = 10 MeV
    :param e_prompt_max: maximal prompt energy from energy cut in MeV, normally e_prompt_max = 105 MeV
    :param e_delayed_min: minimal delayed energy from energy cut in MeV, normally e_delayed_min = 1.9 MeV
    :param e_delayed_max: maximal delayed energy from energy cut in MeV, normally e_delayed_max = 2.5 MeV
    :param time_cut_min: minimal time difference delayed to prompt in ns, normally time_cut_min = 600 ns
    :param time_cut_max: maximal time difference delayed to prompt in ns, normally time_cut_max = 1 ms = 1 000 000 ns
    :param distance_cut: distance cut between prompt and delayed signal in mm, normally distance_cut < 1.5 m = 1500 mm
    :param time_resolution: time in ns, where two prompt signals can not be separated anymore,
                            normally time_resolution =
    :param number_entries_input: number of entries, that the input files should have (integer), normally = 100

    :return:
    """
    # load the ROOT file:
    rfile = ROOT.TFile(rootfile_input)
    # get the "evt"-TTree from the TFile:
    # rtree_evt = rfile.Get("evt")
    # get the "geninfo"-TTree from the TFile:
    rtree_geninfo = rfile.Get("geninfo")
    # get the "prmtrkdep"-TTree from the TFile:
    rtree_prmtrkdep = rfile.Get("prmtrkdep")

    # get the number of events in the geninfo Tree:
    number_events_geninfo = rtree_geninfo.GetEntries()
    # get the number of events in the prmtrkdep Tree:
    number_events_prmtrkdep = rtree_prmtrkdep.GetEntries()
    if number_events_geninfo == number_events_prmtrkdep:
        number_events = number_events_geninfo
    else:
        # number_events = 0
        # print("ERROR: number of events in the Trees are NOT equal!!")
        sys.exit("ERROR: number of events in t Trees are NOT equal!!")

    # check if number_events is equal to number_entries_input (if not, the detector simulation was incorrect!!):
    if number_events != number_entries_input:
        # number_events = 0
        # print("ERROR: number of events are not equal to {0:d}".format(number_entries_input))
        # print("-> Detector Simulation not correct!!")
        sys.exit("ERROR: number of events are not equal to {0:d} -> Detector Simulation not correct!"
                 .format(number_entries_input))

    # preallocate array of visible energy of prompt signal (energy in MeV) (np.array of float):
    e_vis = np.array([])
    # preallocate array, where event ID of the IBD like events is saved (np.array of float):
    evt_id_ibd = np.array([])

    # preallocate number of events, where there is 1 possible prompt and 1 possible delayed (case0):
    number_case0 = 0
    # preallocate number of events, where there is 1 possible prompt and 1 possible delayed (case0), which make a
    # IBD-like signal:
    number_case0_1ibdlike = 0

    # preallocate number of events in case1:
    number_case1 = 0
    # preallocate number of events, where there are more than 1 possible prompt and 1 possible delayed (case1),
    # but only possible 1 IBD-like signal:
    number_case1_1posibdlike = 0
    # preallocate number of events, where there are more than 1 possible prompt and 1 possible delayed (case1),
    # but 2 possible IBD-like signals:
    number_case1_2posibdlike = 0
    # preallocate number of events in case 1 with 2 possible IBD-like signals, which are NOT added to e_vis:
    number_case1_2posibdlike_notadded = 0
    # preallocate number of events in case 1 with 2 possible IBD-like signals, which are added to e_vis:
    number_case1_2posibdlike_added = 0
    # preallocate number of events in case 1 with 3 or more possible IBD-like events:
    number_case1_moreposibdlike = 0

    # preallocate number of events, where there is 1 possible prompt but more than 1 possible delayed (case2):
    number_case2 = 0
    # preallocate number of events in case2 with 1 possible IBD-like signal:
    number_case2_1posibdlike = 0

    # preallocate number of events, where there is more than 1 possible prompt and more than 1 possible delayed (case3):
    number_case3 = 0
    # preallocate number of events in case3 with 0 possible IBD-like signals:
    number_case3_noibdlike = 0
    # preallocate number of events in case 3 with 1 possible IBD-like signal:
    number_case3_1ibdlike = 0
    # preallocate number of events in case3 with 2 possible IBD-like signals:
    number_case3_2ibdlike = 0
    # preallocate number of events in case3 with more than 2 possible IBD-like signals:
    number_case3_moreibdlike = 0

    # preallocate number of events, where there are 2 prompt and only 1 corresponding delayed signal:
    number_check1 = 0
    # preallocate number of events, where there are 2 possible IBD signals:
    number_check2 = 0

    # loop over every event, i.e. every entry, in the TTree:
    for event in range(number_events):

        """ preallocate arrays: """
        # PDG ID of initial particles of geninfo tree of each particle in the event:
        pdgid_init_geninfo = np.array([])
        # initial position in x-direction in millimeter of each particle in the event:
        x_init = np.array([])
        # initial position in y-direction in millimeter of each particle in the event:
        y_init = np.array([])
        # initial position in z-direction in millimeter of each particle in the event:
        z_init = np.array([])
        # initial time in nanoseconds of each particle in the event:
        time_init = np.array([])
        # exit or stopping position in x-direction in millimeter of each particle in the event:
        x_exit = np.array([])
        # exit or stopping position in y-direction in millimeter of each particle in the event:
        y_exit = np.array([])
        # exit or stopping position in z-direction in millimeter of each particle in the event:
        z_exit = np.array([])
        # exit or stopping time in nanoseconds of each particle in the event:
        time_exit = np.array([])
        # deposited energy of each particle in the event in MeV:
        e_dep = np.array([])
        # visible energy (quenched deposited energy) of each particle in the event in MeV:
        e_qdep = np.array([])

        # PDG ID of each particle in the event:
        pdgid = np.array([])

        """ first read the "geninfo" Tree"""
        # get the current event in the TTree:
        rtree_geninfo.GetEntry(event)

        # get the value of the event ID:
        evt_id_geninfo = int(rtree_geninfo.GetBranch('evtID').GetLeaf('evtID').GetValue())

        # get the value of the number of initial particles:
        n_par_geninfo = int(rtree_geninfo.GetBranch('nInitParticles').GetLeaf('nInitParticles').GetValue())

        # loop over the number of particles to get information about every particle in the event:
        for index in range(n_par_geninfo):

            # get the value of the initial PDG ID:
            init_pdgid_geninfo = int(rtree_geninfo.GetBranch('InitPDGID').GetLeaf('InitPDGID').GetValue(index))
            pdgid_init_geninfo = np.append(pdgid_init_geninfo, init_pdgid_geninfo)

            # get initial x position:
            init_x = rtree_geninfo.GetBranch('InitX').GetLeaf('InitX').GetValue(index)
            x_init = np.append(x_init, init_x)

            # get initial y position:
            init_y = rtree_geninfo.GetBranch('InitY').GetLeaf('InitY').GetValue(index)
            y_init = np.append(y_init, init_y)

            # get initial z position:
            init_z = rtree_geninfo.GetBranch('InitZ').GetLeaf('InitZ').GetValue(index)
            z_init = np.append(z_init, init_z)

            # get initial time:
            init_time = rtree_geninfo.GetBranch('InitTime').GetLeaf('InitTime').GetValue(index)
            time_init = np.append(time_init, init_time)

            # get exit/stopping x-position:
            exit_x = rtree_geninfo.GetBranch('ExitX').GetLeaf('ExitX').GetValue(index)
            x_exit = np.append(x_exit, exit_x)

            # get exit/stopping y-position:
            exit_y = rtree_geninfo.GetBranch('ExitY').GetLeaf('ExitY').GetValue(index)
            y_exit = np.append(y_exit, exit_y)

            # get exit/stopping z-position:
            exit_z = rtree_geninfo.GetBranch('ExitZ').GetLeaf('ExitZ').GetValue(index)
            z_exit = np.append(z_exit, exit_z)

            # get the exit/stopping time:
            exit_time = rtree_geninfo.GetBranch('ExitT').GetLeaf('ExitT').GetValue(index)
            time_exit = np.append(time_exit, exit_time)

        """ then read the "prmtrkdep" Tree"""
        # get the current event in the TTree:
        rtree_prmtrkdep.GetEntry(event)

        # get the value of the event ID:
        evt_id_prmtrkdep = int(rtree_prmtrkdep.GetBranch('evtID').GetLeaf('evtID').GetValue())

        # get the value of the number of initial particles:
        n_par_prmtrkdep = int(rtree_prmtrkdep.GetBranch('nInitParticles').GetLeaf('nInitParticles').GetValue())

        # check event ID of the Trees:
        if evt_id_prmtrkdep == evt_id_geninfo:
            evt_id = evt_id_geninfo
        else:
            # evt_id = 0
            # print("ERROR: event ID in the Trees are NOT equal!!")
            sys.exit("ERROR: event ID in the Trees are NOT equal!")

        # check number of initial particles of the Trees:
        if n_par_prmtrkdep == n_par_geninfo:
            n_par = n_par_geninfo
        else:
            # n_par = 0
            # print("ERROR: number of initial particles in the Trees are NOT equal!!")
            sys.exit("ERROR: number of initial particles in the Trees are NOT equal!")

        # loop over the number of particles to get information about every particle in the event:
        for index in range(n_par):

            # get the value of the PDG ID:
            pdgid_prmtrkdep = int(rtree_prmtrkdep.GetBranch('PDGID').GetLeaf('PDGID').GetValue(index))
            # check PDG ID of the Trees:
            if pdgid_prmtrkdep == pdgid_init_geninfo[index]:
                pdgid = np.append(pdgid, pdgid_prmtrkdep)
            else:
                pdgid = np.append(pdgid, 0)
                print("ERROR: PDG ID in the Trees are NOT equal!!")

            # get deposited energy:
            dep_e = rtree_prmtrkdep.GetBranch('edep').GetLeaf('edep').GetValue(index)
            e_dep = np.append(e_dep, dep_e)

            # get visible energy:
            qdep_e = rtree_prmtrkdep.GetBranch('Qedep').GetLeaf('Qedep').GetValue(index)
            e_qdep = np.append(e_qdep, qdep_e)


        """ Does the event mimic an IBD signal? """
        # preallocate flag (array of boolean):
        is_prompt_signal = np.array([])
        is_delayed_signal = np.array([])

        # set flags:
        for index in range(n_par):
            # calculate the distance of the particle to the center of the event:
            r_init = np.sqrt(x_init[index]**2 + y_init[index]**2 + z_init[index]**2)
            r_exit = np.sqrt(x_exit[index]**2 + y_exit[index]**2 + z_exit[index]**2)

            # set is_prompt_signal flag (criteria: 10 MeV <= edep <= 105 MeV AND r_init < 17m AND r_exit < 17m):
            if e_prompt_min <= e_dep[index] <= e_prompt_max and r_init < r_cut and r_exit < r_cut:
                is_prompt_signal = np.append(is_prompt_signal, True)
            else:
                is_prompt_signal = np.append(is_prompt_signal, False)

            # set is_delayed_signal flag (criteria: 1.9 MeV <= edep <= 2.5 MeV AND r_init < 17m AND r_exit < 17m):
            if e_delayed_min <= e_dep[index] <= e_delayed_max and r_init < r_cut and r_exit < r_cut:
                is_delayed_signal = np.append(is_delayed_signal, True)
            else:
                is_delayed_signal = np.append(is_delayed_signal, False)

        # check if there are prompt and delayed signals in the event
        check_prompt = np.count_nonzero(is_prompt_signal)
        check_delayed = np.count_nonzero(is_delayed_signal)
        if check_prompt == 0 or check_delayed == 0:
            # no prompt signal OR no delayed signal in the event -> go to the next event
            continue

        # get the index/particle of the event that can be prompt or delayed signal, respectively (np.array):
        index_prompt = np.where(is_prompt_signal)[0]
        index_delayed = np.where(is_delayed_signal)[0]


        if len(index_prompt) == 1 and len(index_delayed) == 1:
            """ 1 possible prompt AND 1 possible delayed signal (case0): """
            # only one prompt and one delayed signal. Get first entry of the array:
            index_p = index_prompt[0]
            index_d = index_delayed[0]

            # increment number_case0:
            number_case0 = number_case0 + 1

            # check if initial time of possible prompt and delayed signals in the event is 0:
            if time_init[index_p] == 0 and time_init[index_d] == 0:

                # calculate the time difference delta_t between delayed and prompt signal:
                delta_t = time_exit[index_d] - time_exit[index_p]

                # time cut criteria: 600 ns <= delta_t <= 1.0 ms (1.0 ms = 1000000 ns):
                if time_cut_min < delta_t < time_cut_max:

                    # calculate distance from prompt to delayed:
                    distance_p_d = np.sqrt((x_exit[index_p] - x_exit[index_d])**2 +
                                           (y_exit[index_p] - y_exit[index_d])**2 +
                                           (z_exit[index_p] - z_exit[index_d])**2)

                    # prompt - delayed distance cut: R_prompt_delayed < 1.5 m (1.5 m = 1500 mm)
                    if distance_p_d < distance_cut:

                        # increment number_case0_1ibdlike:
                        number_case0_1ibdlike = number_case0_1ibdlike + 1

                        # append evt_id of the IBD like signal to evt_id_ibd array:
                        evt_id_ibd = np.append(evt_id_ibd, evt_id)

                        # append Qedep of the prompt signal to the e_vis array:
                        e_vis = np.append(e_vis, e_qdep[index_p])
                    else:
                        continue
                else:
                    continue
            else:
                print("WARNING: initial time is not 0 for possible prompt and delayed signal in event {0:d}"
                      .format(evt_id))


        elif len(index_prompt) > 1 and len(index_delayed) == 1:
            """ more than 1 possible prompt signals, BUT only 1 possible delayed signal (case1): """
            # preallocate value that represents the number of IBD-like signals in this event:
            number_ibd_evts = 0

            # check:
            number_case1 = number_case1 + 1

            # preallocate array, where the visible energy of the prompt signals of the IBD-like signal is stored
            # (energy in MeV):
            array_e_vis = np.array([])

            # preallocate array, where the index of the prompt signal of the IBD-like signal is stored:
            array_prompt_index = np.array([])

            # loop over the possible prompt signals in the event:
            for index in range(len(index_prompt)):
                # get the index of the prompt signal in the array:
                index_p = index_prompt[index]
                # only one possible delayed signal in the event:
                index_d = index_delayed[0]

                # check if initial time of possible prompt and delayed signals in the event is 0:
                if time_init[index_p] == 0 and time_init[index_d] == 0:

                    # calculate the time difference delta_t between delayed and prompt signal:
                    delta_t = time_exit[index_d] - time_exit[index_p]

                    # time cut criteria: 600 ns <= delta_t <= 1.0 ms (1.0 ms = 1000000 ns):
                    if time_cut_min < delta_t < time_cut_max:

                        # calculate distance from prompt to delayed:
                        distance_p_d = np.sqrt((x_exit[index_p] - x_exit[index_d])**2 +
                                               (y_exit[index_p] - y_exit[index_d])**2 +
                                               (z_exit[index_p] - z_exit[index_d])**2)

                        # prompt - delayed distance cut: R_prompt_delayed < 1.5 m (1.5 m = 1500 mm)
                        if distance_p_d < distance_cut:

                            # increment number of IBD-like signals in this event:
                            number_ibd_evts = number_ibd_evts + 1

                            # append visible energy of the prompt signal to the array (energy in MeV):
                            array_e_vis = np.append(array_e_vis, e_qdep[index_p])

                            # append index to the array, where index of prompt signal is stored:
                            array_prompt_index = np.append(array_prompt_index, index_p)

                        else:
                            continue
                    else:
                        continue
                else:
                    print("WARNING: initial time is not 0 for possible prompt and delayed signal in event {0:d}"
                          .format(evt_id))

            # if there is 1 IBD-like signal in the event, store information in array. If there is NO IBD-like signal
            # go to the next event. Else: not print warning.
            if number_ibd_evts == 1:
                # check, if len(array_e_vis) is also equal to 1:
                if number_ibd_evts != len(array_e_vis):
                    # print("ERROR: Number of IBD-like events different to length of 'array_e_vis' (evt_ID = {0:d})"
                    #       .format(evt_id))
                    sys.exit("ERROR: Number of IBD-like events different to length of 'array_e_vis' (evt_ID = {0:d})"
                             .format(evt_id))

                # check:
                number_case1_1posibdlike = number_case1_1posibdlike + 1

                # append evt_id of the IBD like signal to evt_id_ibd array:
                evt_id_ibd = np.append(evt_id_ibd, evt_id)

                # append Qedep of the prompt signal to the e_vis array:
                e_vis = np.append(e_vis, array_e_vis[0])

            elif number_ibd_evts == 0:
                continue

            elif number_ibd_evts == 2:
                # 2 possible prompt signal and only 1 corresponding delayed signals:
                # calculate the time difference between the prompt signals in ns:
                # TODO-me: How is time_exit defined exactly?
                delta_t_pp = np.absolute(time_exit[int(array_prompt_index[0])] - time_exit[int(array_prompt_index[1])])

                # check:
                number_case1_2posibdlike = number_case1_2posibdlike + 1

                if time_resolution <= delta_t_pp <= time_cut_max:
                    # this means, one prompt signal lies between the other prompt and the delayed signal and can be
                    # separated from the other prompt signal. -> NO IBD-like signal -> go to next event:

                    # check:
                    number_case1_2posibdlike_notadded = number_case1_2posibdlike_notadded + 1
                    continue

                elif delta_t_pp < time_resolution:
                    # this means, that the two prompt signals can not be separated clearly -> treat the two prompt
                    # signals as one prompt signal with energy = E_prompt1 + E_prompt2:

                    # append evt_id of the IBD like signal to evt_id_ibd array:
                    evt_id_ibd = np.append(evt_id_ibd, evt_id)

                    # append sum of Qedep of first prompt signal and Qedep of second prompt signal to e_vis array:
                    e_vis = np.append(e_vis, array_e_vis[0] + array_e_vis[1])

                    # check:
                    number_case1_2posibdlike_added = number_case1_2posibdlike_added + 1
                    # print(evt_id)

                else:
                    # this means, one prompt signal is after the delayed signal (should be not possible):
                    # print("WARNING in len(index_prompt) > 1 and len(index_delayed) == 1 ------ 1 prompt
                    # after delayed signal")
                    sys.exit("WARNING in len(index_prompt) > 1 and len(index_delayed) == 1 ------ 1 prompt after "
                             "delayed signal")

            else:
                # TODO-me: How to deal with events where there are 3 or more IBD-like events (but only one delayed sig.)
                # to estimate the number of such events, increment variable number_case1_moreposibdlike:
                number_case1_moreposibdlike = number_case1_moreposibdlike + 1
                print("WARNING: more than 1 IBD-like signal in event {2:d}: n_IBD_like = {0:.0f}, prompt = {1:.0f}"
                      .format(number_ibd_evts, check_prompt, evt_id))
                print("-----------> not yet included!")


        elif len(index_prompt) == 1 and len(index_delayed) > 1:
            """ 1 possible prompt signal, BUT more than 1 possible delayed signals (case2): """
            # preallocate value that represents the number of IBD-like signals in this event:
            number_ibd_evts = 0

            # check:
            number_case2 = number_case2 + 1

            # loop over the possible delayed signals in the event:
            for index in range(len(index_delayed)):
                # only one index of the prompt signal in the event:
                index_p = index_prompt[0]
                # get index of possible delayed signal in the array:
                index_d = index_delayed[index]

                # check if initial time of possible prompt and delayed signals in the event is 0:
                if time_init[index_p] == 0 and time_init[index_d] == 0:

                    # calculate the time difference delta_t between delayed and prompt signal:
                    delta_t = time_exit[index_d] - time_exit[index_p]

                    # time cut criteria: 600 ns <= delta_t <= 1.0 ms (1.0 ms = 1000000 ns):
                    if time_cut_min < delta_t < time_cut_max:

                        # calculate distance from prompt to delayed:
                        distance_p_d = np.sqrt((x_exit[index_p] - x_exit[index_d])**2 +
                                               (y_exit[index_p] - y_exit[index_d])**2 +
                                               (z_exit[index_p] - z_exit[index_d])**2)

                        # prompt - delayed distance cut: R_prompt_delayed < 1.5 m (1.5 m = 1500 mm):
                        if distance_p_d < distance_cut:

                            # increment number of IBD-like signals in this event:
                            number_ibd_evts = number_ibd_evts + 1

                        else:
                            continue
                    else:
                        continue
                else:
                    print("WARNING: initial time is not 0 for possible prompt and delayed signal in event {0:d}"
                          .format(evt_id))

            # if there is no or more than 1 IBD-like signal in the event, go to next event (neutron multiplicity cut!).
            # If there is only one IBD-like event, store the visible energy.
            if number_ibd_evts == 1:

                # check:
                number_case2_1posibdlike = number_case2_1posibdlike + 1

                # append evt_id of the IBD like signal to evt_id_ibd array:
                evt_id_ibd = np.append(evt_id_ibd, evt_id)

                # append Qedep of the prompt signal to the e_vis array:
                e_vis = np.append(e_vis, e_qdep[index_prompt[0]])

            else:
                continue


        else:
            """ More than 1 possible prompt signal AND more than 1 possible delayed signal 
                (len(index_prompt) > 1 and len(index_delayed) > 1) """
            # preallocate array that represents the number of IBD-like signals in this event
            # (array of length index_prompt):
            array_number_ibd = np.zeros(len(index_prompt))

            # preallocate array, where the indices of the event of the "real" delayed signals (which correspond to
            # this prompt signal) are stored:
            array_delayed_index = np.array([])

            # check:
            number_case3 = number_case3 + 1

            # loop over possible prompt signals in this event:
            for index1 in range(len(index_prompt)):

                # preallocate value that represents the number of IBD-like signals in this event
                # (this value is then included in array_number_ibd):
                number_ibd_prompt = 0

                # get the index of the prompt signal in the array:
                index_p = index_prompt[index1]

                # loop over the possible delayed signals in the event:
                for index2 in range(len(index_delayed)):

                    # get index of possible delayed signal in the event:
                    index_d = index_delayed[index2]

                    # check if initial time of possible prompt and delayed signals in the event is 0:
                    if time_init[index_p] == 0 and time_init[index_d] == 0:

                        # calculate the time difference delta_t between delayed and prompt signal:
                        delta_t = time_exit[index_d] - time_exit[index_p]

                        # time cut criteria: 600 ns <= delta_t <= 1.0 ms (1.0 ms = 1000000 ns):
                        if time_cut_min < delta_t < time_cut_max:

                            # calculate distance from prompt to delayed:
                            distance_p_d = np.sqrt((x_exit[index_p] - x_exit[index_d])**2 +
                                                   (y_exit[index_p] - y_exit[index_d])**2 +
                                                   (z_exit[index_p] - z_exit[index_d])**2)

                            # prompt - delayed distance cut: R_prompt_delayed < 1.5 m (1.5 m = 1500 mm)
                            if distance_p_d < distance_cut:

                                # increment number of IBD-like signals for this possible prompt signal in this event:
                                number_ibd_prompt = number_ibd_prompt + 1

                                # store the index of the delayed signal in the event:
                                index_real_delayed = index_d

                            else:
                                continue
                        else:
                            continue
                    else:
                        print("WARNING: initial time is not 0 for possible prompt and delayed signal in event {0:d}"
                              .format(evt_id))

                # if there is no or more than 1 IBD-like signal in the event for this prompt signal, go to next event
                # (neutron multiplicity cut!).
                # If there is only one IBD-like event, store the visible energy.
                if number_ibd_prompt == 1:

                    # include number of IBD-like signals for this prompt event to the array:
                    array_number_ibd[index1] = number_ibd_prompt

                    # store index of the delayed signal to the array:
                    array_delayed_index = np.append(array_delayed_index, index_real_delayed)

                else:
                    # 0 or more than one IBD-like signals for this prompt signal in this event
                    # (-> neutron multiplicity cut)
                    # -> go to the next prompt signal
                    continue

            # check array_number_ibd and array_delayed_index:
            if len(array_delayed_index) == 0:
                # this means, that there is no IBD-like signal in this event (either 0 or more than 1 delayed signal to
                # 1 prompt signal)

                # check:
                number_case3_noibdlike = number_case3_noibdlike + 1
                continue

            elif len(array_delayed_index) == 1:
                # this means, that there is 1 prompt signal and 1 corresponding delayed signal
                # -> therefore 1 IBD-like signal:

                # get the index of the prompt signal in the event:
                which_index = np.where(array_number_ibd == 1)[0]
                # check if there is just one entry equal to 1 in the array_delayed_index:
                if len(which_index) == 1:

                    # append evt_id of the IBD like signal to evt_id_ibd array:
                    evt_id_ibd = np.append(evt_id_ibd, evt_id)

                    correct_prompt_index = which_index[0]

                    # append Qedep of this prompt signal to the e_vis array:
                    e_vis = np.append(e_vis, e_qdep[index_prompt[correct_prompt_index]])

                    # check:
                    number_case3_1ibdlike = number_case3_1ibdlike + 1

                else:
                    print("ERROR in line 528")
                    continue

            elif len(array_delayed_index) == 2:
                # this means, that there are 2 possible IBD-like signals (2 prompt and 2 corresponding delayed signals)
                # there are 2 possibilities:

                # check:
                number_case3_2ibdlike = number_case3_2ibdlike + 1

                # 1. the 'two' delayed signals could be the same. Then: 2 prompt signal and 1 corresponding delayed
                # signal.
                # 2. two different delayed signals. Therefore 2 possible IBD-like signals.
                if array_delayed_index[0] == array_delayed_index[1]:
                    # possibility 1:
                    # to estimate the number of such events, increment variable number_check1:
                    number_check1 = number_check1 + 1
                else:
                    # possibility 2:
                    # INFO-me: by analyzing detsim files user_atmoNC_0.root to user_atmoNC_599.root (60000 evts) in
                    # INFO-me: '/local/scratch1/pipc51/astro/blum/detsim_output_data', there are 0 events like this
                    # INFO-me: -> number_check2 = 0 for all 60000 events!!!
                    # to estimate the number of such events, increment variable number_check2:
                    number_check2 = number_check2 + 1

            else:
                # TODO-me: How to deal with events where there are 2 or more IBD-like events (but only one delayed sig.)
                # to estimate the number of such events, increment variable number_case3_moreibdlike:
                number_case3_moreibdlike = number_case3_moreibdlike + 1
                print("WARNING: More than 1 possible prompt signal to only 1 possible delayed signal: evt ID = {0:d}"
                      .format(evt_id))
                print("----------> not yet included!!!!")

    return (number_events, evt_id_ibd, e_vis, number_case0, number_case0_1ibdlike, number_case1,
            number_case1_1posibdlike, number_case1_2posibdlike, number_case1_2posibdlike_added,
            number_case1_2posibdlike_notadded, number_case1_moreposibdlike, number_case2, number_case2_1posibdlike,
            number_case3, number_case3_noibdlike, number_case3_1ibdlike, number_case3_2ibdlike, number_check1,
            number_check2, number_case3_moreibdlike)


def preselect_sample_detsim_user(rootfile_input, r_cut, min_prompt_energy, max_prompt_energy, time_cut_min,
                                 time_cut_max, distance_cut, number_entries_input):
    """
    function to read the sample_detsim_user.root file and do a preselection of possible IBD-like signals

    :param rootfile_input: input root file from tut_detsim.py: sample_detsim_user.root
    :param r_cut: specifies fiducial volume cut, radius is mm, normally r < 16 m = 16000 mm
    :param min_prompt_energy:   sets the minimal energy of the prompt signal in MeV -> compare this with total deposit
                                energy of the event
    :param max_prompt_energy:   maximum of total deposit energy of one event, does not contribute to cuts, just for
                                information (in MeV)
    :param time_cut_min: minimal time difference delayed to prompt in ns
    :param time_cut_max: maximal time difference delayed to prompt in ns, normally time_cut_max = 1 ms = 1 000 000 ns
    :param distance_cut: specifies distance cut between prompt and delayed signal, normally distance < 1.5 m = 1500 mm
    :param number_entries_input: number of entries, that the input files should have (integer), normally = 1000

    :return:
    """
    # load the ROOT file:
    rfile = ROOT.TFile(rootfile_input)
    # get the "evt"-TTree from the TFile:
    rtree_evt = rfile.Get("evt")
    # get the "geninfo"-TTree from the TFile:
    rtree_geninfo = rfile.Get("geninfo")
    # get the "prmtrkdep"-TTree from the TFile:
    rtree_prmtrkdep = rfile.Get("prmtrkdep")
    # get the "nCapture"-TTree from TFile:
    rtree_ncapture = rfile.Get("nCapture")
    # get the "secondaries"-TTree from TFile:
    # rtree_secondaries = rfile.Get("secondaries")

    # get number of events in evt tree:
    number_events_evt = rtree_evt.GetEntries()
    # get the number of events in the geninfo Tree:
    number_events_geninfo = rtree_geninfo.GetEntries()
    # get the number of events in the prmtrkdep tree:
    number_events_prmtrkdep = rtree_prmtrkdep.GetEntries()
    # get number of events in nCapture tree:
    number_events_ncapture = rtree_ncapture.GetEntries()
    # get number of events in secondaries tree:
    # number_events_secondaries = rtree_secondaries.GetEntries()
    number_events_secondaries = 100

    if (number_events_evt == number_events_geninfo and number_events_evt == number_events_ncapture and
            number_events_evt == number_events_secondaries and number_events_evt == number_events_prmtrkdep):
        number_events = number_events_geninfo
    else:
        sys.exit("ERROR: number of events in the Trees are NOT equal!!")

    # check if number_events is equal to number_entries_input (if not, the detector simulation was incorrect!!):
    if number_events != number_entries_input:
        sys.exit("ERROR: number of events {0:d} are not equal to {1:d} -> Detector Simulation not correct!"
                 .format(number_events, number_entries_input))

    # preallocate array, where event ID of events are saved, that pass the preselection (np.array of float):
    evt_id_preselected = np.array([])
    # preallocate array, where the total deposit energy of the event is saved in MeV (np.array of float):
    edep_total = np.array([])
    # preallocate array, where reconstructed x-position of the event is saved in mm:
    x_reco_array = np.array([])
    # preallocate array, where reconstructed y-position of the event is saved in mm:
    y_reco_array = np.array([])
    # preallocate array, where reconstructed z-position of the event is saved in mm:
    z_reco_array = np.array([])

    # preallocate number of events, that pass the preselection (volume cut and neutron-multiplicity cut, time cut and
    # distance cut):
    number_preselected = 0
    # preallocate number of events, which are rejected by preselection criteria:
    number_rejected = 0
    # preallocate number of events with position reconstruction that pass volume cut:
    number_vol_pass = 0
    # preallocate number of events, where initial position (without position reconstruction) pass volume cut:
    number_vol_pass_initial = 0
    # preallocate number of events that are rejected by volume cut:
    number_vol_reject = 0
    # preallocate number of events that pass the energy cut:
    number_e_pass = 0
    # preallocate number of events that are rejected of minimal energy cut:
    number_mine_reject = 0
    # number of events that are rejected of maximum energy cut:
    number_maxe_reject = 0
    # preallocate number of events without nCapture:
    number_without_ncap = 0
    # preallocate, number of events, that pass neutron-multiplicity cut:
    number_nmult_pass = 0
    # preallocate, number of events, that are rejected by neutron-multiplicity cut:
    number_nmult_reject = 0
    # preallocate, number of events, that pass time cut:
    number_time_pass = 0
    # preallocate, number of events, that are rejected by time cut:
    number_time_reject = 0
    # preallocate, number of events, that pass distance cut:
    number_dist_pass = 0
    # preallocate, number of events, that are rejected by distance cut:
    number_dist_reject = 0

    # loop over every event, i.e. every entry, in the TTree:
    for event in range(number_events):

        """ preallocate variables for geninfo tree: """
        # PDG ID of initial particles of geninfo tree of each particle in the event:
        pdgid_init_geninfo = np.array([])
        # initial position in x-direction in mm of initial particle in the event:
        x_init_geninfo = np.array([])
        # initial position in y-direction in mm of initial particle in the event:
        y_init_geninfo = np.array([])
        # initial position in z-direction in mm of initial particle in the event:
        z_init_geninfo = np.array([])
        # initial time in nanoseconds of initial particle in the event:
        time_init_geninfo = np.array([])
        # exit or stopping position in x-direction in mm of initial particle in the event:
        x_exit_geninfo = np.array([])
        # exit or stopping position in y-direction in mm of initial particle in the event:
        y_exit_geninfo = np.array([])
        # exit or stopping position in z-direction in mm of initial particle in the event:
        z_exit_geninfo = np.array([])

        """ preallocate variables for prmtrkdep tree: """
        # stop position, where initial particles deposit their energy:
        x_edep_prmtrkdep = np.array([])
        y_edep_prmtrkdep = np.array([])
        z_edep_prmtrkdep = np.array([])
        # quenched deposit energy of the primary particles in MeV (is used to apply the vertex smearing):
        qedep_prmtrkdep = 0

        """ preallocate variables for nCapture tree: """
        # time, when neutron is captured in ns:
        capture_time_ncap = np.array([])
        # start position of neutron capture in mm:
        x_start_ncap = np.array([])
        y_start_ncap = np.array([])
        z_start_ncap = np.array([])
        # stop position of neutron capture in mm:
        x_stop_ncap = np.array([])
        y_stop_ncap = np.array([])
        z_stop_ncap = np.array([])
        # kinetic energy in MeV of the gamma, that is produced by neutron capture:
        kine_gamma_ncap = np.array([])

        """ read 'evt' tree: """
        # get current event:
        rtree_evt.GetEntry(event)
        # get the value of the event ID:
        evt_id_evt = int(rtree_evt.GetBranch('evtID').GetLeaf('evtID').GetValue())
        # get total deposited energy of this event in MeV:
        edep_evt = float(rtree_evt.GetBranch('edep').GetLeaf('edep').GetValue())

        """ read the 'geninfo' Tree """
        # get the current event in the TTree:
        rtree_geninfo.GetEntry(event)
        # get the value of the event ID:
        evt_id_geninfo = int(rtree_geninfo.GetBranch('evtID').GetLeaf('evtID').GetValue())

        # check event ID of the Trees:
        if evt_id_geninfo == evt_id_evt:
            evt_id = evt_id_evt
        else:
            sys.exit("ERROR: event ID in the Trees are NOT equal!")

        # get the value of the number of initial particles:
        n_par_geninfo = int(rtree_geninfo.GetBranch('nInitParticles').GetLeaf('nInitParticles').GetValue())

        # loop over the number of particles to get information about every particle in the event:
        for index in range(n_par_geninfo):
            # get the value of the initial PDG ID:
            init_pdgid_geninfo = int(rtree_geninfo.GetBranch('InitPDGID').GetLeaf('InitPDGID').GetValue(index))
            pdgid_init_geninfo = np.append(pdgid_init_geninfo, init_pdgid_geninfo)

            # get initial x position:
            init_x = float(rtree_geninfo.GetBranch('InitX').GetLeaf('InitX').GetValue(index))
            x_init_geninfo = np.append(x_init_geninfo, init_x)

            # get initial y position:
            init_y = float(rtree_geninfo.GetBranch('InitY').GetLeaf('InitY').GetValue(index))
            y_init_geninfo = np.append(y_init_geninfo, init_y)

            # get initial z position:
            init_z = float(rtree_geninfo.GetBranch('InitZ').GetLeaf('InitZ').GetValue(index))
            z_init_geninfo = np.append(z_init_geninfo, init_z)

            # get initial time:
            init_time = float(rtree_geninfo.GetBranch('InitTime').GetLeaf('InitTime').GetValue(index))
            time_init_geninfo = np.append(time_init_geninfo, init_time)

            # get exit/stopping x-position:
            exit_x = float(rtree_geninfo.GetBranch('ExitX').GetLeaf('ExitX').GetValue(index))
            x_exit_geninfo = np.append(x_exit_geninfo, exit_x)

            # get exit/stopping y-position:
            exit_y = float(rtree_geninfo.GetBranch('ExitY').GetLeaf('ExitY').GetValue(index))
            y_exit_geninfo = np.append(y_exit_geninfo, exit_y)

            # get exit/stopping z-position:
            exit_z = float(rtree_geninfo.GetBranch('ExitZ').GetLeaf('ExitZ').GetValue(index))
            z_exit_geninfo = np.append(z_exit_geninfo, exit_z)

        """ then read the "prmtrkdep" Tree"""
        # get the current event in the TTree:
        rtree_prmtrkdep.GetEntry(event)
        # get the value of the event ID:
        evt_id_prmtrkdep = int(rtree_prmtrkdep.GetBranch('evtID').GetLeaf('evtID').GetValue())

        # check event ID of the Trees:
        if evt_id_prmtrkdep == evt_id_evt:
            evt_id = evt_id_evt
        else:
            sys.exit("ERROR: event ID in the Trees are NOT equal!")

        # get the value of the number of initial particles:
        n_par_prmtrkdep = int(rtree_prmtrkdep.GetBranch('nInitParticles').GetLeaf('nInitParticles').GetValue())

        # check number of initial particles of the Trees:
        if n_par_prmtrkdep == n_par_geninfo:
            n_par = n_par_geninfo
        else:
            sys.exit("ERROR: number of initial particles in the Trees are NOT equal!")

        # loop over the number of particles to get information about every particle in the event:
        for index in range(n_par):
            # get the value of the PDG ID:
            pdgid_prmtrkdep = int(rtree_prmtrkdep.GetBranch('PDGID').GetLeaf('PDGID').GetValue(index))
            # check PDG ID of the Trees:
            if pdgid_prmtrkdep != pdgid_init_geninfo[index]:
                sys.exit("ERROR: PDG ID in the Trees are NOT equal!!")

            # get position deposited energy:
            x_edep = float(rtree_prmtrkdep.GetBranch('edepX').GetLeaf('edepX').GetValue(index))
            x_edep_prmtrkdep = np.append(x_edep_prmtrkdep, x_edep)
            y_edep = float(rtree_prmtrkdep.GetBranch('edepY').GetLeaf('edepY').GetValue(index))
            y_edep_prmtrkdep = np.append(y_edep_prmtrkdep, y_edep)
            z_edep = float(rtree_prmtrkdep.GetBranch('edepZ').GetLeaf('edepZ').GetValue(index))
            z_edep_prmtrkdep = np.append(z_edep_prmtrkdep, z_edep)

            # get quenched deposited energy in MeV:
            qedep = float(rtree_prmtrkdep.GetBranch('Qedep').GetLeaf('Qedep').GetValue(index))
            # add qedep to qedep_prmtrkdep to get the sum of Qedep of the primary particles:
            qedep_prmtrkdep += qedep

        """ read 'nCapture' tree: """
        # get current event:
        rtree_ncapture.GetEntry(event)
        # get value of eventID:
        evt_id_ncap = int(rtree_ncapture.GetBranch('evtID').GetLeaf('evtID').GetValue())

        # check event ID of the Trees:
        if evt_id_ncap == evt_id_evt:
            evt_id = evt_id_evt
        else:
            sys.exit("ERROR: event ID in the Trees are NOT equal!")

        # get number of neutron captures in the event:
        number_ncap = int(rtree_ncapture.GetBranch('NeutronN').GetLeaf('NeutronN').GetValue())

        # get number of particles that are produced by neutron capture in the whole event:
        n_ncap = int(rtree_ncapture.GetBranch('n').GetLeaf('n').GetValue())

        # loop over number of neutron captures:
        for index in range(number_ncap):

            # get track ID of neutron:
            trkid_neutron = int(rtree_ncapture.GetBranch('NeutronTrkid').GetLeaf('NeutronTrkid').GetValue(index))

            # is there a gamma between 2.2 and 2.25 MeV:
            number_gamma = 0

            # loop over particles that are produced by neutron capture:
            for index2 in range(n_ncap):
                # get PDGID of the particle:
                pdgid = int(rtree_ncapture.GetBranch('pdgid').GetLeaf('pdgid').GetValue(index2))
                # get Track id of the particle:
                trkid = int(rtree_ncapture.GetBranch('trkid').GetLeaf('trkid').GetValue(index2))
                # check if particle is gamma and that it corresponds to the neutron:
                if pdgid == 22 and trkid == trkid_neutron:
                    # get kinetic energy of gamma:
                    kine_gamma = float(rtree_ncapture.GetBranch('kine').GetLeaf('kine').GetValue(index2))

                    if 2.2 < kine_gamma < 2.25:
                        kine_gamma_ncap = np.append(kine_gamma_ncap, kine_gamma)
                        # increment number_gamma:
                        number_gamma += 1

                        # if one 2.2 MeV gamma was found, go to the next neutron capture. This will maybe overestimate
                        # the number of possible IBD-events, but the preselected events will be analyzed in detail
                        # later, so this is no problem.
                        break

                else:
                    continue

            # when there is no gamma between 2.2 and 2.25 MeV corresponding to the neutron capture, so if
            # number_gamma = 0, add 0.0 to kine_gamma_ncap array:
            if number_gamma == 0:
                kine_gamma_ncap = np.append(kine_gamma_ncap, 0.0)

            # time, when neutron is captured
            time_ncap = float(rtree_ncapture.GetBranch('NeutronCaptureT').GetLeaf('NeutronCaptureT').GetValue(index))
            capture_time_ncap = np.append(capture_time_ncap, time_ncap)
            # start position of neutron capture:
            x_start = float(rtree_ncapture.GetBranch('NCStartX').GetLeaf('NCStartX').GetValue(index))
            x_start_ncap = np.append(x_start_ncap, x_start)
            y_start = float(rtree_ncapture.GetBranch('NCStartY').GetLeaf('NCStartY').GetValue(index))
            y_start_ncap = np.append(y_start_ncap, y_start)
            z_start = float(rtree_ncapture.GetBranch('NCStartZ').GetLeaf('NCStartZ').GetValue(index))
            z_start_ncap = np.append(z_start_ncap, z_start)
            # stop position of neutron capture:
            x_stp = float(rtree_ncapture.GetBranch('NCStopX').GetLeaf('NCStopX').GetValue(index))
            x_stop_ncap = np.append(x_stop_ncap, x_stp)
            y_stp = float(rtree_ncapture.GetBranch('NCStopY').GetLeaf('NCStopY').GetValue(index))
            y_stop_ncap = np.append(y_stop_ncap, y_stp)
            z_stp = float(rtree_ncapture.GetBranch('NCStopZ').GetLeaf('NCStopZ').GetValue(index))
            z_stop_ncap = np.append(z_stop_ncap, z_stp)

        # set flag for different cuts of preselection:
        flag_vol_pass = True
        flag_vol_pass_initial = True
        flag_e_pass = False
        flag_nmult_pass = False
        flag_time_pass = False
        # flag_dist_pass: 0 means event is rejected, 1 means event passes distance cut,
        # 2 means 0 or more than 1 neutron capture (do not count these cases)
        flag_dist_pass = 2

        # check if there is a particle in the event:
        if len(pdgid_init_geninfo) == 0:
            # this means pdgid_init_geninfo array is empty -> no particle in the event.
            # therefore: go to next event
            continue

        """ check volume cut: """
        # Do vertex smearing of the position of init_geninfo particles:
        # set x-position of first initial particle:
        x_initial = x_init_geninfo[0]
        # set y-position of first initial particle:
        y_initial = y_init_geninfo[0]
        # set z-position of first initial particles:
        z_initial = z_init_geninfo[0]
        # loop over position of initial particles and check if the positions are equal for all initial particles:
        for index in range(n_par):
            if (x_initial != x_init_geninfo[index] or y_initial != y_init_geninfo[index]
                    or z_initial != z_init_geninfo[index]):
                # position of initial particles differ:
                sys.exit("position of initial particles differ!!!! (file = {0}, evtID = {1:d})"
                         .format(rootfile_input, event))

        # all initial particles have same initial position. Smear x,y and z position of this initial position with
        # function position_smearing(). (returns reconstructed position in mm):
        x_reconstructed = position_smearing(x_initial, qedep_prmtrkdep)
        y_reconstructed = position_smearing(y_initial, qedep_prmtrkdep)
        z_reconstructed = position_smearing(z_initial, qedep_prmtrkdep)

        # calculate the distance to detector center for the reconstructed position in mm:
        r_reconstructed = np.sqrt(x_reconstructed**2 + y_reconstructed**2 + z_reconstructed**2)
        if r_reconstructed >= r_cut:
            flag_vol_pass = False

        # as cross-check calculate the distance to detector center for the initial position in mm (not reconstructed):
        r_initial = np.sqrt(x_initial**2 + y_initial**2 + z_initial**2)
        if r_initial >= r_cut:
            flag_vol_pass_initial = False

        if flag_vol_pass:
            number_vol_pass += 1
        else:
            number_vol_reject += 1

        if flag_vol_pass_initial:
            # cross-check for initial position:
            number_vol_pass_initial += 1

        """ check prompt energy cut: """
        if min_prompt_energy <= edep_evt <= max_prompt_energy:
            flag_e_pass = True
            number_e_pass += 1
        elif edep_evt < min_prompt_energy:
            number_mine_reject += 1
        else:
            # edep_evt > max_prompt_energy:
            number_maxe_reject += 1

        """ check time cut: """
        # check if starting time of initial particles is zero:
        for index in range(n_par):
            if time_init_geninfo[index] == 0:
                time_init = 0
            else:
                sys.exit("ERROR: initial time is not zero in evtid = {0:d}".format(evt_id))

        # check if there are neutron captures in the event:
        if number_ncap == 0:
            # no nCaptures in the event:
            number_without_ncap += 1

        # preallocate number of neutron captures on hydrogen within 1 ms:
        number_ncapture_window = 0
        # preallocate array, where the indices of the neutron capture are saved:
        index_ncapture = np.array([])

        if len(capture_time_ncap) > len(kine_gamma_ncap):
            sys.exit("ERROR: len(capture_time_ncap > len(kine_gamma_ncap) in evtID = {0:d}".format(evt_id))
        elif len(capture_time_ncap) < len(kine_gamma_ncap):
            print("WARNING: len(capture_time_ncap) < len(kine_gamma_ncap)! At least 1 neutron capture, where two 2.2 "
                  "MeV gammas are produced!")
            print("evtID = {0:d}".format(evt_id))
            print("len(capture_time_ncap) = {0:d}".format(len(capture_time_ncap)))
            print("capture_time_ncap = {0}".format(capture_time_ncap))

        # calculate time difference between 0 and neutron capture:
        for index in range(number_ncap):
            if time_cut_min < capture_time_ncap[index] < time_cut_max and 2.2 < kine_gamma_ncap[index] < 2.25:
                # neutron is captured on Hydrogen within time window:
                flag_time_pass = True
                # increment number of n-captures within time window:
                number_ncapture_window += 1
                # save index of the neutron capture:
                index_ncapture = np.append(index_ncapture, index)
            else:
                continue

        if flag_time_pass:
            number_time_pass += 1
        else:
            number_time_reject += 1

        """ check neutron multiplicity cut: """
        # check number of neutron captures on hydrogen within time window:
        if number_ncapture_window == 1:
            flag_nmult_pass = True
            number_nmult_pass += 1
        else:
            number_nmult_reject += 1

        """ check distance cut: """
        # check distance only, when there is just 1 neutron capture on hydrogen in the time window:
        if number_ncapture_window == 1 and len(index_ncapture) == 1:
            # check distance: reconstructed position of initial particles -> neutron capture:
            r_dist_start_nc = np.sqrt((x_stop_ncap[int(index_ncapture[0])] - x_reconstructed)**2 +
                                      (y_stop_ncap[int(index_ncapture[0])] - y_reconstructed)**2 +
                                      (z_stop_ncap[int(index_ncapture[0])] - z_reconstructed)**2)

            if r_dist_start_nc >= distance_cut:
                # event is rejected
                flag_dist_pass = 0
            else:
                # event passes cut
                flag_dist_pass = 1
        else:
            # no or more than one neutron capture:
            flag_dist_pass = 2

        if flag_dist_pass == 1:
            number_dist_pass += 1
        elif flag_dist_pass == 0:
            number_dist_reject += 1

        """ check flags: """
        if flag_vol_pass and flag_e_pass and flag_nmult_pass and flag_time_pass and flag_dist_pass == 1:
            # event passes all cuts!
            # add evt ID to array:
            evt_id_preselected = np.append(evt_id_preselected, evt_id)
            # add total deposit energy to array:
            edep_total = np.append(edep_total, edep_evt)
            # add reconstructed x-position to array:
            x_reco_array = np.append(x_reco_array, x_reconstructed)
            # add reconstructed y-position to array:
            y_reco_array = np.append(y_reco_array, y_reconstructed)
            # add reconstructed z-position to array:
            z_reco_array = np.append(z_reco_array, z_reconstructed)

            # increment number of preselected events:
            number_preselected += 1

        else:
            # event is rejected:
            number_rejected += 1
            continue

    return (number_events, evt_id_preselected, edep_total, x_reco_array, y_reco_array, z_reco_array,
            number_preselected, number_rejected,
            number_vol_pass, number_vol_reject, number_vol_pass_initial,
            number_e_pass, number_mine_reject, number_maxe_reject,
            number_nmult_pass, number_nmult_reject, number_without_ncap,
            number_time_pass, number_time_reject,
            number_dist_pass, number_dist_reject)


def number_c12_atoms(radius_cut):
    """
    function to calculate the number of C12 atoms in the JUNO liquid scintillator for a specific volume of the central
    detector.

    :param radius_cut: radius, which defines the fiducial volume in the central detector in meter, normally 17m is used
    as radius of the fiducial volume like in the calculation of the IBD detection efficiency on page 39 of the yellow
    book (float)
    :return: number of C12 atoms (float)
    """
    # there are 20 ktons LS (which mainly contains LAB) in the central detector with R = 17.7 m.
    # mass LS in ktons:
    mass_ls = 20
    # radius of central detector in meter:
    radius_cd = 17.7
    # INFO-me: approximation, that LS consists only of LAB
    # mass of LS for volume cut with 'radius_cut' in tons:
    mass_ls_cut = radius_cut**3 / radius_cd**3 * mass_ls * 1000

    # the number of C12 atoms depends on the structure formula of LAB. LAB is C_6 H_5 C_n H_(2n+1), where n = 10, 11,
    # 12 or 13. Therefore the number of C12 atoms must be calculated for n= 10, 11, 12 and 13.

    " number of C12 in one LAB molecule: "
    num_c12_lab_n10 = 16
    num_c12_lab_n11 = 17
    num_c12_lab_n12 = 18
    num_c12_lab_n13 = 19

    # atomic mass number u in tons:
    u_in_tons = 1.6605 * 10**(-30)

    " mass of one LAB molecule in tons "
    mass_lab_n10 = (num_c12_lab_n10*12 + 26*1) * u_in_tons
    mass_lab_n11 = (num_c12_lab_n11*12 + 28*1) * u_in_tons
    mass_lab_n12 = (num_c12_lab_n12*12 + 30*1) * u_in_tons
    mass_lab_n13 = (num_c12_lab_n13*12 + 32*1) * u_in_tons

    # number of C12 for different n:
    number_c12_n10 = mass_ls_cut / mass_lab_n10 * num_c12_lab_n10
    number_c12_n11 = mass_ls_cut / mass_lab_n11 * num_c12_lab_n11
    number_c12_n12 = mass_ls_cut / mass_lab_n12 * num_c12_lab_n12
    number_c12_n13 = mass_ls_cut / mass_lab_n13 * num_c12_lab_n13

    # to calculate the number of C12 atoms for this fiducial volume, take the average:
    number_c12_cut = (number_c12_n10 + number_c12_n11 + number_c12_n12 + number_c12_n13) / 4

    return number_c12_cut


def read_xml_xsec(path_file, interval_e):
    """
    function to read xml-file, where the cross-sections of GENIE (gxspl-FNALsmall.xml) are save. These cross-section are
    also used in the GENIE simulation of Julia, where the interactions of neutrinos with C12 are simulated.

    :param path_file: path to the xml-file, where cross-sections are saved (gxspl-FNALsmall_nue.xml,
                      gxspl-FNALsmall_nuebar.xml, gxspl-FNALsmall_numu.xml, gxspl-FNALsmall_numubar.xml) (string)
    :param interval_e: energy interval (bin-width) in MeV

    :return: cross_section: sum of cross-section for certain neutrino flavor in cm**2 as function of the energy
    (array of float)
    """

    # read xml file:
    tree = ET.parse(path_file)
    # get root from xml file, root = genie_xsec_spline_list:
    root = tree.getroot()

    # preallocate energy-array in MeV:
    energy = np.arange(0, 10000+interval_e, interval_e)
    # preallocate cross-section array in cm**2:
    cross_section = np.zeros(len(energy))

    # loop over root, child = spline
    for child in root:
        # get attrib of child (dict):
        spline = child.attrib

        # get the name attribute from spline-dictionary
        # (e.g. spline_name = 'genie::AhrensNCELPXSec/Default/nu:14;tgt:1000060120;N:2112;proc:Weak[NC],QES;',
        # type=string)
        spline_name = spline.get('name')
        # print(spline_name)

        # INFO-me: for E=0MeV, set xsec=0cm**2 -> if not, interp sets all values below 10 MeV to value of xsec(10MeV)!
        # preallocate energy-array for this knot (MeV):
        e_spline = np.array([0])
        # preallocate cross-section-array for this knot (in cm**2):
        xsec_spline = np.array([0])

        # loop over child, subchild = knot, subchild is list with 2 values, first value = E, second value = xsec:
        for subchild in child:
            # get the value of E in GeV (string)
            e = subchild[0].text
            # convert to float and MeV:
            e = float(e) * 1000

            # get the value of cross-section in natural units (xsec) (string):
            xsec = subchild[1].text
            # convert to float and cm**2 (natural unit: 1/GeV**2 entspricht 3.89391289*10**(-28) cm**2):
            xsec = float(xsec) * 3.89391289 * 10**(-28)

            # append e to energy array:
            e_spline = np.append(e_spline, e)
            # append xsec to cross-section-array:
            xsec_spline = np.append(xsec_spline, xsec)

        # interpolate the cross-section of this spline to array "energy" to have the same bin-width and energy range for
        # all cross-section (also same bin-width and energy range like the fluxes):
        xsec_interp = np.interp(energy, e_spline, xsec_spline)

        # add cross-section to the total cross-section in cm**2:
        cross_section = cross_section + xsec_interp

    return cross_section


def event_rate(interval_energy, radius_cut, output_path, plot_flux, show_fluxplot, save_fluxplot, plot_evt_rate,
               show_evt_rate, save_evt_rate):
    """
    function to calculate the event rate of the atmospheric NC background in JUNO
    (equation: dN/dT(E) = A * fluxes * cross_sections):

    - 4*pi: factor to compensate that the flux is given in unit 1/sr

    - A: number of C12 atoms in the fiducial volume of the LS in JUNO, given by function number_c12_atoms(radius_cut)

    - fluxes for nu_e, nu_e_bar, nu_mu and nu_mu_bar:
        - fluxes of HONDA for JUNO site from 100 MeV to 10 GeV (4*pi: factor to compensate that the flux is given
        in unit 1/sr)
        - for energies below 100 MeV: take the shape of the simulation of FLUKA for Super-K location and normalize it
        to the HONDA flux
        - the total flux for nu_e, nu_e_bar, nu_mu and nu_mu_bar adding these two fluxes. The flux is then given from
        0 MeV to 10 GeV (bin-width in energy given by interval_energy)

    - neutrino neutral current interaction cross-sections on C12 (cross-section from GENIE:
    /home/astro/blum/juno/GENIE/genie_xsec/v2_12_0/NULL/DefaultPlusMECWithNC/data/gxspl-FNALsmall.xml, (also used in
    GENIE simulation)):
        - use function 'read_xml_xsec()':
        - cross-section of all 4 flavour: nu_e, nu_e_bar, nu_mu, nu_mu_bar
        - cross-section from 0 MeV to 10 GeV with a bin-width specified by 'interval_e'
        - for each flavour: total cross-section is sum of the cross-sections of different production channels
        (QES: quasi-elastic scattering, DIS: deep inelastic scattering, RES: resonant neutrino scattering,
        COH: coherent neutrino scattering)
        - cross-sections of gxspl-FNALsmall.xml are copied to 4 smaller xml-files:
        gxspl-FNALsmall_nue.xml, gxspl-FNALsmall_nuebar.xml, gxspl-FNALsmall_numu.xml, gxspl-FNALsmall_numubar.xml
    
    :param interval_energy: energy interval (bin-width) in MeV (float)
    :param radius_cut: radius, which defines the fiducial volume in the central detector in meter (normally 17m is used
    as radius of the fiducial volume like in the calculation of the IBD detection efficiency on page 39 of the yellow
    book) (float)
    :param output_path: path, where the plots are saved
    :param plot_flux: if True, flux is plotted (boolean)
    :param show_fluxplot: if True, plot of flux is shown (boolean)
    :param save_fluxplot: if True, plot of flux is saved (boolean)
    :param plot_evt_rate: if True, event rate is plotted as function of energy (boolean)
    :param show_evt_rate: if True, plot of event rate is shown (boolean)
    :param save_evt_rate: if True, plot of event rate is saved (boolean)

    :return: event_rate: event rate of atmospheric NC neutrino interactions in JUNO in events/second (float)
    """

    """ calculate the neutrino fluxes for all 4 flavours: """

    """ Results of the FLUKA simulation (from the paper of Battistoni2005 'The atmospheric neutrino fluxes below 
    100 MeV: The FLUKA results'): """
    # Neutrino energy in MeV from table 3 from paper 1-s2.0-S0927650505000526-main (np.array of float):
    energy_fluka = np.array([0, 13, 15, 17, 19, 21, 24, 27, 30, 33, 38, 42, 47, 53, 60, 67, 75, 84, 94, 106, 119, 133,
                             150, 168, 188, 211, 237, 266, 299, 335, 376, 422, 473, 531, 596, 668, 750, 841, 944])

    # differential flux from FLUKA in energy for no oscillation for electron-neutrinos for solar average at the site
    # of Super-Kamiokande, in (MeV**(-1) * cm**(-2) * s**(-1)) (np.array of float).
    # Assumption: for energy = 0 MeV, the flux is also 0!
    flux_nue_fluka = 10 ** (-4) * np.array([0, 69.6, 74.6, 79.7, 87.4, 94.2, 101., 103., 109., 108., 107., 101., 88.5,
                                            69.6, 64.4, 59.3, 54.3, 49.7, 45.1, 40.6, 35.8, 31.7, 27.3, 23.9, 20.4,
                                            17.0, 14.5, 12.0, 9.96, 8.11, 6.62, 5.27, 4.23, 3.37, 2.66, 2.09, 1.62,
                                            1.24, 0.950])

    # differential flux from FLUKA in energy for no oscillation for electron-antineutrinos for solar average at the site
    # of Super-Kamiokande, in (MeV**(-1) * cm**(-2) * s**(-1)) (np.array of float).
    # Assumption: for energy = 0 MeV, the flux is also 0!
    flux_nuebar_fluka = 10 ** (-4) * np.array([0, 63.7, 69.7, 79.5, 84.2, 89.4, 95.0, 99.3, 103., 104., 101., 96.1,
                                               83.5, 65.9, 60.0, 56.4, 51.4, 46.3, 43.0, 37.2, 32.9, 28.8, 24.9, 21.3,
                                               18.3, 15.4, 12.9, 10.6, 8.80, 7.13, 5.75, 4.60, 3.68, 2.88, 2.28,
                                               1.87, 1.37, 1.06, 0.800])

    # differential flux from FLUKA in energy for no oscillation for muon-neutrinos for solar average at the site
    # of Super-Kamiokande, in (MeV**(-1) * cm**(-2) * s**(-1)) (np.array of float).
    # Assumption: for energy = 0 MeV, the flux is also 0!
    flux_numu_fluka = 10 ** (-4) * np.array([0, 114., 124., 138., 146., 155., 159., 164., 181., 174., 179., 178., 176.,
                                             153., 131., 123., 114., 107., 96.3, 84.2, 72.7, 63.5, 55.2, 47.7, 41.2,
                                             34.4, 28.4, 23.6, 19.6, 15.8, 12.8, 10.3, 8.20, 6.49, 5.15, 3.98, 3.13,
                                             2.41, 1.82])

    # differential flux from FLUKA in energy for no oscillation for muon-antineutrinos for solar average at the site of
    # Super-K, in (MeV**(-1) * cm**(-2) * s**(-1)) (np.array of float).
    # Assumption: for energy = 0 MeV, the flux is also 0!
    flux_numubar_fluka = 10 ** (-4) * np.array([0, 116., 128., 136., 150., 158., 162., 170., 196., 177., 182., 183.,
                                                181., 155., 132., 123., 112., 101., 92.1, 82.2, 72.5, 64.0, 55.6,
                                                47.6, 40.8, 34.1, 28.6, 23.5, 19.3, 15.7, 12.6, 10.2, 8.15, 6.48,
                                                5.02, 3.94, 3.03, 2.33, 1.79])

    """ Results of the HONDA simulation (based on the paper of Honda2015: 'Atmospheric neutrino flux calculation using
    the NRLMSISE-00 atmospheric model'), (for solar maximum (HONDA_juno-ally-01-01-solmax.d)): """
    # Neutrino energy in MeV from the table from file HONDA_juno-ally-01-01-solmax.d (np.array of float):
    energy_honda = 10 ** 3 * np.array([1.0000E-01, 1.1220E-01, 1.2589E-01, 1.4125E-01, 1.5849E-01, 1.7783E-01,
                                       1.9953E-01, 2.2387E-01, 2.5119E-01, 2.8184E-01, 3.1623E-01, 3.5481E-01,
                                       3.9811E-01, 4.4668E-01, 5.0119E-01, 5.6234E-01, 6.3096E-01, 7.0795E-01,
                                       7.9433E-01, 8.9125E-01, 1.0000E+00, 1.1220E+00, 1.2589E+00, 1.4125E+00,
                                       1.5849E+00, 1.7783E+00, 1.9953E+00, 2.2387E+00, 2.5119E+00, 2.8184E+00,
                                       3.1623E+00, 3.5481E+00, 3.9811E+00, 4.4668E+00, 5.0119E+00, 5.6234E+00,
                                       6.3096E+00, 7.0795E+00, 7.9433E+00, 8.9125E+00, 1.0000E+01])

    # all-direction averaged flux for no oscillation for electron-neutrinos for solar maximum at the site of JUNO
    # (WITHOUT mountain over the detector), in (MeV**(-1) * cm**(-2) * s**(-1) * sr**(-1)) (np.array of float):
    flux_nue_max_honda = 10 ** (-7) * np.array([2.7743E+03, 2.4562E+03, 2.1530E+03, 1.8723E+03, 1.6181E+03, 1.3882E+03,
                                                1.1819E+03, 9.9837E+02, 8.3703E+02, 6.9614E+02, 5.7337E+02, 4.6903E+02,
                                                3.8128E+02, 3.0772E+02, 2.4680E+02, 1.9673E+02, 1.5600E+02, 1.2295E+02,
                                                9.6275E+01, 7.4975E+01, 5.8069E+01, 4.4733E+01, 3.4262E+01, 2.6103E+01,
                                                1.9770E+01, 1.4881E+01, 1.1137E+01, 8.2775E+00, 6.1088E+00, 4.4822E+00,
                                                3.2629E+00, 2.3653E+00, 1.7104E+00, 1.2266E+00, 8.7045E-01, 6.1557E-01,
                                                4.3368E-01, 3.0448E-01, 2.1286E-01, 1.4843E-01, 1.0281E-01])

    # all-direction averaged flux for no oscillation for electron-antineutrinos for solar maximum at the site of JUNO
    # (WITHOUT mountain over the detector), in (MeV**(-1) * cm**(-2) * s**(-1) * sr**(-1)) (np.array of float):
    flux_nuebar_max_honda = 10 ** (-7) * np.array([2.7733E+03, 2.4332E+03, 2.1124E+03, 1.8187E+03, 1.5545E+03,
                                                   1.3190E+03, 1.1105E+03, 9.2820E+02, 7.7040E+02, 6.3403E+02,
                                                   5.1790E+02, 4.1997E+02, 3.3811E+02, 2.7054E+02, 2.1539E+02,
                                                   1.7049E+02, 1.3418E+02, 1.0499E+02, 8.1651E+01, 6.3166E+01,
                                                   4.8654E+01, 3.7230E+01, 2.8329E+01, 2.1428E+01, 1.6121E+01,
                                                   1.2064E+01, 8.9697E+00, 6.6258E+00, 4.8598E+00, 3.5435E+00,
                                                   2.5685E+00, 1.8478E+00, 1.3252E+00, 9.4491E-01, 6.6836E-01,
                                                   4.7226E-01, 3.3159E-01, 2.3192E-01, 1.6107E-01, 1.1131E-01,
                                                   7.7823E-02])

    # all-direction averaged flux for no oscillation for muon-neutrinos for solar maximum at the site of JUNO
    # (WITHOUT mountain over the detector), in (MeV**(-1) * cm**(-2) * s**(-1) * sr^(-1)) (np.array of float):
    flux_numu_max_honda = 10 ** (-7) * np.array([5.7913E+03, 5.0884E+03, 4.4520E+03, 3.8714E+03, 3.3388E+03, 2.8520E+03,
                                                 2.4128E+03, 2.0226E+03, 1.6807E+03, 1.3858E+03, 1.1351E+03, 9.2472E+02,
                                                 7.4912E+02, 6.0324E+02, 4.8323E+02, 3.8514E+02, 3.0543E+02, 2.4122E+02,
                                                 1.8959E+02, 1.4845E+02, 1.1579E+02, 8.9940E+01, 6.9618E+01, 5.3647E+01,
                                                 4.1114E+01, 3.1343E+01, 2.3751E+01, 1.7914E+01, 1.3453E+01, 1.0049E+01,
                                                 7.4735E+00, 5.5296E+00, 4.0719E+00, 2.9889E+00, 2.1817E+00, 1.5909E+00,
                                                 1.1558E+00, 8.3657E-01, 6.0575E-01, 4.3508E-01, 3.1237E-01])

    # all-direction averaged flux for no oscillation for muon-antineutrinos for solar maximum at the site of JUNO
    # (WITHOUT mountain over the detector), in (MeV**(-1) * cm**(-2) * s**(-1) * sr^(-1)) (np.array of float):
    flux_numubar_max_honda = 10 ** (-7) * np.array([5.8966E+03, 5.1676E+03, 4.5104E+03, 3.9127E+03, 3.3665E+03,
                                                    2.8701E+03, 2.4238E+03, 2.0277E+03, 1.6821E+03, 1.3857E+03,
                                                    1.1333E+03, 9.2144E+02, 7.4476E+02, 5.9875E+02, 4.7865E+02,
                                                    3.8024E+02, 3.0060E+02, 2.3645E+02, 1.8519E+02, 1.4444E+02,
                                                    1.1204E+02, 8.6529E+01, 6.6529E+01, 5.0910E+01, 3.8731E+01,
                                                    2.9299E+01, 2.2048E+01, 1.6504E+01, 1.2291E+01, 9.1084E+00,
                                                    6.7266E+00, 4.9403E+00, 3.6136E+00, 2.6356E+00, 1.9115E+00,
                                                    1.3828E+00, 9.9752E-01, 7.1482E-01, 5.1189E-01, 3.6743E-01,
                                                    2.6256E-01])

    """ to compensate unit 1/sr from HONDA fluxes, multiply fluxes with 4*pi: """
    flux_nue_max_honda = flux_nue_max_honda * 4*np.pi
    flux_nuebar_max_honda = flux_nuebar_max_honda * 4*np.pi
    flux_numu_max_honda = flux_numu_max_honda * 4*np.pi
    flux_numubar_max_honda = flux_numubar_max_honda * 4*np.pi

    """ Extrapolate the HONDA flux to the energies of the FLUKA simulation from 0 MeV to 100 MeV: """
    """ Assumption:
        1. the shape of the FLUKA flux as function of energy do NOT depend on the location
            -> the shape of the flux at Super-K can also be used at JUNO site

        2. the absolute value of the FLUKA flux at Super-K should be normalized to the location of JUNO
            ->  therefore get the normalization factor by comparing the HONDA flux and the FLUKA flux in the energy 
                range from 100 MeV to 10 GeV        
    """
    # define the energy range from 0 MeV to 10000 MeV, corresponding to flux_nue_juno (array of float):
    energy_neutrino = np.arange(0, 10000 + interval_energy, interval_energy)

    # define the energy-array, in which the normalization will be calculated (neutrino energy in MeV)
    # (np.array of float):
    energy_norm = np.arange(min(energy_honda), max(energy_fluka) + 0.1, 0.1)

    """ For electron neutrinos: """
    # Interpolate the flux of FLUKA to get the differential flux in the energy range from 100 MeV to 950 MeV,
    # in 1/(MeV * cm**2 * s) (np.array of float):
    flux_nue_fluka_interpolated = np.interp(energy_norm, energy_fluka, flux_nue_fluka)

    # Interpolate the flux of HONDA to get the differential flux in the energy range from 100 MeV to 950 MeV,
    # in 1/(MeV * cm**2 * s) (np.array of float):
    flux_nue_honda_interpolated = np.interp(energy_norm, energy_honda, flux_nue_max_honda)

    # Calculate the integral of the FLUKA flux in the energy range given by energy_norm (float):
    integral_nue_fluka = np.trapz(flux_nue_fluka_interpolated, energy_norm)

    # Calculate the integral of the HONDA flux in the energy range given by energy_norm (float):
    integral_nue_honda = np.trapz(flux_nue_honda_interpolated, energy_norm)

    # Interpolate the part of the FLUKA flux in the energy range from 0 MeV to 99.9 MeV, in 1/(MeV*s*cm**2)
    # (np.array of float):
    flux_nue_fluka_interesting = np.interp(np.arange(0, 100, interval_energy), energy_fluka, flux_nue_fluka)

    # Normalize flux_nue_fluka_interesting at Super-K to the electron-neutrino flux at JUNO,
    # in 1/(MeV * s * cm**2) (np.array of float):
    flux_nue_fluka_norm = flux_nue_fluka_interesting * integral_nue_honda / integral_nue_fluka

    # interpolate flux_nue_max_honda in the energy range from 100 to 10000 MeV, in 1/(MeV*s*cm^2) (np.array of float):
    flux_nue_honda_interp = np.interp(np.arange(100, 10000+interval_energy, interval_energy), energy_honda,
                                      flux_nue_max_honda)

    # combine the normalized flux of FLUKA with the flux of HONDA (from 0 MeV to 100 MeV: FLUKA, above 100 MeV: Honda)
    # in 1/(MeV * s * cm**2) (np.array of float):
    flux_nue_juno = np.append(flux_nue_fluka_norm, flux_nue_honda_interp)

    """ For electron antineutrinos: """
    # Interpolate the flux of FLUKA to get the differential flux in the energy range from 100 MeV to 950 MeV,
    # in 1/(MeV * cm**2 * s) (np.array of float):
    flux_nuebar_fluka_interpolated = np.interp(energy_norm, energy_fluka, flux_nuebar_fluka)

    # Interpolate the flux of HONDA to get the differential flux in the energy range from 100 MeV to 950 MeV,
    # in 1/(MeV * cm**2 * s) (np.array of float):
    flux_nuebar_honda_interpolated = np.interp(energy_norm, energy_honda, flux_nuebar_max_honda)

    # Calculate the integral of the FLUKA flux in the energy range given by energy_norm (float):
    integral_nuebar_fluka = np.trapz(flux_nuebar_fluka_interpolated, energy_norm)

    # Calculate the integral of the HONDA flux in the energy range given by energy_norm (float):
    integral_nuebar_honda = np.trapz(flux_nuebar_honda_interpolated, energy_norm)

    # Interpolate the part of the FLUKA flux in the energy range from 0 MeV to 99.9 MeV, in 1/(MeV*s*cm**2)
    # (np.array of float):
    flux_nuebar_fluka_interesting = np.interp(np.arange(0, 100, interval_energy), energy_fluka, flux_nuebar_fluka)

    # Normalize flux_nuebar_fluka_interesting at Super-K to the electron-antineutrino flux at JUNO,
    # in 1/(MeV * s * cm**2) (np.array of float):
    flux_nuebar_fluka_norm = flux_nuebar_fluka_interesting * integral_nuebar_honda / integral_nuebar_fluka

    # interpolate flux_nuebar_max_honda in the energy range from 100 to 10000 MeV, in 1/(MeV*s*cm^2)
    # (np.array of float):
    flux_nuebar_honda_interp = np.interp(np.arange(100, 10000+interval_energy, interval_energy), energy_honda,
                                         flux_nuebar_max_honda)

    # combine the normalized flux of FLUKA with the flux of HONDA (from 0 MeV to 100 MeV: FLUKA, above 100 MeV: Honda)
    # in 1/(MeV * s * cm**2) (np.array of float):
    flux_nuebar_juno = np.append(flux_nuebar_fluka_norm, flux_nuebar_honda_interp)

    """ For muon neutrinos: """
    # Interpolate the flux of FLUKA to get the differential flux in the energy range from 100 MeV to 950 MeV,
    # in 1/(MeV * cm**2 * s) (np.array of float):
    flux_numu_fluka_interpolated = np.interp(energy_norm, energy_fluka, flux_numu_fluka)

    # Interpolate the flux of HONDA to get the differential flux in the energy range from 100 MeV to 950 MeV,
    # in 1/(MeV * cm**2 * s) (np.array of float):
    flux_numu_honda_interpolated = np.interp(energy_norm, energy_honda, flux_numu_max_honda)

    # Calculate the integral of the FLUKA flux in the energy range given by energy_norm (float):
    integral_numu_fluka = np.trapz(flux_numu_fluka_interpolated, energy_norm)

    # Calculate the integral of the HONDA flux in the energy range given by energy_norm (float):
    integral_numu_honda = np.trapz(flux_numu_honda_interpolated, energy_norm)

    # Interpolate the part of the FLUKA flux in the energy range from 0 MeV to 99.9 MeV, in 1/(MeV*s*cm**2)
    # (np.array of float):
    flux_numu_fluka_interesting = np.interp(np.arange(0, 100, interval_energy), energy_fluka, flux_numu_fluka)

    # Normalize flux_numu_fluka_interesting at Super-K to the muon-neutrino flux at JUNO,
    # in 1/(MeV * s * cm**2) (np.array of float):
    flux_numu_fluka_norm = flux_numu_fluka_interesting * integral_numu_honda / integral_numu_fluka

    # interpolate flux_numu_max_honda in the energy range from 100 to 10000 MeV, in 1/(MeV*s*cm^2)
    # (np.array of float):
    flux_numu_honda_interp = np.interp(np.arange(100, 10000+interval_energy, interval_energy), energy_honda,
                                       flux_numu_max_honda)

    # combine the normalized flux of FLUKA with the flux of HONDA (from 0 MeV to 100 MeV: FLUKA, above 100 MeV: Honda)
    # in 1/(MeV * s * cm**2) (np.array of float):
    flux_numu_juno = np.append(flux_numu_fluka_norm, flux_numu_honda_interp)

    """ For muon antineutrinos: """
    # Interpolate the flux of FLUKA to get the differential flux in the energy range from 100 MeV to 950 MeV,
    # in 1/(MeV * cm**2 * s) (np.array of float):
    flux_numubar_fluka_interpolated = np.interp(energy_norm, energy_fluka, flux_numubar_fluka)

    # Interpolate the flux of HONDA to get the differential flux in the energy range from 100 MeV to 950 MeV,
    # in 1/(MeV * cm**2 * s) (np.array of float):
    flux_numubar_honda_interpolated = np.interp(energy_norm, energy_honda, flux_numubar_max_honda)

    # Calculate the integral of the FLUKA flux in the energy range given by energy_norm (float):
    integral_numubar_fluka = np.trapz(flux_numubar_fluka_interpolated, energy_norm)

    # Calculate the integral of the HONDA flux in the energy range given by energy_norm (float):
    integral_numubar_honda = np.trapz(flux_numubar_honda_interpolated, energy_norm)

    # Interpolate the part of the FLUKA flux in the energy range from 0 MeV to 99.9 MeV, in 1/(MeV*s*cm**2)
    # (np.array of float):
    flux_numubar_fluka_interesting = np.interp(np.arange(0, 100, interval_energy), energy_fluka, flux_numubar_fluka)

    # Normalize flux_numubar_fluka_interesting at Super-K to the muon-antineutrino flux at JUNO,
    # in 1/(MeV * s * cm**2) (np.array of float):
    flux_numubar_fluka_norm = flux_numubar_fluka_interesting * integral_numubar_honda / integral_numubar_fluka

    # interpolate flux_numubar_max_honda in the energy range from 100 to 10000 MeV, in 1/(MeV*s*cm^2)
    # (np.array of float):
    flux_numubar_honda_interp = np.interp(np.arange(100, 10000+interval_energy, interval_energy), energy_honda,
                                          flux_numubar_max_honda)

    # combine the normalized flux of FLUKA with the flux of HONDA (from 0 MeV to 100 MeV: FLUKA, above 100 MeV: Honda)
    # in 1/(MeV * s * cm**2) (np.array of float):
    flux_numubar_juno = np.append(flux_numubar_fluka_norm, flux_numubar_honda_interp)

    if plot_flux:
        h1 = plt.figure(1, figsize=(15, 8))
        plt.semilogy(energy_neutrino, flux_nue_juno, "r--", label="$\\nu_e$")
        plt.semilogy(energy_neutrino, flux_nuebar_juno, "g--", label="$\\bar{\\nu}_e$")
        plt.semilogy(energy_neutrino, flux_numu_juno, "b--", label="$\\nu_\\mu$")
        plt.semilogy(energy_neutrino, flux_numubar_juno, "m--", label="$\\bar{\\nu}_\\mu$")
        plt.semilogy(energy_neutrino, flux_nue_juno + flux_nuebar_juno + flux_numu_juno + flux_numubar_juno, "k",
                     label="total flux")
        plt.xlabel("Neutrino energy $E_{\\nu}$ in MeV")
        plt.ylabel("Neutrino flux in (MeV $\\cdot$ s $\\cdot$ cm$^2$)$^{(-1)}$")
        plt.title("Neutrino fluxes at JUNO site")
        plt.xlim(xmin=0.0, xmax=10000.0)
        plt.ylim(ymin=1E-8)
        plt.legend()
        plt.grid()

        if show_fluxplot:
            plt.show()

        if save_fluxplot:
            plt.savefig(output_path + "neutrino_flux.png")
            plt.close(h1)

    """ Calculate the interaction cross-sections of neutrinos with C12 for each neutrino flavour: """
    # path, where cross-sections are saved:
    path_xsec = "/home/astro/blum/juno/GENIE/genie_xsec_2.12.0_eventrate/genie_xsec/v2_12_0/NULL/DefaultPlusMECWithNC" \
                "/data/"

    # for electron-neutrino:
    # define path, where cross-sections are saved (string):
    path_xsec_nue = path_xsec + "gxspl-FNALsmall_nue.xml"
    # calculate total cross-section with function 'read_xml_xsec()' (total cross-section in cm**2, array of float):
    xsec_nue = read_xml_xsec(path_xsec_nue, interval_energy)

    # for electron-antineutrino:
    # define path, where cross-sections are saved (string):
    path_xsec_nuebar = path_xsec + "gxspl-FNALsmall_nuebar.xml"
    # calculate total cross-section with function 'read_xml_xsec()' (total cross-section in cm**2, array of float):
    xsec_nuebar = read_xml_xsec(path_xsec_nuebar, interval_energy)

    # for muon-neutrino:
    # define path, where cross-sections are saved (string):
    path_xsec_numu = path_xsec + "gxspl-FNALsmall_numu.xml"
    # calculate total cross-section with function 'read_xml_xsec()' (total cross-section in cm**2, array of float):
    xsec_numu = read_xml_xsec(path_xsec_numu, interval_energy)

    # for muon-antineutrino:
    # define path, where cross-sections are saved (string):
    path_xsec_numubar = path_xsec + "gxspl-FNALsmall_numubar.xml"
    # calculate total cross-section with function 'read_xml_xsec()' (total cross-section in cm**2, array of float):
    xsec_numubar = read_xml_xsec(path_xsec_numubar, interval_energy)

    """ Calculate the number of C12 atoms in the LS for fiducial volume with radius = radius_cut: """
    # number of C12 atoms in fiducial volume (float):
    number_c12 = number_c12_atoms(radius_cut)

    """ Event rate of atmospheric NC neutrino interactions in JUNO as function of energy: """
    # event rate as function of energy in events/(MeV * s) (array of float) (equ. 1 in AtmNeuBkgStudies_DocDB3884.pdf):
    evt_rate_per_energy = number_c12 * (flux_nue_juno * xsec_nue + flux_nuebar_juno * xsec_nuebar +
                                        flux_numu_juno * xsec_numu + flux_numubar_juno * xsec_numubar)

    if plot_evt_rate:
        h2 = plt.figure(2, figsize=(15, 8))
        plt.plot(energy_neutrino, evt_rate_per_energy, "b")
        plt.xlabel("Neutrino energy $E_{\\nu}$ in MeV")
        plt.ylabel("Event rate in evts/(MeV $\\cdot$ s)")
        plt.title("Event rate of atmospheric NC neutrino interaction in JUNO detector")
        plt.xlim(xmin=0.0, xmax=10000.0)
        plt.ylim(ymin=0.0)
        plt.grid()

        if show_evt_rate:
            plt.show()

        if save_evt_rate:
            plt.savefig(output_path + "event_rate.png")
            plt.close(h2)

    """ Event rate of atmospheric NC neutrino interactions in JUNO: """
    # integrate evt_rate_per_energy over the energy to get the total event rate in events/s (float):
    evt_rate = np.trapz(evt_rate_per_energy, energy_neutrino)

    return evt_rate


def convert_genie_file_for_generator(rootfile_input, path_output):
    """
    function to convert the 'original' GENIE root-file of Julia to a root-file, which can be used as input for the
    DSNB-NC.exe generator of JUNO offline.

    INFO:
    The variable isopdg can only be set for the following isotopes: C11, B11, C10, B10, Be10, C9, B9, Be9, Li9,
    B8, Li8, Be7 and Li7.
    For all other isotopes (e.g. C8, Be8, He8, B7, He7, H7 and lighter isotopes), NO isopdg should be set, because there
    are no TALYS files for the de-excitation of these isotopes and therefore you get an error in DSNB-NC.exe.

    :param rootfile_input: path to the original GENIE ROOT-file (for example: gntp.101.gst.root (string)
    :param path_output: path, where the output ROOT file is saved (string)
    :return:
    """
    # load the ROOT file:
    rfile_input = ROOT.TFile(rootfile_input)
    # get the TTree from the TFile:
    rtree_input = rfile_input.Get("gst")

    # Info-me: "gst;13" is a copy of meta data of "gst;14", "gst;14" contains correct data and is read

    # set new ROOT file:
    rfile_output = ROOT.TFile(path_output + "genie_data_NC.root", "recreate")
    # set new ROOT Tree:
    rtree_output = ROOT.TTree("particleT", "particle Tree")

    # get the number of entries in the ROOT-file:
    number_entries = rtree_input.GetEntries()
    # number_entries = 100

    # set the maximal number of final particles for one event:
    max_n_par = 20

    """ preallocate all arrays: """
    # INFO-me: type code 'd' is float-type in python and double-type in C and ROOT
    # event number (integer):
    event_number = array('i', [0])
    # neutrino PDG code (integer):
    p_pdg = array('i', [0])
    # Nuclear target PDG code (integer):
    t_pdg = array('i', [0])
    # channel id of the NC interaction (integer):
    channel_id = array('i', [0])
    # incoming neutrino energy in GeV (double):
    p_en = array('d', [0.])
    # incoming neutrino x-momentum in GeV (double):
    p_px = array('d', [0.])
    # incoming neutrino y-momentum in GeV (double):
    p_py = array('d', [0.])
    # incoming neutrino z-momentum in GeV (double):
    p_pz = array('d', [0.])
    # PDG code of produced isotope (integer):
    m_isopdg = array('i', [0])
    # momentum of produced isotope in GeV (double):
    m_isopx = array('d', [0.])
    # momentum of produced isotope in GeV (double):
    m_isopy = array('d', [0.])
    # momentum of produced isotope in GeV (double):
    m_isopz = array('d', [0.])
    # mass of produced isotope in GeV (double):
    m_isomass = array('d', [0.])
    # number of final state particles in hadronic system (integer):
    n_pars = array('i', [0])
    # PDG code of i-th final state particle in hadronic system (array of integer):
    pdg = array('i', max_n_par*[0])
    # energy of i-th final state particle in hadronic system in GeV (array of double):
    energy = array('d', max_n_par*[0.])
    # Px of i-th final state particle in hadronic system in GeV (array of double):
    px = array('d', max_n_par*[0.])
    # Py of i-th final state particle in hadronic system in GeV (array of double):
    py = array('d', max_n_par*[0.])
    # Pz of i-th final state particle in hadronic system in GeV (array of double):
    pz = array('d', max_n_par*[0.])
    # mass of the i-th final particle in GeV (array of double):
    mass = array('d', max_n_par*[0.])

    # Add the arrays to the TBranch:
    # INFO-me: D indicates that values are 'double'-type like required from the DSNB-NC.exe generator
    rtree_output.Branch('pPdg', p_pdg, 'pPdg/I')
    rtree_output.Branch('tPdg', t_pdg, 'tPdg/I')
    rtree_output.Branch('channelID', channel_id, "channelID/I")
    rtree_output.Branch('pEn', p_en, "pEn/D")
    rtree_output.Branch('pPx', p_px, "pPx/D")
    rtree_output.Branch('pPy', p_py, "pPy/D")
    rtree_output.Branch('pPz', p_pz, "pPz/D")
    rtree_output.Branch('m_isoPdg', m_isopdg, 'm_isoPdg/I')
    rtree_output.Branch('m_isoPx', m_isopx, 'm_isoPx/D')
    rtree_output.Branch('m_isoPy', m_isopy, 'm_isoPy/D')
    rtree_output.Branch('m_isoPz', m_isopz, 'm_isoPz/D')
    rtree_output.Branch('m_isoMass', m_isomass, 'm_isoMass/D')
    rtree_output.Branch('Npars', n_pars, 'Npars/I')
    rtree_output.Branch('pdg', pdg, 'pdg[Npars]/I')
    rtree_output.Branch('px', px, 'px[Npars]/D')
    rtree_output.Branch('py', py, 'py[Npars]/D')
    rtree_output.Branch('pz', pz, 'pz[Npars]/D')
    rtree_output.Branch('energy', energy, 'energy[Npars]/D')
    rtree_output.Branch('mass', mass, 'mass[Npars]/D')

    """ Read the data from the TTree: """
    # loop over every entry, i.e. every event, in the TTree:
    for event in range(number_entries):

        # get the current event in the TTree:
        rtree_input.GetEntry(event)

        # is it a quasi-elastic scattering event? (0 = no QEL event, 1 = QEL event):
        qel = rtree_input.GetBranch('qel').GetLeaf('qel').GetValue()
        qel = int(qel)

        # is it a NC event? (0 = no NC event, 1 = NC event):
        nc = rtree_input.GetBranch('nc').GetLeaf('nc').GetValue()
        nc = int(nc)

        # get the value of target PDG:
        tgt = rtree_input.GetBranch('tgt').GetLeaf('tgt').GetValue()
        tgt = int(tgt)

        # read only NC and QEL events:
        # if qel == 1 and nc == 1:

        # read only NC events and interactions on C12:
        # if nc == 1 and tgt == 1000060120:

        # read only NC, but all targets:
        if nc == 1:

            # set the event number:
            event_number[0] = event

            # get the value of neutrino PDG:
            neu = rtree_input.GetBranch('neu').GetLeaf('neu').GetValue()
            neu = int(neu)
            p_pdg[0] = neu

            # get the value of target PDG:
            # tgt = rtree_input.GetBranch('tgt').GetLeaf('tgt').GetValue()
            # tgt = int(tgt)
            t_pdg[0] = tgt

            # get the value of number of final p:
            nfp = rtree_input.GetBranch('nfp').GetLeaf('nfp').GetValue()
            nfp = int(nfp)

            # get the value of number of final n:
            nfn = rtree_input.GetBranch('nfn').GetLeaf('nfn').GetValue()
            nfn = int(nfn)

            # get the value of number of final pi_minus:
            nfpim = rtree_input.GetBranch('nfpim').GetLeaf('nfpim').GetValue()
            nfpim = int(nfpim)

            # get the value of number of final pi_plus:
            nfpip = rtree_input.GetBranch('nfpip').GetLeaf('nfpip').GetValue()
            nfpip = int(nfpip)

            # get the value of number of final Kaon_minus:
            nfkm = rtree_input.GetBranch('nfkm').GetLeaf('nfkm').GetValue()
            nfkm = int(nfkm)

            # get the value of number of final Kaon_plus:
            nfkp = rtree_input.GetBranch('nfkp').GetLeaf('nfkp').GetValue()
            nfkp = int(nfkp)

            # preallocate the channel ID (integer):
            ch_id = int(0)
            # calculate the channel ID (n_nu, n_p, n_n, n_piminus, n_piplus, n_Kminus, n_Kplus)
            if tgt == 1000060120:
                # target C12:
                if nfkm == 0 and nfkp == 0:
                    # no Kaons:
                    ch_id = int(str(1) + str(nfp) + str(nfn) + str(nfpim) + str(nfpip))
                else:
                    # with Kaons:
                    ch_id = int(str(1) + str(nfp) + str(nfn) + str(nfpim) + str(nfpip) + str(nfkm) + str(nfkp))

            elif tgt == 2212:
                # target proton:
                if nfp == 1 and nfn == 0 and nfpim == 0 and nfpip == 0 and nfkm == 0 and nfkp == 0:
                    # interaction channel: nu + p -> nu + p (elastic scattering):
                    ch_id = int(2)

                elif nfkm == 0 and nfkp == 0:
                    # no Kaons:
                    ch_id = int(str(1) + str(nfp) + str(nfn) + str(nfpim) + str(nfpip))

                else:
                    # with Kaons:
                    ch_id = int(str(1) + str(nfp) + str(nfn) + str(nfpim) + str(nfpip) + str(nfkm) + str(nfkp))

            elif tgt == 11 or tgt == 1000080160 or tgt == 1000070140 or tgt == 1000160320:
                # target either electron, N14, O16 or S32 (interaction: nu + tgt -> nu + tgt (elastic scattering)):
                ch_id = int(3)

            else:
                print("other target PDG as expected (no C12, p, e, N14, O16, S32)")
                print(tgt)

            # add ch_id to the tree (integer):
            channel_id[0] = ch_id

            # preallocate the PDG code and the mass of the isotope:
            isopdg = int(0)
            isomass = float(0.0)
            # calculate the PDG of the isotope:
            if tgt == 1000060120:
                # target C12
                if (nfp + nfn) == 1:
                    # possible isotopes: B11, C11
                    if nfp == 1 and nfn == 0 and (nfpim - nfpip) == 0:
                        # interaction: nu + C12 -> nu + B11 + p + (N*pi_minus + N*pi_plus):
                        isopdg = int(1000050110)
                        isomass = float(get_mass_from_pdg(isopdg))
                    elif nfp == 0 and nfn == 1 and (nfpim - nfpip) == -1:
                        # interaction: nu + C12 -> nu + B11 + n + (N*pi_minus + (N+1)*pi_plus):
                        isopdg = int(1000050110)
                        isomass = float(get_mass_from_pdg(isopdg))
                    elif nfp == 0 and nfn == 1 and (nfpim - nfpip) == 0:
                        # channel: nu + C12 -> nu + C11 + n + (N*pi_minus + N*pi_plus):
                        isopdg = int(1000060110)
                        isomass = float(get_mass_from_pdg(isopdg))
                    elif nfp == 1 and nfn == 0 and (nfpim - nfpip) == 1:
                        # channel: nu + C12 -> nu + C11 + p + ((N+1)*pi_minus + N*pi_plus):
                        isopdg = int(1000060110)
                        isomass = float(get_mass_from_pdg(isopdg))
                    else:
                        print("other possible channels with B11 and C11: nfp={0:d}, nfn={1:d}, nfpim={2:d}, nfpip={3:d}"
                              .format(nfp, nfn, nfpim, nfpip))

                elif (nfp + nfn) == 2:
                    # possible isotopes: B10, C10, Be10
                    if nfp == 1 and nfn == 1 and (nfpim - nfpip) == 0:
                        # channel: nu + C12 -> nu + B10 + p + n + (N*pi_minus + N*pi_plus):
                        isopdg = int(1000050100)
                        isomass = float(get_mass_from_pdg(isopdg))
                    elif nfp == 2 and nfn == 0 and (nfpim - nfpip) == 1:
                        # channel: nu + C12 -> nu + B10 + 2p + ((N+1)*pi_minus + N*pi_plus):
                        isopdg = int(1000050100)
                        isomass = float(get_mass_from_pdg(isopdg))
                    elif nfp == 0 and nfn == 2 and (nfpim - nfpip) == -1:
                        # channel: nu + C12 -> nu + B10 + 2n + (N*pi_minus + (N+1)*pi_plus):
                        isopdg = int(1000050100)
                        isomass = float(get_mass_from_pdg(isopdg))
                    elif nfp == 0 and nfn == 2 and (nfpim - nfpip) == 0:
                        # channel: nu + C12 -> nu + C10 + 2n + (N*pi_minus + N*pi_plus):
                        isopdg = int(1000060100)
                        isomass = float(get_mass_from_pdg(isopdg))
                    elif nfp == 1 and nfn == 1 and (nfpim - nfpip) == 1:
                        # channel: nu + C12 -> nu + C10 + p + n + ((N+1)*pi_minus + N*pi_plus):
                        isopdg = int(1000060100)
                        isomass = float(get_mass_from_pdg(isopdg))
                    elif nfp == 2 and nfn == 0 and (nfpim - nfpip) == 2:
                        # channel: nu + C12 -> nu + C10 + 2p + ((N+2)*pi_minus + N*pi_plus):
                        isopdg = int(1000060100)
                        isomass = float(get_mass_from_pdg(isopdg))
                    elif nfp == 2 and nfn == 0 and (nfpim - nfpip) == 0:
                        # channel: nu + C12 -> nu + Be10 + 2p + (N*pi_minus + N*pi_plus):
                        isopdg = int(1000040100)
                        isomass = float(get_mass_from_pdg(isopdg))
                    elif nfp == 1 and nfn == 1 and (nfpim - nfpip) == -1:
                        # channel: nu + C12 -> nu + Be10 + p + n + (N*pi_minus + (N+1)*pi_plus):
                        isopdg = int(1000040100)
                        isomass = float(get_mass_from_pdg(isopdg))
                    elif nfp == 0 and nfn == 2 and (nfpim - nfpip) == -2:
                        # channel: nu + C12 -> nu + Be10 + 2n + (N*pi_minus + (N+2)*pi_plus):
                        isopdg = int(1000040100)
                        isomass = float(get_mass_from_pdg(isopdg))
                    else:
                        print("other possible channels with C10, B10, Be10: nfp={0:d}, nfn={1:d}, nfpim={2:d}, "
                              "nfpip={3:d}".format(nfp, nfn, nfpim, nfpip))

                elif (nfp + nfn) == 3:
                    # possible isotopes: B9, C9, Be9, Li9
                    if nfp == 1 and nfn == 2 and (nfpim - nfpip) == 0:
                        # channel: nu + C12 -> nu + B9 + p + 2n + (N*pi_minus + N*pi_plus):
                        isopdg = int(1000050090)
                        isomass = float(get_mass_from_pdg(isopdg))
                    elif nfp == 2 and nfn == 1 and (nfpim - nfpip) == 1:
                        # channel: nu + C12 -> nu + B9 + 2p + n + ((N+1)*pi_minus + N*pi_plus):
                        isopdg = int(1000050090)
                        isomass = float(get_mass_from_pdg(isopdg))
                    elif nfp == 0 and nfn == 3 and (nfpim - nfpip) == -1:
                        # channel: nu + C12 -> nu + B9 + 3n + (N*pi_minus + (N+1)*pi_plus):
                        isopdg = int(1000050090)
                        isomass = float(get_mass_from_pdg(isopdg))
                    elif nfp == 3 and nfn == 0 and (nfpim - nfpip) == 2:
                        # channel: nu + C12 -> nu + B9 + 3p + ((N+2)*pi_minus + N*pi_plus):
                        isopdg = int(1000050090)
                        isomass = float(get_mass_from_pdg(isopdg))
                    elif nfp == 0 and nfn == 3 and (nfpim - nfpip) == 0:
                        # channel: nu + C12 -> nu + C9 + 3n + (N*pi_minus + N*pi_plus):
                        isopdg = int(1000060090)
                        isomass = float(get_mass_from_pdg(isopdg))
                    elif nfp == 1 and nfn == 2 and (nfpim - nfpip) == 1:
                        # channel: nu + C12 -> nu + C9 + p + 2n + ((N+1)*pi_minus + N*pi_plus):
                        isopdg = int(1000060090)
                        isomass = float(get_mass_from_pdg(isopdg))
                    elif nfp == 2 and nfn == 1 and (nfpim - nfpip) == 2:
                        # channel: nu + C12 -> nu + C9 + 2p + n + ((N+2)*pi_minus + N*pi_plus):
                        isopdg = int(1000060090)
                        isomass = float(get_mass_from_pdg(isopdg))
                    elif nfp == 2 and nfn == 1 and (nfpim - nfpip) == 0:
                        # channel: nu + C12 -> nu + Be9 + 2p + n + (N*pi_minus + N*pi_plus):
                        isopdg = int(1000040090)
                        isomass = float(get_mass_from_pdg(isopdg))
                    elif nfp == 1 and nfn == 2 and (nfpim - nfpip) == -1:
                        # channel: nu + C12 -> nu + Be9 + p + 2n + (N*pi_minus + (N+1)*pi_plus):
                        isopdg = int(1000040090)
                        isomass = float(get_mass_from_pdg(isopdg))
                    elif nfp == 3 and nfn == 0 and (nfpim - nfpip) == 1:
                        # channel: nu + C12 -> nu + Be9 + 3p + ((N+1)*pi_minus + N*pi_plus):
                        isopdg = int(1000040090)
                        isomass = float(get_mass_from_pdg(isopdg))
                    elif nfp == 0 and nfn == 3 and (nfpim - nfpip) == -2:
                        # channel: nu + C12 -> nu + Be9 + 3n + (N*pi_minus + (N+2)*pi_plus):
                        isopdg = int(1000040090)
                        isomass = float(get_mass_from_pdg(isopdg))
                    elif nfp == 3 and nfn == 0 and (nfpim - nfpip) == 0:
                        # channel: nu + C12 -> nu + Li9 + 3p + (N*pi_minus + N*pi_plus):
                        isopdg = int(1000030090)
                        isomass = float(get_mass_from_pdg(isopdg))
                    elif nfp == 2 and nfn == 1 and (nfpim - nfpip) == -1:
                        # channel: nu + C12 -> nu + Li9 + 2p + n + (N*pi_minus + (N+1)*pi_plus):
                        isopdg = int(1000030090)
                        isomass = float(get_mass_from_pdg(isopdg))
                    elif nfp == 1 and nfn == 2 and (nfpim - nfpip) == -2:
                        # channel: nu + C12 -> nu + Li9 + p + 2n + (N*pi_minus + (N+2)*pi_plus):
                        isopdg = int(1000030090)
                        isomass = float(get_mass_from_pdg(isopdg))
                    else:
                        print("other possible channels with B9, C9, Be9, Li9: nfp={0:d}, nfn={1:d}, nfpim={2:d}, "
                              "nfpip={3:d}".format(nfp, nfn, nfpim, nfpip))

                elif (nfp + nfn) == 4:
                    # possible isotopes: B8, Li8
                    if nfp == 1 and nfn == 3 and (nfpim - nfpip) == 0:
                        # channel: nu + C12 -> nu + B8 + p + 3n + (N*pi_minus + N*pi_plus):
                        isopdg = int(1000050080)
                        isomass = float(get_mass_from_pdg(isopdg))
                    elif nfp == 2 and nfn == 2 and (nfpim - nfpip) == 1:
                        # channel: nu + C12 -> nu + B8 + 2p + 2n + ((N+1)*pi_minus + N*pi_plus):
                        isopdg = int(1000050080)
                        isomass = float(get_mass_from_pdg(isopdg))
                    elif nfp == 0 and nfn == 4 and (nfpim - nfpip) == -1:
                        # channel: nu + C12 -> nu + B8 + 4n + (N*pi_minus + (N+1)*pi_plus):
                        isopdg = int(1000050080)
                        isomass = float(get_mass_from_pdg(isopdg))
                    elif nfp == 2 and nfn == 2 and (nfpim - nfpip) == 0:
                        # channel: nu + C12 -> nu + Be8 + 2p + 2n + (N*pi_minus + N*pi_plus):
                        print("Be8")
                        # isopdg = int(1000040080)
                        # isomass = float(get_mass_from_pdg(isopdg))
                    # elif nfp == 3 and nfn == 1 and (nfpim - nfpip) == 1:
                        # channel: nu + C12 -> nu + Be8 + 3p + n + ((N+1)*pi_minus + N*pi_plus):
                        # isopdg = int(1000040080)
                        # isomass = float(get_mass_from_pdg(isopdg))
                    # elif nfp == 1 and nfn == 3 and (nfpim - nfpip) == -1:
                        # channel: nu + C12 -> nu + Be8 + p + 3n + (N*pi_minus + (N+1)*pi_plus):
                        # isopdg = int(1000040080)
                        # isomass = float(get_mass_from_pdg(isopdg))
                    # elif nfp == 0 and nfn == 4 and (nfpim - nfpip) == -2:
                        # channel: nu + C12 -> nu + Be8 + 4n + (N*pi_minus + (N+2)*pi_plus):
                        # isopdg = int(1000040080)
                        # isomass = float(get_mass_from_pdg(isopdg))
                    # elif nfp == 4 and nfn == 0 and (nfpim - nfpip) == 2:
                        # channel: nu + C12 -> nu + Be8 + 4p + ((N+2)*pi_minus + N*pi_plus):
                        # isopdg = int(1000040080)
                        # isomass = float(get_mass_from_pdg(isopdg))
                    elif nfp == 3 and nfn == 1 and (nfpim - nfpip) == 0:
                        # channel: nu + C12 -> nu + Li8 + 3p + n + (N*pi_minus + N*pi_plus):
                        isopdg = int(1000030080)
                        isomass = float(get_mass_from_pdg(isopdg))
                    elif nfp == 4 and nfn == 0 and (nfpim - nfpip) == 1:
                        # channel: nu + C12 -> nu + Li8 + 4p + ((N+1)*pi_minus + N*pi_plus):
                        isopdg = int(1000030080)
                        isomass = float(get_mass_from_pdg(isopdg))
                    elif nfp == 2 and nfn == 2 and (nfpim - nfpip) == -1:
                        # channel: nu + C12 -> nu + Li8 + 2p + 2n + (N*pi_minus + (N+1)*pi_plus):
                        isopdg = int(1000030080)
                        isomass = float(get_mass_from_pdg(isopdg))
                    elif nfp == 1 and nfn == 3 and (nfpim - nfpip) == -2:
                        # channel: nu + C12 -> nu + Li8 + p + 3n + (N*pi_minus + (N+2)*pi_plus):
                        isopdg = int(1000030080)
                        isomass = float(get_mass_from_pdg(isopdg))
                    elif nfp == 0 and nfn == 4 and (nfpim - nfpip) == 0:
                        # channel: nu + C12 -> nu + C8 + 4n + (N*pi_minus + N*pi_plus):
                        print("C8")
                        # isopdg = int(1000060080)
                        # isomass = float(get_mass_from_pdg(isopdg))
                    elif nfp == 4 and nfn == 0 and (nfpim - nfpip) == 0:
                        # channel: nu + C12 -> nu + He8 + 4p + (N*pi_minus + N*pi_plus):
                        print("He8")
                        # isopdg = int(1000020080)
                        # isomass = float(get_mass_from_pdg(isopdg))
                    else:
                        print("other possible channels with B8, Be8, Li8, C8, He8: nfp={0:d}, nfn={1:d}, nfpim={2:d}, "
                              "nfpip={3:d}".format(nfp, nfn, nfpim, nfpip))


                elif (nfn + nfp) == 5:
                    # possible isotopes: Be7, Li7, B7, He7, H7
                    if nfp == 2 and nfn == 3 and (nfpim - nfpip) == 0:
                        # channel: nu + C12 -> nu + Be7 + 2p + 3n + (N*pi_minus + N*pi_plus):
                        isopdg = int(1000040070)
                        isomass = float(get_mass_from_pdg(isopdg))
                    elif nfp == 1 and nfn == 4 and (nfpim - nfpip) == -1:
                        # channel: nu + C12 -> nu + Be7 + p + 4n + (N*pi_minus + (N+1)*pi_plus):
                        isopdg = int(1000040070)
                        isomass = float(get_mass_from_pdg(isopdg))
                    elif nfp == 3 and nfn == 2 and (nfpim - nfpip) == 1:
                        # channel: nu + C12 -> nu + Be7 + 3p + 2n + ((N+1)*pi_minus + N*pi_plus):
                        isopdg = int(1000040070)
                        isomass = float(get_mass_from_pdg(isopdg))
                    elif nfp == 4 and nfn == 1 and (nfpim - nfpip) == 2:
                        # channel: nun + C12 -> nu + Be7 + 4p + n + ((N+2)*pi_minus + N*pi_plus):
                        isopdg = int(1000040070)
                        isomass = float(get_mass_from_pdg(isopdg))
                    elif nfp == 3 and nfn == 2 and (nfpim - nfpip) == 0:
                        # channel: nu + C12 -> nu + Li7 + 3p + 2n + (N*pi_minus + N*pi_plus):
                        isopdg = int(1000030070)
                        isomass = float(get_mass_from_pdg(isopdg))
                    elif nfp == 2 and nfn == 3 and (nfpim - nfpip) == -1:
                        # channel: nu + C12 -> nu + Li7 + 2p + 3n + (N*pi_minus + (N+1)*pi_plus):
                        isopdg = int(1000030070)
                        isomass = float(get_mass_from_pdg(isopdg))
                    elif nfp == 4 and nfn == 1 and (nfpim - nfpip) == 1:
                        # channel: nu + C12 -> nu + Li7 + 4p + n + ((N+1)*pi_minus + N*pi_plus):
                        isopdg = int(1000030070)
                        isomass = float(get_mass_from_pdg(isopdg))
                    elif nfp == 1 and nfn == 4 and (nfpim - nfpip) == 0:
                        # channel: nu + C12 -> nu + B7 + p + 4n + (N*pi_minus + N*pi_plus):
                        print("B7")
                        # isopdg = int(1000050070)
                        # isomass = float(get_mass_from_pdg(isopdg))
                    elif nfp == 4 and nfn == 1 and (nfpim - nfpip) == 0:
                        # channel: nu + C12 -> nu + He7 + 4p + n + (N*pi_minus + N*pi_plus):
                        print("He7")
                        # isopdg = int(1000020070)
                        # isomass = float(get_mass_from_pdg(isopdg))
                    elif nfp == 5 and nfn == 0 and (nfpim - nfpip) == 0:
                        # channel: nu + C12 -> nu + H7 + 5p + (N*pi_minus + N*pi_plus):
                        print("H7")
                        # isopdg = int(1000010070)
                        # isomass = float(get_mass_from_pdg(isopdg))
                    else:
                        print("other possible channels with Be7, Li7, B7, He7, H7: nfp={0:d}, nfn={1:d}, nfpim={2:d}, "
                              "nfpip={3:d}".format(nfp, nfn, nfpim, nfpip))

                # elif (nfn + nfp) == 6:
                # possible isotopes: Li6, Be6, He6, H6:
                #     if nfp == 3 and nfn == 3 and (nfpim - nfpip) == 0:
                #         # channel: nu + C12 -> nu + Li6 + 3p + 3n + (N*pi_minus + N*pi_plus):
                #         isopdg = int(1000030060)
                #         isomass = float(get_mass_from_pdg(isopdg))
                #     elif nfp == 2 and nfn == 4 and (nfpim - nfpip) == -1:
                #         # channel: nu + C12 -> nu + Li6 + 2p + 4n + (N*pi_minus + (N+1)*pi_plus):
                #         isopdg = int(1000030060)
                #         isomass = float(get_mass_from_pdg(isopdg))
                #     elif nfp == 4 and nfn == 2 and (nfpim - nfpip) == 1:
                #         # channel: nu + C12 -> nu + Li6 + 4p + 2n + ((N+1)*pi_minus + N*pi_plus):
                #         isopdg = int(1000030060)
                #         isomass = float(get_mass_from_pdg(isopdg))
                #     elif nfp == 5 and nfn == 1 and (nfpim - nfpip) == 2:
                #         # channel: nu + C12 -> nu + Li6 + 5p + n + ((N+2)*pi_minus + N*pi_plus):
                #         isopdg = int(1000030060)
                #         isomass = float(get_mass_from_pdg(isopdg))
                #     elif nfp == 1 and nfn == 5 and (nfpim - nfpip) == -2:
                #         # channel: nu + C12 -> nu + Li6 + p + 5n + (N*pi_minus + (N+2)*pi_plus):
                #         isopdg = int(1000030060)
                #         isomass = float(get_mass_from_pdg(isopdg))
                #     elif nfp == 2 and nfn == 4 and (nfpim - nfpip) == 0:
                #         # channel: nu + C12 -> nu + Be6 + 2p + 4n + (N*pi_minus + N*pi_plus):
                #         isopdg = int(1000040060)
                #         isomass = float(get_mass_from_pdg(isopdg))
                #     elif nfp == 4 and nfn == 2 and (nfpim - nfpip) == 0:
                #         # channel: nu + C12 -> nu + He6 + 4p + 2n + (N*pi_minus + N*pi_plus):
                #         isopdg = int(1000020060)
                #         isomass = float(get_mass_from_pdg(isopdg))
                #     elif nfp == 5 and nfn == 1 and (nfpim - nfpip) == 0:
                #         # channel: nu + C12 -> nu + H6 + 5p + n + (N*pi_minus + N*pi_plus):
                #         isopdg = int(1000010060)
                #         isomass = float(get_mass_from_pdg(isopdg))
                #     else:
                #         # other possible channel with Li6 ot interactions, where Be6, He6 or H6 are produced:
                #         # dummy isopdg for this case:
                #         isopdg = int(60000000)

                else:
                    # print("other possible channels lighter isotopes (mass<=7): nfp={0:d}, nfn={1:d}, nfpim={2:d}, "
                    #       "nfpip={3:d}".format(nfp, nfn, nfpim, nfpip))
                    lala = 1

            else:
                # target: p, e, N14, O16 or S32 -> no isotope is produced:
                isopdg = int(0)
                isomass = float(0.0)

            # add isopdg to the tree (integer):
            m_isopdg[0] = isopdg
            # add momentum of isotope to the tree (set momentum=0, because there is no information about the momentum)
            # (float):
            m_isopx[0] = float(0)
            m_isopy[0] = float(0)
            m_isopz[0] = float(0)
            # add mass of isotope to the tree (float):
            m_isomass[0] = isomass

            # get the value of energy of incoming neutrino:
            ev = rtree_input.GetBranch('Ev').GetLeaf('Ev').GetValue()
            p_en[0] = ev

            # get the x momentum of incoming neutrino:
            pxv = rtree_input.GetBranch('pxv').GetLeaf('pxv').GetValue()
            p_px[0] = pxv

            # get the y momentum of incoming neutrino:
            pyv = rtree_input.GetBranch('pyv').GetLeaf('pyv').GetValue()
            p_py[0] = pyv

            # get the x momentum of incoming neutrino:
            pzv = rtree_input.GetBranch('pzv').GetLeaf('pzv').GetValue()
            p_pz[0] = pzv

            # get the value of number of final particles:
            nf = rtree_input.GetBranch('nf').GetLeaf('nf').GetValue()
            nf = int(nf)
            n_pars[0] = nf

            # loop over all i final particles in one event:
            for index in range(nf):
                # get the value of the final PDG and append it to the array:
                pdgf = rtree_input.GetBranch('pdgf').GetLeaf('pdgf').GetValue(index)
                pdgf = int(pdgf)
                pdg[index] = pdgf

                # get the value of the x-momentum in GeV and append it to the array:
                pxf = rtree_input.GetBranch('pxf').GetLeaf('pxf').GetValue(index)
                px[index] = pxf

                # get the value of the x-momentum in GeV and append it to the array:
                pyf = rtree_input.GetBranch('pyf').GetLeaf('pyf').GetValue(index)
                py[index] = pyf

                # get the value of the x-momentum in GeV and append it to the array:
                pzf = rtree_input.GetBranch('pzf').GetLeaf('pzf').GetValue(index)
                pz[index] = pzf

                # get the value of energy of the final particle in GeV and append it to the array:
                ef = rtree_input.GetBranch('Ef').GetLeaf('Ef').GetValue(index)
                energy[index] = ef

                # get the mass of the final particle from its PDG code and append it to the array:
                mass_value = get_mass_from_pdg(pdgf)
                mass[index] = mass_value

            # Fill this one event to the TTree and go to the next event:
            rtree_output.Fill()

    # write TTree to the TFile:
    rfile_output.Write()
    rfile_output.Close()

    return


def get_channels_from_original_genie_file(rootfile_input):
    """
    function to read the 'original' GENIE root-file from Julia check the fractions of the different NC interaction
    channels from variables nfp, nfn, nfpim, nfpip

    :param rootfile_input: path to the original GENIE ROOT-file (for example: gntp.101.gst.root (string)

    :return:
    """
    # load the ROOT file:
    rfile_input = ROOT.TFile(rootfile_input)
    # get the TTree from the TFile:
    rtree_input = rfile_input.Get("gst")

    # Info-me: "gst;13" is a copy of meta data of "gst;14", "gst;14" contains correct data and is read

    # get the number of entries in the ROOT-file:
    number_entries = rtree_input.GetEntries()
    # number_entries = 10000

    """ preallocate all arrays: """
    """ B11 """
    # number of interaction channel: nu + C12 -> B11 + p (integer):
    number_c12_b11_p = 0
    # number of interaction channel: nu + C12 -> B11 + n + pi_plus (integer):
    number_c12_b11_n_piplus = 0
    # number of interaction channel: nu + C12 -> B11 + n + pi_minus + 2*pi_plus (integer):
    number_c12_b11_n_piminus_2piplus = 0
    # number of interaction channel: nu + C12 -> B11 + p + pi_minus + pi_plus (integer):
    number_c12_b11_p_piminus_piplus = 0
    # number of interaction channel: nu + C12 -> B11 + p + 2*pi_minus + 2*pi_plus (integer):
    number_c12_b11_p_2piminus_2piplus = 0
    # number of interaction channel: nu + C12 -> B11 + pi_plus (integer):
    number_c12_b11_piplus = 0

    """ C11 """
    # number of interaction channel: nu + C12 -> C11 + n (integer):
    number_c12_c11_n = 0
    # number of interaction channel: nu + C12 -> C11 + p + pi_minus (integer):
    number_c12_c11_p_piminus = 0
    # number of interaction channel: nu + C12 -> C11 + n + pi_minus + pi_plus (integer):
    number_c12_c11_n_piminus_piplus = 0
    # number of interaction channel: nu + C12 -> C11 + p + 2*pi_minus + pi_plus (integer):
    number_c12_c11_p_2piminus_piplus = 0
    # number of interaction channel: nu + C12 -> C11 + p + 3*pi_minus + 2*pi_plus (integer):
    number_c12_c11_p_3piminus_2piplus = 0
    # number of interaction channel: nu + C12 -> C11 + n + 2*pi_minus + 2*pi_plus (integer):
    number_c12_c11_n_2piminus_2piplus = 0

    """ B10 """
    # number of interaction channel: nu + C12 -> B10 + p + n (integer):
    number_c12_b10_p_n = 0
    # number of interaction channel: nu + C12 -> B10 + 2p + pi_minus (integer):
    number_c12_b10_2p_piminus = 0
    # number of interaction channel: nu + C12 -> B10 + p + n + pi_minus + pi_plus (integer):
    number_c12_b10_p_n_piminus_piplus = 0
    # number of interaction channel: nu + C12 -> B10 + 2n + pi_plus (integer):
    number_c12_b10_2n_piplus = 0
    # number of interaction channel: nu + C12 -> B10 + 2n + pi_minus + 2*pi_plus (integer):
    number_c12_b10_2n_piminus_2piplus = 0
    # number of interaction channel: nu + C12 -> B10 + 2p + 2*pi_minus + pi_plus (integer):
    number_c12_b10_2p_2piminus_piplus = 0
    # number of interaction channel: nu + C12 -> B10 + 2p + 3*pi_minus + 2*pi_plus (integer):
    number_c12_b10_2p_3piminus_2piplus = 0
    # number of interaction channel: nu + C12 -> B10 + p + n + 2*pi_minus + 2*pi_plus (integer):
    number_c12_b10_p_n_2piminus_2piplus = 0

    """ C10 """
    # number of interaction channel: nu + C12 -> C10 + 2n (integer):
    number_c12_c10_2n = 0
    # number of interaction channel: nu + C12 -> C10 + p + n + pi_minus (integer):
    number_c12_c10_p_n_piminus = 0
    # number of interaction channel: nu + C12 -> C10 + p + n + 2*pi_minus + pi_plus (integer):
    number_c12_c10_p_n_2piminus_piplus = 0
    # number of interaction channel: nu + C12 -> C10 + 2n + pi_minus + pi_plus (integer):
    number_c12_c10_2n_piminus_piplus = 0
    # number of interaction channel: nu + C12 -> C10 + 2p + 2*pi_minus (integer):
    number_c12_c10_2p_2piminus = 0

    """ Be10 """
    # number of interaction channel: nu + C12 -> Be10 + 2*p (integer):
    number_c12_be10_2p = 0
    # number of interaction channel: nu + C12 -> Be10 + p + n + pi_plus (integer):
    number_c12_be10_p_n_piplus = 0
    # number of interaction channel: nu + C12 -> Be10 + p + n + pi_minus + 2*pi_plus (integer):
    number_c12_be10_p_n_piminus_2piplus = 0
    # number of interaction channel: nu + C12 -> Be10 + 2*p + pi_minus + pi_plus (integer):
    number_c12_be10_2p_piminus_piplus = 0
    # number of interaction channel: nu + C12 -> Be10 + 2*n + 2*pi_plus (integer):
    number_c12_be10_2n_2piplus = 0
    # number of interaction channel: nu + C12 -> Be10 + p + n + 2*pi_minus + 3*pi_plus (integer):
    number_c12_be10_p_n_2piminus_3piplus = 0
    # number of interaction channel: nu + C12 -> Be10 + 2*p + 2*pi_minus + 2*pi_plus (integer):
    number_c12_be10_2p_2piminus_2piplus = 0
    # number of interaction channel: nu + C12 -> Be10 + 2*p + 3*pi_minus + 3*pi_plus (integer):
    number_c12_be10_2p_3piminus_3piplus = 0

    """ B9 """
    # number of interaction channel: nu + C12 -> B9 + p + 2n (integer):
    number_c12_b9_p_2n = 0
    # number of interaction channel: nu + C12 -> B9 + p + 2n + pi_minus + pi_plus (integer):
    number_c12_b9_p_2n_piminus_piplus = 0
    # number of interaction channel: nu + C12 -> B9 + 2p + n + 3*pi_minus + 2*pi_plus (integer):
    number_c12_b9_2p_n_3piminus_2piplus = 0
    # number of interaction channel: nu + C12 -> B9 + 2p + n + pi_minus (integer):
    number_c12_b9_2p_n_piminus = 0
    # number of interaction channel: nu + C12 -> B9 + 3n + pi_plus (integer):
    number_c12_b9_3n_piplus = 0
    # number of interaction channel: nu + C12 -> B9 + p + 2n + 2*pi_minus + 2*pi_plus (integer):
    number_c12_b9_p_2n_2piminus_2piplus = 0
    # number of interaction channel: nu + C12 -> B9 + 2p + n + 2*pi_minus+ pi_plus (integer):
    number_c12_b9_2p_n_2piminus_piplus = 0

    """ Be9 """
    # number of interaction channel: nu + C12 -> Be9 + 2*p + n (integer):
    number_c12_be9_2p_n = 0
    # number of interaction channel: nu + C12 -> Be9 + p + 2n + pi_plus (integer):
    number_c12_be9_p_2n_piplus = 0
    # number of interaction channel: nu + C12 -> Be9 + 3p + pi_minus (integer):
    number_c12_be9_3p_piminus = 0
    # number of interaction channel: nu + C12 -> Be9 + p + 2n + pi_minus + 2*pi_plus (integer):
    number_c12_be9_p_2n_piminus_2piplus = 0
    # number of interaction channel: nu + C12 -> Be9 + 2p + n + pi_minus + pi_plus (integer):
    number_c12_be9_2p_n_piminus_piplus = 0
    # number of interaction channel: nu + C12 -> Be9 + 2p + n + 3*pi_minus + 3*pi_plus (integer):
    number_c12_be9_2p_n_3piminus_3piplus = 0
    # number of interaction channel: nu + C12 -> Be9 + 2p + n + 2*pi_minus + 2*pi_plus (integer):
    number_c12_be9_2p_n_2piminus_2piplus = 0
    # number of interaction channel: nu + C12 -> Be9 + 3n + 2*pi_plus (integer):
    number_c12_be9_3n_2piplus = 0
    # number of interaction channel: nu + C12 -> Be9 + 3p + 2*pi_minus + pi_plus (integer):
    number_c12_be9_3p_2piminus_piplus = 0

    """ Be8 """
    # number of interaction channel: nu + C12 -> Be8 + 2p + 2n (integer):
    number_c12_be8_2p_2n = 0
    # number of interaction channel: nu + C12 -> Be8 + 3p + n + pi_minus (integer):
    number_c12_be8_3p_n_piminus = 0
    # number of interaction channel: nu + C12 -> Be8 + p + 3n + pi_plus (integer):
    number_c12_be8_p_3n_piplus = 0
    # number of interaction channel: nu + C12 -> Be8 + 2p + 2n + 2*pi_minus + 2*pi_plus (integer):
    number_c12_be8_2p_2n_2piminus_2piplus = 0
    # number of interaction channel: nu + C12 -> Be8 + 4n + 2*pi_plus (integer):
    number_c12_be8_4n_2piplus = 0
    # number of interaction channel: nu + C12 -> Be8 + 2p + 2n + pi_minus * pi_plus (integer):
    number_c12_be8_2p_2n_piminus_piplus = 0
    # number of interaction channel: nu + C12 -> Be8 + 3p + n + 2*pi_minus + pi_plus (integer):
    number_c12_be8_3p_n_2piminus_piplus = 0
    # number of interaction channel: nu + C12 -> Be8 + 4p + 2*pi_minus (integer):
    number_c12_be8_4p_2piminus = 0

    """ C9 """
    # number of interaction channel: nu + C12 -> C9 + p + 2n + pi_minus (integer):
    number_c12_c9_p_2n_piminus = 0
    # number of interaction channel: nu + C12 -> C9 + 3n (integer):
    number_c12_c9_3n = 0
    # number of interaction channel: nu + C12 -> C9 + 2p + n + 2*pi_minus (integer):
    number_c12_c9_2p_n_2piminus = 0
    # number of interaction channel: nu + C12 -> C9 + 3n + 2*pi_minus + 2*pi_plus (integer):
    number_c12_c9_3n_2piminus_2piplus = 0

    """ Be7 """
    # number of interaction channel: nu + C12 -> Be7 + 2p + 3n (integer):
    number_c12_be7_2p_3n = 0
    # number of interaction channel: nu + C12 -> Be7 + p + 4n + pi_plus (integer):
    number_c12_be7_p_4n_piplus = 0
    # number of interaction channel: nu + C12 -> Be7 + 2p + 3n + 2*pi_minus + 2*pi_plus (integer):
    number_c12_be7_2p_3n_2piminus_2piplus = 0
    # number of interaction channel: nu + C12 -> Be7 + 3p + 2n + pi_minus (integer):
    number_c12_be7_3p_2n_piminus = 0
    # number of interaction channel: nu + C12 -> Be7 + 4p + n + 2*pi_minus (integer):
    number_c12_be7_4p_n_2piminus = 0
    # number of interaction channel: nu + C12 -> Be7 + 3p + 2n + 2*pi_minus + pi_plus (integer):
    number_c12_be7_3p_2n_2piminus_piplus = 0

    """ Li6 """
    # number of interaction channel: nu + C12 -> Li6 + 3p + 3n (integer):
    number_c12_li6_3p_3n = 0
    # number of interaction channel: nu + C12 -> Li6 + 2p + 4n + pi_plus (integer):
    number_c12_li6_2p_4n_piplus = 0
    # number of interaction channel: nu + C12 -> Li6 + 5p + n + 2*pi_minus (integer):
    number_c12_li6_5p_n_2piminus = 0
    # number of interaction channel: nu + C12 -> Li6 + 2p + 4n + pi_minus + 2*pi_plus (integer):
    number_c12_li6_2p_4n_piminus_2piplus = 0
    # number of interaction channel: nu + C12 -> Li6 + 4p + 2n + pi_minus (integer):
    number_c12_li6_4p_2n_piminus = 0
    # number of interaction channel: nu + C12 -> Li6 + 3p + 3n + pi_minus + pi_plus (integer):
    number_c12_li6_3p_3n_piminus_piplus = 0

    """ Li8 """
    # number of interaction channel: nu + C12 -> Li8 + 3p + n (integer):
    number_c12_li8_3p_n = 0
    # number of interaction channel: nu + C12 -> Li8 + 4p + pi_minus (integer):
    number_c12_li8_4p_piminus = 0
    # number of interaction channel: nu + C12 -> Li8 + 4p + 2*pi_minus + pi_plus (integer):
    number_c12_li8_4p_2piminus_piplus = 0
    # number of interaction channel: nu + C12 -> Li8 + 2p + 2n + pi_plus (integer):
    number_c12_li8_2p_2n_piplus = 0
    # number of interaction channel: nu + C12 -> Li8 + 3p + n + pi_minus + pi_plus (integer):
    number_c12_li8_3p_n_piminus_piplus = 0

    """ Li7 """
    # number of interaction channel: nu + C12 -> Li7 + 2p + 3n + pi_plus (integer):
    number_c12_li7_2p_3n_piplus = 0
    # number of interaction channel: nu + C12 -> Li7 + 4p + n + pi_minus (integer):
    number_c12_li7_4p_n_piminus = 0
    # number of interaction channel: nu + C12 -> Li7 + 3p + 2n (integer):
    number_c12_li7_3p_2n = 0
    # number of interaction channel: nu + C12 -> Li7 + 3p + 2n + pi_minus + pi_plus (integer):
    number_c12_li7_3p_2n_piminus_piplus = 0
    # number of interaction channel: nu + C12 -> Li7 + 4p + n + 2*pi_minus + pi_plus (integer):
    number_c12_li7_4p_n_2piminus_piplus = 0
    # number of interaction channel: nu + C12 -> Li7 + 2p + 3n + pi_minus + 2*pi_plus (integer):
    number_c12_li7_2p_3n_piminus_2piplus = 0

    """ B8 """
    # number of interaction channel: nu + C12 -> B8 + p + 3n (integer):
    number_c12_b8_p_3n = 0
    # number of interaction channel: nu + C12 -> B8 + p + 3n + pi_minus + pi_plus (integer):
    number_c12_b8_p_3n_piminus_piplus = 0
    # number of interaction channel: nu + C12 -> B8 + 2p + 2n + 2*pi_minus + pi_plus (integer):
    number_c12_b8_2p_2n_2piminus_piplus = 0
    # number of interaction channel: nu + C12 -> B8 + 2p + 2n + pi_minus (integer):
    number_c12_b8_2p_2n_piminus = 0
    # number of interaction channel: nu + C12 -> B8 + 4n + pi_plus (integer):
    number_c12_b8_4n_piplus = 0

    """ Li9 """
    # number of interaction channel: nu + C12 -> Li9 + 2p + n + pi_plus (integer):
    number_c12_li9_2p_n_piplus = 0
    # number of interaction channel: nu + C12 -> Li9 + 3p (integer).
    number_c12_li9_3p = 0
    # number of interaction channel: nu + C12 -> Li9 + 3p + pi_minus + pi_plus (integer):
    number_c12_li9_3p_piminus_piplus = 0
    # number of interaction channel: nu + C12 -> Li9 + 2p + n + pi_minus + 2*pi_plus (integer):
    number_c12_li9_2p_n_piminus_2piplus = 0
    # number of interaction channel: nu + C12 -> Li9 + p + 2n + pi_minus + 3*pi_plus (integer):
    number_c12_li9_p_2n_piminus_3piplus = 0

    """ C8 """
    # number of interaction channel: nu + C12 -> C8 + 4n (integer):
    number_c12_c8_4n = 0

    """ He8 """
    # number of interaction channel: nu + C12 -> He8 + 4p (integer):
    number_c12_he8_4p = 0

    """ B7 """
    # number of interaction channel: nu + C12 -> B7 + p + 4n (integer):
    number_c12_b7_p_4n = 0

    """ He7 """
    # number of interaction channel: nu + C12 -> He7 + 4p + n (integer):
    number_c12_he7_4p_n = 0

    """ H7 """
    # number of interaction channel: nu + C12 -> H7 + 4p + n (integer):
    number_c12_h7_5p = 0

    """ Be6 """
    # number of interaction channel: nu + C12 -> Be6 + 2p + 4n (integer):
    number_c12_be6_2p_4n = 0

    """ Li5 """
    # number of interaction channel: nu + C12 -> Li5 + 3p + 4n (integer):
    number_c12_li5_3p_4n = 0

    """ Li4 """
    # number of interaction channel: nu + C12 -> Li4 + 3p + 5n (integer):
    number_c12_li4_3p_5n = 0

    """ He6 """
    # number of interaction channel: nu + C12 -> He6 + 4p + 2n (integer):
    number_c12_he6_4p_2n = 0

    """ He5 """
    # number of interaction channel: nu + C12 -> He5 + 4p + 3n (integer):
    number_c12_he5_4p_3n = 0

    """ He4 """
    # number of interaction channel: nu + C12 -> He4 + 4p + 4n (integer):
    number_c12_he4_4p_4n = 0

    """ He3 """
    # number of interaction channel: nu + C12 -> He3 + 4p + 5n (integer):
    number_c12_he3_4p_5n = 0

    """ H6 """
    # number of interaction channel: nu + C12 -> H6 + 5p + n (integer):
    number_c12_h6_5p_n = 0

    """ H5 """
    # number of interaction channel: nu + C12 -> H5 + 5p + 2n (integer):
    number_c12_h5_5p_2n = 0

    """ H4 """
    # number of interaction channel: nu + C12 -> H4 + 5p + 3n + ... (integer):
    number_c12_h4_5p_3n = 0

    """ H3 = tritium """
    # number of interaction channel: nu + C12 -> H3 + 5p + 4n + ... (integer):
    number_c12_h3_5p_4n = 0

    """ H2 = deuteron """
    # number of interaction channel: nu + C12 -> H2 + 5p + 5n + ... (integer):
    number_c12_h2_5p_5n = 0

    """ C12 """
    # number of interaction channels: nu + C12 -> nu + C12 + other particles (like pi_minus, pi_plus, kaon_minus,
    # koan_plus and so on):
    number_c12_c12 = 0

    """ no isotope (only protons, neutrons, pions): """
    # number of interaction channels with NO isotope (only proton, neutrons, pions):
    number_c12_noiso = 0

    """ missing interaction channels: not yet implemented channels: """
    # number of interaction channels: nu + C12 -> nu + C12 + ...:
    number_c12_missing = 0

    """ Other targets than C12: """
    # number of channels without C12 as target (integer):
    number_no_c12 = 0
    # number of elastic scattering interactions with protons: nu + p -> nu + p + ... (integer):
    number_es_p = 0
    # number of elastic scattering interactions with electrons: nu + electron -> nu + electron + ... (integer):
    number_es_e = 0
    # number of elastic scattering interactions with O16: nu + O16 -> nu + O16 + ... (integer):
    number_es_o16 = 0
    # number of elastic scattering interactions with N14: nu + N14 -> nu + N14 + ... (integer):
    number_es_n14 = 0
    # number of elastic scattering interactions with S32: nu + S32 -> nu + S32 + ... (integer):
    number_es_s32 = 0

    # number of events for the current interactions (e.g. only NC, or NC + QEL, ...):
    number_events = 0

    """ Read the data from the TTree: """
    # loop over every entry, i.e. every event, in the TTree:
    for event in range(number_entries):

        # get the current event in the TTree:
        rtree_input.GetEntry(event)

        # is it a quasi-elastic scattering event? (0 = no QEL event, 1 = QEL event):
        qel = rtree_input.GetBranch('qel').GetLeaf('qel').GetValue()
        qel = int(qel)

        # is it a NC event? (0 = no NC event, 1 = NC event):
        nc = rtree_input.GetBranch('nc').GetLeaf('nc').GetValue()
        nc = int(nc)

        # get the value of target PDG:
        tgt = rtree_input.GetBranch('tgt').GetLeaf('tgt').GetValue()

        # read only NC and QEL events:
        # if qel == 1 and nc == 1:
        # if nc == 1:
        if nc == 1 and tgt == 1000060120:

            # increase the number of events:
            number_events = number_events + 1

            # get the value of target PDG:
            # tgt = rtree_input.GetBranch('tgt').GetLeaf('tgt').GetValue()
            tgt = int(tgt)

            # get the value of number of final p:
            nfp = rtree_input.GetBranch('nfp').GetLeaf('nfp').GetValue()
            nfp = int(nfp)

            # get the value of number of final n:
            nfn = rtree_input.GetBranch('nfn').GetLeaf('nfn').GetValue()
            nfn = int(nfn)

            # get the value of number of final pi_minus:
            nfpim = rtree_input.GetBranch('nfpim').GetLeaf('nfpim').GetValue()
            nfpim = int(nfpim)

            # get the value of number of final pi_plus:
            nfpip = rtree_input.GetBranch('nfpip').GetLeaf('nfpip').GetValue()
            nfpip = int(nfpip)

            # get the value of number of final Kaon_minus:
            nfkm = rtree_input.GetBranch('nfkm').GetLeaf('nfkm').GetValue()
            nfkm = int(nfkm)

            # get the value of number of final Kaon_plus:
            nfkp = rtree_input.GetBranch('nfkp').GetLeaf('nfkp').GetValue()
            nfkp = int(nfkp)


            # Get the NC interaction channel of the event:
            if tgt == 1000060120:
                # target C12

                # B11:
                if nfp == 1 and nfn == 0 and nfpim == 0 and nfpip == 0:
                    # interaction channel: nu + C12 -> B11 + proton:
                    number_c12_b11_p = number_c12_b11_p + 1

                elif nfp == 0 and nfn == 1 and nfpim == 0 and nfpip == 1:
                    # interaction channel: nu + C12 -> B11 + n + pi_plus:
                    number_c12_b11_n_piplus = number_c12_b11_n_piplus + 1

                elif nfp == 0 and nfn == 1 and nfpim == 1 and nfpip == 2:
                    # interaction channel: nu + C12 -> B11 + n + pi_minus * 2*pi_plus:
                    number_c12_b11_n_piminus_2piplus = number_c12_b11_n_piminus_2piplus + 1

                elif nfp == 1 and nfn == 0 and nfpim == 1 and nfpip == 1:
                    # interaction channel: nu + C12 -> B11 + p + pi_minus + pi_plus:
                    number_c12_b11_p_piminus_piplus = number_c12_b11_p_piminus_piplus + 1

                elif nfp == 1 and nfn == 0 and nfpim == 2 and nfpip == 2:
                    # interaction channel: nu + C12 -> B11 + p + 2*pi_minus + 2*pi_plus:
                    number_c12_b11_p_2piminus_2piplus = number_c12_b11_p_2piminus_2piplus + 1

                elif nfp == 0 and nfn == 0 and nfpim == 0 and nfpip == 1:
                    # interaction channel: nu + C12 -> B11 + pi_plus:
                    number_c12_b11_piplus = number_c12_b11_piplus + 1

                # C11:
                elif nfp == 0 and nfn == 1 and nfpim == 0 and nfpip == 0:
                    # interaction channel: nu + C12 -> C11 + n:
                    number_c12_c11_n = number_c12_c11_n + 1

                elif nfp == 1 and nfn == 0 and nfpim == 1 and nfpip == 0:
                    # interaction channel: nu + C12 -> C11 + p + pi_minus:
                    number_c12_c11_p_piminus = number_c12_c11_p_piminus + 1

                elif nfp == 0 and nfn == 1 and nfpim == 1 and nfpip == 1:
                    # interaction channel: nu + C12 -> C11 + n + pi_minus + pi_plus:
                    number_c12_c11_n_piminus_piplus = number_c12_c11_n_piminus_piplus + 1

                elif nfp == 1 and nfn == 0 and nfpim == 2 and nfpip == 1:
                    # interaction channel: nu + C12 -> C11 + p + 2*pi_minus + pi_plus:
                    number_c12_c11_p_2piminus_piplus = number_c12_c11_p_2piminus_piplus + 1

                elif nfp == 1 and nfn == 0 and nfpim == 3 and nfpip == 2:
                    # interaction channel: nu + C12 -> C11 + p + 3*pi_minus + 2*pi_plus:
                    number_c12_c11_p_3piminus_2piplus = number_c12_c11_p_3piminus_2piplus + 1

                elif nfp == 0 and nfn == 1 and nfpim == 2 and nfpip == 2:
                    # interaction channel: nu + C12 -> C11 + n + 2*pi_minus + 2*pi_plus:
                    number_c12_c11_n_2piminus_2piplus = number_c12_c11_n_2piminus_2piplus + 1

                # B10:
                elif nfp == 1 and nfn == 1 and nfpim == 0 and nfpip == 0:
                    # interaction channel: nu + C12 -> B10 + p + n:
                    number_c12_b10_p_n = number_c12_b10_p_n + 1

                elif nfp == 2 and nfn == 0 and nfpim == 1 and nfpip == 0:
                    # interaction channel: nu + C12 -> B10 + 2*p + pi_minus:
                    number_c12_b10_2p_piminus = number_c12_b10_2p_piminus + 1

                elif nfp == 1 and nfn == 1 and nfpim == 1 and nfpip == 1:
                    # interaction channel: nu + C12 -> B10 + p + n + pi_minus + pi_plus:
                    number_c12_b10_p_n_piminus_piplus = number_c12_b10_p_n_piminus_piplus + 1

                elif nfp == 0 and nfn == 2 and nfpim == 0 and nfpip == 1:
                    # interaction channel: nu + C12 -> B10 + 2*n + pi_plus:
                    number_c12_b10_2n_piplus = number_c12_b10_2n_piplus + 1

                elif nfp == 0 and nfn == 2 and nfpim == 1 and nfpip == 2:
                    # interaction channel: nu + C12 -> B10 + 2*n + pi_minus + 2*pi_plus:
                    number_c12_b10_2n_piminus_2piplus = number_c12_b10_2n_piminus_2piplus + 1

                elif nfp == 2 and nfn == 0 and nfpim == 2 and nfpip == 1:
                    # interaction channel: nu + C12 -> B10 + 2*p + 2*pi_minus + pi_plus:
                    number_c12_b10_2p_2piminus_piplus = number_c12_b10_2p_2piminus_piplus + 1

                elif nfp == 2 and nfn == 0 and nfpim == 3 and nfpip == 2:
                    # interaction channel: nu + C12 -> B10 + 2*p + 3*pi_minus + 2*pi_plus:
                    number_c12_b10_2p_3piminus_2piplus = number_c12_b10_2p_3piminus_2piplus + 1

                elif nfp == 1 and nfn == 1 and nfpim == 2 and nfpip == 2:
                    # interaction channel: nu + C12 -> B10 + p + n + 2*pi_minus + 2*pi_plus:
                    number_c12_b10_p_n_2piminus_2piplus = number_c12_b10_p_n_2piminus_2piplus + 1

                # C10:
                elif nfp == 0 and nfn == 2 and nfpim == 0 and nfpip == 0:
                    # interaction channel: nu + C12 -> C10 + 2n:
                    number_c12_c10_2n = number_c12_c10_2n + 1

                elif nfp == 1 and nfn == 1 and nfpim == 1 and nfpip == 0:
                    # interaction channel: nu + C12 -> C10 + p + n + pi_minus:
                    number_c12_c10_p_n_piminus = number_c12_c10_p_n_piminus + 1

                elif nfp == 1 and nfn == 1 and nfpim == 2 and nfpip == 1:
                    # interaction channel: nu + C12 -> C10 + p + n + 2*pi_minus + pi_plus:
                    number_c12_c10_p_n_2piminus_piplus = number_c12_c10_p_n_2piminus_piplus + 1

                elif nfp == 0 and nfn == 2 and nfpim == 1 and nfpip == 1:
                    # interaction channel: nu + C12 -> C10 + 2*n + pi_minus + pi_plus:
                    number_c12_c10_2n_piminus_piplus = number_c12_c10_2n_piminus_piplus + 1

                elif nfp == 2 and nfn == 0 and nfpim == 2 and nfpip == 0:
                    # interaction channel: nu + C12 -> C10 + 2*p + 2*pi_minus:
                    number_c12_c10_2p_2piminus = number_c12_c10_2p_2piminus + 1

                # Be10:
                elif nfp == 2 and nfn == 0 and nfpim == 0 and nfpip == 0:
                    # interaction channel: nu + C12 -> Be10 + 2*p:
                    number_c12_be10_2p = number_c12_be10_2p + 1

                elif nfp == 1 and nfn == 1 and nfpim == 0 and nfpip == 1:
                    # interaction channel: nu + C12 -> Be10 + p + n + pi_plus:
                    number_c12_be10_p_n_piplus = number_c12_be10_p_n_piplus + 1

                elif nfp == 1 and nfn == 1 and nfpim == 1 and nfpip == 2:
                    # interaction channel: nu + C12 -> Be10 + p + n + pi_minus + 2*pi_plus:
                    number_c12_be10_p_n_piminus_2piplus = number_c12_be10_p_n_piminus_2piplus + 1

                elif nfp == 2 and nfn == 0 and nfpim == 1 and nfpip == 1:
                    # interaction channel: nu + C12 -> Be10 + 2*p + pi_minus + pi_plus:
                    number_c12_be10_2p_piminus_piplus = number_c12_be10_2p_piminus_piplus + 1

                elif nfp == 0 and nfn == 2 and nfpim == 0 and nfpip == 2:
                    # interaction channel: nu + C12 -> Be10 + 2n + 2*pi_plus:
                    number_c12_be10_2n_2piplus = number_c12_be10_2n_2piplus + 1

                elif nfp == 1 and nfn == 1 and nfpim == 2 and nfpip == 3:
                    # interaction channel: nu + C12 -> Be10 + p + n + 2*pi_minus + 3*pi_plus:
                    number_c12_be10_p_n_2piminus_3piplus = number_c12_be10_p_n_2piminus_3piplus + 1

                elif nfp == 2 and nfn == 0 and nfpim == 2 and nfpip == 2:
                    # interaction channel: nu + C12 -> Be10 + 2p + 2*pi_minus + 2*pi_plus:
                    number_c12_be10_2p_2piminus_2piplus = number_c12_be10_2p_2piminus_2piplus + 1

                elif nfp == 2 and nfn == 0 and nfpim == 3 and nfpip == 3:
                    # interaction channel: nu + C12 -> Be10 + 2p + 3*pi_minus + 3*pi_plus:
                    number_c12_be10_2p_3piminus_3piplus = number_c12_be10_2p_3piminus_3piplus + 1

                # B9:
                elif nfp == 1 and nfn == 2 and nfpim == 0 and nfpip == 0:
                    # interaction channel: nu + C12 -> B9 + p + 2*n:
                    number_c12_b9_p_2n = number_c12_b9_p_2n + 1

                elif nfp == 1 and nfn == 2 and nfpim == 1 and nfpip == 1:
                    # interaction channel: nu + C12 -> B9 + p + 2n + pi_minus + pi_plus:
                    number_c12_b9_p_2n_piminus_piplus = number_c12_b9_p_2n_piminus_piplus + 1

                elif nfp == 2 and nfn == 1 and nfpim == 3 and nfpip == 2:
                    # interaction channel: nu + C12 -> B9 + 2p + n + 3*pi_minus + 2*pi_plus:
                    number_c12_b9_2p_n_3piminus_2piplus = number_c12_b9_2p_n_3piminus_2piplus + 1

                elif nfp == 2 and nfn == 1 and nfpim == 1 and nfpip == 0:
                    # interaction channel: nu + C12 -> B9 + 2p + n + pi_minus:
                    number_c12_b9_2p_n_piminus = number_c12_b9_2p_n_piminus + 1

                elif nfp == 0 and nfn == 3 and nfpim == 0 and nfpip == 1:
                    # interaction channel: nu + C12 -> B9 + 3n + pi_plus:
                    number_c12_b9_3n_piplus = number_c12_b9_3n_piplus + 1

                elif nfp == 1 and nfn == 2 and nfpim == 2 and nfpip == 2:
                    # interaction channel: nu + C12 -> B9 + p + 2n + 2*pi_minus + 2*pi_plus:
                    number_c12_b9_p_2n_2piminus_2piplus = number_c12_b9_p_2n_2piminus_2piplus + 1

                elif nfp == 2 and nfn == 1 and nfpim == 2 and nfpip == 1:
                    # interaction channel: nu + C12 -> B9 + 2p + n + 2*pi_minus + pi_plus:
                    number_c12_b9_2p_n_2piminus_piplus = number_c12_b9_2p_n_2piminus_piplus + 1

                # Be9:
                elif nfp == 2 and nfn == 1 and nfpim == 0 and nfpip == 0:
                    # interaction channel: nu + C12 -> Be9 + 2*p + n:
                    number_c12_be9_2p_n = number_c12_be9_2p_n + 1

                elif nfp == 1 and nfn == 2 and nfpim == 0 and nfpip == 1:
                    # interaction channel: nu + C12 -> Be9 + p + 2*n + pi_plus:
                    number_c12_be9_p_2n_piplus = number_c12_be9_p_2n_piplus + 1

                elif nfp == 3 and nfn == 0 and nfpim == 1 and nfpip == 0:
                    # interaction channel: nu + C12 -> Be9 + 3p + pi_minus:
                    number_c12_be9_3p_piminus = number_c12_be9_3p_piminus + 1

                elif nfp == 1 and nfn == 2 and nfpim == 1 and nfpip == 2:
                    # interaction channel: nu + C12 -> Be9 + p + 2*n + pi_minus + 2*pi_plus:
                    number_c12_be9_p_2n_piminus_2piplus = number_c12_be9_p_2n_piminus_2piplus + 1

                elif nfp == 2 and nfn == 1 and nfpim == 1 and nfpip == 1:
                    # interaction channel: nu + C12 -> Be9 + 2*p + n + pi_minus + pi_plus:
                    number_c12_be9_2p_n_piminus_piplus = number_c12_be9_2p_n_piminus_piplus + 1

                elif nfp == 2 and nfn == 1 and nfpim == 3 and nfpip == 3:
                    # interaction channel: nu + C12 -> Be9 + 2*p + n + 3*pi_minus + 3*pi_plus:
                    number_c12_be9_2p_n_3piminus_3piplus = number_c12_be9_2p_n_3piminus_3piplus + 1

                elif nfp == 2 and nfn == 1 and nfpim == 2 and nfpip == 2:
                    # interaction channel: nu + C12 -> Be9 + 2*p + n + 2*pi_minus + 2*pi_plus:
                    number_c12_be9_2p_n_2piminus_2piplus = number_c12_be9_2p_n_2piminus_2piplus + 1

                elif nfp == 0 and nfn == 3 and nfpim == 0 and nfpip == 2:
                    # interaction channel: nu + C12 -> Be9 + 3*n + 2*pi_plus:
                    number_c12_be9_3n_2piplus = number_c12_be9_3n_2piplus + 1

                elif nfp == 3 and nfn == 0 and nfpim == 2 and nfpip == 1:
                    # interaction channel: nu + C12 -> Be9 + 3*p + 2*pi_minus + pi_plus:
                    number_c12_be9_3p_2piminus_piplus = number_c12_be9_3p_2piminus_piplus + 1

                # C9:
                elif nfp == 1 and nfn == 2 and nfpim == 1 and nfpip == 0:
                    # interaction channel: nu + C12 -> C9 + p + 2*n + pi_minus:
                    number_c12_c9_p_2n_piminus = number_c12_c9_p_2n_piminus + 1

                elif nfp == 0 and nfn == 3 and nfpim == 0 and nfpip == 0:
                    # interaction channel: nu + C12 -> C9 + 3*n:
                    number_c12_c9_3n = number_c12_c9_3n + 1

                elif nfp == 2 and nfn == 1 and nfpim == 2 and nfpip == 0:
                    # interaction channel: nu + C12 -> C9 + 2*p + n + 2*pi_minus:
                    number_c12_c9_2p_n_2piminus = number_c12_c9_2p_n_2piminus + 1

                elif nfp == 0 and nfn == 3 and nfpim == 2 and nfpip == 2:
                    # interaction channel: nu + C12 -> C9 + 3*n + 2*pi_minus + 2*pi_plus:
                    number_c12_c9_3n_2piminus_2piplus = number_c12_c9_3n_2piminus_2piplus + 1

                # Li9:
                elif nfp == 2 and nfn == 1 and nfpim == 0 and nfpip == 1:
                    # interaction channel: nu + C12 -> Li9 + 2*p + n + pi_plus:
                    number_c12_li9_2p_n_piplus = number_c12_li9_2p_n_piplus + 1

                elif nfp == 3 and nfn == 0 and nfpim == 0 and nfpip == 0:
                    # interaction channel: nu + C12 -> Li9 + 3*p:
                    number_c12_li9_3p = number_c12_li9_3p + 1

                elif nfp == 3 and nfn == 0 and nfpim == 1 and nfpip == 1:
                    # interaction channel: nu + C12 -> Li9 + 3*p + pi_minus + pi_plus:
                    number_c12_li9_3p_piminus_piplus = number_c12_li9_3p_piminus_piplus + 1

                elif nfp == 2 and nfn == 1 and nfpim == 1 and nfpip == 2:
                    # interaction channel: nu + C12 -> Li9 + 2*p + n + pi_minus + 2*pi_plus:
                    number_c12_li9_2p_n_piminus_2piplus = number_c12_li9_2p_n_piminus_2piplus + 1

                elif nfp == 1 and nfn == 2 and nfpim == 1 and nfpip == 3:
                    # interaction channel: nu + C12 -> Li9 + p + 2*n + pi_minus + 3*pi_plus:
                    number_c12_li9_p_2n_piminus_3piplus = number_c12_li9_p_2n_piminus_3piplus + 1

                # C8:
                elif nfp == 0 and nfn == 4 and nfpim == 0 and nfpip == 0:
                    # interaction channel: nu + C12 -> C8 + 4n:
                    number_c12_c8_4n = number_c12_c8_4n + 1

                # B8:
                elif nfp == 1 and nfn == 3 and nfpim == 0 and nfpip == 0:
                    # interaction channel: nu + C12 -> B8 + p + 3*n:
                    number_c12_b8_p_3n = number_c12_b8_p_3n + 1

                elif nfp == 1 and nfn == 3 and nfpim == 1 and nfpip == 1:
                    # interaction channel: nu + C12 -> B8 + p + 3*n + pi_minus + pi_plus:
                    number_c12_b8_p_3n_piminus_piplus = number_c12_b8_p_3n_piminus_piplus + 1

                elif nfp == 2 and nfn == 2 and nfpim == 2 and nfpip == 1:
                    # interaction channel: nu + C12 -> B8 + 2*p + 2*n + 2*pi_minus + pi_plus:
                    number_c12_b8_2p_2n_2piminus_piplus = number_c12_b8_2p_2n_2piminus_piplus + 1

                elif nfp == 2 and nfn == 2 and nfpim == 1 and nfpip == 0:
                    # interaction channel: nu + C12 -> B8 + 2*p + 2*n + pi_minus:
                    number_c12_b8_2p_2n_piminus = number_c12_b8_2p_2n_piminus + 1

                elif nfp == 0 and nfn == 4 and nfpim == 0 and nfpip == 1:
                    # interaction channel: nu + C12 -> B8 + 4*n + pi_plus:
                    number_c12_b8_4n_piplus = number_c12_b8_4n_piplus + 1

                # Be8:
                elif nfp == 2 and nfn == 2 and nfpim == 0 and nfpip == 0:
                    # interaction channel: nu + C12 -> Be8 + 2*p + 2*n:
                    number_c12_be8_2p_2n = number_c12_be8_2p_2n + 1

                elif nfp == 3 and nfn == 1 and nfpim == 1 and nfpip == 0:
                    # interaction channel: nu + C12 -> Be8 + 3*p + n + pi_minus:
                    number_c12_be8_3p_n_piminus = number_c12_be8_3p_n_piminus + 1

                elif nfp == 1 and nfn == 3 and nfpim == 0 and nfpip == 1:
                    # interaction channel: nu + C12 -> Be8 + p + 3*n + pi_plus:
                    number_c12_be8_p_3n_piplus = number_c12_be8_p_3n_piplus + 1

                elif nfp == 2 and nfn == 2 and nfpim == 2 and nfpip == 2:
                    # interaction channel: nu + C12 -> Be8 + 2*p + 2*n + 2*pi_minus + 2*pi_plus:
                    number_c12_be8_2p_2n_2piminus_2piplus = number_c12_be8_2p_2n_2piminus_2piplus + 1

                elif nfp == 0 and nfn == 4 and nfpim == 0 and nfpip == 2:
                    # interaction channel: nu + C12 -> Be8 + 4*n + 2*pi_plus:
                    number_c12_be8_4n_2piplus = number_c12_be8_4n_2piplus + 1

                elif nfp == 2 and nfn == 2 and nfpim == 1 and nfpip == 1:
                    # interaction channel: nu + C12 -> Be8 + 2*p + 2*n + pi_minus + pi_plus:
                    number_c12_be8_2p_2n_piminus_piplus = number_c12_be8_2p_2n_piminus_piplus + 1

                elif nfp == 3 and nfn == 1 and nfpim == 2 and nfpip == 1:
                    # interaction channel: nu + C12 -> Be8 + 3*p + n + 2*pi_minus + pi_plus:
                    number_c12_be8_3p_n_2piminus_piplus = number_c12_be8_3p_n_2piminus_piplus + 1

                elif nfp == 4 and nfn == 0 and nfpim == 2 and nfpip == 0:
                    # interaction channel: nu + C12 -> Be8 + 4*p + 2*pi_minus:
                    number_c12_be8_4p_2piminus = number_c12_be8_4p_2piminus + 1

                # Li8:
                elif nfp == 3 and nfn == 1 and nfpim == 0 and nfpip == 0:
                    # interaction channel: nu + C12 -> Li8 + 3*p + n:
                    number_c12_li8_3p_n = number_c12_li8_3p_n + 1

                elif nfp == 4 and nfn == 0 and nfpim == 1 and nfpip == 0:
                    # interaction channel: nu + C12 -> Li8 + 4*p + pi_minus:
                    number_c12_li8_4p_piminus = number_c12_li8_4p_piminus + 1

                elif nfp == 4 and nfn == 0 and nfpim == 2 and nfpip == 1:
                    # interaction channel: nu + C12 -> Li8 + 4*p + 2*pi_minus + pi_plus:
                    number_c12_li8_4p_2piminus_piplus = number_c12_li8_4p_2piminus_piplus + 1

                elif nfp == 2 and nfn == 2 and nfpim == 0 and nfpip == 1:
                    # interaction channel: nu + C12 -> Li8 + 2*p + 2*n + pi_plus:
                    number_c12_li8_2p_2n_piplus = number_c12_li8_2p_2n_piplus + 1

                elif nfp == 3 and nfn == 1 and nfpim == 1 and nfpip == 1:
                    # interaction channel: nu + C12 -> Li8 + 3*p + n + pi_minus + pi_plus:
                    number_c12_li8_3p_n_piminus_piplus = number_c12_li8_3p_n_piminus_piplus + 1

                # He8:
                elif nfp == 4 and nfn == 0 and nfpim == 0 and nfpip == 0:
                    # interaction channel: nu + C12 -> He8 + 4p:
                    number_c12_he8_4p = number_c12_he8_4p + 1

                # B7:
                elif nfp == 1 and nfn == 4 and nfpim == 0 and nfpip == 0:
                    # interaction channel: nu + C12 -> B7 + p + 4n:
                    number_c12_b7_p_4n = number_c12_b7_p_4n + 1

                # Be7:
                elif nfp == 2 and nfn == 3 and nfpim == 0 and nfpip == 0:
                    # interaction channel: nu + C12 -> Be7 + 2*p + 3*n:
                    number_c12_be7_2p_3n = number_c12_be7_2p_3n + 1

                elif nfp == 1 and nfn == 4 and nfpim == 0 and nfpip == 1:
                    # interaction channel: nu + C12 -> Be7 + p + 4*n + pi_plus:
                    number_c12_be7_p_4n_piplus = number_c12_be7_p_4n_piplus + 1

                elif nfp == 2 and nfn == 3 and nfpim == 2 and nfpip == 2:
                    # interaction channel: nu + C12 -> Be7 + 2*p + 3*n + 2*pi_minus + 2*pi_plus:
                    number_c12_be7_2p_3n_2piminus_2piplus = number_c12_be7_2p_3n_2piminus_2piplus + 1

                elif nfp == 3 and nfn == 2 and nfpim == 1 and nfpip == 0:
                    # interaction channel: nu + C12 -> Be7 + 3*p + 2*n + pi_minus:
                    number_c12_be7_3p_2n_piminus = number_c12_be7_3p_2n_piminus + 1

                elif nfp == 4 and nfn == 1 and nfpim == 2 and nfpip == 0:
                    # interaction channel: nu + C12 -> Be7 + 4*p + n + 2*pi_minus:
                    number_c12_be7_4p_n_2piminus = number_c12_be7_4p_n_2piminus + 1

                elif nfp == 3 and nfn == 2 and nfpim == 2 and nfpip == 1:
                    # interaction channel: nu + C12 -> Be7 + 3*p + 2*n + 2*pi_minus + pi_plus:
                    number_c12_be7_3p_2n_2piminus_piplus = number_c12_be7_3p_2n_2piminus_piplus + 1

                # Li7:
                elif nfp == 2 and nfn == 3 and nfpim == 0 and nfpip == 1:
                    # interaction channel: nu + C12 -> Li7 + 2*p + 3*n + pi_plus:
                    number_c12_li7_2p_3n_piplus = number_c12_li7_2p_3n_piplus + 1

                elif nfp == 4 and nfn == 1 and nfpim == 1 and nfpip == 0:
                    # interaction channel: nu + C12 -> Li7 + 4*p + n + pi_minus:
                    number_c12_li7_4p_n_piminus = number_c12_li7_4p_n_piminus + 1

                elif nfp == 3 and nfn == 2 and nfpim == 0 and nfpip == 0:
                    # interaction channel: nu + C12 -> Li7 + 3*p + 2*n:
                    number_c12_li7_3p_2n = number_c12_li7_3p_2n + 1

                elif nfp == 3 and nfn == 2 and nfpim == 1 and nfpip == 1:
                    # interaction channel: nu + C12 -> Li7 + 3*p + 2*n + pi_minus + pi_plus:
                    number_c12_li7_3p_2n_piminus_piplus = number_c12_li7_3p_2n_piminus_piplus + 1

                elif nfp == 4 and nfn == 1 and nfpim == 2 and nfpip == 1:
                    # interaction channel: nu + C12 -> Li7 + 4*p + n + 2*pi_minus + pi_plus:
                    number_c12_li7_4p_n_2piminus_piplus = number_c12_li7_4p_n_2piminus_piplus + 1

                elif nfp == 2 and nfn == 3 and nfpim == 1 and nfpip == 2:
                    # interaction channel: nu + C12 -> Li7 + 2*p + 3*n + pi_minus + 2*pi_plus:
                    number_c12_li7_2p_3n_piminus_2piplus = number_c12_li7_2p_3n_piminus_2piplus + 1

                # He7:
                elif nfp == 4 and nfn == 1 and nfpim == 0 and nfpip == 0:
                    # interaction channel: nu + C12 -> He7 + 4p + n:
                    number_c12_he7_4p_n = number_c12_he7_4p_n + 1

                # H7:
                elif nfp == 5 and nfn == 0 and nfpim == 0 and nfpip == 0:
                    # interaction channel: nu + C12 -> H7 + 5p:
                    number_c12_h7_5p = number_c12_h7_5p + 1

                # Be6:
                elif nfp == 2 and nfn == 4 and nfpim == 0 and nfpip == 0:
                    # interaction channel: nu + C12 -> Be6 + 2p + 4n:
                    number_c12_be6_2p_4n = number_c12_be6_2p_4n + 1

                # Li6:
                elif nfp == 3 and nfn == 3 and nfpim == 0 and nfpip == 0:
                    # interaction channel: nu + C12 -> Li6 + 3*p + 3*n:
                    number_c12_li6_3p_3n = number_c12_li6_3p_3n + 1

                elif nfp == 2 and nfn == 4 and nfpim == 0 and nfpip == 1:
                    # interaction channel: nu + C12 -> Li6 + 2*p + 4*n + pi_plus:
                    number_c12_li6_2p_4n_piplus = number_c12_li6_2p_4n_piplus + 1

                elif nfp == 5 and nfn == 1 and nfpim == 2 and nfpip == 0:
                    # interaction channel: nu + C12 -> Li6 + 5*p + n + 2*pi_minus:
                    number_c12_li6_5p_n_2piminus = number_c12_li6_5p_n_2piminus + 1

                elif nfp == 2 and nfn == 4 and nfpim == 1 and nfpip == 2:
                    # interaction channel: nu + C12 -> Li6 + 2*p + 4*n + pi_minus + 2*pi_plus:
                    number_c12_li6_2p_4n_piminus_2piplus = number_c12_li6_2p_4n_piminus_2piplus + 1

                elif nfp == 4 and nfn == 2 and nfpim == 1 and nfpip == 0:
                    # interaction channel: nu + C12 -> Li6 + 4*p + 2*n + pi_minus:
                    number_c12_li6_4p_2n_piminus = number_c12_li6_4p_2n_piminus + 1

                elif nfp == 3 and nfn == 3 and nfpim == 1 and nfpip == 1:
                    # interaction channel: nu + C12 -> Li6 + 3*p + 3*n + pi_minus + pi_plus:
                    number_c12_li6_3p_3n_piminus_piplus = number_c12_li6_3p_3n_piminus_piplus + 1

                # Li5:
                elif nfp == 3 and nfn == 4 and nfpim == 0 and nfpip == 0:
                    # interaction channel: nu + C12 -> Li5 + 3p + 4n:
                    number_c12_li5_3p_4n = number_c12_li5_3p_4n + 1

                # Li4:
                elif nfp == 3 and nfn == 5 and nfpim == 0 and nfpip == 0:
                    # interaction channel: nu + C12 -> Li4 + 3p + 5n:
                    number_c12_li4_3p_5n = number_c12_li4_3p_5n + 1

                # He6:
                elif nfp == 4 and nfn == 2 and nfpim == 0 and nfpip == 0:
                    # interaction channel: nu + C12 -> He6 + 4p + 2n:
                    number_c12_he6_4p_2n = number_c12_he6_4p_2n + 1

                # He5:
                elif nfp == 4 and nfn == 3 and nfpim == 0 and nfpip == 0:
                    # interaction channel: nu + C12 -> He5 + 4p + 3n:
                    number_c12_he5_4p_3n = number_c12_he5_4p_3n + 1

                # He4:
                elif nfp == 4 and nfn == 4 and nfpim == 0 and nfpip == 0:
                    # interaction channel: nu + C12 -> He4 + 4p + 4n:
                    number_c12_he4_4p_4n = number_c12_he4_4p_4n + 1

                # He3:
                elif nfp == 4 and nfn == 5 and nfpim == 0 and nfpip == 0:
                    # interaction channel: nu + C12 -> He3 + 4p + 5n:
                    number_c12_he3_4p_5n = number_c12_he3_4p_5n + 1

                # H6:
                elif nfp == 5 and nfn == 1 and nfpim == 0 and nfpip == 0:
                    # interaction channel: nu + C12 -> H6 + 5p + n:
                    number_c12_h6_5p_n = number_c12_h6_5p_n + 1

                # H5:
                elif nfp == 5 and nfn == 2 and nfpim == 0 and nfpip == 0:
                    # interaction channel: nu + C12 -> H5 + 5p + 2n:
                    number_c12_h5_5p_2n = number_c12_h5_5p_2n + 1

                # H4:
                elif nfp == 5 and nfn == 3:
                    # interaction channel: nu + C12 -> H4 + 5p + 3n:
                    number_c12_h4_5p_3n = number_c12_h4_5p_3n + 1

                # H3:
                elif nfp == 5 and nfn == 4:
                    # interaction channel: nu + C12 -> H3 + 5p + 4n:
                    number_c12_h3_5p_4n = number_c12_h3_5p_4n + 1

                # H2:
                elif nfp == 5 and nfn == 5:
                    # interaction channel: nu + C12 -> H2 + 5p + 5n:
                    number_c12_h2_5p_5n = number_c12_h2_5p_5n + 1

                # C12:
                elif nfp == 0 and nfn == 0 and nfpim == 0 and nfpip == 0:
                    # interaction channel: nu + C12 -> nu + C12:
                    number_c12_c12 = number_c12_c12 + 1

                # no isotope:
                elif (nfp == 6 and nfn == 6 and nfpim == 0 and nfpip == 0) or \
                        (nfp == 5 and nfn == 6 and nfpim == 0 and nfpip == 0) or \
                        (nfp == 6 and nfn == 5 and nfpim == 0 and nfpip == 0) or \
                        (nfp == 6 and nfn == 4 and nfpim == 0 and nfpip == 0) or \
                        (nfp == 4 and nfn == 6 and nfpim == 0 and nfpip == 0):
                    # interaction channel: nu + C12 -> nu + 6p + 6n:
                    number_c12_noiso = number_c12_noiso + 1

                else:
                    # interaction channels, that are not covered above, and channels with Kaon or Sigma-Baryon:
                    number_c12_missing = number_c12_missing + 1
                    print("new interaction channel with nu + C12 ->: nfp={0:d}, nfn={1:d}, nfpim={2:d}, nfpip={3:d}, "
                          "nfkm={4:d}, nfkp={5:d}".format(nfp, nfn, nfpim, nfpip, nfkm, nfkp))

            else:
                # NC channel with other target than C12:
                number_no_c12 = number_no_c12 + 1

                if tgt == 2212:
                    # proton as target: Es interaction: nu + proton -> nu + proton
                    number_es_p = number_es_p + 1

                elif tgt == 11:
                    # electron as target: ES interaction: nu + electron -> nu + electron (maybe also pi_zero or gammas):
                    number_es_e = number_es_e + 1

                elif tgt == 1000080160:
                    # O16 as target: ES interaction: nu + O16 -> nu + O16 (maybe also pi_zero or gammas):
                    number_es_o16 = number_es_o16 + 1

                elif tgt == 1000070140:
                    # N14 as target: ES interaction: nu + N14 -> nu + N14 (maybe also pi_zero or gammas):
                    number_es_n14 = number_es_n14 + 1

                elif tgt == 1000160320:
                    # S32 as target: ES interaction: nu + S32 -> nu + S32 (maybe also pi_zero or gammas):
                    number_es_s32 = number_es_s32 + 1

                else:
                    print("other target than C12, p, e, N14, O16, S32: tgt = {0:d}".format(tgt))


    """ calculate the fraction of the different NC interaction channels in PERCENT (float): """
    # B11:
    frac_c12_b11_p = float(number_c12_b11_p) / float(number_events) * 100
    frac_c12_b11_n_piplus = float(number_c12_b11_n_piplus) / float(number_events) * 100
    frac_c12_b11_n_piminus_2piplus = float(number_c12_b11_n_piminus_2piplus) / float(number_events) * 100
    frac_c12_b11_p_piminus_piplus = float(number_c12_b11_p_piminus_piplus) / float(number_events) * 100
    frac_c12_b11_p_2piminus_2piplus = float(number_c12_b11_p_2piminus_2piplus) / float(number_events) * 100
    frac_c12_b11_piplus = float(number_c12_b11_piplus) / float(number_events) * 100

    # C11:
    frac_c12_c11_n = float(number_c12_c11_n) / float(number_events) * 100
    frac_c12_c11_p_piminus = float(number_c12_c11_p_piminus) / float(number_events) * 100
    frac_c12_c11_n_piminus_piplus = float(number_c12_c11_n_piminus_piplus) / float(number_events) * 100
    frac_c12_c11_p_2piminus_piplus = float(number_c12_c11_p_2piminus_piplus) / float(number_events) * 100
    frac_c12_c11_p_3piminus_2piplus = float(number_c12_c11_p_3piminus_2piplus) / float(number_events) * 100
    frac_c12_c11_n_2piminus_2piplus = float(number_c12_c11_n_2piminus_2piplus) / float(number_events) * 100

    # B10:
    frac_c12_b10_p_n = float(number_c12_b10_p_n) / float(number_events) * 100
    frac_c12_b10_2p_piminus = float(number_c12_b10_2p_piminus) / float(number_events) * 100
    frac_c12_b10_p_n_piminus_piplus = float(number_c12_b10_p_n_piminus_piplus) / float(number_events) * 100
    frac_c12_b10_2n_piplus = float(number_c12_b10_2n_piplus) / float(number_events) * 100
    frac_c12_b10_2n_piminus_2piplus = float(number_c12_b10_2n_piminus_2piplus) / float(number_events) * 100
    frac_c12_b10_2p_2piminus_piplus = float(number_c12_b10_2p_2piminus_piplus) / float(number_events) * 100
    frac_c12_b10_2p_3piminus_2piplus = float(number_c12_b10_2p_3piminus_2piplus) / float(number_events) * 100
    frac_c12_b10_p_n_2piminus_2piplus = float(number_c12_b10_p_n_2piminus_2piplus) / float(number_events) * 100

    # C10:
    frac_c12_c10_2n = float(number_c12_c10_2n) / float(number_events) * 100
    frac_c12_c10_p_n_piminus = float(number_c12_c10_p_n_piminus) / float(number_events) * 100
    frac_c12_c10_p_n_2piminus_piplus = float(number_c12_c10_p_n_2piminus_piplus) / float(number_events) * 100
    frac_c12_c10_2n_piminus_piplus = float(number_c12_c10_2n_piminus_piplus) / float(number_events) * 100
    frac_c12_c10_2p_2piminus = float(number_c12_c10_2p_2piminus) / float(number_events) * 100

    # Be10:
    frac_c12_be10_2p = float(number_c12_be10_2p) / float(number_events) * 100
    frac_c12_be10_p_n_piplus = float(number_c12_be10_p_n_piplus) / float(number_events) * 100
    frac_c12_be10_p_n_piminus_2piplus = float(number_c12_be10_p_n_piminus_2piplus) / float(number_events) * 100
    frac_c12_be10_2p_piminus_piplus = float(number_c12_be10_2p_piminus_piplus) / float(number_events) * 100
    frac_c12_be10_2n_2piplus = float(number_c12_be10_2n_2piplus) / float(number_events) * 100
    frac_c12_be10_p_n_2piminus_3piplus = float(number_c12_be10_p_n_2piminus_3piplus) / float(number_events) * 100
    frac_c12_be10_2p_2piminus_2piplus = float(number_c12_be10_2p_2piminus_2piplus) / float(number_events) * 100
    frac_c12_be10_2p_3piminus_3piplus = float(number_c12_be10_2p_3piminus_3piplus) / float(number_events) * 100

    # B9:
    frac_c12_b9_p_2n = float(number_c12_b9_p_2n) / float(number_events) * 100
    frac_c12_b9_p_2n_piminus_piplus = float(number_c12_b9_p_2n_piminus_piplus) / float(number_events) * 100
    frac_c12_b9_2p_n_3piminus_2piplus = float(number_c12_b9_2p_n_3piminus_2piplus) / float(number_events) * 100
    frac_c12_b9_2p_n_piminus = float(number_c12_b9_2p_n_piminus) / float(number_events) * 100
    frac_c12_b9_3n_piplus = float(number_c12_b9_3n_piplus) / float(number_events) * 100
    frac_c12_b9_p_2n_2piminus_2piplus = float(number_c12_b9_p_2n_2piminus_2piplus) / float(number_events) * 100
    frac_c12_b9_2p_n_2piminus_piplus = float(number_c12_b9_2p_n_2piminus_piplus) / float(number_events) * 100

    # Be9:
    frac_c12_be9_2p_n = float(number_c12_be9_2p_n) / float(number_events) * 100
    frac_c12_be9_p_2n_piplus = float(number_c12_be9_p_2n_piplus) / float(number_events) * 100
    frac_c12_be9_3p_piminus = float(number_c12_be9_3p_piminus) / float(number_events) * 100
    frac_c12_be9_p_2n_piminus_2piplus = float(number_c12_be9_p_2n_piminus_2piplus) / float(number_events) * 100
    frac_c12_be9_2p_n_piminus_piplus = float(number_c12_be9_2p_n_piminus_piplus) / float(number_events) * 100
    frac_c12_be9_2p_n_3piminus_3piplus = float(number_c12_be9_2p_n_3piminus_3piplus) / float(number_events) * 100
    frac_c12_be9_2p_n_2piminus_2piplus = float(number_c12_be9_2p_n_2piminus_2piplus) / float(number_events) * 100
    frac_c12_be9_3n_2piplus = float(number_c12_be9_3n_2piplus) / float(number_events) * 100
    frac_c12_be9_3p_2piminus_piplus = float(number_c12_be9_3p_2piminus_piplus) / float(number_events) * 100

    # Be8:
    frac_c12_be8_2p_2n = float(number_c12_be8_2p_2n) / float(number_events) * 100
    frac_c12_be8_3p_n_piminus = float(number_c12_be8_3p_n_piminus) / float(number_events) * 100
    frac_c12_be8_p_3n_piplus = float(number_c12_be8_p_3n_piplus) / float(number_events) * 100
    frac_c12_be8_2p_2n_2piminus_2piplus = float(number_c12_be8_2p_2n_2piminus_2piplus) / float(number_events) * 100
    frac_c12_be8_4n_2piplus = float(number_c12_be8_4n_2piplus) / float(number_events) * 100
    frac_c12_be8_2p_2n_piminus_piplus = float(number_c12_be8_2p_2n_piminus_piplus) / float(number_events) * 100
    frac_c12_be8_3p_n_2piminus_piplus = float(number_c12_be8_3p_n_2piminus_piplus) / float(number_events) * 100
    frac_c12_be8_4p_2piminus = float(number_c12_be8_4p_2piminus) / float(number_events) * 100

    # C9:
    frac_c12_c9_p_2n_piminus = float(number_c12_c9_p_2n_piminus) / float(number_events) * 100
    frac_c12_c9_3n = float(number_c12_c9_3n) / float(number_events) * 100
    frac_c12_c9_2p_n_2piminus = float(number_c12_c9_2p_n_2piminus) / float(number_events) * 100
    frac_c12_c9_3n_2piminus_2piplus = float(number_c12_c9_3n_2piminus_2piplus) / float(number_events) * 100

    # Be7:
    frac_c12_be7_2p_3n = float(number_c12_be7_2p_3n) / float(number_events) * 100
    frac_c12_be7_p_4n_piplus = float(number_c12_be7_p_4n_piplus) / float(number_events) * 100
    frac_c12_be7_2p_3n_2piminus_2piplus = float(number_c12_be7_2p_3n_2piminus_2piplus) / float(number_events) * 100
    frac_c12_be7_3p_2n_piminus = float(number_c12_be7_3p_2n_piminus) / float(number_events) * 100
    frac_c12_be7_4p_n_2piminus = float(number_c12_be7_4p_n_2piminus) / float(number_events) * 100
    frac_c12_be7_3p_2n_2piminus_piplus = float(number_c12_be7_3p_2n_2piminus_piplus) / float(number_events) * 100

    # Li6:
    frac_c12_li6_3p_3n = float(number_c12_li6_3p_3n) / float(number_events) * 100
    frac_c12_li6_2p_4n_piplus = float(number_c12_li6_2p_4n_piplus) / float(number_events) * 100
    frac_c12_li6_5p_n_2piminus = float(number_c12_li6_5p_n_2piminus) / float(number_events) * 100
    frac_c12_li6_2p_4n_piminus_2piplus = float(number_c12_li6_2p_4n_piminus_2piplus) / float(number_events) * 100
    frac_c12_li6_4p_2n_piminus = float(number_c12_li6_4p_2n_piminus) / float(number_events) * 100
    frac_c12_li6_3p_3n_piminus_piplus = float(number_c12_li6_3p_3n_piminus_piplus) / float(number_events) * 100

    # Li8:
    frac_c12_li8_3p_n = float(number_c12_li8_3p_n) / float(number_events) * 100
    frac_c12_li8_4p_piminus = float(number_c12_li8_4p_piminus) / float(number_events) * 100
    frac_c12_li8_4p_2piminus_piplus = float(number_c12_li8_4p_2piminus_piplus) / float(number_events) * 100
    frac_c12_li8_2p_2n_piplus = float(number_c12_li8_2p_2n_piplus) / float(number_events) * 100
    frac_c12_li8_3p_n_piminus_piplus = float(number_c12_li8_3p_n_piminus_piplus) / float(number_events) * 100

    # Li7:
    frac_c12_li7_2p_3n_piplus = float(number_c12_li7_2p_3n_piplus) / float(number_events) * 100
    frac_c12_li7_4p_n_piminus = float(number_c12_li7_4p_n_piminus) / float(number_events) * 100
    frac_c12_li7_3p_2n = float(number_c12_li7_3p_2n) / float(number_events) * 100
    frac_c12_li7_3p_2n_piminus_piplus = float(number_c12_li7_3p_2n_piminus_piplus) / float(number_events) * 100
    frac_c12_li7_4p_n_2piminus_piplus = float(number_c12_li7_4p_n_2piminus_piplus) / float(number_events) * 100
    frac_c12_li7_2p_3n_piminus_2piplus = float(number_c12_li7_2p_3n_piminus_2piplus) / float(number_events) * 100

    # B8:
    frac_c12_b8_p_3n = float(number_c12_b8_p_3n) / float(number_events) * 100
    frac_c12_b8_p_3n_piminus_piplus = float(number_c12_b8_p_3n_piminus_piplus) / float(number_events) * 100
    frac_c12_b8_2p_2n_2piminus_piplus = float(number_c12_b8_2p_2n_2piminus_piplus) / float(number_events) * 100
    frac_c12_b8_2p_2n_piminus = float(number_c12_b8_2p_2n_piminus) / float(number_events) * 100
    frac_c12_b8_4n_piplus = float(number_c12_b8_4n_piplus) / float(number_events) * 100

    # Li9:
    frac_c12_li9_2p_n_piplus = float(number_c12_li9_2p_n_piplus) / float(number_events) * 100
    frac_c12_li9_3p = float(number_c12_li9_3p) / float(number_events) * 100
    frac_c12_li9_3p_piminus_piplus = float(number_c12_li9_3p_piminus_piplus) / float(number_events) * 100
    frac_c12_li9_2p_n_piminus_2piplus = float(number_c12_li9_2p_n_piminus_2piplus) / float(number_events) * 100
    frac_c12_li9_p_2n_piminus_3piplus = float(number_c12_li9_p_2n_piminus_3piplus) / float(number_events) * 100

    # C8:
    frac_c12_c8_4n = float(number_c12_c8_4n) / float(number_events) * 100

    # He8:
    frac_c12_he8_4p = float(number_c12_he8_4p) / float(number_events) * 100

    # B7:
    frac_c12_b7_p_4n = float(number_c12_b7_p_4n) / float(number_events) * 100

    # He7:
    frac_c12_he7_4p_n = float(number_c12_he7_4p_n) / float(number_events) * 100

    # H7:
    frac_c12_h7_5p = float(number_c12_h7_5p) / float(number_events) * 100

    # Be6:
    frac_c12_be6_2p_4n = float(number_c12_be6_2p_4n) / float(number_events) * 100

    # Li5:
    frac_c12_li5_3p_4n = float(number_c12_li5_3p_4n) / float(number_events) * 100

    # Li4:
    frac_c12_li4_3p_5n = float(number_c12_li4_3p_5n) / float(number_events) * 100

    # He6:
    frac_c12_he6_4p_2n = float(number_c12_he6_4p_2n) / float(number_events) * 100

    # He5:
    frac_c12_he5_4p_3n = float(number_c12_he5_4p_3n) / float(number_events) * 100

    # He4:
    frac_c12_he4_4p_4n = float(number_c12_he4_4p_4n) / float(number_events) * 100

    # He3:
    frac_c12_he3_4p_5n = float(number_c12_he3_4p_5n) / float(number_events) * 100

    # H6:
    frac_c12_h6_5p_n = float(number_c12_h6_5p_n) / float(number_events) * 100

    # H5:
    frac_c12_h5_5p_2n = float(number_c12_h5_5p_2n) / float(number_events) * 100

    # H4:
    frac_c12_h4_5p_3n = float(number_c12_h4_5p_3n) / float(number_events) * 100

    # H3:
    frac_c12_h3_5p_4n = float(number_c12_h3_5p_4n) / float(number_events) * 100

    # H2:
    frac_c12_h2_5p_5n = float(number_c12_h2_5p_5n) / float(number_events) * 100

    # C12:
    frac_c12_c12 = float(number_c12_c12) / float(number_events) * 100

    # no isotope (only protons, neutrons, pions):
    frac_c12_noiso = float(number_c12_noiso) / float(number_events) * 100

    # missing interaction channels:
    frac_c12_missing = float(number_c12_missing) / float(number_events) * 100

    # Other targets than C12:
    frac_no_c12 = float(number_no_c12) / float(number_events) * 100
    frac_es_p = float(number_es_p) / float(number_events) * 100
    frac_es_e = float(number_es_e) / float(number_events) * 100
    frac_es_o16 = float(number_es_o16) / float(number_events) * 100
    frac_es_n14 = float(number_es_n14) / float(number_events) * 100
    frac_es_s32 = float(number_es_s32) / float(number_events) * 100


    return (number_events,
            frac_c12_b11_p, frac_c12_b11_n_piplus, frac_c12_b11_n_piminus_2piplus, frac_c12_b11_p_piminus_piplus,
            frac_c12_b11_p_2piminus_2piplus, frac_c12_b11_piplus,
            frac_c12_c11_n, frac_c12_c11_p_piminus, frac_c12_c11_n_piminus_piplus, frac_c12_c11_p_2piminus_piplus,
            frac_c12_c11_p_3piminus_2piplus, frac_c12_c11_n_2piminus_2piplus,
            frac_c12_b10_p_n, frac_c12_b10_2p_piminus, frac_c12_b10_p_n_piminus_piplus, frac_c12_b10_2n_piplus,
            frac_c12_b10_2n_piminus_2piplus, frac_c12_b10_2p_2piminus_piplus, frac_c12_b10_2p_3piminus_2piplus,
            frac_c12_b10_p_n_2piminus_2piplus,
            frac_c12_c10_2n, frac_c12_c10_p_n_piminus, frac_c12_c10_p_n_2piminus_piplus, frac_c12_c10_2n_piminus_piplus,
            frac_c12_c10_2p_2piminus,
            frac_c12_be10_2p, frac_c12_be10_p_n_piplus, frac_c12_be10_p_n_piminus_2piplus,
            frac_c12_be10_2p_piminus_piplus, frac_c12_be10_2n_2piplus, frac_c12_be10_p_n_2piminus_3piplus,
            frac_c12_be10_2p_2piminus_2piplus, frac_c12_be10_2p_3piminus_3piplus,
            frac_c12_b9_p_2n, frac_c12_b9_p_2n_piminus_piplus, frac_c12_b9_2p_n_3piminus_2piplus,
            frac_c12_b9_2p_n_piminus, frac_c12_b9_3n_piplus, frac_c12_b9_p_2n_2piminus_2piplus,
            frac_c12_b9_2p_n_2piminus_piplus,
            frac_c12_be9_2p_n, frac_c12_be9_p_2n_piplus, frac_c12_be9_3p_piminus, frac_c12_be9_p_2n_piminus_2piplus,
            frac_c12_be9_2p_n_piminus_piplus, frac_c12_be9_2p_n_3piminus_3piplus, frac_c12_be9_2p_n_2piminus_2piplus,
            frac_c12_be9_3n_2piplus, frac_c12_be9_3p_2piminus_piplus,
            frac_c12_be8_2p_2n, frac_c12_be8_3p_n_piminus, frac_c12_be8_p_3n_piplus,
            frac_c12_be8_2p_2n_2piminus_2piplus, frac_c12_be8_4n_2piplus, frac_c12_be8_2p_2n_piminus_piplus,
            frac_c12_be8_3p_n_2piminus_piplus, frac_c12_be8_4p_2piminus,
            frac_c12_c9_p_2n_piminus, frac_c12_c9_3n, frac_c12_c9_2p_n_2piminus, frac_c12_c9_3n_2piminus_2piplus,
            frac_c12_be7_2p_3n, frac_c12_be7_p_4n_piplus, frac_c12_be7_2p_3n_2piminus_2piplus,
            frac_c12_be7_3p_2n_piminus, frac_c12_be7_4p_n_2piminus, frac_c12_be7_3p_2n_2piminus_piplus,
            frac_c12_li6_3p_3n, frac_c12_li6_2p_4n_piplus, frac_c12_li6_5p_n_2piminus,
            frac_c12_li6_2p_4n_piminus_2piplus, frac_c12_li6_4p_2n_piminus, frac_c12_li6_3p_3n_piminus_piplus,
            frac_c12_li8_3p_n, frac_c12_li8_4p_piminus, frac_c12_li8_4p_2piminus_piplus, frac_c12_li8_2p_2n_piplus,
            frac_c12_li8_3p_n_piminus_piplus,
            frac_c12_li7_2p_3n_piplus, frac_c12_li7_4p_n_piminus, frac_c12_li7_3p_2n, frac_c12_li7_3p_2n_piminus_piplus,
            frac_c12_li7_4p_n_2piminus_piplus, frac_c12_li7_2p_3n_piminus_2piplus,
            frac_c12_b8_p_3n, frac_c12_b8_p_3n_piminus_piplus, frac_c12_b8_2p_2n_2piminus_piplus,
            frac_c12_b8_2p_2n_piminus, frac_c12_b8_4n_piplus,
            frac_c12_li9_2p_n_piplus, frac_c12_li9_3p, frac_c12_li9_3p_piminus_piplus,
            frac_c12_li9_2p_n_piminus_2piplus, frac_c12_li9_p_2n_piminus_3piplus,
            frac_c12_c8_4n, frac_c12_he8_4p, frac_c12_b7_p_4n, frac_c12_he7_4p_n, frac_c12_h7_5p, frac_c12_be6_2p_4n,
            frac_c12_li5_3p_4n, frac_c12_li4_3p_5n,
            frac_c12_he6_4p_2n, frac_c12_he5_4p_3n, frac_c12_he4_4p_4n, frac_c12_he3_4p_5n, frac_c12_h6_5p_n,
            frac_c12_h5_5p_2n, frac_c12_h4_5p_3n, frac_c12_h3_5p_4n,
            frac_c12_h2_5p_5n, frac_c12_c12,
            frac_c12_noiso,
            frac_no_c12, frac_es_p, frac_es_e, frac_es_o16, frac_es_n14, frac_es_s32, frac_c12_missing)


def get_mass_from_pdg(pdg):
    """
    function to get the mass in GeV of a particle from its PDG ID (Monte Carlo Particle Number Scheme)
    (masses taken from NCGenerator.cc, )

    :param pdg: PDG ID of the particle (integer)
    :return: mass[pdg]: mass of the particle in GeV (float)
    """
    mass = dict()
    # mass of gamma:
    mass[22] = 0
    # mass of electron:
    mass[11] = 0.000511
    # mass of positron:
    mass[-11] = 0.000511
    # mass of electron-neutrino:
    mass[12] = 0
    # mass of electron-antineutrino:
    mass[-12] = 0
    # mass of muon (https://de.wikipedia.org/wiki/Myon):
    mass[13] = 0.105658
    # mass of anti-muon (https://de.wikipedia.org/wiki/Myon):
    mass[-13] = 0.105658
    # mass of muon-neutrino:
    mass[14] = 0
    # mass of muon-antineutrino:
    mass[-14] = 0
    # mass of pion_0:
    mass[111] = 0.13957
    # mass of pion_plus:
    mass[211] = 0.13957
    # mass of pion_minus:
    mass[-211] = 0.13957
    # mass of Kaon plus (http://pdg.lbl.gov/2018/reviews/rpp2018-rev-charged-kaon-mass.pdf):
    mass[321] = 0.493677
    # mass of Kaon minus (http://pdg.lbl.gov/2018/reviews/rpp2018-rev-charged-kaon-mass.pdf):
    mass[-321] = 0.493677
    # mass of Kaon 0 (http://pdg.lbl.gov/2015/tables/rpp2015-tab-mesons-strange.pdf):
    mass[311] = 0.497611
    # mass of anti Kaon 0 (http://pdg.lbl.gov/2015/tables/rpp2015-tab-mesons-strange.pdf):
    mass[-311] = 0.497611
    # mass of Lambda (http://pdg.lbl.gov/2017/tables/rpp2017-tab-baryons-Lambda.pdf):
    mass[3122] = 1.115683
    # mass of anti-Lambda (http://pdg.lbl.gov/2017/tables/rpp2017-tab-baryons-Lambda.pdf):
    mass[-3122] = 1.115683
    # mass of sigma plus (http://pdg.lbl.gov/2017/tables/rpp2017-tab-baryons-Sigma.pdf):
    mass[3222] = 1.18937
    # mass of sigma minus (http://pdg.lbl.gov/2017/tables/rpp2017-tab-baryons-Sigma.pdf):
    mass[3112] = 1.197449
    # mass of sigma_0 (http://pdg.lbl.gov/2017/tables/rpp2017-tab-baryons-Sigma.pdf):
    mass[3212] = 1.192642
    # mass of anti_sigma_0 (http://pdg.lbl.gov/2017/tables/rpp2017-tab-baryons-Sigma.pdf):
    mass[-3212] = 1.192642
    # mass of neutron:
    mass[2112] = 0.93957
    # mass of anti-neutron:
    mass[-2112] = 0.93957
    # proton:
    mass[2212] = 0.93827
    # mass of anti-proton:
    mass[-2212] = 0.93827
    # deuterium H2 (stable):
    mass[1000010020] = 1.8756
    # tritium H3:
    mass[1000010030] = 2.8089
    # H6 (https://en.wikipedia.org/wiki/Isotopes_of_hydrogen) (6.044u * 0.93149 GeV/u):
    mass[1000010060] = 6.044 * 0.93149
    # H7 (https://en.wikipedia.org/wiki/Isotopes_of_hydrogen):
    mass[1000010070] = 7.052 * 0.93149
    # He3 (stable):
    mass[1000020030] = 2.8084
    # He4 or alpha (stable):
    mass[1000020040] = 3.7274
    # He6 (https://en.wikipedia.org/wiki/Isotopes_of_helium) (6.018u * 0.93149 GeV/u):
    mass[1000020060] = 6.019 * 0.93149
    # He7 (https://en.wikipedia.org/wiki/Isotopes_of_helium):
    mass[1000020070] = 7.028 * 0.93149
    # He8 (https://en.wikipedia.org/wiki/Isotopes_of_helium):
    mass[1000020080] = 8.034 * 0.93149
    # Li6 (stable):
    mass[1000030060] = 5.6015
    # Li7 (stable):
    mass[1000030070] = 6.5335
    # Li8:
    mass[1000030080] = 7.4708
    # Li9:
    mass[1000030090] = 8.4061
    # Be6 (https://en.wikipedia.org/wiki/Isotopes_of_beryllium) (6.020u * 0.93149 GeV/u):
    mass[1000040060] = 6.020 * 0.93149
    # Be7:
    mass[1000040070] = 6.5344
    # Be8:
    mass[1000040080] = 7.4548
    # Be9 (stable):
    mass[1000040090] = 8.3925
    # Be10:
    mass[1000040100] = 9.3249
    # B7 (https://en.wikipedia.org/wiki/Isotopes_of_boron) (7.030u * 0.93149 GeV/u):
    mass[1000050070] = 7.030 * 0.93149
    # B8:
    mass[1000050080] = 7.4728
    # B9:
    mass[1000050090] = 8.3935
    # B10 (stable):
    mass[1000050100] = 9.3244
    # B11 (stable):
    mass[1000050110] = 10.2522
    # C8 (https://en.wikipedia.org/wiki/Isotopes_of_carbon) (8.038u * 0.93149 GeV/u):
    mass[1000060080] = 8.038 * 0.93149
    # C9:
    mass[1000060090] = 8.4100
    # C10:
    mass[1000060100] = 9.3280
    # C11:
    mass[1000060110] = 10.2542
    # C12 (stable):
    mass[1000060120] = 11.1748
    # O16 (stable) (15.999u * 0.93149 GeV/u):
    mass[1000080160] = 15.999 * 0.93149
    # N14 (stable):
    mass[1000070140] = 14.0067 * 0.93149
    # S32 (stable):
    mass[1000160320] = 32.065 * 0.93149

    """ set dummy value to very rare particles: """
    # Lambda_c+:
    mass[4122] = 0
    # Sigma_c+:
    mass[4212] = 0
    # Sigma_c++:
    mass[4222] = 0
    # D_0:
    mass[421] = 0
    # anti-D_0:
    mass[-421] = 0
    # D_plus:
    mass[411] = 0
    # D_minus:
    mass[-411] = 0
    # Kaon long:
    mass[130] = 0
    # Kaon short:
    mass[310] = 0
    # D_s plus:
    mass[431] = 0
    # anti D_s minus:
    mass[-431] = 0


    return mass[pdg]


def get_number_of_p_and_n_of_isotope(pdg_id):
    """
    function to get the number of protons and neutrons of a nuclei from its PDG ID

    :param pdg_id: PDG ID of the nuclei (float)
    :return:
    """
    # get the number of protons (integer)
    number_p = int((pdg_id - 1000000000) / 10000)

    # get the sum of protons and neutrons (integer):
    number_p_and_n = int((pdg_id - 1000000000 - number_p*10000) / 10)

    # calculate the number of neutrons (integer):
    number_n = number_p_and_n - number_p

    return number_p, number_n


def get_number_of_particles_of_channelid(channel_id):
    """
    Function to calculate the number of neutrinos, protons, neutron, pion_minus and pion_plus from the channel ID

    :param channel_id: Channel ID of the NC interaction, represents which particles are produced via the NC
    interaction (float)
    :return:
    """
    # preallocate variables:
    number_p = 0
    number_n = 0
    number_pion_minus = 0
    number_pion_plus = 0

    if channel_id > 10000:
        # channel_id > 10000 -> at least one pion is created:
        # number of protons from channel ID (integer):
        number_p = int((channel_id - 10000) / 1000)
        # number of neutrons from channel ID (integer):
        number_n = int((channel_id - 10000 - number_p*1000) / 100)
        # number of pion_minus from channel ID (integer):
        number_pion_minus = int((channel_id - 10000 - number_p*1000 - number_n*100) / 10)
        # number of pion_plus from channel ID (integer):
        number_pion_plus = int(channel_id - 10000 - number_p*1000 - number_n*100 - number_pion_minus*10)

    elif channel_id > 3:
        # channel_id > 3 -> only protons and neutrons are created:
        # number of protons from channel ID (integer):
        number_p = int((channel_id - 100) / 10)
        # number of neutrons from channel ID (integer):
        number_n = int(channel_id - 100 - number_p*10)

    elif channel_id == 2:
        print("channel ID = 2 in get_number_of_particles_of_channelid()")

    elif channel_id == 3:
        print("channel ID = 3 in get_number_of_particles_of_channelid()")

    return number_p, number_n, number_pion_minus, number_pion_plus


def get_number_of_particles_of_deexid(deex_id):
    """
    Function to calculate the number of different particles which are produced by deexcitation of isotopes.

    :param deex_id: deexcitation channel ID from DSNB-NC generator (more information in ~/juno/test_output_DSNB_gen/)

    :return:
    """
    # preallocate variables:
    number_n = 0
    number_p = 0
    number_deuterium = 0
    number_tritium = 0
    number_he3 = 0
    number_alpha = 0

    if deex_id > 0:
        # deex_id > 0 -> nucleus is de-excited:
        # number of neutrons:
        number_n = int((deex_id - 1000000) / 100000)
        # number of protons:
        number_p = int((deex_id - 1000000 - number_n*100000) / 10000)
        # number of deuterium:
        number_deuterium = int((deex_id - 1000000 - number_n*100000 - number_p*10000) / 1000)
        # number of tritium:
        number_tritium = int((deex_id - 1000000 - number_n*100000 - number_p*10000 - number_deuterium*1000) / 100)
        # number of He3:
        number_he3 = int((deex_id - 1000000 - number_n*100000 - number_p*10000 - number_deuterium*1000 -
                          number_tritium*100) / 10)
        # number of alpha/He4:
        number_alpha = int(deex_id - 1000000 - number_n*100000 - number_p*10000 - number_deuterium*1000 -
                           number_tritium*100 - number_he3*10)

    elif deex_id == 0:
        # set all numbers to 0:
        number_n = 0
        number_p = 0
        number_deuterium = 0
        number_tritium = 0
        number_he3 = 0
        number_alpha = 0

    else:
        print("ERROR in get_number_of_particles_of deexid: deex_id is negative: deex_id = {0:d}".format(deex_id))

    return number_n, number_p, number_deuterium, number_tritium, number_he3, number_alpha


def read_nc_data(rootfile):
    """
    function reads a ROOT-file and saves the values from the root-tree to numpy arrays.

    :param rootfile: path to the input root file (string)

    :return:
    """
    # load the ROOT file:
    rfile = ROOT.TFile(rootfile)
    # get the TTree from the TFile:
    rtree = rfile.Get("genEvt")

    # get the number of entries, i.e. events, in the ROOT-file:
    number_entries = rtree.GetEntries()

    "preallocate all arrays: "
    # event ID (starts from 0) (1d array of integer):
    event_id = np.array([])
    # PDG ID of the projectile (i.e. which neutrino is interacting) (1d array integer):
    projectile_pdg = np.array([])
    # energy of the incoming neutrino in GeV (1d array of float):
    projectile_energy = np.array([])
    # PDG ID of the target particle (either C12 or proton) (1d array integer):
    target_pdg = np.array([])
    # Channel ID of the NC interaction, represents which particles are produced via the NC interaction
    # (1d array integer):
    nc_interaction_ch_id = np.array([])
    # Channel ID of the deexcitation, represents which particles are produced via the deexication of the produced
    # excited isotope (1d array integer):
    deexcitation_id = np.array([])
    # PDG ID of the isotope after the NC interaction, BUT before the deexcitation (1d array integer):
    isotope_pdg = np.array([])
    # number of final particles after NC interactions and deexcitation (1d array integer):
    n_particles = np.array([])

    # PDG ID of the final particles (list of np.arrays of integers):
    final_pdg = []
    # momentum in x-direction of the final particles (list of np.arrays of floats):
    final_px = []
    # momentum in y-direction of the final particles (list of np.arrays of floats):
    final_py = []
    # momentum in z-direction of the final particles (list of np.arrays of floats):
    final_pz = []

    # loop over every entry, i.e. every event, in the TTree:
    for event in range(number_entries):

        # get the current event in the TTree:
        rtree.GetEntry(event)

        # get the value of the event ID and append it to the array:
        evt_id = rtree.GetBranch('t_evtID').GetLeaf('t_evtID').GetValue()
        event_id = np.append(event_id, evt_id)

        # get the value of the projectile PDG and append it to the array:
        pjt_pdg = rtree.GetBranch('t_pPdg').GetLeaf('t_pPdg').GetValue()
        projectile_pdg = np.append(projectile_pdg, pjt_pdg)

        # get the value of the neutrino energy and append it to the array:
        pjt_e = rtree.GetBranch('t_pEn').GetLeaf('t_pEn').GetValue()
        projectile_energy = np.append(projectile_energy, pjt_e)

        # get the value of the target PDG and append it to the array:
        trt_pdg = rtree.GetBranch('t_tPdg').GetLeaf('t_tPdg').GetValue()
        target_pdg = np.append(target_pdg, trt_pdg)

        # get the value of the channel ID and append it to the array:
        ch_id = rtree.GetBranch('t_channelID').GetLeaf('t_channelID').GetValue()
        nc_interaction_ch_id = np.append(nc_interaction_ch_id, ch_id)

        # get the value of the deexcitation ID  and append it to the array:
        deex_id = rtree.GetBranch('t_deexID').GetLeaf('t_deexID').GetValue()
        deexcitation_id = np.append(deexcitation_id, deex_id)

        # get the value of the PDG of the produced isotope and append it to the array:
        iso_pdg = rtree.GetBranch('t_isoPdg').GetLeaf('t_isoPdg').GetValue()
        isotope_pdg = np.append(isotope_pdg, iso_pdg)

        # get the value of the number of particles and append it to the array:
        n_par = rtree.GetBranch('t_Npars').GetLeaf('t_Npars').GetValue()
        n_particles = np.append(n_particles, n_par)

        # get final PDGs of all final particles
        # preallocate an array, where all "n_par" values are stored:
        f_pdg_array = np.array([])

        # loop over the number of particles, get the final PDG and append it to the array:
        for index in range(int(n_par)):
            # get the value of the final PDG and append it to the array:
            f_pdg = rtree.GetBranch('t_pdg').GetLeaf('t_pdg').GetValue(index)
            f_pdg_array = np.append(f_pdg_array, f_pdg)

        # append the np.array to the list:
        final_pdg.append(f_pdg_array)

        # get final momentum Px of all final particles
        # preallocate an array, where all "n_par" values are stored:
        f_px_array = np.array([])

        # loop over the number of particles, get the final momentum and append it to the array:
        for index in range(int(n_par)):
            # get the value of the final momentum Px and append it to the array:
            f_px = rtree.GetBranch('t_px').GetLeaf('t_px').GetValue(index)
            f_px_array = np.append(f_px_array, f_px)

        # append the np.array to the list:
        final_px.append(f_px_array)

        # get final momentum Py of all final particles
        # preallocate an array, where all "n_par" values are stored:
        f_py_array = np.array([])

        # loop over the number of particles, get the final momentum and append it to the array:
        for index in range(int(n_par)):
            # get the value of the final momentum Py and append it to the array:
            f_py = rtree.GetBranch('t_py').GetLeaf('t_py').GetValue(index)
            f_py_array = np.append(f_py_array, f_py)

        # append the np.array to the list:
        final_py.append(f_py_array)

        # get final momentum Pz of all final particles
        # preallocate an array, where all "n_par" values are stored:
        f_pz_array = np.array([])

        # loop over the number of particles, get the final momentum and append it to the array:
        for index in range(int(n_par)):
            # get the value of the final momentum Pz and append it to the array:
            f_pz = rtree.GetBranch('t_pz').GetLeaf('t_pz').GetValue(index)
            f_pz_array = np.append(f_pz_array, f_pz)

        # append the np.array to the list:
        final_pz.append(f_py_array)

    return (event_id, projectile_pdg, projectile_energy, target_pdg, nc_interaction_ch_id, deexcitation_id,
            isotope_pdg, n_particles, final_pdg, final_px, final_py, final_pz)


def read_nc_data_ibdlike_signal(rootfile, event_number):
    """
    function reads a ROOT-file and saves the values from the root-tree to numpy arrays.

    This is the same function like read_nc_data(), BUT not all events in 'rootfile' are read, but only some specific
    ones (only atmospheric NC events that cause an IBD-like signal in JUNO detector)

    :param rootfile: path to the input root file (string)
    :param event_number: list of the event numbers in ascending order (array/list)

    :return:
    """
    # load the ROOT file:
    rfile = ROOT.TFile(rootfile)
    # get the TTree from the TFile:
    rtree = rfile.Get("genEvt")

    # get the number of entries, i.e. events, in the ROOT-file:
    number_entries = rtree.GetEntries()

    "preallocate all arrays: "
    # event ID (starts from 0) (1d array of integer):
    event_id = np.array([])
    # PDG ID of the projectile (i.e. which neutrino is interacting) (1d array integer):
    projectile_pdg = np.array([])
    # energy of the incoming neutrino in GeV (1d array of float):
    projectile_energy = np.array([])
    # PDG ID of the target particle (either C12 or proton) (1d array integer):
    target_pdg = np.array([])
    # Channel ID of the NC interaction, represents which particles are produced via the NC interaction
    # (1d array integer):
    nc_interaction_ch_id = np.array([])
    # Channel ID of the deexcitation, represents which particles are produced via the deexication of the produced
    # excited isotope (1d array integer):
    deexcitation_id = np.array([])
    # PDG ID of the isotope after the NC interaction, BUT before the deexcitation (1d array integer):
    isotope_pdg = np.array([])
    # number of final particles after NC interactions and deexcitation (1d array integer):
    n_particles = np.array([])

    # PDG ID of the final particles (list of np.arrays of integers):
    final_pdg = []
    # momentum in x-direction of the final particles (list of np.arrays of floats):
    final_px = []
    # momentum in y-direction of the final particles (list of np.arrays of floats):
    final_py = []
    # momentum in z-direction of the final particles (list of np.arrays of floats):
    final_pz = []

    # loop over the entries that are specified by event_number (only IBDlike events),
    # index is counter starting from 0 to len(event_number),
    # event is the value of event_number[index]:
    for index, event in enumerate(event_number):

        # get the current event in the TTree:
        rtree.GetEntry(event)

        # get the value of the event ID and append it to the array:
        evt_id = rtree.GetBranch('t_evtID').GetLeaf('t_evtID').GetValue()
        event_id = np.append(event_id, evt_id)

        # get the value of the projectile PDG and append it to the array:
        pjt_pdg = rtree.GetBranch('t_pPdg').GetLeaf('t_pPdg').GetValue()
        projectile_pdg = np.append(projectile_pdg, pjt_pdg)

        # get the value of the neutrino energy and append it to the array:
        pjt_e = rtree.GetBranch('t_pEn').GetLeaf('t_pEn').GetValue()
        projectile_energy = np.append(projectile_energy, pjt_e)

        # get the value of the target PDG and append it to the array:
        trt_pdg = rtree.GetBranch('t_tPdg').GetLeaf('t_tPdg').GetValue()
        target_pdg = np.append(target_pdg, trt_pdg)

        # get the value of the channel ID and append it to the array:
        ch_id = rtree.GetBranch('t_channelID').GetLeaf('t_channelID').GetValue()
        nc_interaction_ch_id = np.append(nc_interaction_ch_id, ch_id)

        # get the value of the deexcitation ID  and append it to the array:
        deex_id = rtree.GetBranch('t_deexID').GetLeaf('t_deexID').GetValue()
        deexcitation_id = np.append(deexcitation_id, deex_id)

        # get the value of the PDG of the produced isotope and append it to the array:
        iso_pdg = rtree.GetBranch('t_isoPdg').GetLeaf('t_isoPdg').GetValue()
        isotope_pdg = np.append(isotope_pdg, iso_pdg)

        # get the value of the number of particles and append it to the array:
        n_par = rtree.GetBranch('t_Npars').GetLeaf('t_Npars').GetValue()
        n_particles = np.append(n_particles, n_par)

        # get final PDGs of all final particles
        # preallocate an array, where all "n_par" values are stored:
        f_pdg_array = np.array([])

        # loop over the number of particles, get the final PDG and append it to the array:
        for index in range(int(n_par)):
            # get the value of the final PDG and append it to the array:
            f_pdg = rtree.GetBranch('t_pdg').GetLeaf('t_pdg').GetValue(index)
            f_pdg_array = np.append(f_pdg_array, f_pdg)

        # append the np.array to the list:
        final_pdg.append(f_pdg_array)

        # get final momentum Px of all final particles
        # preallocate an array, where all "n_par" values are stored:
        f_px_array = np.array([])

        # loop over the number of particles, get the final momentum and append it to the array:
        for index in range(int(n_par)):
            # get the value of the final momentum Px and append it to the array:
            f_px = rtree.GetBranch('t_px').GetLeaf('t_px').GetValue(index)
            f_px_array = np.append(f_px_array, f_px)

        # append the np.array to the list:
        final_px.append(f_px_array)

        # get final momentum Py of all final particles
        # preallocate an array, where all "n_par" values are stored:
        f_py_array = np.array([])

        # loop over the number of particles, get the final momentum and append it to the array:
        for index in range(int(n_par)):
            # get the value of the final momentum Py and append it to the array:
            f_py = rtree.GetBranch('t_py').GetLeaf('t_py').GetValue(index)
            f_py_array = np.append(f_py_array, f_py)

        # append the np.array to the list:
        final_py.append(f_py_array)

        # get final momentum Pz of all final particles
        # preallocate an array, where all "n_par" values are stored:
        f_pz_array = np.array([])

        # loop over the number of particles, get the final momentum and append it to the array:
        for index in range(int(n_par)):
            # get the value of the final momentum Pz and append it to the array:
            f_pz = rtree.GetBranch('t_pz').GetLeaf('t_pz').GetValue(index)
            f_pz_array = np.append(f_pz_array, f_pz)

        # append the np.array to the list:
        final_pz.append(f_py_array)

    return (event_id, projectile_pdg, projectile_energy, target_pdg, nc_interaction_ch_id, deexcitation_id,
            isotope_pdg, n_particles, final_pdg, final_px, final_py, final_pz)


def get_neutrino_energy(projectile_pdg, projectile_energy, bin_width, event_rate, time):
    """
    function to get the number of events (number of NC interactions of atmospheric neutrinos) as function of the
    neutrino energy of the different types of neutrinos (electron-neutrino, electron-antineutrino, muon-neutrino,
    muon-antineutrino, tau-neutrino, tau-antineutrino).
    Also the total number of events for each neutrino type and the fraction of each neutrino type to the total number
    of events is calculated.

    :param projectile_pdg: PDG ID of the projectile, i.e. of the incoming neutrino (array float)
    :param projectile_energy: energy of the projectile, i.e. of the incoming neutrinos, in GeV (array of float)
    :param bin_width: bin width of the array, which represents the energy of incoming neutrinos in GeV (float)
    :param event_rate: NC interaction event rate in units of 1/(sec * 20 kton) (float)
    :param time: total exposure time in seconds (float)

    :return:
    """

    """ Preallocate all arrays: """
    # energy of incoming electron neutrinos in GeV (np.array of float):
    energy_nu_e_incoming = np.array([])
    # energy of incoming electron antineutrinos in GeV (np.array of float):
    energy_nu_e_bar_incoming = np.array([])
    # energy of incoming muon neutrinos in GeV (np.array of float):
    energy_nu_mu_incoming = np.array([])
    # energy of incoming muon antineutrinos in GeV (np.array of float):
    energy_nu_mu_bar_incoming = np.array([])
    # energy of incoming tau neutrinos in GeV (np.array of float):
    energy_nu_tau_incoming = np.array([])
    # energy of incoming tau antineutrinos in GeV (np.array of float):
    energy_nu_tau_bar_incoming = np.array([])

    # check, if input arrays have same length:
    if len(projectile_pdg) != len(projectile_energy):
        print("WARNING (in function get_neutrino_energy()): input arrays don't have same length!!")

    # get number of entries in the arrays:
    number_entries = len(projectile_pdg)

    # loop over all entries of the arrays:
    for index in range(number_entries):

        # check the PDG ID:
        if projectile_pdg[index] == 12:
            # for electron-neutrino PDG = 12:
            energy = projectile_energy[index]
            energy_nu_e_incoming = np.append(energy_nu_e_incoming, energy)

        elif projectile_pdg[index] == -12:
            # for electron-antineutrino PDG = -12:
            energy = projectile_energy[index]
            energy_nu_e_bar_incoming = np.append(energy_nu_e_bar_incoming, energy)

        elif projectile_pdg[index] == 14:
            # for muon-neutrino PDG = 14:
            energy = projectile_energy[index]
            energy_nu_mu_incoming = np.append(energy_nu_mu_incoming, energy)

        elif projectile_pdg[index] == -14:
            # for muon-antineutrino PDG == -14:
            energy = projectile_energy[index]
            energy_nu_mu_bar_incoming = np.append(energy_nu_mu_bar_incoming, energy)

        elif projectile_pdg[index] == 16:
            # for tau-neutrino PDG = 16:
            energy = projectile_energy[index]
            energy_nu_tau_incoming = np.append(energy_nu_tau_incoming, energy)

        elif projectile_pdg[index] == -16:
            # for tau-antineutrino PDG = -16:
            energy = projectile_energy[index]
            energy_nu_tau_bar_incoming = np.append(energy_nu_tau_bar_incoming, energy)

        else:
            print("WARNING (in function get_neutrino_energy(): NOT only neutrinos as projectile!")

    """ get number of events from NC interactions of different neutrino types: """
    # calculate the number of events from NC interaction of electron-neutrinos (integer):
    n_nu_e = len(energy_nu_e_incoming)
    # calculate the number of events from NC interaction of electron-antineutrinos (integer):
    n_nu_e_bar = len(energy_nu_e_bar_incoming)
    # calculate the number of events from NC interaction of muon-neutrinos (integer):
    n_nu_mu = len(energy_nu_mu_incoming)
    # calculate the number of events from NC interaction of muon-antineutrinos (integer):
    n_nu_mu_bar = len(energy_nu_mu_bar_incoming)
    # calculate the number of events from NC interaction of tau-neutrinos (integer):
    n_nu_tau = len(energy_nu_tau_incoming)
    # calculate the number of events from NC interaction of tau-antineutrinos (integer):
    n_nu_tau_bar = len(energy_nu_tau_bar_incoming)

    """ get fraction of events from NC interactions of different neutrino types (IN PERCENT): """
    # calculate the fraction of events from electron-neutrinos of the all events in % (float):
    fraction_nu_e = float(n_nu_e)/float(number_entries)*100
    # calculate the fraction of events from electron-antineutrinos of the all events in % (float):
    fraction_nu_e_bar = float(n_nu_e_bar)/float(number_entries)*100
    # calculate the fraction of events from muon-neutrinos of the all events in % (float):
    fraction_nu_mu = float(n_nu_mu)/float(number_entries)*100
    # calculate the fraction of events from muon-antineutrinos of the all events in % (float):
    fraction_nu_mu_bar = float(n_nu_mu_bar)/float(number_entries)*100
    # calculate the fraction of events from tau-neutrinos of the all events in % (float):
    fraction_nu_tau = float(n_nu_tau)/float(number_entries)*100
    # calculate the fraction of events from tau-antineutrinos of the all events in % (float):
    fraction_nu_tau_bar = float(n_nu_tau_bar)/float(number_entries)*100

    """ create histograms with the energy arrays from above to get the number of events per bin: """
    energy_range = np.arange(0, np.max(projectile_energy)+2, bin_width)
    events_nu_e_in, bins1 = np.histogram(energy_nu_e_incoming, energy_range)
    events_nu_e_bar_in, bins1 = np.histogram(energy_nu_e_bar_incoming, energy_range)
    events_nu_mu_in, bins1 = np.histogram(energy_nu_mu_incoming, energy_range)
    events_nu_mu_bar_in, bins1 = np.histogram(energy_nu_mu_bar_incoming, energy_range)
    events_nu_tau_in, bins1 = np.histogram(energy_nu_tau_incoming, energy_range)
    events_nu_tau_bar_in, bins1 = np.histogram(energy_nu_tau_bar_incoming, energy_range)

    """ calculate the number of neutrino NC interactions as function of energy per "time" per 20 ktons: """
    # TODO-me: event rate and exposure time are NOT included correctly!!!
    # INFO-me: Calculation is only correct, when event_rate=1 and time=1!!!
    event_nu_e_incoming = events_nu_e_in * event_rate * time
    event_nu_e_bar_incoming = events_nu_e_bar_in * event_rate * time
    event_nu_mu_incoming = events_nu_mu_in * event_rate * time
    event_nu_mu_bar_incoming = events_nu_mu_bar_in * event_rate * time
    event_nu_tau_incoming = events_nu_tau_in * event_rate * time
    event_nu_tau_bar_incoming = events_nu_tau_bar_in * event_rate * time

    """ calculate the total number of neutrino NC interactions per "time" and 20 ktons: """
    number_nu_e_incoming = np.sum(event_nu_e_incoming)
    number_nu_e_bar_incoming = np.sum(event_nu_e_bar_incoming)
    number_nu_mu_incoming = np.sum(event_nu_mu_incoming)
    number_nu_mu_bar_incoming = np.sum(event_nu_mu_bar_incoming)
    number_nu_tau_incoming = np.sum(event_nu_tau_incoming)
    number_nu_tau_bar_incoming = np.sum(event_nu_tau_bar_incoming)

    return (energy_range, event_nu_e_incoming, event_nu_e_bar_incoming, event_nu_mu_incoming, event_nu_mu_bar_incoming,
            event_nu_tau_incoming, event_nu_tau_bar_incoming, number_nu_e_incoming, number_nu_e_bar_incoming,
            number_nu_mu_incoming, number_nu_mu_bar_incoming, number_nu_tau_incoming, number_nu_tau_bar_incoming,
            fraction_nu_e, fraction_nu_e_bar, fraction_nu_mu, fraction_nu_mu_bar, fraction_nu_tau,
            fraction_nu_tau_bar)


def get_target_ratio(projectile_energy, target_pdg, bin_width):
    """
    function to get the fraction of NC interaction events on target particles (C12, N14, O16, S32) and of elastic
    scattering events on free protons and electrons, and to get the number of events of the different target particles
    as function of the energy of the incoming neutrinos.

    :param projectile_energy: energy of the projectile, i.e. of the incoming neutrinos, in GeV (array of float)
    :param target_pdg: PDG ID of the target particle (array of float)
    :param bin_width: bin width of the array, which represents the energy of incoming neutrinos in GeV (float)

    :return:
    """

    """ preallocate the arrays: """
    # energy of the incoming neutrino interacting via NC with C12:
    energy_nu_c12 = np.array([])
    # energy of the incoming neutrino interacting via ES with free protons:
    energy_nu_proton = np.array([])
    # energy of the incoming neutrino interacting via NC with N14:
    energy_nu_n14 = np.array([])
    # energy of the incoming neutrino interacting via NC with O16:
    energy_nu_o16 = np.array([])
    # energy of the incoming neutrino interacting via ES with an electron:
    energy_nu_electron = np.array([])
    # energy of the incoming neutrino interacting via NC with S32:
    energy_nu_s32 = np.array([])

    # check, if input arrays have same length:
    if len(projectile_energy) != len(target_pdg):
        print("WARNING (in function get_target_ratio()): input arrays don't have same length!!")

    # get number of entries in the arrays:
    number_entries = len(projectile_energy)

    # loop over all entries in the arrays:
    for index in range(number_entries):

        # check the PDG ID:
        if target_pdg[index] == 1000060120:
            # PDG ID of C12 = 1000060120:
            energy = projectile_energy[index]
            energy_nu_c12 = np.append(energy_nu_c12, energy)

        elif target_pdg[index] == 2212:
            # PDG ID of proton = 2212:
            energy = projectile_energy[index]
            energy_nu_proton = np.append(energy_nu_proton, energy)

        elif target_pdg[index] == 1000070140:
            # PDG ID of N14 = 1000070140:
            energy = projectile_energy[index]
            energy_nu_n14 = np.append(energy_nu_n14, energy)

        elif target_pdg[index] == 1000080160:
            # PDG ID of O16 = 1000080160:
            energy = projectile_energy[index]
            energy_nu_o16 = np.append(energy_nu_o16, energy)

        elif target_pdg[index] == 11:
            # PDG ID of electron = 11:
            energy = projectile_energy[index]
            energy_nu_electron = np.append(energy_nu_electron, energy)

        elif target_pdg[index] == 1000160320:
            # PDG ID of S32 = 1000160320:
            energy = projectile_energy[index]
            energy_nu_s32 = np.append(energy_nu_s32, energy)

        else:
            print("WARNING (in function get_target_ratio(): new target particle!")
            print(target_pdg[index])


    """ get number of events for the two different targets: """
    # calculate the number of NC interaction events on C12 (integer):
    n_c12 = len(energy_nu_c12)
    # calculate the number of elastic scattering events on free protons (integer):
    n_proton = len(energy_nu_proton)
    # calculate the number of NC interaction events on N14 (integer):
    n_n14 = len(energy_nu_n14)
    # calculate the number of NC interaction events on O16 (integer):
    n_o16 = len(energy_nu_o16)
    # calculate the number of elastic scattering events on electrons (integer):
    n_electron = len(energy_nu_electron)
    # calculate the number of NC interaction events on S32 (integer):
    n_s32 = len(energy_nu_s32)

    """ get fraction of events for the two different targets (IN PERCENT): """
    # calculate the fraction of NC interaction events on C12 in % (float):
    fraction_c12 = float(n_c12)/float(number_entries)*100
    # calculate the fraction of ES events on free protons in % (float):
    fraction_proton = float(n_proton)/float(number_entries)*100
    # calculate the fraction of NC interaction events on N14 in % (float):
    fraction_n14 = float(n_n14)/float(number_entries)*100
    # calculate the fraction of NC interaction events on O16 in % (float):
    fraction_o16 = float(n_o16)/float(number_entries)*100
    # calculate the fraction of ES events on electrons in % (float):
    fraction_electron = float(n_electron)/float(number_entries)*100
    # calculate the fraction of NC interaction events on S32 in % (float):
    fraction_s32 = float(n_s32)/float(number_entries)*100

    """ create histograms with the energy arrays from above to get the number of events per bin: """
    energy_range = np.arange(0, np.max(projectile_energy)+2, bin_width)
    events_c12, bins1 = np.histogram(energy_nu_c12, energy_range)
    events_proton, bins1 = np.histogram(energy_nu_proton, energy_range)
    events_n14, bins1 = np.histogram(energy_nu_n14, energy_range)
    events_o16, bins1 = np.histogram(energy_nu_o16, energy_range)
    events_electron, bins1 = np.histogram(energy_nu_electron, energy_range)
    events_s32, bins1 = np.histogram(energy_nu_s32, energy_range)

    # TODO-me: the event rate and exposure time is NOT included yet!!!

    return (energy_range, events_c12, events_proton, events_n14, events_o16, events_electron, events_s32,
            n_c12, n_proton, n_n14, n_o16, n_electron, n_s32,
            fraction_c12, fraction_proton, fraction_n14, fraction_o16, fraction_electron, fraction_s32)


def get_combined_channel(n_particles, particle_pdg, target_pdg):
    """
    function to get the combined interaction channels (combined = NC interaction channel + deexcitation channel) from
    the interaction of neutrinos with the target particles

    :param n_particles: array with number of final particles per event
    :param particle_pdg: list of array with the PDG ID of the final particles per event
    :param target_pdg: array with the PDG ID of the target
    :return:
    """
    # get the number of entries of the array (integer):
    number_entries = len(target_pdg)

    """ preallocate array and variables: """
    # INFO-me: gamma's are not considered in specification of the channels -> all channels can be with additional gammas
    # B11:
    number_b11 = 0
    number_b11_n_piplus = 0
    number_b11_p = 0
    # C11:
    number_c11 = 0
    number_c11_n = 0
    number_c11_p_piminus = 0

    # C10:
    number_c10 = 0
    number_c10_2n = 0
    number_c10_n_p_piminus = 0
    # B10:
    number_b10 = 0
    number_b10_n_p = 0
    # Be10:
    number_be10 = 0
    number_be10_2p = 0
    number_be10_n_p_piplus = 0

    # C9:
    number_c9 = 0
    number_c9_3n = 0
    number_c9_2n_p_piminus = 0
    # B9:
    number_b9 = 0
    number_b9_n_d = 0
    number_b9_2n_p = 0
    number_b9_n_2p_piminus = 0
    number_b9_3n_piplus = 0
    # Be9:
    number_be9 = 0
    number_be9_n_2p = 0
    number_be9_p_d = 0
    number_be9_2n_p_piplus = 0
    # Li9:
    number_li9 = 0
    number_li9_3p = 0
    number_li9_n_2p_piplus = 0

    # B8:
    number_b8 = 0
    number_b8_3n_p = 0
    number_b8_2n_d = 0
    # Be8:
    number_be8 = 0
    number_be8_n_p_d = 0
    number_be8_2n_2p = 0
    number_be8_n_3p_piminus = 0
    number_be8_n_he3 = 0
    number_be8_p_t = 0
    # Li8:
    number_li8 = 0
    number_li8_n_3p = 0
    # He8:
    number_he8 = 0
    number_he8_4p = 0

    # B7:
    number_b7 = 0
    # Be7:
    number_be7 = 0
    number_be7_n_p_t = 0
    number_be7_3n_2p = 0
    number_be7_n_alpha = 0
    number_be7_2n_p_d = 0
    # Li7:
    number_li7 = 0
    number_li7_n_p_he3 = 0
    number_li7_2n_3p = 0
    number_li7_n_2p_d = 0
    number_li7_p_alpha = 0
    number_li7_n_alpha_piplus = 0
    # He7:
    number_he7 = 0
    number_he7_n_4p = 0

    # Be6:
    number_be6 = 0
    number_be6_3n_p_d = 0
    number_be6_2n_2d = 0
    number_be6_4n_2p = 0
    number_be6_2n_p_t = 0
    # Li6:
    number_li6 = 0
    number_li6_n_p_alpha = 0
    number_li6_2n_2p_d = 0
    number_li6_n_2p_t = 0
    number_li6_3n_3p = 0
    number_li6_2n_p_he3 = 0
    number_li6_2n_alpha_piplus = 0
    number_li6_2n_4p_piminus = 0
    number_li6_n_2p_he3_piminus = 0
    # He6:
    number_he6 = 0
    number_he6_2n_4p = 0
    number_he6_n_3p_d = 0
    number_he6_n_2p_he3 = 0
    # H6:
    number_h6 = 0
    number_h6_n_5p = 0

    # Li5:
    number_li5 = 0
    number_li5_2n_p_2d = 0
    number_li5_2n_p_alpha = 0
    number_li5_n_alpha_d = 0
    number_li5_3n_p_he3 = 0
    number_li5_3n_2p_d = 0
    number_li5_2n_he3_d = 0
    # He5:
    number_he5 = 0
    number_he5_n_3p_t = 0
    number_he5_n_2p_alpha = 0
    number_he5_n_2p_2d = 0
    number_he5_3n_4p = 0
    # H5:
    number_h5 = 0
    number_h5_2n_5p = 0

    # Li4:
    number_li4 = 0
    number_li4_3n_p_alpha = 0
    # He4:
    number_he4 = 0
    number_he4_n_2p_t_d = 0
    number_he4_3n_2p_he3 = 0
    number_he4_n_alpha_he3 = 0
    number_he4_2n_p_he3_d = 0
    number_he4_4n_4p = 0
    number_he4_2n_2he3 = 0
    # H4:
    number_h4 = 0
    number_h4_n_3p_alpha = 0
    number_h4_3n_5p = 0

    # FALSE channel:
    number_false_channel = 0

    """ target """
    # number of events with C12 as target:
    number_c12 = 0
    # number of events without C12 as target (integer):
    number_no_c12 = 0

    # loop over all entries of the array:
    for index in range(number_entries):

        # check, if target is C12 (PDG ID = 1000060120):
        if target_pdg[index] == 1000060120:
            # increment number_c12:
            number_c12 += 1

            # get the final PDG IDs of this event (final_pdg_event is array):
            final_pdg_event = particle_pdg[index]

            # check if n_particles[index] is equal to len(final_pdg_event:
            if n_particles[index] != len(final_pdg_event):
                sys.exit("ERROR: number of particles is NOT equal to length of particle_pdg")

            # preallocate number of single particles:
            num_b11 = 0
            num_c11 = 0
            num_c10 = 0
            num_b10 = 0
            num_be10 = 0
            num_c9 = 0
            num_b9 = 0
            num_be9 = 0
            num_li9 = 0
            num_c8 = 0
            num_b8 = 0
            num_be8 = 0
            num_li8 = 0
            num_he8 = 0
            num_b7 = 0
            num_be7 = 0
            num_li7 = 0
            num_he7 = 0
            num_be6 = 0
            num_li6 = 0
            num_he6 = 0
            num_h6 = 0
            num_li5 = 0
            num_he5 = 0
            num_alpha = 0
            num_he3 = 0
            num_triton = 0
            num_deuteron = 0
            num_proton = 0
            num_neutron = 0
            num_gamma = 0
            num_piplus = 0
            num_piminus = 0

            # loop over final_pdg_event to check the PDG of the particles:
            for index1 in range(int(n_particles[index])):
                if final_pdg_event[index1] == 1000050110:
                    num_b11 += 1
                elif final_pdg_event[index1] == 1000060110:
                    num_c11 += 1
                elif final_pdg_event[index1] == 1000060100:
                    num_c10 += 1
                elif final_pdg_event[index1] == 1000050100:
                    num_b10 += 1
                elif final_pdg_event[index1] == 1000040100:
                    num_be10 += 1
                elif final_pdg_event[index1] == 1000060090:
                    num_c9 += 1
                elif final_pdg_event[index1] == 1000050090:
                    num_b9 += 1
                elif final_pdg_event[index1] == 1000040090:
                    num_be9 += 1
                elif final_pdg_event[index1] == 1000030090:
                    num_li9 += 1
                elif final_pdg_event[index1] == 1000060080:
                    num_c8 += 1
                elif final_pdg_event[index1] == 1000050080:
                    num_b8 += 1
                elif final_pdg_event[index1] == 1000040080:
                    num_be8 += 1
                elif final_pdg_event[index1] == 1000030080:
                    num_li8 += 1
                elif final_pdg_event[index1] == 1000020080:
                    num_he8 += 1
                elif final_pdg_event[index1] == 1000050070:
                    num_b7 += 1
                elif final_pdg_event[index1] == 1000040070:
                    num_be7 += 1
                elif final_pdg_event[index1] == 1000030070:
                    num_li7 += 1
                elif final_pdg_event[index1] == 1000020070:
                    num_he7 += 1
                elif final_pdg_event[index1] == 1000040060:
                    num_be6 += 1
                elif final_pdg_event[index1] == 1000030060:
                    num_li6 += 1
                elif final_pdg_event[index1] == 1000020060:
                    num_he6 += 1
                elif final_pdg_event[index1] == 1000010060:
                    num_h6 += 1
                elif final_pdg_event[index1] == 1000030050:
                    num_li5 += 1
                elif final_pdg_event[index1] == 1000020050:
                    num_he5 += 1
                elif final_pdg_event[index1] == 1000020040:
                    num_alpha += 1
                elif final_pdg_event[index1] == 1000020030:
                    num_he3 += 1
                elif final_pdg_event[index1] == 1000010030:
                    num_triton += 1
                elif final_pdg_event[index1] == 1000010020:
                    num_deuteron += 1
                elif final_pdg_event[index1] == 2212:
                    num_proton += 1
                elif final_pdg_event[index1] == 2112:
                    num_neutron += 1
                elif final_pdg_event[index1] == 22:
                    num_gamma += 1
                elif final_pdg_event[index1] == 211:
                    num_piplus += 1
                elif final_pdg_event[index1] == -211:
                    num_piminus += 1
                else:
                    print("particle not yet included: PDG = {0:.0f}".format(final_pdg_event[index1]))

            if num_b11 == 1:
                # event with B11:
                number_b11 += 1
                if (num_neutron == num_piplus == 1 and num_proton == num_alpha == num_he3 == num_triton
                        == num_deuteron == num_piminus == 0):
                    number_b11_n_piplus += 1
                elif (num_proton == 1 and num_neutron == num_alpha == num_he3 == num_triton == num_deuteron
                      == num_piplus == num_piminus == 0):
                    number_b11_p += 1
                else:
                    print("B11 channel")

            elif num_c11 == 1:
                # event with C11:
                number_c11 += 1
                if (num_neutron == 1 and num_proton == num_alpha == num_he3 == num_triton == num_deuteron == num_piplus
                        == num_piminus == 0):
                    number_c11_n += 1
                elif (num_proton == num_piminus == 1 and num_neutron == num_alpha == num_he3 == num_triton
                      == num_deuteron == num_piplus == 0):
                    number_c11_p_piminus += 1
                else:
                    print("C11 channel")

            elif num_c10 == 1:
                # event with C10:
                number_c10 += 1
                if (num_neutron == 2 and num_proton == num_alpha == num_he3 == num_triton == num_deuteron == num_piplus
                        == num_piminus == 0):
                    number_c10_2n += 1
                elif (num_neutron == num_proton == num_piminus == 1 and num_alpha == num_he3 == num_triton
                      == num_deuteron == num_piplus == 0):
                    number_c10_n_p_piminus += 1
                else:
                    print("C10 channel")

            elif num_b10 == 1:
                # event with B10:
                number_b10 += 1
                if (num_neutron == num_proton == 1 and num_alpha == num_he3 == num_triton == num_deuteron == num_piplus
                        == num_piminus == 0):
                    number_b10_n_p += 1
                else:
                    print("B10 channel")

            elif num_be10 == 1:
                # event with Be10:
                number_be10 += 1
                if (num_proton == 2 and num_neutron == num_alpha == num_he3 == num_triton == num_deuteron == num_piplus
                        == num_piminus == 0):
                    number_be10_2p += 1
                elif (num_neutron == num_proton == num_piplus == 1 and num_alpha == num_he3 == num_triton
                      == num_deuteron == num_piminus == 0):
                    number_be10_n_p_piplus += 1
                else:
                    print("Be10 channel")

            elif num_c9 == 1:
                # event with C9:
                number_c9 += 1
                if (num_neutron == 3 and num_proton == num_alpha == num_he3 == num_triton == num_deuteron == num_piplus
                        == num_piminus == 0):
                    number_c9_3n += 1
                elif (num_neutron == 2 and num_proton == num_piminus == 1 and num_alpha == num_he3 == num_triton
                      == num_deuteron == num_piplus == 0):
                    number_c9_2n_p_piminus += 1
                else:
                    print("C9 channel")

            elif num_b9 == 1:
                # event with B9:
                number_b9 += 1
                if (num_neutron == 2 and num_proton == 1 and num_alpha == num_he3 == num_triton == num_deuteron
                        == num_piplus == num_piminus == 0):
                    number_b9_2n_p += 1
                elif (num_neutron == num_deuteron == 1 and num_proton == num_alpha == num_he3 == num_triton
                      == num_piplus == num_piminus == 0):
                    number_b9_n_d += 1
                elif (num_proton == 2 and num_neutron == num_piminus == 1 and num_alpha == num_he3 == num_triton
                      == num_deuteron == num_piplus == 0):
                    number_b9_n_2p_piminus += 1
                elif (num_neutron == 3 and num_piplus == 1 and num_proton == num_alpha == num_he3 == num_triton
                      == num_deuteron == num_piminus == 0):
                    number_b9_3n_piplus += 1
                else:
                    print("B9 channel")

            elif num_be9 == 1:
                # event with Be9:
                number_be9 += 1
                if (num_neutron == 1 and num_proton == 2 and num_alpha == num_he3 == num_triton == num_deuteron
                        == num_piplus == num_piminus == 0):
                    number_be9_n_2p += 1
                elif (num_proton == num_deuteron == 1 and num_neutron == num_alpha == num_he3 == num_triton
                      == num_piplus == num_piminus == 0):
                    number_be9_p_d += 1
                elif (num_neutron == 2 and num_proton == num_piplus == 1 and num_alpha == num_he3 == num_triton
                      == num_deuteron == num_piminus == 0):
                    number_be9_2n_p_piplus += 1
                else:
                    print("Be9 channel")

            elif num_li9 == 1:
                # event with Li9:
                number_li9 += 1
                if (num_proton == 3 and num_neutron == num_alpha == num_he3 == num_triton == num_deuteron == num_piplus
                        == num_piminus == 0):
                    number_li9_3p += 1
                elif (num_proton == 2 and num_neutron == num_piplus == 1 and num_alpha == num_he3 == num_triton
                        == num_deuteron == num_piminus == 0):
                    number_li9_n_2p_piplus += 1
                else:
                    print("Li9 channel")

            elif num_b8 == 1:
                # event with B8:
                number_b8 += 1
                if (num_neutron == 3 and num_proton == 1 and num_alpha == num_he3 == num_triton == num_deuteron
                        == num_piplus == num_piminus == 0):
                    number_b8_3n_p += 1
                elif (num_neutron == 2 and num_deuteron == 1 and num_proton == num_alpha == num_he3 == num_triton
                        == num_piplus == num_piminus == 0):
                    number_b8_2n_d += 1
                else:
                    print("B8 channel")

            elif num_be8 == 1:
                # event with Be8:
                number_be8 += 1
                if (num_neutron == num_proton == num_deuteron == 1 and num_alpha == num_he3 == num_triton
                        == num_piplus == num_piminus == 0):
                    number_be8_n_p_d += 1
                elif (num_neutron == num_proton == 2 and num_alpha == num_he3 == num_triton == num_deuteron
                      == num_piplus == num_piminus == 0):
                    number_be8_2n_2p += 1
                    print("index = {0:d}".format(index))
                    interesting_index = index
                elif (num_proton == 3 and num_neutron == num_piminus == 1 and num_alpha == num_he3 == num_triton
                      == num_deuteron == num_piplus == 0):
                    number_be8_n_3p_piminus += 1
                elif (num_neutron == num_he3 == 1 and num_proton == num_alpha == num_triton == num_deuteron
                      == num_piplus == num_piminus == 0):
                    number_be8_n_he3 += 1
                elif (num_proton == num_triton == 1 and num_neutron == num_alpha == num_he3 == num_deuteron
                        == num_piplus == num_piminus == 0):
                    number_be8_p_t += 1
                else:
                    print("Be8 channel")

            elif num_li8 == 1:
                # event with Li8:
                number_li8 += 1
                if (num_neutron == 1 and num_proton == 3 and num_alpha == num_he3 == num_triton == num_deuteron
                        == num_piplus == num_piminus == 0):
                    number_li8_n_3p += 1
                else:
                    print("Li8 channel")

            elif num_be7 == 1:
                # event with Be7:
                number_be7 += 1
                if (num_neutron == num_proton == num_triton == 1 and num_alpha == num_he3 == num_deuteron
                        == num_piplus == num_piminus == 0):
                    number_be7_n_p_t += 1
                elif (num_neutron == 3 and num_proton == 2 and num_alpha == num_he3 == num_triton == num_deuteron
                      == num_piplus == num_piminus == 0):
                    number_be7_3n_2p += 1
                elif (num_neutron == num_alpha == 1 and num_proton == num_he3 == num_triton == num_deuteron
                      == num_piplus == num_piminus == 0):
                    number_be7_n_alpha += 1
                elif (num_neutron == 2 and num_proton == num_deuteron == 1 and num_alpha == num_he3 == num_triton
                      == num_piplus == num_piminus == 0):
                    number_be7_2n_p_d += 1
                else:
                    print("Be7 channel")

            elif num_li7 == 1:
                # event with Li7:
                number_li7 += 1
                if (num_neutron == num_proton == num_he3 == 1 and num_alpha == num_triton == num_deuteron
                        == num_piplus == num_piminus == 0):
                    number_li7_n_p_he3 += 1
                elif (num_neutron == 2 and num_proton == 3 and num_alpha == num_he3 == num_triton == num_deuteron
                      == num_piplus == num_piminus == 0):
                    number_li7_2n_3p += 1
                elif (num_proton == 2 and num_neutron == num_deuteron == 1 and num_alpha == num_he3 == num_triton
                        == num_piplus == num_piminus == 0):
                    number_li7_n_2p_d += 1
                elif (num_proton == num_alpha == 1 and num_neutron == num_he3 == num_triton == num_deuteron
                        == num_piplus == num_piminus == 0):
                    number_li7_p_alpha += 1
                elif (num_neutron == num_alpha == num_piplus == 1 and num_proton == num_he3 == num_triton
                      == num_deuteron == num_piminus == 0):
                    number_li7_n_alpha_piplus += 1
                else:
                    print("Li7 channel")

            elif num_li6 == 1:
                # event with Li6:
                number_li6 += 1
                if (num_neutron == num_proton == num_alpha == 1 and num_he3 == num_triton == num_deuteron
                        == num_piplus == num_piminus == 0):
                    number_li6_n_p_alpha += 1
                elif (num_neutron == num_proton == 2 and num_deuteron == 1 and num_alpha == num_he3 == num_triton
                        == num_piplus == num_piminus == 0):
                    number_li6_2n_2p_d += 1
                elif (num_proton == 2 and num_neutron == num_triton == 1 and num_alpha == num_he3 == num_deuteron
                        == num_piplus == num_piminus == 0):
                    number_li6_n_2p_t += 1
                elif (num_neutron == num_proton == 3 and num_alpha == num_he3 == num_triton == num_deuteron
                        == num_piplus == num_piminus == 0):
                    number_li6_3n_3p += 1
                elif (num_neutron == 2 and num_proton == num_he3 == 1 and num_alpha == num_triton == num_deuteron
                        == num_piplus == num_piminus == 0):
                    number_li6_2n_p_he3 += 1
                elif (num_neutron == 2 and num_alpha == num_piplus == 1 and num_proton == num_he3 == num_triton
                      == num_deuteron == num_piminus == 0):
                    number_li6_2n_alpha_piplus += 1
                elif (num_proton == 2 and num_neutron == num_he3 == num_piminus == 1 and num_alpha == num_triton
                      == num_deuteron == num_piplus == 0):
                    number_li6_n_2p_he3_piminus += 1
                else:
                    print("Li6 channel")

            else:
                # channels, where no TALYS simulation of the deexcitation of the isotope exist:
                if (num_neutron == num_proton == 2 and num_alpha == num_he3 == num_triton == num_deuteron
                        == num_piplus == num_piminus == 0):
                    # isotope Be8:
                    number_be8 += 1
                    number_be8_2n_2p += 1

                elif (num_proton == 4 and num_neutron == num_alpha == num_he3 == num_triton == num_deuteron
                        == num_piplus == num_piminus == 0):
                    # isotope He8:
                    number_he8 += 1
                    number_he8_4p += 1

                elif (num_neutron == 1 and num_proton == 4 and num_alpha == num_he3 == num_triton == num_deuteron
                        == num_piplus == num_piminus == 0):
                    # isotope He7:
                    number_he7 += 1
                    number_he7_n_4p += 1

                elif (num_neutron == 3 and num_proton == num_deuteron == 1 and num_alpha == num_he3 == num_triton
                        == num_piplus == num_piminus == 0):
                    # isotope Be6:
                    number_be6 += 1
                    number_be6_3n_p_d += 1
                elif (num_neutron == num_deuteron == 2 and num_proton == num_alpha == num_he3 == num_triton
                        == num_piplus == num_piminus == 0):
                    # isotope Be6:
                    number_be6 += 1
                    number_be6_2n_2d += 1
                elif (num_neutron == 4 and num_proton == 2 and num_alpha == num_he3 == num_triton == num_deuteron
                        == num_piplus == num_piminus == 0):
                    # isotope Be6:
                    number_be6 += 1
                    number_be6_4n_2p += 1
                elif (num_neutron == 2 and num_proton == num_triton == 1 and num_alpha == num_he3 == num_deuteron
                        == num_piplus == num_piminus == 0):
                    # isotope Be6:
                    number_be6 += 1
                    number_be6_2n_p_t += 1

                elif (num_neutron == 2 and num_proton == 4 and num_piminus == 1 and num_alpha == num_he3 == num_triton
                        == num_deuteron == num_piplus == 0):
                    # isotope Be6:
                    number_li6 += 1
                    number_li6_2n_4p_piminus += 1
                elif (num_neutron == num_proton == 3 and num_alpha == num_he3 == num_triton == num_deuteron
                      == num_piplus == num_piminus == 0):
                    # isotope Be6:
                    number_li6 += 1
                    number_li6_3n_3p += 1

                elif (num_neutron == 2 and num_proton == 4 and num_alpha == num_he3 == num_triton == num_deuteron
                        == num_piplus == num_piminus == 0):
                    # isotope He6:
                    number_he6 += 1
                    number_he6_2n_4p += 1
                elif (num_proton == 3 and num_neutron == num_deuteron == 1 and num_alpha == num_he3 == num_triton
                        == num_piplus == num_piminus == 0):
                    # isotope He6:
                    number_he6 += 1
                    number_he6_n_3p_d += 1
                elif (num_proton == 2 and num_neutron == num_he3 == 1 and num_alpha == num_triton == num_deuteron
                        == num_piplus == num_piminus == 0):
                    # isotope He6:
                    number_he6 += 1
                    number_he6_n_2p_he3 += 1

                elif (num_neutron == 1 and num_proton == 5 and num_alpha == num_he3 == num_triton == num_deuteron
                      == num_piplus == num_piminus == 0):
                    # isotope H6:
                    number_h6 += 1
                    number_h6_n_5p += 1

                elif (num_neutron == num_deuteron == 2 and num_proton == 1 and num_alpha == num_he3 == num_triton
                      == num_piplus == num_piminus == 0):
                    # isotope Li5:
                    number_li5 += 1
                    number_li5_2n_p_2d += 1
                elif (num_neutron == 2 and num_proton == num_alpha == 1 and num_he3 == num_triton == num_deuteron
                      == num_piplus == num_piminus == 0):
                    # isotope Li5:
                    number_li5 += 1
                    number_li5_2n_p_alpha += 1
                elif (num_neutron == num_alpha == num_deuteron == 1 and num_proton == num_he3 == num_triton
                      == num_piplus == num_piminus == 0):
                    # isotope Li5:
                    number_li5 += 1
                    number_li5_n_alpha_d += 1
                elif (num_neutron == 3 and num_proton == num_he3 == 1 and num_alpha == num_triton == num_deuteron
                      == num_piplus == num_piminus == 0):
                    # isotope Li5:
                    number_li5 += 1
                    number_li5_3n_p_he3 += 1
                elif (num_neutron == 3 and num_proton == 2 and num_deuteron == 1 and num_alpha == num_he3 == num_triton
                      == num_piplus == num_piminus == 0):
                    # isotope Li5:
                    number_li5 += 1
                    number_li5_3n_2p_d += 1
                elif (num_neutron == 2 and num_he3 == num_deuteron == 1 and num_proton == num_alpha == num_triton
                      == num_piplus == num_piminus == 0):
                    # isotope Li5:
                    number_li5 += 1
                    number_li5_2n_he3_d += 1

                elif (num_proton == 3 and num_neutron == num_triton == 1 and num_alpha == num_he3 == num_deuteron
                      == num_piplus == num_piminus == 0):
                    # isotope He5:
                    number_he5 += 1
                    number_he5_n_3p_t += 1
                elif (num_proton == 2 and num_neutron == num_alpha == 1 and num_he3 == num_triton == num_deuteron
                      == num_piplus == num_piminus == 0):
                    # isotope He5:
                    number_he5 += 1
                    number_he5_n_2p_alpha += 1
                elif (num_proton == num_deuteron == 2 and num_neutron == 1 and num_alpha == num_he3 == num_triton
                      == num_piplus == num_piminus == 0):
                    # isotope He5:
                    number_he5 += 1
                    number_he5_n_2p_2d += 1
                elif (num_neutron == 3 and num_proton == 4 and num_alpha == num_he3 == num_triton == num_deuteron
                      == num_piplus == num_piminus == 0):
                    # isotope He5:
                    number_he5 += 1
                    number_he5_3n_4p += 1

                elif (num_neutron == 2 and num_proton == 5 and num_alpha == num_he3 == num_triton == num_deuteron
                      == num_piplus == num_piminus == 0):
                    # isotope H5:
                    number_h5 += 1
                    number_h5_2n_5p += 1

                elif (num_neutron == 3 and num_proton == num_alpha == 1 and num_he3 == num_triton == num_deuteron
                      == num_piplus == num_piminus == 0):
                    # isotope Li4:
                    number_li4 += 1
                    number_li4_3n_p_alpha += 1

                elif (num_proton == 2 and num_neutron == num_triton == num_deuteron == 1 and num_alpha == num_he3
                      == num_piplus == num_piminus == 0):
                    # isotope He4:
                    number_he4 += 1
                    number_he4_n_2p_t_d += 1
                elif (num_neutron == 3 and num_proton == 2 and num_he3 == 1 and num_alpha == num_triton == num_deuteron
                      == num_piplus == num_piminus == 0):
                    # isotope He4:
                    number_he4 += 1
                    number_he4_3n_2p_he3 += 1
                elif (num_neutron == 1 and num_alpha == 2 and num_proton == num_he3 == num_triton == num_deuteron
                      == num_piplus == num_piminus == 0):
                    # isotope He4:
                    number_he4 += 1
                    number_he4_n_alpha_he3 += 1
                elif (num_proton == 2 and num_neutron == num_alpha == num_deuteron == 1 and num_he3 == num_triton
                      == num_piplus == num_piminus == 0):
                    # isotope He4:
                    number_he4 += 1
                    number_he4_n_2p_t_d += 1
                elif (num_neutron == 2 and num_proton == num_alpha == num_deuteron == 1 and num_he3 == num_triton
                      == num_piplus == num_piminus == 0):
                    # isotope He4 (actually He3 is the isotope, but take He4 because it is heavier):
                    number_he4 += 1
                    number_he4_2n_p_he3_d += 1
                elif (num_neutron == 2 and num_proton == num_he3 == num_deuteron == 1 and num_alpha == num_triton
                      == num_piplus == num_piminus == 0):
                    # isotope He4:
                    number_he4 += 1
                    number_he4_2n_p_he3_d += 1
                elif (num_neutron == num_proton == 4 and num_alpha == num_he3 == num_triton == num_deuteron
                      == num_piplus == num_piminus == 0):
                    # isotope He4:
                    number_he4 += 1
                    number_he4_4n_4p += 1
                elif (num_neutron == 2 and num_alpha == num_he3 == 1 and num_proton == num_triton == num_deuteron
                      == num_piplus == num_piminus == 0):
                    # isotope He4:
                    number_he4 += 1
                    number_he4_2n_2he3 += 1

                elif (num_proton == 3 and num_neutron == num_alpha == 1 and num_he3 == num_triton == num_deuteron
                      == num_piplus == num_piminus == 0):
                    # isotope H4:
                    number_h4 += 1
                    number_h4_n_3p_alpha += 1
                elif (num_proton == 5 and num_neutron == 3 and num_alpha == num_he3 == num_triton == num_deuteron
                      == num_piplus == num_piminus == 0):
                    # isotope H4:
                    number_h4 += 1
                    number_h4_3n_5p += 1

                elif (num_neutron == 1 and num_proton == 6 and num_alpha == num_he3 == num_triton == num_deuteron
                      == num_piplus == num_piminus == 0):
                    # WRONG channel: 5 neutrons are left -> no isotope:
                    number_false_channel += 1
                elif (num_neutron == 2 and num_proton == 6 and num_alpha == num_he3 == num_triton == num_deuteron
                      == num_piplus == num_piminus == 0):
                    # WRONG channel: 4 neutrons are left -> no isotope:
                    number_false_channel += 1
                elif (num_neutron == 3 and num_proton == 6 and num_alpha == num_he3 == num_triton == num_deuteron
                      == num_piplus == num_piminus == 0):
                    # WRONG channel: 3 neutrons are left -> no isotope:
                    number_false_channel += 1

                else:
                    print("new channel")

        else:
            # target is not C12:
            number_no_c12 += 1

    return (number_c12, number_no_c12,
            number_b11, number_b11_n_piplus, number_b11_p,
            number_c11, number_c11_n, number_c11_p_piminus,
            number_c10, number_c10_2n, number_c10_n_p_piminus,
            number_b10, number_b10_n_p,
            number_be10, number_be10_2p, number_be10_n_p_piplus,
            number_c9, number_c9_3n, number_c9_2n_p_piminus,
            number_b9, number_b9_n_d, number_b9_2n_p, number_b9_n_2p_piminus, number_b9_3n_piplus,
            number_be9, number_be9_n_2p, number_be9_p_d, number_be9_2n_p_piplus,
            number_li9, number_li9_3p, number_li9_n_2p_piplus,
            number_b8, number_b8_3n_p, number_b8_2n_d,
            number_be8, number_be8_n_p_d, number_be8_2n_2p, number_be8_n_3p_piminus, number_be8_n_he3, number_be8_p_t,
            number_li8, number_li8_n_3p,
            number_he8, number_he8_4p,
            number_b7,
            number_be7, number_be7_n_p_t, number_be7_3n_2p, number_be7_n_alpha, number_be7_2n_p_d,
            number_li7, number_li7_n_p_he3, number_li7_2n_3p, number_li7_n_2p_d, number_li7_p_alpha,
            number_li7_n_alpha_piplus,
            number_he7, number_he7_n_4p,
            number_be6, number_be6_3n_p_d, number_be6_2n_2d, number_be6_4n_2p, number_be6_2n_p_t,
            number_li6, number_li6_n_p_alpha, number_li6_2n_2p_d, number_li6_n_2p_t, number_li6_3n_3p,
            number_li6_2n_p_he3, number_li6_2n_alpha_piplus, number_li6_2n_4p_piminus, number_li6_n_2p_he3_piminus,
            number_he6, number_he6_2n_4p, number_he6_n_3p_d, number_he6_n_2p_he3,
            number_h6, number_h6_n_5p,
            number_li5, number_li5_2n_p_2d, number_li5_2n_p_alpha, number_li5_n_alpha_d, number_li5_3n_p_he3,
            number_li5_3n_2p_d, number_li5_2n_he3_d,
            number_he5, number_he5_n_3p_t, number_he5_n_2p_alpha, number_he5_n_2p_2d, number_he5_3n_4p,
            number_h5, number_h5_2n_5p,
            number_li4, number_li4_3n_p_alpha,
            number_he4, number_he4_n_2p_t_d, number_he4_3n_2p_he3, number_he4_n_alpha_he3, number_he4_2n_p_he3_d,
            number_he4_4n_4p, number_he4_2n_2he3,
            number_h4, number_h4_n_3p_alpha, number_h4_3n_5p,
            number_false_channel, interesting_index)


def get_interaction_channel(channel_id, isotope_pdg, target_pdg):
    """
    function to get the different NC interaction channels from the interaction of neutrinos with the target particles

    :param channel_id: ID of the NC interaction channel, defines the product particles of the NC interactions
    (array of float)
    :param isotope_pdg: PDG ID of the isotope, that is created via the NC interaction (before deexcitation)
    (array of float)
    :param target_pdg: PDG ID of the target particle (array of float)

    :return:
    """

    # get the number of entries of the array (integer):
    number_entries = len(channel_id)

    """ preallocate arrays and variables: """
    """ B11 """
    # number of interaction channel: nu + C12 -> B11 + p (integer):
    number_c12_b11_p = 0
    # number of interaction channel: nu + C12 -> B11 + n + pi_plus (integer):
    number_c12_b11_n_piplus = 0
    # number of interaction channel: nu + C12 -> B11 + n + pi_minus + 2*pi_plus (integer):
    number_c12_b11_n_piminus_2piplus = 0
    # number of interaction channel: nu + C12 -> B11 + p + pi_minus + pi_plus (integer):
    number_c12_b11_p_piminus_piplus = 0
    # number of interaction channel: nu + C12 -> B11 + p + 2*pi_minus + 2*pi_plus (integer):
    number_c12_b11_p_2piminus_2piplus = 0
    # number of interaction channel: nu + C12 -> B11 + pi_plus (integer):
    number_c12_b11_piplus = 0
    # number of OTHER interaction channels: nu + C12 -> B11 + ...:
    number_c12_b11_other = 0

    """ C11 """
    # number of interaction channel: nu + C12 -> C11 + n (integer):
    number_c12_c11_n = 0
    # number of interaction channel: nu + C12 -> C11 + p + pi_minus (integer):
    number_c12_c11_p_piminus = 0
    # number of interaction channel: nu + C12 -> C11 + n + pi_minus + pi_plus (integer):
    number_c12_c11_n_piminus_piplus = 0
    # number of interaction channel: nu + C12 -> C11 + p + 2*pi_minus + pi_plus (integer):
    number_c12_c11_p_2piminus_piplus = 0
    # number of interaction channel: nu + C12 -> C11 + p + 3*pi_minus + 2*pi_plus (integer):
    number_c12_c11_p_3piminus_2piplus = 0
    # number of interaction channel: nu + C12 -> C11 + n + 2*pi_minus + 2*pi_plus (integer):
    number_c12_c11_n_2piminus_2piplus = 0
    # number of OTHER interaction channels: nu + C12 -> C11 + ...:
    number_c12_c11_other = 0

    """ B10 """
    # number of interaction channel: nu + C12 -> B10 + p + n (integer):
    number_c12_b10_p_n = 0
    # number of interaction channel: nu + C12 -> B10 + 2p + pi_minus (integer):
    number_c12_b10_2p_piminus = 0
    # number of interaction channel: nu + C12 -> B10 + p + n + pi_minus + pi_plus (integer):
    number_c12_b10_p_n_piminus_piplus = 0
    # number of interaction channel: nu + C12 -> B10 + 2n + pi_plus (integer):
    number_c12_b10_2n_piplus = 0
    # number of interaction channel: nu + C12 -> B10 + 2n + pi_minus + 2*pi_plus (integer):
    number_c12_b10_2n_piminus_2piplus = 0
    # number of interaction channel: nu + C12 -> B10 + 2p + 2*pi_minus + pi_plus (integer):
    number_c12_b10_2p_2piminus_piplus = 0
    # number of interaction channel: nu + C12 -> B10 + 2p + 3*pi_minus + 2*pi_plus (integer):
    number_c12_b10_2p_3piminus_2piplus = 0
    # number of interaction channel: nu + C12 -> B10 + p + n + 2*pi_minus + 2*pi_plus (integer):
    number_c12_b10_p_n_2piminus_2piplus = 0
    # number of OTHER interaction channels: nu + C12 -> B10 + ...:
    number_c12_b10_other = 0

    """ C10 """
    # number of interaction channel: nu + C12 -> C10 + 2n (integer):
    number_c12_c10_2n = 0
    # number of interaction channel: nu + C12 -> C10 + p + n + pi_minus (integer):
    number_c12_c10_p_n_piminus = 0
    # number of interaction channel: nu + C12 -> C10 + p + n + 2*pi_minus + pi_plus (integer):
    number_c12_c10_p_n_2piminus_piplus = 0
    # number of interaction channel: nu + C12 -> C10 + 2n + pi_minus + pi_plus (integer):
    number_c12_c10_2n_piminus_piplus = 0
    # number of interaction channel: nu + C12 -> C10 + 2p + 2*pi_minus (integer):
    number_c12_c10_2p_2piminus = 0
    # number of OTHER interaction channels: nu + C12 -> C10 + ...:
    number_c12_c10_other = 0

    """ Be10 """
    # number of interaction channel: nu + C12 -> Be10 + 2*p (integer):
    number_c12_be10_2p = 0
    # number of interaction channel: nu + C12 -> Be10 + p + n + pi_plus (integer):
    number_c12_be10_p_n_piplus = 0
    # number of interaction channel: nu + C12 -> Be10 + p + n + pi_minus + 2*pi_plus (integer):
    number_c12_be10_p_n_piminus_2piplus = 0
    # number of interaction channel: nu + C12 -> Be10 + 2*p + pi_minus + pi_plus (integer):
    number_c12_be10_2p_piminus_piplus = 0
    # number of interaction channel: nu + C12 -> Be10 + 2*n + 2*pi_plus (integer):
    number_c12_be10_2n_2piplus = 0
    # number of interaction channel: nu + C12 -> Be10 + p + n + 2*pi_minus + 3*pi_plus (integer):
    number_c12_be10_p_n_2piminus_3piplus = 0
    # number of interaction channel: nu + C12 -> Be10 + 2*p + 2*pi_minus + 2*pi_plus (integer):
    number_c12_be10_2p_2piminus_2piplus = 0
    # number of interaction channel: nu + C12 -> Be10 + 2*p + 3*pi_minus + 3*pi_plus (integer):
    number_c12_be10_2p_3piminus_3piplus = 0
    # number of OTHER interaction channels: nu + C12 -> Be10 + ...:
    number_c12_be10_other = 0

    """ B9 """
    # number of interaction channel: nu + C12 -> B9 + p + 2n (integer):
    number_c12_b9_p_2n = 0
    # number of interaction channel: nu + C12 -> B9 + p + 2n + pi_minus + pi_plus (integer):
    number_c12_b9_p_2n_piminus_piplus = 0
    # number of interaction channel: nu + C12 -> B9 + 2p + n + 3*pi_minus + 2*pi_plus (integer):
    number_c12_b9_2p_n_3piminus_2piplus = 0
    # number of interaction channel: nu + C12 -> B9 + 2p + n + pi_minus (integer):
    number_c12_b9_2p_n_piminus = 0
    # number of interaction channel: nu + C12 -> B9 + 3n + pi_plus (integer):
    number_c12_b9_3n_piplus = 0
    # number of interaction channel: nu + C12 -> B9 + p + 2n + 2*pi_minus + 2*pi_plus (integer):
    number_c12_b9_p_2n_2piminus_2piplus = 0
    # number of interaction channel: nu + C12 -> B9 + 2p + n + 2*pi_minus+ pi_plus (integer):
    number_c12_b9_2p_n_2piminus_piplus = 0
    # number of OTHER interaction channels: nu + C12 -> B9 + ...
    number_c12_b9_other = 0

    """ Be9 """
    # number of interaction channel: nu + C12 -> Be9 + 2*p + n (integer):
    number_c12_be9_2p_n = 0
    # number of interaction channel: nu + C12 -> Be9 + p + 2n + pi_plus (integer):
    number_c12_be9_p_2n_piplus = 0
    # number of interaction channel: nu + C12 -> Be9 + 3p + pi_minus (integer):
    number_c12_be9_3p_piminus = 0
    # number of interaction channel: nu + C12 -> Be9 + p + 2n + pi_minus + 2*pi_plus (integer):
    number_c12_be9_p_2n_piminus_2piplus = 0
    # number of interaction channel: nu + C12 -> Be9 + 2p + n + pi_minus + pi_plus (integer):
    number_c12_be9_2p_n_piminus_piplus = 0
    # number of interaction channel: nu + C12 -> Be9 + 2p + n + 3*pi_minus + 3*pi_plus (integer):
    number_c12_be9_2p_n_3piminus_3piplus = 0
    # number of interaction channel: nu + C12 -> Be9 + 2p + n + 2*pi_minus + 2*pi_plus (integer):
    number_c12_be9_2p_n_2piminus_2piplus = 0
    # number of interaction channel: nu + C12 -> Be9 + 3n + 2*pi_plus (integer):
    number_c12_be9_3n_2piplus = 0
    # number of interaction channel: nu + C12 -> Be9 + 3p + 2*pi_minus + pi_plus (integer):
    number_c12_be9_3p_2piminus_piplus = 0
    # number of OTHER interaction channels: nu + C12 -> Be9 + ...:
    number_c12_be9_other = 0

    """ Be8 """
    # number of interaction channel: nu + C12 -> Be8 + 2p + 2n (integer):
    number_c12_be8_2p_2n = 0
    # number of interaction channel: nu + C12 -> Be8 + 3p + n + pi_minus (integer):
    number_c12_be8_3p_n_piminus = 0
    # number of interaction channel: nu + C12 -> Be8 + p + 3n + pi_plus (integer):
    number_c12_be8_p_3n_piplus = 0
    # number of interaction channel: nu + C12 -> Be8 + 2p + 2n + 2*pi_minus + 2*pi_plus (integer):
    number_c12_be8_2p_2n_2piminus_2piplus = 0
    # number of interaction channel: nu + C12 -> Be8 + 4n + 2*pi_plus (integer):
    number_c12_be8_4n_2piplus = 0
    # number of interaction channel: nu + C12 -> Be8 + 2p + 2n + pi_minus * pi_plus (integer):
    number_c12_be8_2p_2n_piminus_piplus = 0
    # number of interaction channel: nu + C12 -> Be8 + 3p + n + 2*pi_minus + pi_plus (integer):
    number_c12_be8_3p_n_2piminus_piplus = 0
    # number of interaction channel: nu + C12 -> Be8 + 4p + 2*pi_minus (integer):
    number_c12_be8_4p_2piminus = 0
    # number of OTHER interaction channels: nu + C12 -> Be8 + ...:
    number_c12_be8_other = 0

    """ C9 """
    # number of interaction channel: nu + C12 -> C9 + p + 2n + pi_minus (integer):
    number_c12_c9_p_2n_piminus = 0
    # number of interaction channel: nu + C12 -> C9 + 3n (integer):
    number_c12_c9_3n = 0
    # number of interaction channel: nu + C12 -> C9 + 2p + n + 2*pi_minus (integer):
    number_c12_c9_2p_n_2piminus = 0
    # number of interaction channel: nu + C12 -> C9 + 3n + 2*pi_minus + 2*pi_plus (integer):
    number_c12_c9_3n_2piminus_2piplus = 0
    # number of OTHER interaction channels: nu + C12 -> C9 + ...:
    number_c12_c9_other = 0

    """ Be7 """
    # number of interaction channel: nu + C12 -> Be7 + 2p + 3n (integer):
    number_c12_be7_2p_3n = 0
    # number of interaction channel: nu + C12 -> Be7 + p + 4n + pi_plus (integer):
    number_c12_be7_p_4n_piplus = 0
    # number of interaction channel: nu + C12 -> Be7 + 2p + 3n + 2*pi_minus + 2*pi_plus (integer):
    number_c12_be7_2p_3n_2piminus_2piplus = 0
    # number of interaction channel: nu + C12 -> Be7 + 3p + 2n + pi_minus (integer):
    number_c12_be7_3p_2n_piminus = 0
    # number of interaction channel: nu + C12 -> Be7 + 4p + n + 2*pi_minus (integer):
    number_c12_be7_4p_n_2piminus = 0
    # number of interaction channel: nu + C12 -> Be7 + 3p + 2n + 2*pi_minus + pi_plus (integer):
    number_c12_be7_3p_2n_2piminus_piplus = 0
    # number of OTHER interaction channels: nu + C12 -> Be7 + ...:
    number_c12_be7_other = 0

    """ Li6 """
    # number of interaction channel: nu + C12 -> Li6 + 3p + 3n (integer):
    number_c12_li6_3p_3n = 0
    # number of interaction channel: nu + C12 -> Li6 + 2p + 4n + pi_plus (integer):
    number_c12_li6_2p_4n_piplus = 0
    # number of interaction channel: nu + C12 -> Li6 + 5p + n + 2*pi_minus (integer):
    number_c12_li6_5p_n_2piminus = 0
    # number of interaction channel: nu + C12 -> Li6 + 2p + 4n + pi_minus + 2*pi_plus (integer):
    number_c12_li6_2p_4n_piminus_2piplus = 0
    # number of interaction channel: nu + C12 -> Li6 + 4p + 2n + pi_minus (integer):
    number_c12_li6_4p_2n_piminus = 0
    # number of interaction channel: nu + C12 -> Li6 + 3p + 3n + pi_minus + pi_plus (integer):
    number_c12_li6_3p_3n_piminus_piplus = 0
    # number of OTHER interaction channels: nu + C12 -> Li6 + ...:
    number_c12_li6_other = 0

    """ Li8 """
    # number of interaction channel: nu + C12 -> Li8 + 3p + n (integer):
    number_c12_li8_3p_n = 0
    # number of interaction channel: nu + C12 -> Li8 + 4p + pi_minus (integer):
    number_c12_li8_4p_piminus = 0
    # number of interaction channel: nu + C12 -> Li8 + 4p + 2*pi_minus + pi_plus (integer):
    number_c12_li8_4p_2piminus_piplus = 0
    # number of interaction channel: nu + C12 -> Li8 + 2p + 2n + pi_plus (integer):
    number_c12_li8_2p_2n_piplus = 0
    # number of interaction channel: nu + C12 -> Li8 + 3p + n + pi_minus + pi_plus (integer):
    number_c12_li8_3p_n_piminus_piplus = 0
    # number of OTHER interaction channels: nu + C12 -> Li8 + ...:
    number_c12_li8_other = 0

    """ Li7 """
    # number of interaction channel: nu + C12 -> Li7 + 2p + 3n + pi_plus (integer):
    number_c12_li7_2p_3n_piplus = 0
    # number of interaction channel: nu + C12 -> Li7 + 4p + n + pi_minus (integer):
    number_c12_li7_4p_n_piminus = 0
    # number of interaction channel: nu + C12 -> Li7 + 3p + 2n (integer):
    number_c12_li7_3p_2n = 0
    # number of interaction channel: nu + C12 -> Li7 + 3p + 2n + pi_minus + pi_plus (integer):
    number_c12_li7_3p_2n_piminus_piplus = 0
    # number of interaction channel: nu + C12 -> Li7 + 4p + n + 2*pi_minus + pi_plus (integer):
    number_c12_li7_4p_n_2piminus_piplus = 0
    # number of interaction channel: nu + C12 -> Li7 + 2p + 3n + pi_minus + 2*pi_plus (integer):
    number_c12_li7_2p_3n_piminus_2piplus = 0
    # number of OTHER interaction channels: nu + C12 -> Li7 + ...:
    number_c12_li7_other = 0

    """ B8 """
    # number of interaction channel: nu + C12 -> B8 + p + 3n (integer):
    number_c12_b8_p_3n = 0
    # number of interaction channel: nu + C12 -> B8 + p + 3n + pi_minus + pi_plus (integer):
    number_c12_b8_p_3n_piminus_piplus = 0
    # number of interaction channel: nu + C12 -> B8 + 2p + 2n + 2*pi_minus + pi_plus (integer):
    number_c12_b8_2p_2n_2piminus_piplus = 0
    # number of interaction channel: nu + C12 -> B8 + 2p + 2n + pi_minus (integer):
    number_c12_b8_2p_2n_piminus = 0
    # number of interaction channel: nu + C12 -> B8 + 4n + pi_plus (integer):
    number_c12_b8_4n_piplus = 0
    # number of OTHER interaction channels: nu + C12 -> B8 + ...:
    number_c12_b8_other = 0

    """ Li9 """
    # number of interaction channel: nu + C12 -> Li9 + 2p + n + pi_plus (integer):
    number_c12_li9_2p_n_piplus = 0
    # number of interaction channel: nu + C12 -> Li9 + 3p (integer).
    number_c12_li9_3p = 0
    # number of interaction channel: nu + C12 -> Li9 + 3p + pi_minus + pi_plus (integer):
    number_c12_li9_3p_piminus_piplus = 0
    # number of interaction channel: nu + C12 -> Li9 + 2p + n + pi_minus + 2*pi_plus (integer):
    number_c12_li9_2p_n_piminus_2piplus = 0
    # number of interaction channel: nu + C12 -> Li9 + p + 2n + pi_minus + 3*pi_plus (integer):
    number_c12_li9_p_2n_piminus_3piplus = 0
    # number of interaction channel: nu + C12 -> Li9 + ...:
    number_c12_li9_other = 0

    """ C8 """
    # number of interaction channel: nu + C12 -> C8 + 4n (integer):
    number_c12_c8_4n = 0
    # number of interaction channel: nu + C12 -> C8 + 4n + N*pi_minus + N*pi_minus (integer):
    number_c12_c8_4n_other = 0

    """ He8 """
    # number of interaction channel: nu + C12 -> He8 + 4p (integer):
    number_c12_he8_4p = 0
    # number of interaction channel: nu + C12 -> He8 + 4p + N*pi_minus + N*pi_minus (integer):
    number_c12_he8_4p_other = 0

    """ B7 """
    # number of interaction channel: nu + C12 -> B7 + p + 4n (integer):
    number_c12_b7_p_4n = 0
    # number of interaction channel: nu + C12 -> B7 + p + 4n + N*pi_minus + N*pi_minus (integer):
    number_c12_b7_p_4n_other = 0

    """ He7 """
    # number of interaction channel: nu + C12 -> He7 + 4p + n (integer):
    number_c12_he7_4p_n = 0
    # number of interaction channel: nu + C12 -> He7 + 4p + n + N*pi_minus + N*pi_minus (integer):
    number_c12_he7_4p_n_other = 0

    """ H7 """
    # number of interaction channel: nu + C12 -> H7 + 4p + n (integer):
    number_c12_h7_5p = 0
    # number of interaction channel: nu + C12 -> H7 + 4p + n + N*pi_minus + N*pi_minus (integer):
    number_c12_h7_5p_other = 0

    """ Be6 """
    # number of interaction channel: nu + C12 -> Be6 + 2p + 4n (integer):
    number_c12_be6_2p_4n = 0
    # number of interaction channel: nu + C12 -> Be6 + 2p + 4n + N*pi_minus + N*pi_minus (integer):
    number_c12_be6_2p_4n_other = 0

    """ He6 """
    # number of interaction channel: nu + C12 -> He6 + 4p + 2n (integer):
    number_c12_he6_4p_2n = 0
    # number of interaction channel: nu + C12 -> He6 + 4p + 2n + N*pi_minus + N*pi_minus (integer):
    number_c12_he6_4p_2n_other = 0

    """ H6 """
    # number of interaction channel: nu + C12 -> H6 + 5p + n (integer):
    number_c12_h6_5p_n = 0
    # number of interaction channel: nu + C12 -> H6 + 5p + n + N*pi_minus + N*pi_minus (integer):
    number_c12_h6_5p_n_other = 0

    """ C12 """
    # number of interaction channels: nu + C12 -> nu + C12 + other particles (like pi_minus, pi_plus, kaon_minus,
    # koan_plus and so on):
    number_c12_c12 = 0

    """ no isotope (only protons, neutrons, pions): """
    # number of interaction channels with NO isotope (only proton, neutrons, pions):
    number_c12_noiso = 0
    # number of interaction channels with no isotope (nu + C12 -> nu + 5p + 6n + pi_plus):
    number_c12_noiso_5p_6n = 0

    """ other channels, where C11 or B11 is produced: """
    number_c12_mass11u = 0

    """ other channels, where C10, B10 or Be10 is produced: """
    number_c12_mass10u = 0

    """ other channels, where C9, B9, Be9 or Li9 is produced:"""
    number_c12_mass9u = 0

    """ other channels, where C8, B8, Be8, Li8 or He8 is produced: """
    number_c12_mass8u = 0

    """ other channels, where B7, Be7, Li7, He7 or H7 is produced: """
    number_c12_mass7u = 0

    """ other channels, where Be6, Li6, He6 or H6 is produced: """
    number_c12_mass6u = 0

    """ other channels, where isotopes have mass <= 5u: """
    number_c12_mass5orless = 0

    """ faulty interaction channels: isotopes are missing (no isotopes with Z<3 and (N-Z)<3 are considered): """
    # number of interaction channels: nu + C12 -> nu + C12 + ...:
    number_c12_faulty = 0

    """ Other targets than C12: """
    # number of channels without C12 as target (integer):
    number_no_c12 = 0
    # number of elastic scattering interactions with protons: nu + p -> nu + p + ... (integer):
    number_es_p = 0
    # number of elastic scattering interactions with electrons: nu + electron -> nu + electron + ... (integer):
    number_es_e = 0
    # number of elastic scattering interactions with O16: nu + O16 -> nu + O16 + ... (integer):
    number_es_o16 = 0
    # number of elastic scattering interactions with N14: nu + N14 -> nu + N14 + ... (integer):
    number_es_n14 = 0
    # number of elastic scattering interactions with S32: nu + S32 -> nu + S32 + ... (integer):
    number_es_s32 = 0

    # loop over all entries of the array:
    for index in range(number_entries):

        # check, if target is C12 (PDG ID = 1000060120):
        if target_pdg[index] == 1000060120:

            # check the PDG ID of the created isotopes:
            if isotope_pdg[index] == 1000050110:
                # for B11:
                num_p, num_n, num_pi_minus, num_pi_plus = get_number_of_particles_of_channelid(channel_id[index])

                if num_p == 1 and num_n == 0 and num_pi_minus == 0 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> B11 + proton:
                    number_c12_b11_p = number_c12_b11_p + 1

                elif num_p == 0 and num_n == 1 and num_pi_minus == 0 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> B11 + n + pi_plus:
                    number_c12_b11_n_piplus = number_c12_b11_n_piplus + 1

                elif num_p == 0 and num_n == 1 and num_pi_minus == 1 and num_pi_plus == 2:
                    # interaction channel: nu + C12 -> B11 + n + pi_minus * 2*pi_plus:
                    number_c12_b11_n_piminus_2piplus = number_c12_b11_n_piminus_2piplus + 1

                elif num_p == 1 and num_n == 0 and num_pi_minus == 1 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> B11 + p + pi_minus + pi_plus:
                    number_c12_b11_p_piminus_piplus = number_c12_b11_p_piminus_piplus + 1

                elif num_p == 1 and num_n == 0 and num_pi_minus == 2 and num_pi_plus == 2:
                    # interaction channel: nu + C12 -> B11 + p + 2*pi_minus + 2*pi_plus:
                    number_c12_b11_p_2piminus_2piplus = number_c12_b11_p_2piminus_2piplus + 1

                elif num_p == 0 and num_n == 0 and num_pi_minus == 0 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> B11 + pi_plus:
                    number_c12_b11_piplus = number_c12_b11_piplus + 1

                else:
                    # interaction channels, that are not covered above, and channels with Kaon or Sigma-Baryon:
                    number_c12_b11_other = number_c12_b11_other + 1
                    # print("new interaction channel with nu + C12 -> B11: channel ID = {0:.0f}"
                    #       .format(channel_id[index]))


            elif isotope_pdg[index] == 1000060110:
                # for C11:
                num_p, num_n, num_pi_minus, num_pi_plus = get_number_of_particles_of_channelid(channel_id[index])

                if num_p == 0 and num_n == 1 and num_pi_minus == 0 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> C11 + n:
                    number_c12_c11_n = number_c12_c11_n + 1

                elif num_p == 1 and num_n == 0 and num_pi_minus == 1 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> C11 + p + pi_minus:
                    number_c12_c11_p_piminus = number_c12_c11_p_piminus + 1

                elif num_p == 0 and num_n == 1 and num_pi_minus == 1 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> C11 + n + pi_minus + pi_plus:
                    number_c12_c11_n_piminus_piplus = number_c12_c11_n_piminus_piplus + 1

                elif num_p == 1 and num_n == 0 and num_pi_minus == 2 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> C11 + p + 2*pi_minus + pi_plus:
                    number_c12_c11_p_2piminus_piplus = number_c12_c11_p_2piminus_piplus + 1

                elif num_p == 1 and num_n == 0 and num_pi_minus == 3 and num_pi_plus == 2:
                    # interaction channel: nu + C12 -> C11 + p + 3*pi_minus + 2*pi_plus:
                    number_c12_c11_p_3piminus_2piplus = number_c12_c11_p_3piminus_2piplus + 1

                elif num_p == 0 and num_n == 1 and num_pi_minus == 2 and num_pi_plus == 2:
                    # interaction channel: nu + C12 -> C11 + n + 2*pi_minus + 2*pi_plus:
                    number_c12_c11_n_2piminus_2piplus = number_c12_c11_n_2piminus_2piplus + 1

                else:
                    # interaction channels, that are not covered above, and channels with Kaon or Sigma-Baryon:
                    number_c12_c11_other = number_c12_c11_other + 1
                    # print("new interaction channel with nu + C12 -> C11: channel ID = {0:.0f}"
                    #       .format(channel_id[index]))


            elif isotope_pdg[index] == 1000050100:
                # for B10:
                num_p, num_n, num_pi_minus, num_pi_plus = get_number_of_particles_of_channelid(channel_id[index])

                if num_p == 1 and num_n == 1 and num_pi_minus == 0 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> B10 + p + n:
                    number_c12_b10_p_n = number_c12_b10_p_n + 1

                elif num_p == 2 and num_n == 0 and num_pi_minus == 1 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> B10 + 2*p + pi_minus:
                    number_c12_b10_2p_piminus = number_c12_b10_2p_piminus + 1

                elif num_p == 1 and num_n == 1 and num_pi_minus == 1 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> B10 + p + n + pi_minus + pi_plus:
                    number_c12_b10_p_n_piminus_piplus = number_c12_b10_p_n_piminus_piplus + 1

                elif num_p == 0 and num_n == 2 and num_pi_minus == 0 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> B10 + 2*n + pi_plus:
                    number_c12_b10_2n_piplus = number_c12_b10_2n_piplus + 1

                elif num_p == 0 and num_n == 2 and num_pi_minus == 1 and num_pi_plus == 2:
                    # interaction channel: nu + C12 -> B10 + 2*n + pi_minus + 2*pi_plus:
                    number_c12_b10_2n_piminus_2piplus = number_c12_b10_2n_piminus_2piplus + 1

                elif num_p == 2 and num_n == 0 and num_pi_minus == 2 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> B10 + 2*p + 2*pi_minus + pi_plus:
                    number_c12_b10_2p_2piminus_piplus = number_c12_b10_2p_2piminus_piplus + 1

                elif num_p == 2 and num_n == 0 and num_pi_minus == 3 and num_pi_plus == 2:
                    # interaction channel: nu + C12 -> B10 + 2*p + 3*pi_minus + 2*pi_plus:
                    number_c12_b10_2p_3piminus_2piplus = number_c12_b10_2p_3piminus_2piplus + 1

                elif num_p == 1 and num_n == 1 and num_pi_minus == 2 and num_pi_plus == 2:
                    # interaction channel: nu + C12 -> B10 + p + n + 2*pi_minus + 2*pi_plus:
                    number_c12_b10_p_n_2piminus_2piplus = number_c12_b10_p_n_2piminus_2piplus + 1

                else:
                    # interaction channels, that are not covered above, and channels with Kaon or Sigma-Baryon:
                    number_c12_b10_other = number_c12_b10_other + 1
                    # print("new interaction channel with nu + C12 -> B10: channel ID = {0:.0f}"
                    #       .format(channel_id[index]))


            elif isotope_pdg[index] == 1000060100:
                # for C10:
                num_p, num_n, num_pi_minus, num_pi_plus = get_number_of_particles_of_channelid(channel_id[index])

                if num_p == 0 and num_n == 2 and num_pi_minus == 0 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> C10 + 2n:
                    number_c12_c10_2n = number_c12_c10_2n + 1

                elif num_p == 1 and num_n == 1 and num_pi_minus == 1 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> C10 + p + n + pi_minus:
                    number_c12_c10_p_n_piminus = number_c12_c10_p_n_piminus + 1

                elif num_p == 1 and num_n == 1 and num_pi_minus == 2 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> C10 + p + n + 2*pi_minus + pi_plus:
                    number_c12_c10_p_n_2piminus_piplus = number_c12_c10_p_n_2piminus_piplus + 1

                elif num_p == 0 and num_n == 2 and num_pi_minus == 1 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> C10 + 2*n + pi_minus + pi_plus:
                    number_c12_c10_2n_piminus_piplus = number_c12_c10_2n_piminus_piplus + 1

                elif num_p == 2 and num_n == 0 and num_pi_minus == 2 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> C10 + 2*p + 2*pi_minus:
                    number_c12_c10_2p_2piminus = number_c12_c10_2p_2piminus + 1

                else:
                    # interaction channels, that are not covered above, and channels with Kaon or Sigma-Baryon:
                    number_c12_c10_other = number_c12_c10_other + 1
                    # print("new interaction channel with nu + C12 -> C10: channel ID = {0:.0f}"
                    #       .format(channel_id[index]))


            elif isotope_pdg[index] == 1000040100:
                # for Be10:
                num_p, num_n, num_pi_minus, num_pi_plus = get_number_of_particles_of_channelid(channel_id[index])

                if num_p == 2 and num_n == 0 and num_pi_minus == 0 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> Be10 + 2*p:
                    number_c12_be10_2p = number_c12_be10_2p + 1

                elif num_p == 1 and num_n == 1 and num_pi_minus == 0 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> Be10 + p + n + pi_plus:
                    number_c12_be10_p_n_piplus = number_c12_be10_p_n_piplus + 1

                elif num_p == 1 and num_n == 1 and num_pi_minus == 1 and num_pi_plus == 2:
                    # interaction channel: nu + C12 -> Be10 + p + n + pi_minus + 2*pi_plus:
                    number_c12_be10_p_n_piminus_2piplus = number_c12_be10_p_n_piminus_2piplus + 1

                elif num_p == 2 and num_n == 0 and num_pi_minus == 1 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> Be10 + 2*p + pi_minus + pi_plus:
                    number_c12_be10_2p_piminus_piplus = number_c12_be10_2p_piminus_piplus + 1

                elif num_p == 0 and num_n == 2 and num_pi_minus == 0 and num_pi_plus == 2:
                    # interaction channel: nu + C12 -> Be10 + 2n + 2*pi_plus:
                    number_c12_be10_2n_2piplus = number_c12_be10_2n_2piplus + 1

                elif num_p == 1 and num_n == 1 and num_pi_minus == 2 and num_pi_plus == 3:
                    # interaction channel: nu + C12 -> Be10 + p + n + 2*pi_minus + 3*pi_plus:
                    number_c12_be10_p_n_2piminus_3piplus = number_c12_be10_p_n_2piminus_3piplus + 1

                elif num_p == 2 and num_n == 0 and num_pi_minus == 2 and num_pi_plus == 2:
                    # interaction channel: nu + C12 -> Be10 + 2p + 2*pi_minus + 2*pi_plus:
                    number_c12_be10_2p_2piminus_2piplus = number_c12_be10_2p_2piminus_2piplus + 1

                elif num_p == 2 and num_n == 0 and num_pi_minus == 3 and num_pi_plus == 3:
                    # interaction channel: nu + C12 -> Be10 + 2p + 3*pi_minus + 3*pi_plus:
                    number_c12_be10_2p_3piminus_3piplus = number_c12_be10_2p_3piminus_3piplus + 1

                else:
                    # interaction channels, that are not covered above, and channels with Kaon or Sigma-Baryon:
                    number_c12_be10_other = number_c12_be10_other + 1
                    # print("new interaction channel with nu + C12 -> Be10: channel ID = {0:.0f}"
                    #       .format(channel_id[index]))


            elif isotope_pdg[index] == 1000050090:
                # for B9:
                num_p, num_n, num_pi_minus, num_pi_plus = get_number_of_particles_of_channelid(channel_id[index])

                if num_p == 1 and num_n == 2 and num_pi_minus == 0 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> B9 + p + 2*n:
                    number_c12_b9_p_2n = number_c12_b9_p_2n + 1

                elif num_p == 1 and num_n == 2 and num_pi_minus == 1 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> B9 + p + 2n + pi_minus + pi_plus:
                    number_c12_b9_p_2n_piminus_piplus = number_c12_b9_p_2n_piminus_piplus + 1

                elif num_p == 2 and num_n == 1 and num_pi_minus == 3 and num_pi_plus == 2:
                    # interaction channel: nu + C12 -> B9 + 2p + n + 3*pi_minus + 2*pi_plus:
                    number_c12_b9_2p_n_3piminus_2piplus = number_c12_b9_2p_n_3piminus_2piplus + 1

                elif num_p == 2 and num_n == 1 and num_pi_minus == 1 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> B9 + 2p + n + pi_minus:
                    number_c12_b9_2p_n_piminus = number_c12_b9_2p_n_piminus + 1

                elif num_p == 0 and num_n == 3 and num_pi_minus == 0 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> B9 + 3n + pi_plus:
                    number_c12_b9_3n_piplus = number_c12_b9_3n_piplus + 1

                elif num_p == 1 and num_n == 2 and num_pi_minus == 2 and num_pi_plus == 2:
                    # interaction channel: nu + C12 -> B9 + p + 2n + 2*pi_minus + 2*pi_plus:
                    number_c12_b9_p_2n_2piminus_2piplus = number_c12_b9_p_2n_2piminus_2piplus + 1

                elif num_p == 2 and num_n == 1 and num_pi_minus == 2 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> B9 + 2p + n + 2*pi_minus + pi_plus:
                    number_c12_b9_2p_n_2piminus_piplus = number_c12_b9_2p_n_2piminus_piplus + 1

                else:
                    # interaction channels, that are not covered above, and channels with Kaon or Sigma-Baryon:
                    number_c12_b9_other = number_c12_b9_other + 1
                    # print("new interaction channel with nu + C12 -> B9: channel ID = {0:.0f}"
                    # .format(channel_id[index]))


            elif isotope_pdg[index] == 1000040090:
                # for Be9:
                num_p, num_n, num_pi_minus, num_pi_plus = get_number_of_particles_of_channelid(channel_id[index])

                if num_p == 2 and num_n == 1 and num_pi_minus == 0 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> Be9 + 2*p + n:
                    number_c12_be9_2p_n = number_c12_be9_2p_n + 1

                elif num_p == 1 and num_n == 2 and num_pi_minus == 0 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> Be9 + p + 2*n + pi_plus:
                    number_c12_be9_p_2n_piplus = number_c12_be9_p_2n_piplus + 1

                elif num_p == 3 and num_n == 0 and num_pi_minus == 1 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> Be9 + 3p + pi_minus:
                    number_c12_be9_3p_piminus = number_c12_be9_3p_piminus + 1

                elif num_p == 1 and num_n == 2 and num_pi_minus == 1 and num_pi_plus == 2:
                    # interaction channel: nu + C12 -> Be9 + p + 2*n + pi_minus + 2*pi_plus:
                    number_c12_be9_p_2n_piminus_2piplus = number_c12_be9_p_2n_piminus_2piplus + 1

                elif num_p == 2 and num_n == 1 and num_pi_minus == 1 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> Be9 + 2*p + n + pi_minus + pi_plus:
                    number_c12_be9_2p_n_piminus_piplus = number_c12_be9_2p_n_piminus_piplus + 1

                elif num_p == 2 and num_n == 1 and num_pi_minus == 3 and num_pi_plus == 3:
                    # interaction channel: nu + C12 -> Be9 + 2*p + n + 3*pi_minus + 3*pi_plus:
                    number_c12_be9_2p_n_3piminus_3piplus = number_c12_be9_2p_n_3piminus_3piplus + 1

                elif num_p == 2 and num_n == 1 and num_pi_minus == 2 and num_pi_plus == 2:
                    # interaction channel: nu + C12 -> Be9 + 2*p + n + 2*pi_minus + 2*pi_plus:
                    number_c12_be9_2p_n_2piminus_2piplus = number_c12_be9_2p_n_2piminus_2piplus + 1

                elif num_p == 0 and num_n == 3 and num_pi_minus == 0 and num_pi_plus == 2:
                    # interaction channel: nu + C12 -> Be9 + 3*n + 2*pi_plus:
                    number_c12_be9_3n_2piplus = number_c12_be9_3n_2piplus + 1

                elif num_p == 3 and num_n == 0 and num_pi_minus == 2 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> Be9 + 3*p + 2*pi_minus + pi_plus:
                    number_c12_be9_3p_2piminus_piplus = number_c12_be9_3p_2piminus_piplus + 1

                else:
                    # interaction channels, that are not covered above, and channels with Kaon or Sigma-Baryon:
                    number_c12_be9_other = number_c12_be9_other + 1
                    # print("new interaction channel with nu + C12 -> Be9: channel ID = {0:.0f}"
                    #       .format(channel_id[index]))


            elif isotope_pdg[index] == 1000040080:
                # for Be8:
                num_p, num_n, num_pi_minus, num_pi_plus = get_number_of_particles_of_channelid(channel_id[index])

                if num_p == 2 and num_n == 2 and num_pi_minus == 0 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> Be8 + 2*p + 2*n:
                    number_c12_be8_2p_2n = number_c12_be8_2p_2n + 1

                elif num_p == 3 and num_n == 1 and num_pi_minus == 1 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> Be8 + 3*p + n + pi_minus:
                    number_c12_be8_3p_n_piminus = number_c12_be8_3p_n_piminus + 1

                elif num_p == 1 and num_n == 3 and num_pi_minus == 0 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> Be8 + p + 3*n + pi_plus:
                    number_c12_be8_p_3n_piplus = number_c12_be8_p_3n_piplus + 1

                elif num_p == 2 and num_n == 2 and num_pi_minus == 2 and num_pi_plus == 2:
                    # interaction channel: nu + C12 -> Be8 + 2*p + 2*n + 2*pi_minus + 2*pi_plus:
                    number_c12_be8_2p_2n_2piminus_2piplus = number_c12_be8_2p_2n_2piminus_2piplus + 1

                elif num_p == 0 and num_n == 4 and num_pi_minus == 0 and num_pi_plus == 2:
                    # interaction channel: nu + C12 -> Be8 + 4*n + 2*pi_plus:
                    number_c12_be8_4n_2piplus = number_c12_be8_4n_2piplus + 1

                elif num_p == 2 and num_n == 2 and num_pi_minus == 1 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> Be8 + 2*p + 2*n + pi_minus + pi_plus:
                    number_c12_be8_2p_2n_piminus_piplus = number_c12_be8_2p_2n_piminus_piplus + 1

                elif num_p == 3 and num_n == 1 and num_pi_minus == 2 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> Be8 + 3*p + n + 2*pi_minus + pi_plus:
                    number_c12_be8_3p_n_2piminus_piplus = number_c12_be8_3p_n_2piminus_piplus + 1

                elif num_p == 4 and num_n == 0 and num_pi_minus == 2 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> Be8 + 4*p + 2*pi_minus:
                    number_c12_be8_4p_2piminus = number_c12_be8_4p_2piminus + 1

                else:
                    # interaction channels, that are not covered above, and channels with Kaon or Sigma-Baryon:
                    number_c12_be8_other = number_c12_be8_other + 1
                    # print("new interaction channel with nu + C12 -> Be8: channel ID = {0:.0f}"
                    #       .format(channel_id[index]))


            elif isotope_pdg[index] == 1000060090:
                # for C9:
                num_p, num_n, num_pi_minus, num_pi_plus = get_number_of_particles_of_channelid(channel_id[index])

                if num_p == 1 and num_n == 2 and num_pi_minus == 1 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> C9 + p + 2*n + pi_minus:
                    number_c12_c9_p_2n_piminus = number_c12_c9_p_2n_piminus + 1

                elif num_p == 0 and num_n == 3 and num_pi_minus == 0 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> C9 + 3*n:
                    number_c12_c9_3n = number_c12_c9_3n + 1

                elif num_p == 2 and num_n == 1 and num_pi_minus == 2 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> C9 + 2*p + n + 2*pi_minus:
                    number_c12_c9_2p_n_2piminus = number_c12_c9_2p_n_2piminus + 1

                elif num_p == 0 and num_n == 3 and num_pi_minus == 2 and num_pi_plus == 2:
                    # interaction channel: nu + C12 -> C9 + 3*n + 2*pi_minus + 2*pi_plus:
                    number_c12_c9_3n_2piminus_2piplus = number_c12_c9_3n_2piminus_2piplus + 1

                else:
                    # interaction channels, that are not covered above, and channels with Kaon or Sigma-Baryon:
                    number_c12_c9_other = number_c12_c9_other + 1
                    # print("new interaction channel with nu + C12 -> C9: channel ID = {0:.0f}"
                    # .format(channel_id[index]))


            elif isotope_pdg[index] == 1000040070:
                # for Be7:
                num_p, num_n, num_pi_minus, num_pi_plus = get_number_of_particles_of_channelid(channel_id[index])

                if num_p == 2 and num_n == 3 and num_pi_minus == 0 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> Be7 + 2*p + 3*n:
                    number_c12_be7_2p_3n = number_c12_be7_2p_3n + 1

                elif num_p == 1 and num_n == 4 and num_pi_minus == 0 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> Be7 + p + 4*n + pi_plus:
                    number_c12_be7_p_4n_piplus = number_c12_be7_p_4n_piplus + 1

                elif num_p == 2 and num_n == 3 and num_pi_minus == 2 and num_pi_plus == 2:
                    # interaction channel: nu + C12 -> Be7 + 2*p + 3*n + 2*pi_minus + 2*pi_plus:
                    number_c12_be7_2p_3n_2piminus_2piplus = number_c12_be7_2p_3n_2piminus_2piplus + 1

                elif num_p == 3 and num_n == 2 and num_pi_minus == 1 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> Be7 + 3*p + 2*n + pi_minus:
                    number_c12_be7_3p_2n_piminus = number_c12_be7_3p_2n_piminus + 1

                elif num_p == 4 and num_n == 1 and num_pi_minus == 2 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> Be7 + 4*p + n + 2*pi_minus:
                    number_c12_be7_4p_n_2piminus = number_c12_be7_4p_n_2piminus + 1

                elif num_p == 3 and num_n == 2 and num_pi_minus == 2 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> Be7 + 3*p + 2*n + 2*pi_minus + pi_plus:
                    number_c12_be7_3p_2n_2piminus_piplus = number_c12_be7_3p_2n_2piminus_piplus + 1

                else:
                    # interaction channels, that are not covered above, and channels with Kaon or Sigma-Baryon:
                    number_c12_be7_other = number_c12_be7_other + 1
                    # print("new interaction channel with nu + C12 -> Be7: channel ID = {0:.0f}"
                    #       .format(channel_id[index]))


            elif isotope_pdg[index] == 1000030060:
                # for Li6:
                num_p, num_n, num_pi_minus, num_pi_plus = get_number_of_particles_of_channelid(channel_id[index])

                if num_p == 3 and num_n == 3 and num_pi_minus == 0 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> Li6 + 3*p + 3*n:
                    number_c12_li6_3p_3n = number_c12_li6_3p_3n + 1

                elif num_p == 2 and num_n == 4 and num_pi_minus == 0 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> Li6 + 2*p + 4*n + pi_plus:
                    number_c12_li6_2p_4n_piplus = number_c12_li6_2p_4n_piplus + 1

                elif num_p == 5 and num_n == 1 and num_pi_minus == 2 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> Li6 + 5*p + n + 2*pi_minus:
                    number_c12_li6_5p_n_2piminus = number_c12_li6_5p_n_2piminus + 1

                elif num_p == 2 and num_n == 4 and num_pi_minus == 1 and num_pi_plus == 2:
                    # interaction channel: nu + C12 -> Li6 + 2*p + 4*n + pi_minus + 2*pi_plus:
                    number_c12_li6_2p_4n_piminus_2piplus = number_c12_li6_2p_4n_piminus_2piplus + 1

                elif num_p == 4 and num_n == 2 and num_pi_minus == 1 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> Li6 + 4*p + 2*n + pi_minus:
                    number_c12_li6_4p_2n_piminus = number_c12_li6_4p_2n_piminus + 1

                elif num_p == 3 and num_n == 3 and num_pi_minus == 1 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> Li6 + 3*p + 3*n + pi_minus + pi_plus:
                    number_c12_li6_3p_3n_piminus_piplus = number_c12_li6_3p_3n_piminus_piplus + 1

                else:
                    # interaction channels, that are not covered above, and channels with Kaon or Sigma-Baryon:
                    number_c12_li6_other = number_c12_li6_other + 1
                    # print("new interaction channel with nu + C12 -> Li6: channel ID = {0:.0f}"
                    #       .format(channel_id[index]))


            elif isotope_pdg[index] == 1000030080:
                # for Li8:
                num_p, num_n, num_pi_minus, num_pi_plus = get_number_of_particles_of_channelid(channel_id[index])

                if num_p == 3 and num_n == 1 and num_pi_minus == 0 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> Li8 + 3*p + n:
                    number_c12_li8_3p_n = number_c12_li8_3p_n + 1

                elif num_p == 4 and num_n == 0 and num_pi_minus == 1 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> Li8 + 4*p + pi_minus:
                    number_c12_li8_4p_piminus = number_c12_li8_4p_piminus + 1

                elif num_p == 4 and num_n == 0 and num_pi_minus == 2 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> Li8 + 4*p + 2*pi_minus + pi_plus:
                    number_c12_li8_4p_2piminus_piplus = number_c12_li8_4p_2piminus_piplus + 1

                elif num_p == 2 and num_n == 2 and num_pi_minus == 0 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> Li8 + 2*p + 2*n + pi_plus:
                    number_c12_li8_2p_2n_piplus = number_c12_li8_2p_2n_piplus + 1

                elif num_p == 3 and num_n == 1 and num_pi_minus == 1 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> Li8 + 3*p + n + pi_minus + pi_plus:
                    number_c12_li8_3p_n_piminus_piplus = number_c12_li8_3p_n_piminus_piplus + 1

                else:
                    # interaction channels, that are not covered above, and channels with Kaon or Sigma-Baryon:
                    number_c12_li8_other = number_c12_li8_other + 1
                    # print("new interaction channel with nu + C12 -> Li8: channel ID = {0:.0f}"
                    #       .format(channel_id[index]))


            elif isotope_pdg[index] == 1000030070:
                # for Li7:
                num_p, num_n, num_pi_minus, num_pi_plus = get_number_of_particles_of_channelid(channel_id[index])

                if num_p == 2 and num_n == 3 and num_pi_minus == 0 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> Li7 + 2*p + 3*n + pi_plus:
                    number_c12_li7_2p_3n_piplus = number_c12_li7_2p_3n_piplus + 1

                elif num_p == 4 and num_n == 1 and num_pi_minus == 1 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> Li7 + 4*p + n + pi_minus:
                    number_c12_li7_4p_n_piminus = number_c12_li7_4p_n_piminus + 1

                elif num_p == 3 and num_n == 2 and num_pi_minus == 0 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> Li7 + 3*p + 2*n:
                    number_c12_li7_3p_2n = number_c12_li7_3p_2n + 1

                elif num_p == 3 and num_n == 2 and num_pi_minus == 1 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> Li7 + 3*p + 2*n + pi_minus + pi_plus:
                    number_c12_li7_3p_2n_piminus_piplus = number_c12_li7_3p_2n_piminus_piplus + 1

                elif num_p == 4 and num_n == 1 and num_pi_minus == 2 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> Li7 + 4*p + n + 2*pi_minus + pi_plus:
                    number_c12_li7_4p_n_2piminus_piplus = number_c12_li7_4p_n_2piminus_piplus + 1

                elif num_p == 2 and num_n == 3 and num_pi_minus == 1 and num_pi_plus == 2:
                    # interaction channel: nu + C12 -> Li7 + 2*p + 3*n + pi_minus + 2*pi_plus:
                    number_c12_li7_2p_3n_piminus_2piplus = number_c12_li7_2p_3n_piminus_2piplus + 1

                else:
                    # interaction channels, that are not covered above, and channels with Kaon or Sigma-Baryon:
                    number_c12_li7_other = number_c12_li7_other + 1
                    # print("new interaction channel with nu + C12 -> Li7: channel ID = {0:.0f}"
                    #       .format(channel_id[index]))


            elif isotope_pdg[index] == 1000050080:
                # for B8:
                num_p, num_n, num_pi_minus, num_pi_plus = get_number_of_particles_of_channelid(channel_id[index])

                if num_p == 1 and num_n == 3 and num_pi_minus == 0 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> B8 + p + 3*n:
                    number_c12_b8_p_3n = number_c12_b8_p_3n + 1

                elif num_p == 1 and num_n == 3 and num_pi_minus == 1 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> B8 + p + 3*n + pi_minus + pi_plus:
                    number_c12_b8_p_3n_piminus_piplus = number_c12_b8_p_3n_piminus_piplus + 1

                elif num_p == 2 and num_n == 2 and num_pi_minus == 2 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> B8 + 2*p + 2*n + 2*pi_minus + pi_plus:
                    number_c12_b8_2p_2n_2piminus_piplus = number_c12_b8_2p_2n_2piminus_piplus + 1

                elif num_p == 2 and num_n == 2 and num_pi_minus == 1 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> B8 + 2*p + 2*n + pi_minus:
                    number_c12_b8_2p_2n_piminus = number_c12_b8_2p_2n_piminus + 1

                elif num_p == 0 and num_n == 4 and num_pi_minus == 0 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> B8 + 4*n + pi_plus:
                    number_c12_b8_4n_piplus = number_c12_b8_4n_piplus + 1

                else:
                    # interaction channels, that are not covered above, and channels with Kaon or Sigma-Baryon:
                    number_c12_b8_other = number_c12_b8_other + 1
                    # print("new interaction channel with nu + C12 -> B8: channel ID = {0:.0f}"
                    #       .format(channel_id[index]))


            elif isotope_pdg[index] == 1000030090:
                # for Li9:
                num_p, num_n, num_pi_minus, num_pi_plus = get_number_of_particles_of_channelid(channel_id[index])

                if num_p == 2 and num_n == 1 and num_pi_minus == 0 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> Li9 + 2*p + n + pi_plus:
                    number_c12_li9_2p_n_piplus = number_c12_li9_2p_n_piplus + 1

                elif num_p == 3 and num_n == 0 and num_pi_minus == 0 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> Li9 + 3*p:
                    number_c12_li9_3p = number_c12_li9_3p + 1

                elif num_p == 3 and num_n == 0 and num_pi_minus == 1 and num_pi_plus == 1:
                    # interaction channel: nu + C12 -> Li9 + 3*p + pi_minus + pi_plus:
                    number_c12_li9_3p_piminus_piplus = number_c12_li9_3p_piminus_piplus + 1

                elif num_p == 2 and num_n == 1 and num_pi_minus == 1 and num_pi_plus == 2:
                    # interaction channel: nu + C12 -> Li9 + 2*p + n + pi_minus + 2*pi_plus:
                    number_c12_li9_2p_n_piminus_2piplus = number_c12_li9_2p_n_piminus_2piplus + 1

                elif num_p == 1 and num_n == 2 and num_pi_minus == 1 and num_pi_plus == 3:
                    # interaction channel: nu + C12 -> Li9 + p + 2*n + pi_minus + 3*pi_plus:
                    number_c12_li9_p_2n_piminus_3piplus = number_c12_li9_p_2n_piminus_3piplus + 1

                else:
                    # interaction channels, that are not covered above, and channels with Kaon or Sigma-Baryon:
                    number_c12_li9_other = number_c12_li9_other + 1
                    # print("new interaction channel with nu + C12 -> Li9: channel ID = {0:.0f}"
                    #       .format(channel_id[index]))


            elif isotope_pdg[index] == 1000060080:
                # for C8:
                num_p, num_n, num_pi_minus, num_pi_plus = get_number_of_particles_of_channelid(channel_id[index])

                if num_p == 0 and num_n == 4 and num_pi_minus == 0 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> C8 + 4n:
                    number_c12_c8_4n = number_c12_c8_4n + 1

                else:
                    # interaction channels like above, BUT with N pi_minus and pi_plus
                    # (C8 + 4n + N*pi_minus + N*pi_plus):
                    number_c12_c8_4n_other = number_c12_c8_4n_other + 1


            elif isotope_pdg[index] == 1000020080:
                # for He8:
                num_p, num_n, num_pi_minus, num_pi_plus = get_number_of_particles_of_channelid(channel_id[index])

                if num_p == 4 and num_n == 0 and num_pi_minus == 0 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> He8 + 4p:
                    number_c12_he8_4p = number_c12_he8_4p + 1

                else:
                    # interaction channels like above, BUT with N pi_minus and pi_plus
                    # (He8 + 4p + N*pi_minus + N*pi_plus):
                    number_c12_he8_4p_other = number_c12_he8_4p_other + 1


            elif isotope_pdg[index] == 1000050070:
                # for B7:
                num_p, num_n, num_pi_minus, num_pi_plus = get_number_of_particles_of_channelid(channel_id[index])

                if num_p == 1 and num_n == 4 and num_pi_minus == 0 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> B7 + p + 4n:
                    number_c12_b7_p_4n = number_c12_b7_p_4n + 1

                else:
                    # interaction channels like above, BUT with N pi_minus and pi_plus
                    # (B7 + p + 4n + N*pi_minus + N*pi_plus):
                    number_c12_b7_p_4n_other = number_c12_b7_p_4n_other + 1


            elif isotope_pdg[index] == 1000020070:
                # for He7:
                num_p, num_n, num_pi_minus, num_pi_plus = get_number_of_particles_of_channelid(channel_id[index])

                if num_p == 4 and num_n == 1 and num_pi_minus == 0 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> He7 + 4p + n:
                    number_c12_he7_4p_n = number_c12_he7_4p_n + 1

                else:
                    # interaction channels like above, BUT with N pi_minus and pi_plus
                    # (He7 + 4p + n + N*pi_minus + N*pi_plus):
                    number_c12_he7_4p_n_other = number_c12_he7_4p_n_other + 1


            elif isotope_pdg[index] == 1000010070:
                # for H7:
                num_p, num_n, num_pi_minus, num_pi_plus = get_number_of_particles_of_channelid(channel_id[index])

                if num_p == 5 and num_n == 0 and num_pi_minus == 0 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> H7 + 5p:
                    number_c12_h7_5p = number_c12_h7_5p + 1

                else:
                    # interaction channels like above, BUT with N pi_minus and pi_plus
                    # (H7 + 5p + N*pi_minus + N*pi_plus):
                    number_c12_h7_5p_other = number_c12_h7_5p_other + 1


            elif isotope_pdg[index] == 1000040060:
                # for Be6:
                num_p, num_n, num_pi_minus, num_pi_plus = get_number_of_particles_of_channelid(channel_id[index])

                if num_p == 2 and num_n == 4 and num_pi_minus == 0 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> Be6 + 2p + 4n:
                    number_c12_be6_2p_4n = number_c12_be6_2p_4n + 1

                else:
                    # interaction channels like above, BUT with N pi_minus and pi_plus
                    # (Be6 + 2p + 4n + N*pi_minus + N*pi_plus):
                    number_c12_be6_2p_4n_other = number_c12_be6_2p_4n_other + 1


            elif isotope_pdg[index] == 1000020060:
                # for He6:
                num_p, num_n, num_pi_minus, num_pi_plus = get_number_of_particles_of_channelid(channel_id[index])

                if num_p == 4 and num_n == 2 and num_pi_minus == 0 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> He6 + 4p + 2n:
                    number_c12_he6_4p_2n = number_c12_he6_4p_2n + 1

                else:
                    # interaction channels like above, BUT with N pi_minus and pi_plus
                    # (He6 + 4p + 2n + N*pi_minus + N*pi_plus):
                    number_c12_he6_4p_2n_other = number_c12_he6_4p_2n_other + 1


            elif isotope_pdg[index] == 1000010060:
                # for H6:
                num_p, num_n, num_pi_minus, num_pi_plus = get_number_of_particles_of_channelid(channel_id[index])

                if num_p == 5 and num_n == 1 and num_pi_minus == 0 and num_pi_plus == 0:
                    # interaction channel: nu + C12 -> H6 + 5p + n:
                    number_c12_h6_5p_n = number_c12_h6_5p_n + 1

                else:
                    # interaction channels like above, BUT with N pi_minus and pi_plus
                    # (H6 + 5p + n + N*pi_minus + N*pi_plus):
                    number_c12_h6_5p_n_other = number_c12_h6_5p_n_other + 1


            elif isotope_pdg[index] == 110000000:
                # other channels: nu + C12 -> nu + C11/B11 + ...:
                number_c12_mass11u = number_c12_mass11u + 1


            elif isotope_pdg[index] == 100000000:
                # other channels: nu + C12 -> nu + C10/B10/Be10 + ...:
                number_c12_mass10u = number_c12_mass10u + 1


            elif isotope_pdg[index] == 90000000:
                # other channels: nu + C12 -> nu + C9/B9/Be9/Li9 + ...:
                number_c12_mass9u = number_c12_mass9u + 1


            elif isotope_pdg[index] == 80000000:
                # other channels: nu + C12 -> nu + C8/B8/Be8/Li8/He8 + ...:
                number_c12_mass8u = number_c12_mass8u + 1


            elif isotope_pdg[index] == 70000000:
                # other channels: nu + C12 -> nu + B7/Be7/Li7/He7/H7 + ...:
                number_c12_mass7u = number_c12_mass7u + 1


            elif isotope_pdg[index] == 60000000:
                # other channels: nu + C12 -> nu + Be6/Li6/He6/H6 + ...:
                number_c12_mass6u = number_c12_mass6u + 1


            elif isotope_pdg[index] == 50000000:
                # other channels with isotopes with mass <= 5u:
                number_c12_mass5orless = number_c12_mass5orless + 1


            elif isotope_pdg[index] == 1000060120:
                # for C12 as product:
                number_c12_c12 = number_c12_c12 + 1


            elif isotope_pdg[index] == 0:
                # No isotope as product (only proton, neutron, pion, ...):
                # number of p, n, pi_minus, pi_plus after NC interaction:
                num_p, num_n, num_pi_minus, num_pi_plus = get_number_of_particles_of_channelid(channel_id[index])
                # number of p and n of the target particle:
                num_p_target, num_n_target = get_number_of_p_and_n_of_isotope(target_pdg[index])
                # mass number (sum of p and n) of the target particle:
                mass_number_target = num_p_target + num_n_target

                if mass_number_target == (num_p + num_n):
                    # possible interaction channel: nu + C12 -> nu + X*p + Y*n + ...:
                    number_c12_noiso = number_c12_noiso + 1

                elif num_p == 5 and num_n == 6 and (num_pi_plus - num_pi_minus) == 1:
                    # possible interaction channel: nu + C12 -> nu + 5p + 6n + (X*pi_plus - Y*pi_minus)):
                    number_c12_noiso_5p_6n = number_c12_noiso_5p_6n + 1

                else:
                    print("new interaction channel without isotope (nu + C12 -> nu + N*p + M*n + ...), isopdg = {0:.0f}"
                          ", channel ID = {1:.0f}".format(isotope_pdg[index], channel_id[index]))


            else:
                print("other isotope than expected: {0:.0f}, corresponding channel ID = {1:.0f}"
                      .format(isotope_pdg[index], channel_id[index]))


        else:
            # other target than C12

            # number of events with NO C12 as target:
            number_no_c12 = number_no_c12 + 1

            # check the channel ID of the interaction (either channel_id == 2 or channel_id == 3):
            if channel_id[index] == 2:
                # channel_id == 2: ES interaction: nu + p -> nu + p (maybe also pi_zero or gammas):
                number_es_p = number_es_p + 1

            elif channel_id[index] == 3:

                if target_pdg[index] == 11:
                    # electron as target: ES interaction: nu + electron -> nu + electron (maybe also pi_zero or gammas):
                    number_es_e = number_es_e + 1

                elif target_pdg[index] == 1000080160:
                    # O16 as target: ES interaction: nu + O16 -> nu + O16 (maybe also pi_zero or gammas):
                    number_es_o16 = number_es_o16 + 1

                elif target_pdg[index] == 1000070140:
                    # N14 as target: ES interaction: nu + N14 -> nu + N14 (maybe also pi_zero or gammas):
                    number_es_n14 = number_es_n14 + 1

                elif target_pdg[index] == 1000160320:
                    # S32 as target: ES interaction: nu + S32 -> nu + S32 (maybe also pi_zero or gammas):
                    number_es_s32 = number_es_s32 + 1

                else:
                    print("new target PDG for channel ID = 3: target PDG = {0:.0f}".format(target_pdg[index]))

            else:
                print("new channel ID with NO C12 target: target PDG = {0:.0f}, channel ID = {1:.0f}, "
                      "isotope PDG = {2:.0f}".format(target_pdg[index], channel_id[index], isotope_pdg[index]))


    """ calculate the fraction of the different NC interaction channels in PERCENT (float): """
    # B11:
    frac_c12_b11_p = float(number_c12_b11_p) / float(number_entries) * 100
    frac_c12_b11_n_piplus = float(number_c12_b11_n_piplus) / float(number_entries) * 100
    frac_c12_b11_n_piminus_2piplus = float(number_c12_b11_n_piminus_2piplus) / float(number_entries) * 100
    frac_c12_b11_p_piminus_piplus = float(number_c12_b11_p_piminus_piplus) / float(number_entries) * 100
    frac_c12_b11_p_2piminus_2piplus = float(number_c12_b11_p_2piminus_2piplus) / float(number_entries) * 100
    frac_c12_b11_piplus = float(number_c12_b11_piplus) / float(number_entries) * 100
    frac_c12_b11_other = float(number_c12_b11_other) / float(number_entries) * 100

    # C11:
    frac_c12_c11_n = float(number_c12_c11_n) / float(number_entries) * 100
    frac_c12_c11_p_piminus = float(number_c12_c11_p_piminus) / float(number_entries) * 100
    frac_c12_c11_n_piminus_piplus = float(number_c12_c11_n_piminus_piplus) / float(number_entries) * 100
    frac_c12_c11_p_2piminus_piplus = float(number_c12_c11_p_2piminus_piplus) / float(number_entries) * 100
    frac_c12_c11_p_3piminus_2piplus = float(number_c12_c11_p_3piminus_2piplus) / float(number_entries) * 100
    frac_c12_c11_n_2piminus_2piplus = float(number_c12_c11_n_2piminus_2piplus) / float(number_entries) * 100
    frac_c12_c11_other = float(number_c12_c11_other) / float(number_entries) * 100

    # B10:
    frac_c12_b10_p_n = float(number_c12_b10_p_n) / float(number_entries) * 100
    frac_c12_b10_2p_piminus = float(number_c12_b10_2p_piminus) / float(number_entries) * 100
    frac_c12_b10_p_n_piminus_piplus = float(number_c12_b10_p_n_piminus_piplus) / float(number_entries) * 100
    frac_c12_b10_2n_piplus = float(number_c12_b10_2n_piplus) / float(number_entries) * 100
    frac_c12_b10_2n_piminus_2piplus = float(number_c12_b10_2n_piminus_2piplus) / float(number_entries) * 100
    frac_c12_b10_2p_2piminus_piplus = float(number_c12_b10_2p_2piminus_piplus) / float(number_entries) * 100
    frac_c12_b10_2p_3piminus_2piplus = float(number_c12_b10_2p_3piminus_2piplus) / float(number_entries) * 100
    frac_c12_b10_p_n_2piminus_2piplus = float(number_c12_b10_p_n_2piminus_2piplus) / float(number_entries) * 100
    frac_c12_b10_other = float(number_c12_b10_other) / float(number_entries) * 100

    # C10:
    frac_c12_c10_2n = float(number_c12_c10_2n) / float(number_entries) * 100
    frac_c12_c10_p_n_piminus = float(number_c12_c10_p_n_piminus) / float(number_entries) * 100
    frac_c12_c10_p_n_2piminus_piplus = float(number_c12_c10_p_n_2piminus_piplus) / float(number_entries) * 100
    frac_c12_c10_2n_piminus_piplus = float(number_c12_c10_2n_piminus_piplus) / float(number_entries) * 100
    frac_c12_c10_2p_2piminus = float(number_c12_c10_2p_2piminus) / float(number_entries) * 100
    frac_c12_c10_other = float(number_c12_c10_other) / float(number_entries) * 100

    # Be10:
    frac_c12_be10_2p = float(number_c12_be10_2p) / float(number_entries) * 100
    frac_c12_be10_p_n_piplus = float(number_c12_be10_p_n_piplus) / float(number_entries) * 100
    frac_c12_be10_p_n_piminus_2piplus = float(number_c12_be10_p_n_piminus_2piplus) / float(number_entries) * 100
    frac_c12_be10_2p_piminus_piplus = float(number_c12_be10_2p_piminus_piplus) / float(number_entries) * 100
    frac_c12_be10_2n_2piplus = float(number_c12_be10_2n_2piplus) / float(number_entries) * 100
    frac_c12_be10_p_n_2piminus_3piplus = float(number_c12_be10_p_n_2piminus_3piplus) / float(number_entries) * 100
    frac_c12_be10_2p_2piminus_2piplus = float(number_c12_be10_2p_2piminus_2piplus) / float(number_entries) * 100
    frac_c12_be10_2p_3piminus_3piplus = float(number_c12_be10_2p_3piminus_3piplus) / float(number_entries) * 100
    frac_c12_be10_other = float(number_c12_be10_other) / float(number_entries) * 100

    # B9:
    frac_c12_b9_p_2n = float(number_c12_b9_p_2n) / float(number_entries) * 100
    frac_c12_b9_p_2n_piminus_piplus = float(number_c12_b9_p_2n_piminus_piplus) / float(number_entries) * 100
    frac_c12_b9_2p_n_3piminus_2piplus = float(number_c12_b9_2p_n_3piminus_2piplus) / float(number_entries) * 100
    frac_c12_b9_2p_n_piminus = float(number_c12_b9_2p_n_piminus) / float(number_entries) * 100
    frac_c12_b9_3n_piplus = float(number_c12_b9_3n_piplus) / float(number_entries) * 100
    frac_c12_b9_p_2n_2piminus_2piplus = float(number_c12_b9_p_2n_2piminus_2piplus) / float(number_entries) * 100
    frac_c12_b9_2p_n_2piminus_piplus = float(number_c12_b9_2p_n_2piminus_piplus) / float(number_entries) * 100
    frac_c12_b9_other = float(number_c12_b9_other) / float(number_entries) * 100

    # Be9:
    frac_c12_be9_2p_n = float(number_c12_be9_2p_n) / float(number_entries) * 100
    frac_c12_be9_p_2n_piplus = float(number_c12_be9_p_2n_piplus) / float(number_entries) * 100
    frac_c12_be9_3p_piminus = float(number_c12_be9_3p_piminus) / float(number_entries) * 100
    frac_c12_be9_p_2n_piminus_2piplus = float(number_c12_be9_p_2n_piminus_2piplus) / float(number_entries) * 100
    frac_c12_be9_2p_n_piminus_piplus = float(number_c12_be9_2p_n_piminus_piplus) / float(number_entries) * 100
    frac_c12_be9_2p_n_3piminus_3piplus = float(number_c12_be9_2p_n_3piminus_3piplus) / float(number_entries) * 100
    frac_c12_be9_2p_n_2piminus_2piplus = float(number_c12_be9_2p_n_2piminus_2piplus) / float(number_entries) * 100
    frac_c12_be9_3n_2piplus = float(number_c12_be9_3n_2piplus) / float(number_entries) * 100
    frac_c12_be9_3p_2piminus_piplus = float(number_c12_be9_3p_2piminus_piplus) / float(number_entries) * 100
    frac_c12_be9_other = float(number_c12_be9_other) / float(number_entries) * 100

    # Be8:
    frac_c12_be8_2p_2n = float(number_c12_be8_2p_2n) / float(number_entries) * 100
    frac_c12_be8_3p_n_piminus = float(number_c12_be8_3p_n_piminus) / float(number_entries) * 100
    frac_c12_be8_p_3n_piplus = float(number_c12_be8_p_3n_piplus) / float(number_entries) * 100
    frac_c12_be8_2p_2n_2piminus_2piplus = float(number_c12_be8_2p_2n_2piminus_2piplus) / float(number_entries) * 100
    frac_c12_be8_4n_2piplus = float(number_c12_be8_4n_2piplus) / float(number_entries) * 100
    frac_c12_be8_2p_2n_piminus_piplus = float(number_c12_be8_2p_2n_piminus_piplus) / float(number_entries) * 100
    frac_c12_be8_3p_n_2piminus_piplus = float(number_c12_be8_3p_n_2piminus_piplus) / float(number_entries) * 100
    frac_c12_be8_4p_2piminus = float(number_c12_be8_4p_2piminus) / float(number_entries) * 100
    frac_c12_be8_other = float(number_c12_be8_other) / float(number_entries) * 100

    # C9:
    frac_c12_c9_p_2n_piminus = float(number_c12_c9_p_2n_piminus) / float(number_entries) * 100
    frac_c12_c9_3n = float(number_c12_c9_3n) / float(number_entries) * 100
    frac_c12_c9_2p_n_2piminus = float(number_c12_c9_2p_n_2piminus) / float(number_entries) * 100
    frac_c12_c9_3n_2piminus_2piplus = float(number_c12_c9_3n_2piminus_2piplus) / float(number_entries) * 100
    frac_c12_c9_other = float(number_c12_c9_other) / float(number_entries) * 100

    # Be7:
    frac_c12_be7_2p_3n = float(number_c12_be7_2p_3n) / float(number_entries) * 100
    frac_c12_be7_p_4n_piplus = float(number_c12_be7_p_4n_piplus) / float(number_entries) * 100
    frac_c12_be7_2p_3n_2piminus_2piplus = float(number_c12_be7_2p_3n_2piminus_2piplus) / float(number_entries) * 100
    frac_c12_be7_3p_2n_piminus = float(number_c12_be7_3p_2n_piminus) / float(number_entries) * 100
    frac_c12_be7_4p_n_2piminus = float(number_c12_be7_4p_n_2piminus) / float(number_entries) * 100
    frac_c12_be7_3p_2n_2piminus_piplus = float(number_c12_be7_3p_2n_2piminus_piplus) / float(number_entries) * 100
    frac_c12_be7_other = float(number_c12_be7_other) / float(number_entries) * 100

    # Li6:
    frac_c12_li6_3p_3n = float(number_c12_li6_3p_3n) / float(number_entries) * 100
    frac_c12_li6_2p_4n_piplus = float(number_c12_li6_2p_4n_piplus) / float(number_entries) * 100
    frac_c12_li6_5p_n_2piminus = float(number_c12_li6_5p_n_2piminus) / float(number_entries) * 100
    frac_c12_li6_2p_4n_piminus_2piplus = float(number_c12_li6_2p_4n_piminus_2piplus) / float(number_entries) * 100
    frac_c12_li6_4p_2n_piminus = float(number_c12_li6_4p_2n_piminus) / float(number_entries) * 100
    frac_c12_li6_3p_3n_piminus_piplus = float(number_c12_li6_3p_3n_piminus_piplus) / float(number_entries) * 100
    frac_c12_li6_other = float(number_c12_li6_other) / float(number_entries) * 100

    # Li8:
    frac_c12_li8_3p_n = float(number_c12_li8_3p_n) / float(number_entries) * 100
    frac_c12_li8_4p_piminus = float(number_c12_li8_4p_piminus) / float(number_entries) * 100
    frac_c12_li8_4p_2piminus_piplus = float(number_c12_li8_4p_2piminus_piplus) / float(number_entries) * 100
    frac_c12_li8_2p_2n_piplus = float(number_c12_li8_2p_2n_piplus) / float(number_entries) * 100
    frac_c12_li8_3p_n_piminus_piplus = float(number_c12_li8_3p_n_piminus_piplus) / float(number_entries) * 100
    frac_c12_li8_other = float(number_c12_li8_other) / float(number_entries) * 100

    # Li7:
    frac_c12_li7_2p_3n_piplus = float(number_c12_li7_2p_3n_piplus) / float(number_entries) * 100
    frac_c12_li7_4p_n_piminus = float(number_c12_li7_4p_n_piminus) / float(number_entries) * 100
    frac_c12_li7_3p_2n = float(number_c12_li7_3p_2n) / float(number_entries) * 100
    frac_c12_li7_3p_2n_piminus_piplus = float(number_c12_li7_3p_2n_piminus_piplus) / float(number_entries) * 100
    frac_c12_li7_4p_n_2piminus_piplus = float(number_c12_li7_4p_n_2piminus_piplus) / float(number_entries) * 100
    frac_c12_li7_2p_3n_piminus_2piplus = float(number_c12_li7_2p_3n_piminus_2piplus) / float(number_entries) * 100
    frac_c12_li7_other = float(number_c12_li7_other) / float(number_entries) * 100

    # B8:
    frac_c12_b8_p_3n = float(number_c12_b8_p_3n) / float(number_entries) * 100
    frac_c12_b8_p_3n_piminus_piplus = float(number_c12_b8_p_3n_piminus_piplus) / float(number_entries) * 100
    frac_c12_b8_2p_2n_2piminus_piplus = float(number_c12_b8_2p_2n_2piminus_piplus) / float(number_entries) * 100
    frac_c12_b8_2p_2n_piminus = float(number_c12_b8_2p_2n_piminus) / float(number_entries) * 100
    frac_c12_b8_4n_piplus = float(number_c12_b8_4n_piplus) / float(number_entries) * 100
    frac_c12_b8_other = float(number_c12_b8_other) / float(number_entries) * 100

    # Li9:
    frac_c12_li9_2p_n_piplus = float(number_c12_li9_2p_n_piplus) / float(number_entries) * 100
    frac_c12_li9_3p = float(number_c12_li9_3p) / float(number_entries) * 100
    frac_c12_li9_3p_piminus_piplus = float(number_c12_li9_3p_piminus_piplus) / float(number_entries) * 100
    frac_c12_li9_2p_n_piminus_2piplus = float(number_c12_li9_2p_n_piminus_2piplus) / float(number_entries) * 100
    frac_c12_li9_p_2n_piminus_3piplus = float(number_c12_li9_p_2n_piminus_3piplus) / float(number_entries) * 100
    frac_c12_li9_other = float(number_c12_li9_other) / float(number_entries) * 100

    # C8:
    frac_c12_c8_4n = float(number_c12_c8_4n) / float(number_entries) * 100
    frac_c12_c8_4n_other = float(number_c12_c8_4n_other) / float(number_entries) * 100

    # He8:
    frac_c12_he8_4p = float(number_c12_he8_4p) / float(number_entries) * 100
    frac_c12_he8_4p_other = float(number_c12_he8_4p_other) / float(number_entries) * 100

    # B7:
    frac_c12_b7_p_4n = float(number_c12_b7_p_4n) / float(number_entries) * 100
    frac_c12_b7_p_4n_other = float(number_c12_b7_p_4n_other) / float(number_entries) * 100

    # He7:
    frac_c12_he7_4p_n = float(number_c12_he7_4p_n) / float(number_entries) * 100
    frac_c12_he7_4p_n_other = float(number_c12_he7_4p_n_other) / float(number_entries) * 100

    # H7:
    frac_c12_h7_5p = float(number_c12_h7_5p) / float(number_entries) * 100
    frac_c12_h7_5p_other = float(number_c12_h7_5p_other) / float(number_entries) * 100

    # Be6:
    frac_c12_be6_2p_4n = float(number_c12_be6_2p_4n) / float(number_entries) * 100
    frac_c12_be6_2p_4n_other = float(number_c12_be6_2p_4n_other) / float(number_entries) * 100

    # He6:
    frac_c12_he6_4p_2n = float(number_c12_he6_4p_2n) / float(number_entries) * 100
    frac_c12_he6_4p_2n_other = float(number_c12_he6_4p_2n_other) / float(number_entries) * 100

    # H6:
    frac_c12_h6_5p_n = float(number_c12_h6_5p_n) / float(number_entries) * 100
    frac_c12_h6_5p_n_other = float(number_c12_h6_5p_n_other) / float(number_entries) * 100

    # missing channels:
    frac_c12_mass11u = float(number_c12_mass11u) / float(number_entries) * 100
    frac_c12_mass10u = float(number_c12_mass10u) / float(number_entries) * 100
    frac_c12_mass9u = float(number_c12_mass9u) / float(number_entries) * 100
    frac_c12_mass8u = float(number_c12_mass8u) / float(number_entries) * 100
    frac_c12_mass7u = float(number_c12_mass7u) / float(number_entries) * 100
    frac_c12_mass6u = float(number_c12_mass6u) / float(number_entries) * 100
    frac_c12_mass5orless = float(number_c12_mass5orless) / float(number_entries) * 100

    # C12:
    frac_c12_c12 = float(number_c12_c12) / float(number_entries) * 100

    # no isotope (only protons, neutrons, pions):
    frac_c12_noiso = float(number_c12_noiso) / float(number_entries) * 100
    frac_c12_noiso_5p_6n = float(number_c12_noiso_5p_6n) / float(number_entries) * 100

    # faulty interaction channels: isotopes are missing (no isotopes with Z<3 and (N-Z)<3 are considered):
    frac_c12_faulty = float(number_c12_faulty) / float(number_entries) * 100

    # Other targets than C12:
    frac_no_c12 = float(number_no_c12) / float(number_entries) * 100
    frac_es_p = float(number_es_p) / float(number_entries) * 100
    frac_es_e = float(number_es_e) / float(number_entries) * 100
    frac_es_o16 = float(number_es_o16) / float(number_entries) * 100
    frac_es_n14 = float(number_es_n14) / float(number_entries) * 100
    frac_es_s32 = float(number_es_s32) / float(number_entries) * 100


    return (frac_c12_b11_p, frac_c12_b11_n_piplus, frac_c12_b11_n_piminus_2piplus, frac_c12_b11_p_piminus_piplus,
            frac_c12_b11_p_2piminus_2piplus, frac_c12_b11_piplus, frac_c12_b11_other,
            frac_c12_c11_n, frac_c12_c11_p_piminus, frac_c12_c11_n_piminus_piplus, frac_c12_c11_p_2piminus_piplus,
            frac_c12_c11_p_3piminus_2piplus, frac_c12_c11_n_2piminus_2piplus, frac_c12_c11_other,
            frac_c12_b10_p_n, frac_c12_b10_2p_piminus, frac_c12_b10_p_n_piminus_piplus, frac_c12_b10_2n_piplus,
            frac_c12_b10_2n_piminus_2piplus, frac_c12_b10_2p_2piminus_piplus, frac_c12_b10_2p_3piminus_2piplus,
            frac_c12_b10_p_n_2piminus_2piplus, frac_c12_b10_other,
            frac_c12_c10_2n, frac_c12_c10_p_n_piminus, frac_c12_c10_p_n_2piminus_piplus, frac_c12_c10_2n_piminus_piplus,
            frac_c12_c10_2p_2piminus, frac_c12_c10_other,
            frac_c12_be10_2p, frac_c12_be10_p_n_piplus, frac_c12_be10_p_n_piminus_2piplus,
            frac_c12_be10_2p_piminus_piplus, frac_c12_be10_2n_2piplus, frac_c12_be10_p_n_2piminus_3piplus,
            frac_c12_be10_2p_2piminus_2piplus, frac_c12_be10_2p_3piminus_3piplus, frac_c12_be10_other,
            frac_c12_b9_p_2n, frac_c12_b9_p_2n_piminus_piplus, frac_c12_b9_2p_n_3piminus_2piplus,
            frac_c12_b9_2p_n_piminus, frac_c12_b9_3n_piplus, frac_c12_b9_p_2n_2piminus_2piplus,
            frac_c12_b9_2p_n_2piminus_piplus, frac_c12_b9_other,
            frac_c12_be9_2p_n, frac_c12_be9_p_2n_piplus, frac_c12_be9_3p_piminus, frac_c12_be9_p_2n_piminus_2piplus,
            frac_c12_be9_2p_n_piminus_piplus, frac_c12_be9_2p_n_3piminus_3piplus, frac_c12_be9_2p_n_2piminus_2piplus,
            frac_c12_be9_3n_2piplus, frac_c12_be9_3p_2piminus_piplus, frac_c12_be9_other,
            frac_c12_be8_2p_2n, frac_c12_be8_3p_n_piminus, frac_c12_be8_p_3n_piplus,
            frac_c12_be8_2p_2n_2piminus_2piplus, frac_c12_be8_4n_2piplus, frac_c12_be8_2p_2n_piminus_piplus,
            frac_c12_be8_3p_n_2piminus_piplus, frac_c12_be8_4p_2piminus, frac_c12_be8_other,
            frac_c12_c9_p_2n_piminus, frac_c12_c9_3n, frac_c12_c9_2p_n_2piminus, frac_c12_c9_3n_2piminus_2piplus,
            frac_c12_c9_other,
            frac_c12_be7_2p_3n, frac_c12_be7_p_4n_piplus, frac_c12_be7_2p_3n_2piminus_2piplus,
            frac_c12_be7_3p_2n_piminus, frac_c12_be7_4p_n_2piminus, frac_c12_be7_3p_2n_2piminus_piplus,
            frac_c12_be7_other,
            frac_c12_li6_3p_3n, frac_c12_li6_2p_4n_piplus, frac_c12_li6_5p_n_2piminus,
            frac_c12_li6_2p_4n_piminus_2piplus, frac_c12_li6_4p_2n_piminus, frac_c12_li6_3p_3n_piminus_piplus,
            frac_c12_li6_other,
            frac_c12_li8_3p_n, frac_c12_li8_4p_piminus, frac_c12_li8_4p_2piminus_piplus, frac_c12_li8_2p_2n_piplus,
            frac_c12_li8_3p_n_piminus_piplus, frac_c12_li8_other,
            frac_c12_li7_2p_3n_piplus, frac_c12_li7_4p_n_piminus, frac_c12_li7_3p_2n, frac_c12_li7_3p_2n_piminus_piplus,
            frac_c12_li7_4p_n_2piminus_piplus, frac_c12_li7_2p_3n_piminus_2piplus, frac_c12_li7_other,
            frac_c12_b8_p_3n, frac_c12_b8_p_3n_piminus_piplus, frac_c12_b8_2p_2n_2piminus_piplus,
            frac_c12_b8_2p_2n_piminus, frac_c12_b8_4n_piplus, frac_c12_b8_other,
            frac_c12_li9_2p_n_piplus, frac_c12_li9_3p, frac_c12_li9_3p_piminus_piplus,
            frac_c12_li9_2p_n_piminus_2piplus, frac_c12_li9_p_2n_piminus_3piplus, frac_c12_li9_other,
            frac_c12_c8_4n, frac_c12_c8_4n_other, frac_c12_he8_4p, frac_c12_he8_4p_other,
            frac_c12_b7_p_4n, frac_c12_b7_p_4n_other, frac_c12_he7_4p_n, frac_c12_he7_4p_n_other,
            frac_c12_h7_5p, frac_c12_h7_5p_other, frac_c12_be6_2p_4n, frac_c12_be6_2p_4n_other,
            frac_c12_he6_4p_2n, frac_c12_he6_4p_2n_other, frac_c12_h6_5p_n, frac_c12_h6_5p_n_other,
            frac_c12_mass11u, frac_c12_mass10u, frac_c12_mass9u, frac_c12_mass8u, frac_c12_mass7u, frac_c12_mass6u,
            frac_c12_mass5orless,
            frac_c12_c12, frac_c12_noiso, frac_c12_noiso_5p_6n,
            frac_no_c12, frac_es_p, frac_es_e, frac_es_o16, frac_es_n14, frac_es_s32, frac_c12_faulty)


def get_deex_channel(deex_id, isotope_pdg, target_pdg):
    """
    function to calculate the different deexcitation channels of the different isotopes, which were produced by
    neutral current interaction of atmospheric neutrinos (deex_id, isotope_pdg, target_pdg from root file from output
    of DSNB-NC generator)

    :param deex_id: ID of the deexcitation channel: defines, which particles are produced in the deexitation (array)
    :param isotope_pdg: PDG ID of the isotope produced through NC interaction (array)
    :param target_pdg: PDG ID of the target particle (array)

    :return:
    """

    # get the number of entries of the array (integer):
    number_entries = len(target_pdg)
    # number_entries = 1000

    """ preallocate the variables: """
    # preallocate number of events with C12 as target:
    number_target_c12 = 0
    # number of NC interaction without C12 as target:
    number_no_c12 = 0
    # number of NC interaction with 'light' isotopes (isotopes, where no TALYS deexcitation root file exists),
    # where deex_id = 0:
    number_light_iso = 0

    """ C11 """
    # number of events, where C11 is not excited:
    number_c11_notex = 0
    # number of events, where C11 de-excites:
    number_c11_deex = 0
    # C11* -> p + alpha + Li6:
    number_c11_li6_p_alpha = 0
    # C11* -> alpha + Be7:
    number_c11_be7_alpha = 0
    # C11* -> p + B10:
    number_c11_b10_p = 0
    # C11* -> n + p + B9:
    number_c11_b9_n_p = 0
    # C11* -> p + d + Be8:
    number_c11_be8_p_d = 0
    # C11* -> 2p + Be9:
    number_c11_be9_2p = 0
    # C11* -> d + B9:
    number_c11_b9_d = 0
    # C11* -> He3 + Be8:
    number_c11_be8_he3 = 0
    # C11* -> n + C10:
    number_c11_c10_n = 0
    # C11* -> d + alpha + Li5:
    number_c11_li5_d_alpha = 0
    # C11* -> n + p + alpha + Li5:
    number_c11_li5_n_p_alpha = 0
    # deexcitations of C11 not yet included:
    number_c11_missing = 0

    """ B11 """
    # number of events, where B11 is not excited:
    number_b11_notex = 0
    # number of events, where B11 de-excites:
    number_b11_deex = 0
    # B11* -> n + alpha + Li6:
    number_b11_li6_n_alpha = 0
    # B11* -> 2n + B9:
    number_b11_b9_2n = 0
    # B11* -> n + d + Be8:
    number_b11_be8_n_d = 0
    # B11* -> d + Be9:
    number_b11_be9_d = 0
    # B11* -> p + Be10:
    number_b11_be10_p = 0
    # B11* -> n + B10:
    number_b11_b10_n = 0
    # B11* -> n + p + Be9:
    number_b11_be9_n_p = 0
    # B11* -> alpha + Li7:
    number_b11_li7_alpha = 0
    # B11* -> t + Be8:
    number_b11_be8_t = 0
    # B11* -> d + alpha + He5:
    number_b11_he5_d_alpha = 0
    # B11* -> p + alpha + He6:
    number_b11_he6_p_alpha = 0
    # B11* -> 2n + p + Be8:
    number_b11_be8_2n_p = 0
    # deexcitations of B11 not yet included:
    number_b11_missing = 0

    """ C10 """
    # number of events, where C10 is not excited:
    number_c10_notex = 0
    # number of events, where C10 de-excites:
    number_c10_deex = 0
    # C10* -> p + B9:
    number_c10_b9_p = 0
    # C10* -> p + d + Be7:
    number_c10_be7_p_d = 0
    # C10* -> p + He3 + Li6:
    number_c10_li6_p_he3 = 0
    # C10* -> p + d + He3 + He4:
    number_c10_he4_p_d_he3 = 0
    # C10* -> 2p + d + Li6:
    number_c10_li6_2p_d = 0
    # C10* -> 2p + Be8:
    number_c10_be8_2p = 0
    # C10* -> n + 2p + Be7:
    number_c10_be7_n_2p = 0
    # C10* -> n + 3p + Li6:
    number_c10_li6_n_3p = 0
    # C10* -> n + p + d + Be6:
    number_c10_be6_n_p_d = 0
    # C10* -> n + 2p + d + Li5:
    number_c10_li5_n_2p_d = 0
    # C10* -> p + d + alpha + He3:
    number_c10_he3_p_d_alpha = 0
    # C10* -> d + He3 + Li5:
    number_c10_li5_d_he3 = 0
    # C10* -> p + 2d + Li5:
    number_c10_li5_p_2d = 0
    # C10* -> n + 2p + alpha + He3:
    number_c10_he3_n_2p_alpha = 0
    # C10* -> n + p + alpha + Li4:
    number_c10_li4_n_p_alpha = 0
    # C10* -> n + p + B8:
    number_c10_b8_n_p = 0
    # C10* -> d + B8:
    number_c10_b8_d = 0
    # C10* -> p + t + Be6
    number_c10_be6_p_t = 0
    # C10* -> n + 2p + He3 + He4
    number_c10_he4_n_2p_he3 = 0
    # C10* -> n + p + He3 + Li5
    number_c10_li5_n_p_he3 = 0
    # deexcitations of C10 not yet included:
    number_c10_missing = 0

    """ B10 """
    # number of events, where B10 is not excited:
    number_b10_notex = 0
    # number of events, where B10 de-excites:
    number_b10_deex = 0
    # B10* -> p + Be9:
    number_b10_be9_p = 0
    # B10* -> d + Be8:
    number_b10_be8_d = 0
    # B10* -> n + B9:
    number_b10_b9_n = 0
    # B10* -> t + Be7:
    number_b10_be7_t = 0
    # B10* -> n + p + Be8:
    number_b10_be8_n_p = 0
    # B10* -> He3 + Li7:
    number_b10_li7_he3 = 0
    # B10* -> p + alpha + He5:
    number_b10_he5_p_alpha = 0
    # B10* -> alpha + Li6:
    number_b10_li6_alpha = 0
    # B10* -> n + alpha + Li5:
    number_b10_li5_n_alpha = 0
    # B10* -> p + d + Li7:
    number_b10_li7_p_d = 0
    # deexcitations of B10 not yet included:
    number_b10_missing = 0

    """ Be10 """
    # number of events, where Be10 is not excited:
    number_be10_notex = 0
    # number of events, where Be10 de-excites:
    number_be10_deex = 0
    # Be10* -> 2n + d + Li6:
    number_be10_li6_2n_d = 0
    # Be10* -> 2n + p + Li7:
    number_be10_li7_2n_p = 0
    # Be10* -> n + 2p + He7:
    number_be10_he7_n_2p = 0
    # Be10* -> n + t + Li6:
    number_be10_li6_n_t = 0
    # Be10* -> 3n + p + Li6:
    number_be10_li6_3n_p = 0
    # Be10* -> n + 2d + He5:
    number_be10_he5_n_2d = 0
    # Be10* -> 2n + Be8:
    number_be10_be8_2n = 0
    # Be10* -> n + alpha + He5:
    number_be10_he5_n_alpha = 0
    # Be10* -> n + d + alpha + tritium:
    number_be10_t_n_d_alpha = 0
    # Be10* -> n + p + Li8:
    number_be10_li8_n_p = 0
    # Be10* -> n + d + t + He4:
    number_be10_he4_n_d_t = 0
    # Be10* -> n + p + t + He5:
    number_be10_he5_n_p_t = 0
    # Be10* -> n + d + Li7:
    number_be10_li7_n_d = 0
    # Be10* -> n + p + alpha + H4:
    number_be10_h4_n_p_alpha = 0
    # Be10* -> 2n + p + d + He5:
    number_be10_he5_2n_p_d = 0
    # Be10* -> n + Be9:
    number_be10_be9_n = 0
    # Be10* -> n + p + d + He6:
    number_be10_he6_n_p_d = 0
    # Be10* -> d + t + He5:
    number_be10_he5_d_t = 0
    # Be10* -> d + alpha + H4:
    number_be10_h4_d_alpha = 0
    # Be10* -> 2n + p + t + He4:
    number_be10_he4_2n_p_t = 0
    # deexcitations of Be10 not yet included:
    number_be10_missing = 0

    """ C9 """
    # number of events, where C9 is not excited:
    number_c9_notex = 0
    # number of events, where C9 de-excites:
    number_c9_deex = 0
    # C9* -> 2p + Be7:
    number_c9_be7_2p = 0
    # deexcitations of C9 not yet included:
    number_c9_missing = 0

    """ B9 """
    # number of events, where B9 is not excited:
    number_b9_notex = 0
    # number of events, where B9 de-excites:
    number_b9_deex = 0
    # B9* -> p + Be8:
    number_b9_be8_p = 0
    # deexcitations of B9 not yet included:
    number_b9_missing = 0

    """ Be9 """
    # number of events, where Be9 is not excited:
    number_be9_notex = 0
    # number of events, where Be9 de-excites:
    number_be9_deex = 0
    # Be9* -> p + Li8:
    number_be9_li8_p = 0
    # deexcitations of Be9 not yet included:
    number_be9_missing = 0

    """ Li9 """
    # number of events, where Li9 is not excited:
    number_li9_notex = 0
    # number of events, where Li9 de-excites:
    number_li9_deex = 0
    # Li9* -> n + alpha + H4:
    number_li9_h4_n_alpha = 0
    # Li9* -> d + He7:
    number_li9_he7_d = 0
    # Li9* -> n + Li8:
    number_li9_li8_n = 0
    # deexcitations of Li9 not yet included:
    number_li9_missing = 0

    """ B8 """
    # number of events, where B8 is not excited:
    number_b8_notex = 0
    # number of events, where B8 de-excites:
    number_b8_deex = 0
    # B8* -> 2p + Li6:
    number_b8_li6_2p = 0
    # deexcitations of B8 not yet included:
    number_b8_missing = 0

    """ Li8 """
    # number of events, where Li8 is not excited:
    number_li8_notex = 0
    # number of events, where Li8 de-excites:
    number_li8_deex = 0
    # Li8* -> n + Li7:
    number_li8_li7_n = 0
    # Li8* -> 2n + Li6:
    number_li8_li6_2n = 0
    # deexcitations of Li8 not yet included:
    number_li8_missing = 0

    """ Be7 """
    # number of events, where Be7 is not excited:
    number_be7_notex = 0
    # number of events, where Be7 de-excites:
    number_be7_deex = 0
    # Be7* -> d + Li5:
    number_be7_li5_d = 0
    # Be7* -> p + Li6:
    number_be7_li6_p = 0
    # deexcitations of Be7 not yet included:
    number_be7_missing = 0

    """ Li7 """
    # number of events, where Li7 is not excited:
    number_li7_notex = 0
    # number of events, where Li7 de-excites:
    number_li7_deex = 0
    # Li7* -> n + Li6:
    number_li7_li6_n = 0
    # deexcitations of Li7 not yet included:
    number_li7_missing = 0


    # loop over all entries of the array:
    for index in range(number_entries):

        # check, if target is C12 (PDG ID = 1000060120):
        if target_pdg[index] == 1000060120:

            number_target_c12 = number_target_c12 + 1

            # check the PDG ID of the created isotopes:
            if isotope_pdg[index] == 1000060110:
                # C11:
                if deex_id[index] == 0:
                    # C11 is not excited:
                    number_c11_notex = number_c11_notex + 1

                else:
                    # C11 deexcitation:
                    number_c11_deex = number_c11_deex + 1
                    # get number of particles:
                    num_n, num_p, num_d, num_t, num_he3, num_alpha = get_number_of_particles_of_deexid(deex_id[index])

                    if num_n == 0 and num_p == 1 and num_d == 0 and num_t == 0 and num_he3 == 0 and num_alpha == 1:
                        # C11* -> p + alpha + Li6:
                        number_c11_li6_p_alpha = number_c11_li6_p_alpha + 1
                    elif num_n == 0 and num_p == 0 and num_d == 0 and num_t == 0 and num_he3 == 0 and num_alpha == 1:
                        # C11* -> alpha + Be7:
                        number_c11_be7_alpha = number_c11_be7_alpha + 1
                    elif num_n == 0 and num_p == 1 and num_d == 0 and num_t == 0 and num_he3 == 0 and num_alpha == 0:
                        # C11* -> p + B10:
                        number_c11_b10_p = number_c11_b10_p + 1
                    elif num_n == 1 and num_p == 1 and num_d == 0 and num_t == 0 and num_he3 == 0 and num_alpha == 0:
                        # C11* -> n + p + B9:
                        number_c11_b9_n_p = number_c11_b9_n_p + 1
                    elif num_n == 0 and num_p == 1 and num_d == 1 and num_t == 0 and num_he3 == 0 and num_alpha == 0:
                        # C11* -> p + d + Be8:
                        number_c11_be8_p_d = number_c11_be8_p_d + 1
                    elif num_n == 0 and num_p == 2 and num_d == 0 and num_t == 0 and num_he3 == 0 and num_alpha == 0:
                        # C11* -> 2p + Be9:
                        number_c11_be9_2p = number_c11_be9_2p + 1
                    elif num_n == 0 and num_p == 0 and num_d == 1 and num_t == 0 and num_he3 == 0 and num_alpha == 0:
                        # C11* -> d + B9:
                        number_c11_b9_d = number_c11_b9_d + 1
                    elif num_n == 0 and num_p == 0 and num_d == 0 and num_t == 0 and num_he3 == 1 and num_alpha == 0:
                        # C11* -> He3 + Be8:
                        number_c11_be8_he3 = number_c11_be8_he3 + 1
                    elif num_n == 1 and num_p == 0 and num_d == 0 and num_t == 0 and num_he3 == 0 and num_alpha == 0:
                        # C11* -> n + C10:
                        number_c11_c10_n = number_c11_c10_n + 1
                    elif num_n == 0 and num_p == 0 and num_d == 1 and num_t == 0 and num_he3 == 0 and num_alpha == 1:
                        # C11* -> d + alpha + Li5:
                        number_c11_li5_d_alpha = number_c11_li5_d_alpha + 1
                    elif num_n == 1 and num_p == 1 and num_d == 0 and num_t == 0 and num_he3 == 0 and num_alpha == 1:
                        # C11* -> n + p + alpha + Li5:
                        number_c11_li5_n_p_alpha = number_c11_li5_n_p_alpha + 1
                    else:
                        number_c11_missing = number_c11_missing + 1
                        print("----------C11-------")
                        print(deex_id[index])


            elif isotope_pdg[index] == 1000050110:
                # B11:
                if deex_id[index] == 0:
                    # B11 is not excited:
                    number_b11_notex = number_b11_notex + 1

                else:
                    # B11 deexcitation:
                    number_b11_deex = number_b11_deex + 1
                    # get number of particles:
                    num_n, num_p, num_d, num_t, num_he3, num_alpha = get_number_of_particles_of_deexid(deex_id[index])

                    if num_n == 1 and num_p == 0 and num_d == 0 and num_t == 0 and num_he3 == 0 and num_alpha == 1:
                        # B11* -> n + alpha + Li6:
                        number_b11_li6_n_alpha = number_b11_li6_n_alpha + 1
                    elif num_n == 2 and num_p == 0 and num_d == 0 and num_t == 0 and num_he3 == 0 and num_alpha == 0:
                        # B11* -> 2n + B9:
                        number_b11_b9_2n = number_b11_b9_2n + 1
                    elif num_n == 1 and num_p == 0 and num_d == 1 and num_t == 0 and num_he3 == 0 and num_alpha == 0:
                        # B11* -> n + d + Be8:
                        number_b11_be8_n_d = number_b11_be8_n_d + 1
                    elif num_n == 0 and num_p == 0 and num_d == 1 and num_t == 0 and num_he3 == 0 and num_alpha == 0:
                        # B11* -> d + Be9:
                        number_b11_be9_d = number_b11_be9_d + 1
                    elif num_n == 0 and num_p == 1 and num_d == 0 and num_t == 0 and num_he3 == 0 and num_alpha == 0:
                        # B11* -> p + Be10:
                        number_b11_be10_p = number_b11_be10_p + 1
                    elif num_n == 1 and num_p == 0 and num_d == 0 and num_t == 0 and num_he3 == 0 and num_alpha == 0:
                        # B11* -> n + B10:
                        number_b11_b10_n = number_b11_b10_n + 1
                    elif num_n == 1 and num_p == 1 and num_d == 0 and num_t == 0 and num_he3 == 0 and num_alpha == 0:
                        # B11* -> n + p + Be9:
                        number_b11_be9_n_p = number_b11_be9_n_p + 1
                    elif num_n == 0 and num_p == 0 and num_d == 0 and num_t == 0 and num_he3 == 0 and num_alpha == 1:
                        # B11* -> alpha + Li7:
                        number_b11_li7_alpha = number_b11_li7_alpha + 1
                    elif num_n == 0 and num_p == 0 and num_d == 0 and num_t == 1 and num_he3 == 0 and num_alpha == 0:
                        # B11* -> t + Be8:
                        number_b11_be8_t = number_b11_be8_t + 1
                    elif num_n == 0 and num_p == 0 and num_d == 1 and num_t == 0 and num_he3 == 0 and num_alpha == 1:
                        # B11* -> d + alpha + He5:
                        number_b11_he5_d_alpha = number_b11_he5_d_alpha + 1
                    elif num_n == 0 and num_p == 1 and num_d == 0 and num_t == 0 and num_he3 == 0 and num_alpha == 1:
                        # B11* -> p + alpha + He6:
                        number_b11_he6_p_alpha = number_b11_he6_p_alpha + 1
                    elif num_n == 2 and num_p == 1 and num_d == 0 and num_t == 0 and num_he3 == 0 and num_alpha == 0:
                        # B11* -> 2n + p + Be8:
                        number_b11_be8_2n_p = number_b11_be8_2n_p + 1
                    else:
                        number_b11_missing = number_b11_missing + 1
                        print("----------B11-------")
                        print(deex_id[index])


            elif isotope_pdg[index] == 1000060100:
                # C10:
                if deex_id[index] == 0:
                    # C10 is not excited:
                    number_c10_notex = number_c10_notex + 1

                else:
                    # C10 deexcitation:
                    number_c10_deex = number_c10_deex + 1
                    # get number of particles:
                    num_n, num_p, num_d, num_t, num_he3, num_alpha = get_number_of_particles_of_deexid(deex_id[index])

                    if num_n == 0 and num_p == 1 and num_d == 0 and num_t == 0 and num_he3 == 0 and num_alpha == 0:
                        # C10* -> p + B9:
                        number_c10_b9_p = number_c10_b9_p + 1
                    elif num_n == 0 and num_p == 1 and num_d == 1 and num_t == 0 and num_he3 == 0 and num_alpha == 0:
                        # C10* -> p + d + Be7:
                        number_c10_be7_p_d = number_c10_be7_p_d + 1
                    elif num_n == 0 and num_p == 1 and num_d == 0 and num_t == 0 and num_he3 == 1 and num_alpha == 0:
                        # C10* -> p + He3 + Li6:
                        number_c10_li6_p_he3 = number_c10_li6_p_he3 + 1
                    elif num_n == 0 and num_p == 1 and num_d == 1 and num_t == 0 and num_he3 == 1 and num_alpha == 0:
                        # C10* -> p + d + He3 + He4:
                        number_c10_he4_p_d_he3 = number_c10_he4_p_d_he3 + 1
                    elif num_n == 0 and num_p == 2 and num_d == 1 and num_t == 0 and num_he3 == 0 and num_alpha == 0:
                        # C10* -> 2p + d + Li6:
                        number_c10_li6_2p_d = number_c10_li6_2p_d + 1
                    elif num_n == 0 and num_p == 2 and num_d == 0 and num_t == 0 and num_he3 == 0 and num_alpha == 0:
                        # C10* -> 2p + Be8:
                        number_c10_be8_2p = number_c10_be8_2p + 1
                    elif num_n == 1 and num_p == 2 and num_d == 0 and num_t == 0 and num_he3 == 0 and num_alpha == 0:
                        # C10* -> n + 2p + Be7:
                        number_c10_be7_n_2p = number_c10_be7_n_2p + 1
                    elif num_n == 1 and num_p == 3 and num_d == 0 and num_t == 0 and num_he3 == 0 and num_alpha == 0:
                        # C10* -> n + 3p + Li6:
                        number_c10_li6_n_3p = number_c10_li6_n_3p + 1
                    elif num_n == 1 and num_p == 1 and num_d == 1 and num_t == 0 and num_he3 == 0 and num_alpha == 0:
                        # C10* -> n + p + d + Be6:
                        number_c10_be6_n_p_d = number_c10_be6_n_p_d + 1
                    elif num_n == 1 and num_p == 2 and num_d == 1 and num_t == 0 and num_he3 == 0 and num_alpha == 0:
                        # C10* -> n + 2p + d + Li5:
                        number_c10_li5_n_2p_d = number_c10_li5_n_2p_d + 1
                    elif num_n == 0 and num_p == 1 and num_d == 1 and num_t == 0 and num_he3 == 0 and num_alpha == 1:
                        # C10* -> p + d + alpha + He3:
                        number_c10_he3_p_d_alpha = number_c10_he3_p_d_alpha + 1
                    elif num_n == 0 and num_p == 0 and num_d == 1 and num_t == 0 and num_he3 == 1 and num_alpha == 0:
                        # C10* -> d + He3 + Li5:
                        number_c10_li5_d_he3 = number_c10_li5_d_he3 + 1
                    elif num_n == 0 and num_p == 1 and num_d == 2 and num_t == 0 and num_he3 == 0 and num_alpha == 0:
                        # C10* -> p + 2d + Li5:
                        number_c10_li5_p_2d = number_c10_li5_p_2d + 1
                    elif num_n == 1 and num_p == 2 and num_d == 0 and num_t == 0 and num_he3 == 0 and num_alpha == 1:
                        # C10* -> n + 2p + alpha + He3:
                        number_c10_he3_n_2p_alpha = number_c10_he3_n_2p_alpha + 1
                    elif num_n == 1 and num_p == 1 and num_d == 0 and num_t == 0 and num_he3 == 0 and num_alpha == 1:
                        # C10* -> n + p + alpha + Li4:
                        number_c10_li4_n_p_alpha = number_c10_li4_n_p_alpha + 1
                    elif num_n == 1 and num_p == 1 and num_d == 0 and num_t == 0 and num_he3 == 0 and num_alpha == 0:
                        # C10* -> n + p + B8:
                        number_c10_b8_n_p = number_c10_b8_n_p + 1
                    elif num_n == 0 and num_p == 0 and num_d == 1 and num_t == 0 and num_he3 == 0 and num_alpha == 0:
                        # C10* -> d + B8:
                        number_c10_b8_d = number_c10_b8_d + 1
                    elif num_n == 0 and num_p == 1 and num_d == 0 and num_t == 1 and num_he3 == 0 and num_alpha == 0:
                        # C10* -> p + t + Be6
                        number_c10_be6_p_t = number_c10_be6_p_t + 1
                    elif num_n == 1 and num_p == 2 and num_d == 0 and num_t == 0 and num_he3 == 1 and num_alpha == 0:
                        # C10* -> n + 2p + He3 + He4
                        number_c10_he4_n_2p_he3 = number_c10_he4_n_2p_he3 + 1
                    elif num_n == 1 and num_p == 1 and num_d == 0 and num_t == 0 and num_he3 == 1 and num_alpha == 0:
                        # C10* -> n + p + He3 + Li5
                        number_c10_li5_n_p_he3 = number_c10_li5_n_p_he3 + 1
                    else:
                        number_c10_missing = number_c10_missing + 1
                        print("----------C10-------")
                        print(deex_id[index])


            elif isotope_pdg[index] == 1000050100:
                # B10:
                if deex_id[index] == 0:
                    # B10 is not excited:
                    number_b10_notex = number_b10_notex + 1

                else:
                    # B10 deexcitation:
                    number_b10_deex = number_b10_deex + 1
                    # get number of particles:
                    num_n, num_p, num_d, num_t, num_he3, num_alpha = get_number_of_particles_of_deexid(deex_id[index])

                    if num_n == 0 and num_p == 1 and num_d == 0 and num_t == 0 and num_he3 == 0 and num_alpha == 0:
                        # B10* -> p + Be9:
                        number_b10_be9_p = number_b10_be9_p + 1
                    elif num_n == 0 and num_p == 0 and num_d == 1 and num_t == 0 and num_he3 == 0 and num_alpha == 0:
                        # B10* -> d + Be8:
                        number_b10_be8_d = number_b10_be8_d + 1
                    elif num_n == 1 and num_p == 0 and num_d == 0 and num_t == 0 and num_he3 == 0 and num_alpha == 0:
                        # B10* -> n + B9:
                        number_b10_b9_n = number_b10_b9_n + 1
                    elif num_n == 0 and num_p == 0 and num_d == 0 and num_t == 1 and num_he3 == 0 and num_alpha == 0:
                        # B10* -> t + Be7:
                        number_b10_be7_t = number_b10_be7_t + 1
                    elif num_n == 1 and num_p == 1 and num_d == 0 and num_t == 0 and num_he3 == 0 and num_alpha == 0:
                        # B10* -> n + p + Be8:
                        number_b10_be8_n_p = number_b10_be8_n_p + 1
                    elif num_n == 0 and num_p == 0 and num_d == 0 and num_t == 0 and num_he3 == 1 and num_alpha == 0:
                        # B10* -> He3 + Li7:
                        number_b10_li7_he3 = number_b10_li7_he3 + 1
                    elif num_n == 0 and num_p == 1 and num_d == 0 and num_t == 0 and num_he3 == 0 and num_alpha == 1:
                        # B10* -> p + alpha + He5:
                        number_b10_he5_p_alpha = number_b10_he5_p_alpha + 1
                    elif num_n == 0 and num_p == 0 and num_d == 0 and num_t == 0 and num_he3 == 0 and num_alpha == 1:
                        # B10* -> alpha + Li6:
                        number_b10_li6_alpha = number_b10_li6_alpha + 1
                    elif num_n == 1 and num_p == 0 and num_d == 0 and num_t == 0 and num_he3 == 0 and num_alpha == 1:
                        # B10* -> n + alpha + Li5:
                        number_b10_li5_n_alpha = number_b10_li5_n_alpha + 1
                    elif num_n == 0 and num_p == 1 and num_d == 1 and num_t == 0 and num_he3 == 0 and num_alpha == 0:
                        # B10* -> p + d + Li7:
                        number_b10_li7_p_d = number_b10_li7_p_d + 1
                    else:
                        number_b10_missing = number_b10_missing + 1
                        print("----------B10-------")
                        print(deex_id[index])


            elif isotope_pdg[index] == 1000040100:
                # Be10:
                if deex_id[index] == 0:
                    # Be10 is not excited:
                    number_be10_notex = number_be10_notex + 1

                else:
                    # Be10 deexcitation:
                    number_be10_deex = number_be10_deex + 1
                    # get number of particles:
                    num_n, num_p, num_d, num_t, num_he3, num_alpha = get_number_of_particles_of_deexid(deex_id[index])

                    if num_n == 2 and num_p == 0 and num_d == 1 and num_t == 0 and num_he3 == 0 and num_alpha == 0:
                        # Be10* -> 2n + d + Li6:
                        number_be10_li6_2n_d = number_be10_li6_2n_d + 1
                    elif num_n == 2 and num_p == 1 and num_d == 0 and num_t == 0 and num_he3 == 0 and num_alpha == 0:
                        # Be10* -> 2n + p + Li7:
                        number_be10_li7_2n_p = number_be10_li7_2n_p + 1
                    elif num_n == 1 and num_p == 2 and num_d == 0 and num_t == 0 and num_he3 == 0 and num_alpha == 0:
                        # Be10* -> n + 2p + He7:
                        number_be10_he7_n_2p = number_be10_he7_n_2p + 1
                    elif num_n == 1 and num_p == 0 and num_d == 0 and num_t == 1 and num_he3 == 0 and num_alpha == 0:
                        # Be10* -> n + t + Li6:
                        number_be10_li6_n_t = number_be10_li6_n_t + 1
                    elif num_n == 3 and num_p == 1 and num_d == 0 and num_t == 0 and num_he3 == 0 and num_alpha == 0:
                        # Be10* -> 3n + p + Li6:
                        number_be10_li6_3n_p = number_be10_li6_3n_p + 1
                    elif num_n == 1 and num_p == 0 and num_d == 2 and num_t == 0 and num_he3 == 0 and num_alpha == 0:
                        # Be10* -> n + 2d + He5:
                        number_be10_he5_n_2d = number_be10_he5_n_2d + 1
                    elif num_n == 2 and num_p == 0 and num_d == 0 and num_t == 0 and num_he3 == 0 and num_alpha == 0:
                        # Be10* -> 2n + Be8:
                        number_be10_be8_2n = number_be10_be8_2n + 1
                    elif num_n == 1 and num_p == 0 and num_d == 0 and num_t == 0 and num_he3 == 0 and num_alpha == 1:
                        # Be10* -> n + alpha + He5:
                        number_be10_he5_n_alpha = number_be10_he5_n_alpha + 1
                    elif num_n == 1 and num_p == 0 and num_d == 1 and num_t == 0 and num_he3 == 0 and num_alpha == 1:
                        # Be10* -> n + d + alpha + tritium:
                        number_be10_t_n_d_alpha = number_be10_t_n_d_alpha + 1
                    elif num_n == 1 and num_p == 1 and num_d == 0 and num_t == 0 and num_he3 == 0 and num_alpha == 0:
                        # Be10* -> n + p + Li8:
                        number_be10_li8_n_p = number_be10_li8_n_p + 1
                    elif num_n == 1 and num_p == 0 and num_d == 1 and num_t == 1 and num_he3 == 0 and num_alpha == 0:
                        # Be10* -> n + d + t + He4:
                        number_be10_he4_n_d_t = number_be10_he4_n_d_t + 1
                    elif num_n == 1 and num_p == 1 and num_d == 0 and num_t == 1 and num_he3 == 0 and num_alpha == 0:
                        # Be10* -> n + p + t + He5:
                        number_be10_he5_n_p_t = number_be10_he5_n_p_t + 1
                    elif num_n == 1 and num_p == 0 and num_d == 1 and num_t == 0 and num_he3 == 0 and num_alpha == 0:
                        # Be10* -> n + d + Li7:
                        number_be10_li7_n_d = number_be10_li7_n_d + 1
                    elif num_n == 1 and num_p == 1 and num_d == 0 and num_t == 0 and num_he3 == 0 and num_alpha == 1:
                        # Be10* -> n + p + alpha + H4:
                        number_be10_h4_n_p_alpha = number_be10_h4_n_p_alpha + 1
                    elif num_n == 2 and num_p == 1 and num_d == 1 and num_t == 0 and num_he3 == 0 and num_alpha == 0:
                        # Be10* -> 2n + p + d + He5:
                        number_be10_he5_2n_p_d = number_be10_he5_2n_p_d + 1
                    elif num_n == 1 and num_p == 0 and num_d == 0 and num_t == 0 and num_he3 == 0 and num_alpha == 0:
                        # Be10* -> n + Be9:
                        number_be10_be9_n = number_be10_be9_n + 1
                    elif num_n == 1 and num_p == 1 and num_d == 1 and num_t == 0 and num_he3 == 0 and num_alpha == 0:
                        # Be10* -> n + p + d + He6:
                        number_be10_he6_n_p_d = number_be10_he6_n_p_d + 1
                    elif num_n == 0 and num_p == 0 and num_d == 1 and num_t == 1 and num_he3 == 0 and num_alpha == 0:
                        # Be10* -> d + t + He5:
                        number_be10_he5_d_t = number_be10_he5_d_t + 1
                    elif num_n == 0 and num_p == 0 and num_d == 1 and num_t == 0 and num_he3 == 0 and num_alpha == 1:
                        # Be10* -> d + alpha + H4:
                        number_be10_h4_d_alpha = number_be10_h4_d_alpha + 1
                    elif num_n == 2 and num_p == 1 and num_d == 0 and num_t == 1 and num_he3 == 0 and num_alpha == 0:
                        # Be10* -> 2n + p + t + He4:
                        number_be10_he4_2n_p_t = number_be10_he4_2n_p_t + 1
                    else:
                        number_be10_missing = number_be10_missing + 1
                        print("----------Be10-------")
                        print(deex_id[index])


            elif isotope_pdg[index] == 1000060090:
                # C9:
                if deex_id[index] == 0:
                    # C9 is not excited:
                    number_c9_notex = number_c9_notex + 1

                else:
                    # C9 deexcitation:
                    number_c9_deex = number_c9_deex + 1
                    # get number of particles:
                    num_n, num_p, num_d, num_t, num_he3, num_alpha = get_number_of_particles_of_deexid(deex_id[index])

                    if num_n == 0 and num_p == 2 and num_d == 0 and num_t == 0 and num_he3 == 0 and num_alpha == 0:
                        # C9* -> 2p + Be7:
                        number_c9_be7_2p = number_c9_be7_2p + 1
                    else:
                        number_c9_missing = number_c9_missing + 1
                        print("----------C9-------")
                        print(deex_id[index])


            elif isotope_pdg[index] == 1000050090:
                # B9:
                if deex_id[index] == 0:
                    # B9 is not excited:
                    number_b9_notex = number_b9_notex + 1

                else:
                    # B9 deexcitation:
                    number_b9_deex = number_b9_deex + 1
                    # get number of particles:
                    num_n, num_p, num_d, num_t, num_he3, num_alpha = get_number_of_particles_of_deexid(deex_id[index])

                    if num_n == 0 and num_p == 1 and num_d == 0 and num_t == 0 and num_he3 == 0 and num_alpha == 0:
                        # B9* -> p + Be8:
                        number_b9_be8_p = number_b9_be8_p + 1
                    else:
                        number_b9_missing = number_b9_missing + 1
                        print("----------B9-------")
                        print(deex_id[index])


            elif isotope_pdg[index] == 1000040090:
                # Be9:
                if deex_id[index] == 0:
                    # Be9 is not excited:
                    number_be9_notex = number_be9_notex + 1

                else:
                    # Be9 deexcitation:
                    number_be9_deex = number_be9_deex + 1
                    # get number of particles:
                    num_n, num_p, num_d, num_t, num_he3, num_alpha = get_number_of_particles_of_deexid(deex_id[index])

                    if num_n == 0 and num_p == 1 and num_d == 0 and num_t == 0 and num_he3 == 0 and num_alpha == 0:
                        # Be9* -> p + Li8:
                        number_be9_li8_p = number_be9_li8_p + 1
                    else:
                        number_be9_missing = number_be9_missing + 1
                        print("----------Be9-------")
                        print(deex_id[index])


            elif isotope_pdg[index] == 1000030090:
                # Li9:
                if deex_id[index] == 0:
                    # Li9 is not excited:
                    number_li9_notex = number_li9_notex + 1

                else:
                    # Li9 deexcitation:
                    number_li9_deex = number_li9_deex + 1
                    # get number of particles:
                    num_n, num_p, num_d, num_t, num_he3, num_alpha = get_number_of_particles_of_deexid(deex_id[index])

                    if num_n == 1 and num_p == 0 and num_d == 0 and num_t == 0 and num_he3 == 0 and num_alpha == 1:
                        # Li9* -> n + alpha + H4:
                        number_li9_h4_n_alpha = number_li9_h4_n_alpha + 1
                    elif num_n == 0 and num_p == 0 and num_d == 1 and num_t == 0 and num_he3 == 0 and num_alpha == 0:
                        # Li9* -> d + He7:
                        number_li9_he7_d = number_li9_he7_d + 1
                    elif num_n == 1 and num_p == 0 and num_d == 0 and num_t == 0 and num_he3 == 0 and num_alpha == 0:
                        # Li9* -> n + Li8:
                        number_li9_li8_n = number_li9_li8_n + 1
                    else:
                        number_li9_missing = number_li9_missing + 1
                        print("----------Li9-------")
                        print(deex_id[index])


            elif isotope_pdg[index] == 1000050080:
                # B8:
                if deex_id[index] == 0:
                    # B8 is not excited:
                    number_b8_notex = number_b8_notex + 1

                else:
                    # B8 deexcitation:
                    number_b8_deex = number_b8_deex + 1
                    # get number of particles:
                    num_n, num_p, num_d, num_t, num_he3, num_alpha = get_number_of_particles_of_deexid(deex_id[index])

                    if num_n == 0 and num_p == 2 and num_d == 0 and num_t == 0 and num_he3 == 0 and num_alpha == 0:
                        # B8* -> 2p + Li6:
                        number_b8_li6_2p = number_b8_li6_2p + 1
                    else:
                        number_b8_missing = number_b8_missing + 1
                        print("----------B8-------")
                        print(deex_id[index])


            elif isotope_pdg[index] == 1000030080:
                # Li8:
                if deex_id[index] == 0:
                    # Li8 is not excited:
                    number_li8_notex = number_li8_notex + 1

                else:
                    # Li8 deexcitation:
                    number_li8_deex = number_li8_deex + 1
                    # get number of particles:
                    num_n, num_p, num_d, num_t, num_he3, num_alpha = get_number_of_particles_of_deexid(deex_id[index])

                    if num_n == 1 and num_p == 0 and num_d == 0 and num_t == 0 and num_he3 == 0 and num_alpha == 0:
                        # Li8* -> n + Li7:
                        number_li8_li7_n = number_li8_li7_n + 1
                    elif num_n == 2 and num_p == 0 and num_d == 0 and num_t == 0 and num_he3 == 0 and num_alpha == 0:
                        # Li8* -> 2n + Li6:
                        number_li8_li6_2n = number_li8_li6_2n + 1
                    else:
                        number_li8_missing = number_li8_missing + 1
                        print("----------Li8-------")
                        print(deex_id[index])


            elif isotope_pdg[index] == 1000040070:
                # Be7:
                if deex_id[index] == 0:
                    # Be7 is not excited:
                    number_be7_notex = number_be7_notex + 1

                else:
                    # Be7 deexcitation:
                    number_be7_deex = number_be7_deex + 1
                    # get number of particles:
                    num_n, num_p, num_d, num_t, num_he3, num_alpha = get_number_of_particles_of_deexid(deex_id[index])

                    if num_n == 0 and num_p == 0 and num_d == 1 and num_t == 0 and num_he3 == 0 and num_alpha == 0:
                        # Be7* -> d + Li5:
                        number_be7_li5_d = number_be7_li5_d + 1
                    elif num_n == 0 and num_p == 1 and num_d == 0 and num_t == 0 and num_he3 == 0 and num_alpha == 0:
                        # Be7* -> p + Li6:
                        number_be7_li6_p = number_be7_li6_p + 1
                    else:
                        number_be7_missing = number_be7_missing + 1
                        print("----------Be7-------")
                        print(deex_id[index])


            elif isotope_pdg[index] == 1000030070:
                # Li7:
                if deex_id[index] == 0:
                    # Li7 is not excited:
                    number_li7_notex = number_li7_notex + 1

                else:
                    # Li7 deexcitation:
                    number_li7_deex = number_li7_deex + 1
                    # get number of particles:
                    num_n, num_p, num_d, num_t, num_he3, num_alpha = get_number_of_particles_of_deexid(deex_id[index])

                    if num_n == 1 and num_p == 0 and num_d == 0 and num_t == 0 and num_he3 == 0 and num_alpha == 0:
                        # Li7* -> n + Li6:
                        number_li7_li6_n = number_li7_li6_n + 1
                    else:
                        number_li7_missing = number_li7_missing + 1
                        print("----------Li7-------")
                        print(deex_id[index])


            else:
                # check if deex_if = 0 for all other isotopes (C8, Be8, He8, B7, He7, H7, Be6, Li6, He6, H6, Li5, He5,
                # H5, Li4, He4, H4, He3, H3, H2):
                if deex_id[index] == 0:
                    number_light_iso = number_light_iso + 1
                else:
                    print("WARNING: deex_id = {0:d}, BUT isotope has no deexcitation root file!!")

        else:
            # other target than C12:
            number_no_c12 = number_no_c12 + 1

    return (number_entries, number_target_c12, number_no_c12, number_light_iso,
            number_c11_notex, number_c11_deex, number_c11_li6_p_alpha, number_c11_be7_alpha, number_c11_b10_p,
            number_c11_b9_n_p, number_c11_be8_p_d, number_c11_be9_2p, number_c11_b9_d, number_c11_be8_he3,
            number_c11_c10_n, number_c11_li5_d_alpha, number_c11_li5_n_p_alpha, number_c11_missing,
            number_b11_notex, number_b11_deex, number_b11_li6_n_alpha, number_b11_b9_2n, number_b11_be8_n_d,
            number_b11_be9_d, number_b11_be10_p, number_b11_b10_n, number_b11_be9_n_p, number_b11_li7_alpha,
            number_b11_be8_t, number_b11_he5_d_alpha, number_b11_he6_p_alpha, number_b11_be8_2n_p, number_b11_missing,
            number_c10_notex, number_c10_deex, number_c10_b9_p, number_c10_be7_p_d, number_c10_li6_p_he3,
            number_c10_he4_p_d_he3, number_c10_li6_2p_d, number_c10_be8_2p, number_c10_be7_n_2p, number_c10_li6_n_3p,
            number_c10_be6_n_p_d, number_c10_li5_n_2p_d, number_c10_he3_p_d_alpha, number_c10_li5_d_he3,
            number_c10_li5_p_2d, number_c10_he3_n_2p_alpha, number_c10_li4_n_p_alpha, number_c10_b8_n_p,
            number_c10_b8_d, number_c10_be6_p_t, number_c10_he4_n_2p_he3, number_c10_li5_n_p_he3, number_c10_missing,
            number_b10_notex, number_b10_deex, number_b10_be9_p, number_b10_be8_d, number_b10_b9_n, number_b10_be7_t,
            number_b10_be8_n_p, number_b10_li7_he3, number_b10_he5_p_alpha, number_b10_li6_alpha,
            number_b10_li5_n_alpha, number_b10_li7_p_d, number_b10_missing,
            number_be10_notex, number_be10_deex, number_be10_li6_2n_d, number_be10_li7_2n_p, number_be10_he7_n_2p,
            number_be10_li6_n_t, number_be10_li6_3n_p, number_be10_he5_n_2d, number_be10_be8_2n,
            number_be10_he5_n_alpha, number_be10_t_n_d_alpha, number_be10_li8_n_p, number_be10_he4_n_d_t,
            number_be10_he5_n_p_t, number_be10_li7_n_d, number_be10_h4_n_p_alpha, number_be10_he5_2n_p_d,
            number_be10_be9_n, number_be10_he6_n_p_d, number_be10_he5_d_t, number_be10_h4_d_alpha,
            number_be10_he4_2n_p_t, number_be10_missing,
            number_c9_notex, number_c9_deex, number_c9_be7_2p, number_c9_missing,
            number_b9_notex, number_b9_deex, number_b9_be8_p, number_b9_missing,
            number_be9_notex, number_be9_deex, number_be9_li8_p, number_be9_missing,
            number_li9_notex, number_li9_deex, number_li9_h4_n_alpha, number_li9_he7_d, number_li9_li8_n,
            number_li9_missing,
            number_b8_notex, number_b8_deex, number_b8_li6_2p, number_b8_missing,
            number_li8_notex, number_li8_deex, number_li8_li7_n, number_li8_li6_2n, number_li8_missing,
            number_be7_notex, number_be7_deex, number_be7_li5_d, number_be7_li6_p, number_be7_missing,
            number_li7_notex, number_li7_deex, number_li7_li6_n, number_li7_missing)


def get_residual_isotopes_before_deex(projectile_energy, isotope_pdg, target_pdg, bin_width):
    """
    function to get the number of events as function of the energy of the incoming neutrino for the different produced
    residual isotopes (isotopes after NC interaction, but BEFORE deexcitation)

    :param projectile_energy: energy of the projectile, i.e. of the incoming neutrinos, in GeV (array of float)
    :param isotope_pdg: PDG ID of the isotope, that is created via the NC interaction (before deexcitation)
    (array of float)
    :param target_pdg: PDG ID of the target particle (array of float)
    :param bin_width: bin width of the array, which represents the energy of incoming neutrinos in GeV (float)
    :return:
    """
    # get the number of entries of the array (integer):
    number_entries = len(projectile_energy)

    # preallocate arrays, where projectile energy for the different isotopes is saved and number of events corresponding
    # to this isotope:
    # C12:
    energy_c12 = np.array([])
    number_c12 = 0
    # C11:
    energy_c11 = np.array([])
    number_c11 = 0
    # B11:
    energy_b11 = np.array([])
    number_b11 = 0
    # C10:
    energy_c10 = np.array([])
    number_c10 = 0
    # B10:
    energy_b10 = np.array([])
    number_b10 = 0
    # Be10:
    energy_be10 = np.array([])
    number_be10 = 0
    # C9:
    energy_c9 = np.array([])
    number_c9 = 0
    # B9:
    energy_b9 = np.array([])
    number_b9 = 0
    # Be9:
    energy_be9 = np.array([])
    number_be9 = 0
    # Li9:
    energy_li9 = np.array([])
    number_li9 = 0
    # C8:
    energy_c8 = np.array([])
    number_c8 = 0
    # B8:
    energy_b8 = np.array([])
    number_b8 = 0
    # Be8:
    energy_be8 = np.array([])
    number_be8 = 0
    # Li8:
    energy_li8 = np.array([])
    number_li8 = 0
    # He8:
    energy_he8 = np.array([])
    number_he8 = 0
    # B7:
    energy_b7 = np.array([])
    number_b7 = 0
    # Be7:
    energy_be7 = np.array([])
    number_be7 = 0
    # Li7:
    energy_li7 = np.array([])
    number_li7 = 0
    # He7:
    energy_he7 = np.array([])
    number_he7 = 0
    # H7:
    energy_h7 = np.array([])
    number_h7 = 0
    # Be6:
    energy_be6 = np.array([])
    number_be6 = 0
    # Li6:
    energy_li6 = np.array([])
    number_li6 = 0
    # He6:
    energy_he6 = np.array([])
    number_he6 = 0
    # H6:
    energy_h6 = np.array([])
    number_h6 = 0
    # Rest (isotopes not yet included into channel ID):
    number_rest = 0

    # loop over all entries (e.g. events):
    for index in range(number_entries):
        # check if target is C12 (1000060120):
        if target_pdg[index] == 1000060120:
            # check the isotope PDG:
            if isotope_pdg[index] == 1000060120:
                # C12:
                energy = projectile_energy[index]
                energy_c12 = np.append(energy_c12, energy)
                number_c12 += 1
            elif isotope_pdg[index] == 1000060110:
                # C11:
                energy = projectile_energy[index]
                energy_c11 = np.append(energy_c11, energy)
                number_c11 += 1
            elif isotope_pdg[index] == 1000050110:
                # B11:
                energy = projectile_energy[index]
                energy_b11 = np.append(energy_b11, energy)
                number_b11 += 1
            elif isotope_pdg[index] == 1000060100:
                # C10:
                energy = projectile_energy[index]
                energy_c10 = np.append(energy_c10, energy)
                number_c10 += 1
            elif isotope_pdg[index] == 1000050100:
                # B10:
                energy = projectile_energy[index]
                energy_b10 = np.append(energy_b10, energy)
                number_b10 += 1
            elif isotope_pdg[index] == 1000040100:
                # Be10:
                energy = projectile_energy[index]
                energy_be10 = np.append(energy_be10, energy)
                number_be10 += 1
            elif isotope_pdg[index] == 1000060090:
                # C9:
                energy = projectile_energy[index]
                energy_c9 = np.append(energy_c9, energy)
                number_c9 += 1
            elif isotope_pdg[index] == 1000050090:
                # B9:
                energy = projectile_energy[index]
                energy_b9 = np.append(energy_b9, energy)
                number_b9 += 1
            elif isotope_pdg[index] == 1000040090:
                # Be9:
                energy = projectile_energy[index]
                energy_be9 = np.append(energy_be9, energy)
                number_be9 += 1
            elif isotope_pdg[index] == 1000030090:
                # Li9:
                energy = projectile_energy[index]
                energy_li9 = np.append(energy_li9, energy)
                number_li9 += 1
            elif isotope_pdg[index] == 1000060080:
                # C8:
                energy = projectile_energy[index]
                energy_c8 = np.append(energy_c8, energy)
                number_c8 += 1
            elif isotope_pdg[index] == 1000050080:
                # B8:
                energy = projectile_energy[index]
                energy_b8 = np.append(energy_b8, energy)
                number_b8 += 1
            elif isotope_pdg[index] == 1000040080:
                # Be8:
                energy = projectile_energy[index]
                energy_be8 = np.append(energy_be8, energy)
                number_be8 += 1
            elif isotope_pdg[index] == 1000030080:
                # Li8:
                energy = projectile_energy[index]
                energy_li8 = np.append(energy_li8, energy)
                number_li8 += 1
            elif isotope_pdg[index] == 1000020080:
                # He8:
                energy = projectile_energy[index]
                energy_he8 = np.append(energy_he8, energy)
                number_he8 += 1
            elif isotope_pdg[index] == 1000050070:
                # B7:
                energy = projectile_energy[index]
                energy_b7 = np.append(energy_b7, energy)
                number_b7 += 1
            elif isotope_pdg[index] == 1000040070:
                # Be7:
                energy = projectile_energy[index]
                energy_be7 = np.append(energy_be7, energy)
                number_be7 += 1
            elif isotope_pdg[index] == 1000030070:
                # Li7:
                energy = projectile_energy[index]
                energy_li7 = np.append(energy_li7, energy)
                number_li7 += 1
            elif isotope_pdg[index] == 1000020070:
                # He7:
                energy = projectile_energy[index]
                energy_he7 = np.append(energy_he7, energy)
                number_he7 += 1
            elif isotope_pdg[index] == 1000010070:
                # H7:
                energy = projectile_energy[index]
                energy_h7 = np.append(energy_h7, energy)
                number_h7 += 1
            elif isotope_pdg[index] == 1000040060:
                # Be6:
                energy = projectile_energy[index]
                energy_be6 = np.append(energy_be6, energy)
                number_be6 += 1
            elif isotope_pdg[index] == 1000030060:
                # Li6:
                energy = projectile_energy[index]
                energy_li6 = np.append(energy_li6, energy)
                number_li6 += 1
            elif isotope_pdg[index] == 1000020060:
                # He6:
                energy = projectile_energy[index]
                energy_he6 = np.append(energy_he6, energy)
                number_he6 += 1
            elif isotope_pdg[index] == 1000010060:
                # H6:
                energy = projectile_energy[index]
                energy_h6 = np.append(energy_h6, energy)
                number_h6 += 1
            else:
                # isotopes not yet included into channel ID:
                number_rest += 1

        else:
            print("-----WARNING: other target particle than C12!")

    """ get fraction of events for the different residual isotopes BEFORE de-exciation (IN PERCENT): """
    # fraction of NC interaction events with residual isotopes in % (float):
    fraction_c12 = float(number_c12)/float(number_entries)*100
    fraction_c11 = float(number_c11)/float(number_entries)*100
    fraction_b11 = float(number_b11)/float(number_entries)*100
    fraction_c10 = float(number_c10)/float(number_entries)*100
    fraction_b10 = float(number_b10)/float(number_entries)*100
    fraction_be10 = float(number_be10)/float(number_entries)*100
    fraction_c9 = float(number_c9)/float(number_entries)*100
    fraction_b9 = float(number_b9)/float(number_entries)*100
    fraction_be9 = float(number_be9)/float(number_entries)*100
    fraction_li9 = float(number_li9)/float(number_entries)*100
    fraction_c8 = float(number_c8)/float(number_entries)*100
    fraction_b8 = float(number_b8)/float(number_entries)*100
    fraction_be8 = float(number_be8)/float(number_entries)*100
    fraction_li8 = float(number_li8)/float(number_entries)*100
    fraction_he8 = float(number_he8)/float(number_entries)*100
    fraction_b7 = float(number_b7)/float(number_entries)*100
    fraction_be7 = float(number_be7)/float(number_entries)*100
    fraction_li7 = float(number_li7)/float(number_entries)*100
    fraction_he7 = float(number_he7)/float(number_entries)*100
    fraction_h7 = float(number_h7)/float(number_entries)*100
    fraction_be6 = float(number_be6)/float(number_entries)*100
    fraction_li6 = float(number_li6)/float(number_entries)*100
    fraction_he6 = float(number_he6)/float(number_entries)*100
    fraction_h6 = float(number_h6)/float(number_entries)*100
    fraction_rest = float(number_rest)/float(number_entries)*100

    """ create histograms with the energy arrays from above to get the number of events per bin: """
    energy_range = np.arange(0, np.max(projectile_energy)+2, bin_width)
    events_c12, bins1 = np.histogram(energy_c12, energy_range)
    events_c11, bins1 = np.histogram(energy_c11, energy_range)
    events_b11, bins1 = np.histogram(energy_b11, energy_range)
    events_c10, bins1 = np.histogram(energy_c10, energy_range)
    events_b10, bins1 = np.histogram(energy_b10, energy_range)
    events_be10, bins1 = np.histogram(energy_be10, energy_range)
    events_c9, bins1 = np.histogram(energy_c9, energy_range)
    events_b9, bins1 = np.histogram(energy_b9, energy_range)
    events_be9, bins1 = np.histogram(energy_be9, energy_range)
    events_li9, bins1 = np.histogram(energy_li9, energy_range)
    events_c8, bins1 = np.histogram(energy_c8, energy_range)
    events_b8, bins1 = np.histogram(energy_b8, energy_range)
    events_be8, bins1 = np.histogram(energy_be8, energy_range)
    events_li8, bins1 = np.histogram(energy_li8, energy_range)
    events_he8, bins1 = np.histogram(energy_he8, energy_range)
    events_b7, bins1 = np.histogram(energy_b7, energy_range)
    events_be7, bins1 = np.histogram(energy_be7, energy_range)
    events_li7, bins1 = np.histogram(energy_li7, energy_range)
    events_he7, bins1 = np.histogram(energy_he7, energy_range)
    events_h7, bins1 = np.histogram(energy_h7, energy_range)
    events_be6, bins1 = np.histogram(energy_be6, energy_range)
    events_li6, bins1 = np.histogram(energy_li6, energy_range)
    events_he6, bins1 = np.histogram(energy_he6, energy_range)
    events_h6, bins1 = np.histogram(energy_h6, energy_range)

    return (energy_range, events_c12, events_c11, events_b11, events_c10, events_b10, events_be10, events_c9, events_b9,
            events_be9, events_li9, events_c8, events_b8, events_be8, events_li8, events_he8, events_b7, events_be7,
            events_li7, events_he7, events_h7, events_be6, events_li6, events_he6, events_h6,
            number_c12, number_c11, number_b11, number_c10, number_b10, number_be10, number_c9, number_b9, number_be9,
            number_li9, number_c8, number_b8, number_be8, number_li8, number_he8, number_b7, number_be7, number_li7,
            number_he7, number_h7, number_be6, number_li6, number_he6, number_h6, number_rest,
            fraction_c12, fraction_c11, fraction_b11, fraction_c10, fraction_b10, fraction_be10, fraction_c9,
            fraction_b9, fraction_be9, fraction_li9, fraction_c8, fraction_b8, fraction_be8, fraction_li8, fraction_he8,
            fraction_b7, fraction_be7, fraction_li7, fraction_he7, fraction_h7, fraction_be6, fraction_li6,
            fraction_he6, fraction_h6, fraction_rest)


def check_high_channelid(channel_id, final_pdg):
    """
    function to check the channel IDs with large value and to see which particle are produced for this channels.

    :param channel_id: ID of the NC interaction channel, defines the product particles of the NC interactions
    (array of float)
    :param final_pdg: PDG ID of all final particles for each event (list of arrays of floats)
    :return:
    """

    # get the number of entries of the array (integer):
    number_entries = len(channel_id)

    # loop over all entries:
    for index in range(number_entries):

        # check, if channel ID is larger than 5 digits:
        if channel_id[index] > 100000:

            print("----------------------------------")
            print("channel ID: {0:.0f}".format(channel_id[index]))
            print("final particles:")
            print(final_pdg[index][:])
            print("----------------------------------")

        else:
            continue

    return


def check_nc_qel_from_original_genie_file(rootfile_input):
    """
    function to check the variables 'nc' and 'qel' of the 'original' GENIE root-file from Julia.

    :param rootfile_input: path to the original GENIE ROOT-file (for example: gntp.101.gst.root (string)

    :return:
    """
    # load the ROOT file:
    rfile_input = ROOT.TFile(rootfile_input)
    # get the TTree from the TFile:
    rtree_input = rfile_input.Get("gst")

    # Info-me: "gst;13" is a copy of meta data of "gst;14", "gst;14" contains correct data and is read

    # get the number of entries in the ROOT-file:
    # number_entries = rtree_input.GetEntries()
    number_entries = 10

    """ Read the data from the TTree: """
    # loop over every entry, i.e. every event, in the TTree:
    for event in range(number_entries):

        # get the current event in the TTree:
        rtree_input.GetEntry(event)

        # is it a quasi-elastic scattering event? (0 = no QEL event, 1 = QEL event):
        qel = rtree_input.GetBranch('qel').GetLeaf('qel').GetValue()
        qel = int(qel)

        # is it a NC event? (0 = no NC event, 1 = NC event):
        nc = rtree_input.GetBranch('nc').GetLeaf('nc').GetValue()
        nc = int(nc)

        # get the value of target PDG:
        tgt = rtree_input.GetBranch('tgt').GetLeaf('tgt').GetValue()
        tgt = int(tgt)

        # final particles:
        # get the value of number of final p:
        nfp = rtree_input.GetBranch('nfp').GetLeaf('nfp').GetValue()
        nfp = int(nfp)

        # get the value of number of final n:
        nfn = rtree_input.GetBranch('nfn').GetLeaf('nfn').GetValue()
        nfn = int(nfn)

        # get the value of number of final pi_minus:
        nfpim = rtree_input.GetBranch('nfpim').GetLeaf('nfpim').GetValue()
        nfpim = int(nfpim)

        # get the value of number of final pi_plus:
        nfpip = rtree_input.GetBranch('nfpip').GetLeaf('nfpip').GetValue()
        nfpip = int(nfpip)

        # get value of number of final pi_zero:
        nfpi0 = rtree_input.GetBranch('nfpi0').Getleaf('nfpi0').GetValue()
        nfpi0 = int(nfpi0)

        # get the value of number of final Kaon_minus:
        nfkm = rtree_input.GetBranch('nfkm').GetLeaf('nfkm').GetValue()
        nfkm = int(nfkm)

        # get the value of number of final Kaon_plus:
        nfkp = rtree_input.GetBranch('nfkp').GetLeaf('nfkp').GetValue()
        nfkp = int(nfkp)

        # get value of number of final K_zero:
        nfk0 = rtree_input.GetBranch('nfk0').Getleaf('nfk0').GetValue()
        nfk0 = int(nfk0)

        # get value of number of final gamma, electron, positron:
        nfem = rtree_input.GetBranch('nfem').Getleaf('nfem').GetValue()
        nfem = int(nfem)

        # get value of number of final other hadrons:
        nfother = rtree_input.GetBranch('nfother').Getleaf('nfother').GetValue()
        nfother = int(nfother)

        # primary particles:
        # get the value of number of primary p:
        nip = rtree_input.GetBranch('nip').GetLeaf('nip').GetValue()
        nip = int(nip)

        # get the value of number of primary n:
        nin = rtree_input.GetBranch('nin').GetLeaf('nin').GetValue()
        nin = int(nin)

        # get the value of number of primary pi_minus:
        nipim = rtree_input.GetBranch('nipim').GetLeaf('nipim').GetValue()
        nipim = int(nipim)

        # get the value of number of primary pi_plus:
        nipip = rtree_input.GetBranch('nipip').GetLeaf('nipip').GetValue()
        nipip = int(nipip)

        # get value of number of primary pi_zero:
        nipi0 = rtree_input.GetBranch('nipi0').Getleaf('nipi0').GetValue()
        nipi0 = int(nipi0)

        # get the value of number of primary Kaon_minus:
        nikm = rtree_input.GetBranch('nikm').GetLeaf('nikm').GetValue()
        nikm = int(nikm)

        # get the value of number of primary Kaon_plus:
        nikp = rtree_input.GetBranch('nikp').GetLeaf('nikp').GetValue()
        nikp = int(nikp)

        # get value of number of primary K_zero:
        nik0 = rtree_input.GetBranch('nik0').Getleaf('nik0').GetValue()
        nik0 = int(nik0)

        # get value of number of primary gamma, electron, positron:
        niem = rtree_input.GetBranch('niem').Getleaf('niem').GetValue()
        niem = int(niem)

        # get value of number of primary other hadrons:
        niother = rtree_input.GetBranch('niother').Getleaf('niother').GetValue()
        niother = int(niother)

        # check primary and final particles:
        if tgt == 1000060120:
            # target C12

            # B11:
            if nfp == 1 and nfn == 0 and nfpim == 0 and nfpip == 0 and nfkm == 0 and nfkp == 0:
                # interaction channel: nu + C12 -> B11 + proton:
                print("nu + C12 -> nu + B11 + p: nc={0}, qel={1},\n "
                      "nfpi0={2:d}, nfk0={3:d}, nfem={4:d}, nfother={5:d}\n"
                      "nip={6:d}, nin={7:d}, nipip={8:d}, nipim={9:d}, nipi0={10:d}, nikp={11:d}, nikm={12:d}, "
                      "nik0={13:d}, niem={14:d}, niother={15:d}".format(nc, qel, nfpi0, nfk0, nfem, nfother, nip, nin,
                                                                        nipip, nipim, nipi0, nikp, nikm, nik0, niem,
                                                                        niother))

            elif nfp == 0 and nfn == 1 and nfpim == 0 and nfpip == 1 and nfkm == 0 and nfkp == 0:
                # interaction channel: nu + C12 -> B11 + n + pi_plus:
                print("nu + C12 -> nu + B11 + n + pi_plus: nc={0}, qel={1},\n "
                      "nfpi0={2:d}, nfk0={3:d}, nfem={4:d}, nfother={5:d}\n"
                      "nip={6:d}, nin={7:d}, nipip={8:d}, nipim={9:d}, nipi0={10:d}, nikp={11:d}, nikm={12:d}, "
                      "nik0={13:d}, niem={14:d}, niother={15:d}".format(nc, qel, nfpi0, nfk0, nfem, nfother, nip, nin,
                                                                        nipip, nipim, nipi0, nikp, nikm, nik0, niem,
                                                                        niother))
            # C11:
            elif nfp == 0 and nfn == 1 and nfpim == 0 and nfpip == 0 and nfkm == 0 and nfkp == 0:
                # interaction channel: nu + C12 -> C11 + n:
                print("nu + C12 -> nu + C11 + n: nc={0}, qel={1},\n "
                      "nfpi0={2:d}, nfk0={3:d}, nfem={4:d}, nfother={5:d}\n"
                      "nip={6:d}, nin={7:d}, nipip={8:d}, nipim={9:d}, nipi0={10:d}, nikp={11:d}, nikm={12:d}, "
                      "nik0={13:d}, niem={14:d}, niother={15:d}".format(nc, qel, nfpi0, nfk0, nfem, nfother, nip, nin,
                                                                        nipip, nipim, nipi0, nikp, nikm, nik0, niem,
                                                                        niother))

            elif nfp == 1 and nfn == 0 and nfpim == 1 and nfpip == 0 and nfkm == 0 and nfkp == 0:
                # interaction channel: nu + C12 -> C11 + p + pi_minus:
                print("nu + C12 -> nu + C11 + p + pi_minus: nc={0}, qel={1},\n "
                      "nfpi0={2:d}, nfk0={3:d}, nfem={4:d}, nfother={5:d}\n"
                      "nip={6:d}, nin={7:d}, nipip={8:d}, nipim={9:d}, nipi0={10:d}, nikp={11:d}, nikm={12:d}, "
                      "nik0={13:d}, niem={14:d}, niother={15:d}".format(nc, qel, nfpi0, nfk0, nfem, nfother, nip, nin,
                                                                        nipip, nipim, nipi0, nikp, nikm, nik0, niem,
                                                                        niother))

            # B10:
            elif nfp == 1 and nfn == 1 and nfpim == 0 and nfpip == 0 and nfkm == 0 and nfkp == 0:
                # interaction channel: nu + C12 -> B10 + p + n:
                print("nu + C12 -> nu + B10 + p + n: nc={0}, qel={1},\n "
                      "nfpi0={2:d}, nfk0={3:d}, nfem={4:d}, nfother={5:d}\n"
                      "nip={6:d}, nin={7:d}, nipip={8:d}, nipim={9:d}, nipi0={10:d}, nikp={11:d}, nikm={12:d}, "
                      "nik0={13:d}, niem={14:d}, niother={15:d}".format(nc, qel, nfpi0, nfk0, nfem, nfother, nip, nin,
                                                                        nipip, nipim, nipi0, nikp, nikm, nik0, niem,
                                                                        niother))

            elif nfp == 2 and nfn == 0 and nfpim == 1 and nfpip == 0 and nfkm == 0 and nfkp == 0:
                # interaction channel: nu + C12 -> B10 + 2*p + pi_minus:
                print("nu + C12 -> nu + B10 + 2p + pi_minus: nc={0}, qel={1},\n "
                      "nfpi0={2:d}, nfk0={3:d}, nfem={4:d}, nfother={5:d}\n"
                      "nip={6:d}, nin={7:d}, nipip={8:d}, nipim={9:d}, nipi0={10:d}, nikp={11:d}, nikm={12:d}, "
                      "nik0={13:d}, niem={14:d}, niother={15:d}".format(nc, qel, nfpi0, nfk0, nfem, nfother, nip, nin,
                                                                        nipip, nipim, nipi0, nikp, nikm, nik0, niem,
                                                                        niother))

            else:
                continue

    return
