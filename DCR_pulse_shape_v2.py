""" function to load pulse shape of atmoNC, IBD and fast neutron events from file and add dark counts to the
    pulse shape.
    These pulse shapes with dark counts are then also saved and can be analyzed with analyze_PSD_cut_v2.py.

"""
import datetime
import os
import re
import numpy as np


def get_numbers_from_filename(filename):
    """
    function to get number as integer out of a string. For example: filename='file235' -> num = 235 of type integer

    :param filename: string of one part of the filename 'file{}_evt{}_prompt_signal.txt'
    :return:
    """
    # get the number out of filename and convert it into integer:
    num = int(re.search(r'\d+', filename).group(0))

    return num


def pulse_shape_dcr(pathname, filename, dcr_20inch, n_20inch, dcr_3inch, n_3inch, evt_type):
    """
    function to add dark counts to pulse shape

    :param pathname: name of the path, where hittimes are saved
    :param filename: name of the pulse shape file
    :param dcr_20inch: DCR of 20inch PMTs
    :param n_20inch: number of 20inch PMTs
    :param dcr_3inch: DCR of 3inch PMTs
    :param n_3inch: number of 3inch PMTs
    :param evt_type: event type: 'atmoNC', 'IBD' or 'fastN'
    :return:
    """
    # split string file_ibdlike into two parts: 1. part: 'file{}', 2. part: 'vt{}_prompt_signal.txt' or
    # 'vt{}_pulse_shape_R16m.txt' for fast neutron events
    x = filename.split("_e")

    # x[0] is string 'file{}':
    file_string = x[0]
    # x[1] is string 'vt{}_prompt_signal.txt' or 'vt{}_pulse_shape_R16m.txt':
    end_string = x[1]
    # split end_string into two parts: 1. part 'vt{}_prompt'/'vt{}_pulse, 2.part: 'ignal.txt' or 'hape_R16m.txt':
    y = end_string.split("_s")
    # y[0] is string 'vt{}_prompt':
    event_string = y[0]
    # y[1] is string 'nal.txt' or 'nal_R16m.txt':
    rest_string = y[1]

    # get file_number of file_string:
    file_number = get_numbers_from_filename(file_string)
    # get evtID of event_string:
    evtid = get_numbers_from_filename(event_string)

    # for fast neutron events: get radius from the filename:
    if evt_type == 'fastN':
        cut_radius = get_numbers_from_filename(rest_string)

    # read file:
    pulse_shape = np.loadtxt(pathname + filename)

    # get reconstructed position in mm:
    x_reco = pulse_shape[0]
    y_reco = pulse_shape[1]
    z_reco = pulse_shape[2]

    # get start time of pulse shape in ns:
    time_start = pulse_shape[3]
    # get end time of pulse shape in ns:
    time_end = pulse_shape[4]
    # get bin-width of pulse shape in ns:
    bin_width = pulse_shape[5]

    # the rest is the pulse shape:
    pulse_shape = pulse_shape[6:]

    # time-window of the pulse shape in ns:
    time_window = time_end - time_start

    # calculate number of dark counts in the time window:
    number_dc = int(dcr_20inch * time_window * 10**(-9) * n_20inch +
                    dcr_3inch * time_window * 10**(-9) * n_3inch)

    # generate the time of the dark count with uniformly distributed random number in time_window
    # (for all number_dc):
    time_dc = np.random.uniform(time_start, time_end, size=number_dc)

    # build histogram with time_DC (must have same shape like pulse_shape):
    bins_hittime = np.arange(time_start, time_end+bin_width, bin_width)
    npe_dc, bin_edges_dc = np.histogram(time_dc, bins_hittime)

    # add npe_dc to pulse_shape to get pulse shape with DCR considered:
    pulse_shape_dc = pulse_shape + npe_dc

    # save new pulse shape to file:
    pulse_shape_dc_save = [x_reco, y_reco, z_reco]
    pulse_shape_dc_save.extend([time_start, time_end, bin_width])
    pulse_shape_dc_save.extend(pulse_shape_dc)
    if evt_type == 'fastN':
        np.savetxt(pathname + "/file{0:d}_evt{1:d}_pulse_shape_R{2:d}_DCR.txt".format(file_number, evtid, cut_radius),
                   pulse_shape_dc_save, fmt='%1.2f',
                   header="Pulse shape of prompt signal with DCR: Number of pe as function of the time "
                          "(time-of-flight correction, TTS smearing, DCR considered) of file user_{6}_{0:d}.root,"
                          "\nevent {1:d}, {2}:"
                          "\ntime window of pulse shape: from {3:.3f} ns to {4:.3f} ns with bin-width = {5:0.3f} "
                          "ns,"
                   .format(file_number, evtid, now, time_start, time_end, bin_width, evt_type))

    else:
        np.savetxt(pathname + "/file{0:d}_evt{1:d}_prompt_signal_DCR.txt".format(file_number, evtid),
                   pulse_shape_dc_save, fmt='%1.2f',
                   header="Pulse shape of prompt signal with DCR: Number of pe as function of the time "
                          "(time-of-flight correction, TTS smearing, DCR considered) of file user_{6}_{0:d}.root,"
                          "\nevent {1:d}, {2}:"
                          "\ntime window of pulse shape: from {3:.3f} ns to {4:.3f} ns with bin-width = {5:0.3f} "
                          "ns,"
                   .format(file_number, evtid, now, time_start, time_end, bin_width, evt_type))

    return


# get the date and time, when the script was run:
date = datetime.datetime.now()
now = date.strftime("%Y-%m-%d %H:%M")

# path, where the pulse shapes of NC events are stored:
path_NC = "/home/astro/blum/juno/atmoNC/data_NC/output_detsim_v2/hittimes/"
# path, where the pulse shapes of IBD events are stored:
path_IBD = "/home/astro/blum/juno/IBD_events/hittimes/"
# path, where the pulse shapes of Fast Neutron events are stored:
path_FN = "/home/astro/blum/PhD/work/MeVDM_JUNO/fast_neutrons/hittimes/"

""" Dark count rate parameters: """
# from file PmtData.root and PMT_position.root in folder /home/astro/blum/juno/atmoNC/PMT_information/:
# Dark count rate of small PMTs can be neglected
# total number of 20 inch PMTs:
number_20inchPMT = 17739
# number of Hamamatsu PMTs:
number_Hama = 4998
# number of MCP PMTs:
number_MCP = 12741
# total number of 3inch Pmts:
number_3inchPMT = 36572

# mean of Dark count rate of 20inch PMTs in Hz (/home/astro/blum/PhD/paper/PMT/20190114The progress of PMT test.pdf):
DCR_20inch = 31500.0
# Dark count rate of Hamamatsu PMTs in Hz:
DCR_Hama = 15500.0
# Dark count rate of MCP PMTs in Hz:
DCR_MCP = 46100.0
# Dark count rate of 3inch PMTs (/home/astro/blum/PhD/paper/PMT/hem_181115_spmt_review.pdf):
DCR_3inch = 550.0

""" loop over pulse shapes of prompt signals of atmo. NC events: """
# print("atmo NC events...")
# for file_NC in os.listdir(path_NC):
#     # read only files that start with 'file' and end with 'prompt_signal.txt'
#     if file_NC.startswith("file") and file_NC.endswith("prompt_signal.txt"):
#
#         pulse_shape_dcr(path_NC, file_NC, DCR_20inch, number_20inchPMT, DCR_3inch, number_3inchPMT, "atmoNC")

""" loop over pulse shapes of prompt signals of IBD events: """
# print("IBD events...")
# for file_IBD in os.listdir(path_IBD):
#     # read only files that start with 'file' and end with 'prompt_signal.txt'
#     if file_IBD.startswith("file") and file_IBD.endswith("prompt_signal.txt"):
#
#         pulse_shape_dcr(path_IBD, file_IBD, DCR_20inch, number_20inchPMT, DCR_3inch, number_3inchPMT, "IBD")

""" loop over pulse shapes of prompt signals of fast neutron events: """
print("fastN events...")
for file_FN in os.listdir(path_FN):
    # read only files that start with 'file' and end with 'prompt_signal.txt'
    if file_FN.startswith("file") and \
            (file_FN.endswith("pulse_shape_R16.txt") or file_FN.endswith("pulse_shape_R17.txt")):

        pulse_shape_dcr(path_FN, file_FN, DCR_20inch, number_20inchPMT, DCR_3inch, number_3inchPMT, "fastN")

