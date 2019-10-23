import ROOT
import numpy as np
from matplotlib import pyplot as plt
import NC_background_functions

# # path, where the information about the cuts on the IBD events is stored:
# input_path_IBD = "/home/astro/blum/juno/IBD_events/"
#
# # path, where the numbers_.txt files of the delayed cut for IBD events are stored:
# input_path_IBD_del = (input_path_IBD + "delayed_cut_{0:.0f}nsto{1:.0f}ms_mult{2:d}_{3:.0f}PEto{4:.0f}PE_dist{5:.0f}mm_"
#                                        "R{6:.0f}mm/".format(500.0, 1.0, 1, 2400.0, 3400.0, 500.0, 16000.0))
#
# # load files, where filenumber and evtID of events that pass the cut are stored (real/reconstructed data):
# array_pass_volume_cut_real_IBD = np.loadtxt(input_path_IBD + "filenumber_evtID_volume_cut_IBD_{0:.0f}mm.txt"
#                                             .format(16000.0))
#
# array_pass_delayed_cut_real_IBD = np.loadtxt(input_path_IBD_del + "filenumber_evtID_delayed_cut_IBD_{0:.0f}ns_to_"
#                                                                   "{1:.0f}ms_mult{2:d}_{3:.0f}PE_{4:.0f}PE_"
#                                                                   "dist{5:.0f}mm_R{6:.0f}mm_0.txt"
#                                              .format(500.0, 1.0, 1, 2400.0, 3400.0, 500.0, 16000.0))
#
# number = 0
#
# for index in range(len(array_pass_volume_cut_real_IBD)):
#
#     filenumber_vol = array_pass_volume_cut_real_IBD[index][0]
#     evtID_vol = array_pass_volume_cut_real_IBD[index][1]
#
#     for index1 in range(len(array_pass_delayed_cut_real_IBD)):
#
#         filenumber_del = array_pass_delayed_cut_real_IBD[index1][0]
#         evtID_del = array_pass_delayed_cut_real_IBD[index1][1]
#
#         if filenumber_del == filenumber_vol and evtID_del == evtID_vol:
#             number += 1
#             continue
#
# print(number)



# filenumber = 0
# event = 65
# filenumber2 = 200
# event2 = 5
# tail_start = 275
# tail_end = 600
# ttr_cut = 0.01660
#
# hittime = np.loadtxt("/home/astro/blum/juno/IBD_events/hittimes/file{0:d}_evt{1:d}_prompt_signal_DCR.txt"
#                      .format(filenumber, event))
# t_min = hittime[0]
# t_max = hittime[1]
# binwidth = hittime[2]
# nPE = hittime[3:(int((1000.0 + binwidth + np.abs(t_min)) / binwidth)+3)]
# nPE = nPE / np.sum(nPE)
# bins = np.arange(t_min, 1000.0+binwidth, binwidth)
#
# for index in range(len(bins)):
#     if bins[index] == 275:
#         start_index = index
#     elif bins[index] == 600:
#         stop_index = index
#
# ttr = np.sum(nPE[start_index:stop_index]) / np.sum(nPE)
#
#
# # hittime1 = np.loadtxt("/home/astro/blum/juno/IBD_events/hittimes/file{0:d}_evt{1:d}_prompt_signal.txt"
# #                       .format(filenumber, event))
# hittime1 = np.loadtxt("/home/astro/blum/juno/atmoNC/data_NC/output_detsim_v2/hittimes/"
#                       "file{0:d}_evt{1:d}_prompt_signal_DCR.txt".format(filenumber, event))
# t_min1 = hittime1[0]
# t_max1 = hittime1[1]
# binwidth1 = hittime1[2]
# nPE1 = hittime1[3:(int((1000.0 + binwidth1 + np.abs(t_min1)) / binwidth1)+3)]
# nPE1 = nPE1 / np.sum(nPE1)
# bins1 = np.arange(t_min1, 1000.0+binwidth1, binwidth1)
#
# for index in range(len(bins1)):
#     if bins1[index] == 275:
#         start_index1 = index
#     elif bins1[index] == 600:
#         stop_index1 = index
#
# ttr1 = np.sum(nPE1[start_index1:stop_index1]) / np.sum(nPE1)
#
# hittime2 = np.loadtxt("/home/astro/blum/PhD/work/MeVDM_JUNO/fast_neutrons/hittimes/"
#                       "file{0:d}_evt{1:d}_pulse_shape_R16_DCR.txt".format(filenumber2, event2))
# t_min2 = hittime2[0]
# t_max2 = hittime2[1]
# binwidth2 = hittime2[2]
# nPE2 = hittime2[3:(int((1000.0 + binwidth2 + np.abs(t_min2)) / binwidth2)+3)]
# nPE2 = nPE2 / np.sum(nPE2)
# bins2 = np.arange(t_min2, 1000.0+binwidth2, binwidth2)
#
# for index in range(len(bins2)):
#     if bins2[index] == 275:
#         start_index2 = index
#     elif bins2[index] == 600:
#         stop_index2 = index
#
# ttr2 = np.sum(nPE2[start_index2:stop_index2]) / np.sum(nPE2)
#
# plt.semilogy(bins, nPE, "r", label="IBD w/ DCR\nttr = {0:.5f}".format(ttr), drawstyle='steps')
# plt.semilogy(bins1, nPE1, "b", label="NC w/ DCR\nttr = {0:.5f}".format(ttr1), drawstyle='steps')
# plt.semilogy(bins2, nPE2, "g", label="Fast neutron w/ DCR\nttr = {0:.5f}".format(ttr2), drawstyle='steps')
# plt.vlines(tail_start, 0, np.max(nPE), "k", label="ttr cut value = {0:.5f}".format(ttr_cut))
# plt.vlines(tail_end, 0, np.max(nPE), "k")
# plt.xlabel("time in ns")
# plt.legend()
# plt.show()




# ttr_o = np.loadtxt("/home/astro/blum/juno/atmoNC/data_NC/output_PSD_v2/test_neff_1_55/TTR_IBD_350ns_1000ns_0.txt")
# ttr_n = np.loadtxt("/home/astro/blum/juno/atmoNC/data_NC/output_PSD_v2/TTR_IBD_350ns_1000ns_0.txt")
#
# ttr_old = []
# ttr_new = []
#
# for index in range(len(ttr_o)):
#     ttr_old.append(ttr_o[index][2])
# for index in range(len(ttr_n)):
#     ttr_new.append(ttr_n[index][2])
#
# mean_old = np.mean(ttr_old)
# mean_new = np.mean(ttr_new)
#
# print(mean_old)
# print(mean_new)
#
# Bins = np.arange(0, 0.1, 0.001)
# plt.hist(ttr_old, bins=Bins, histtype="step", align='mid', color="r", linewidth=1.5, label='c_eff = 194.67 mm/ns')
# plt.hist(ttr_new, bins=Bins, histtype="step", align='mid', color="b", linewidth=1.5, label="c_eff = 193.41 mm/ns")
# plt.xlim(xmin=0.0, xmax=0.1)
# plt.xlabel("tail-to-total ratio")
# plt.ylabel("events")
# plt.legend()
# plt.grid()
# plt.show()




