""" script to get the NC interaction channels (including gammas) of genie_data.root file and the deexcitation channels
    (including gammas) of gen_NC_onlyC12_250000evts_seed1.root

    Difference to checkout_NCgen.py:
        - the channels are calculated from the final PDG ID and not from channelID or deexID.
        -> also gammas are considered in this script!
"""
import ROOT
import datetime
import numpy as np
import NC_background_functions
from matplotlib import pyplot as plt

# get the date and time, when the script was run:
date = datetime.datetime.now()
now = date.strftime("%Y-%m-%d %H:%M")


def get_number_of_particles(pdg_arr):
    """
    function to get the number of neutron, proton, gamma, ... from array, where PDGs are stored
    :param pdg_arr: array with PDG
    :return: number of neutron, proton, piminus, piplus, pi0, gamma, ...
    """
    n = 0
    p = 0
    piminus = 0
    piplus = 0
    pi0 = 0
    gamma = 0
    deuteron = 0
    triton = 0
    he3 = 0
    alpha = 0
    li6 = 0
    li7 = 0
    li8 = 0
    li9 = 0
    be7 = 0
    be8 = 0
    be9 = 0
    be10 = 0
    b8 = 0
    b9 = 0
    b10 = 0
    b11 = 0
    c9 = 0
    c10 = 0
    c11 = 0
    c12 = 0

    # loop over pdg_array:
    for ind in range(len(pdg_arr)):

        if pdg_arr[ind] == 2112:
            n += 1
        elif pdg_arr[ind] == 2212:
            p += 1
        elif pdg_arr[ind] == -211:
            piminus += 1
        elif pdg_arr[ind] == 211:
            piplus += 1
        elif pdg_arr[ind] == 111:
            pi0 += 1
        elif pdg_arr[ind] == 22:
            gamma += 1
        elif pdg_arr[ind] == 1000010020:
            deuteron += 1
        elif pdg_arr[ind] == 1000010030:
            triton += 1
        elif pdg_arr[ind] == 1000020030:
            he3 += 1
        elif pdg_arr[ind] == 1000020040:
            alpha += 1
        elif pdg_arr[ind] == 1000030060:
            li6 += 1
        elif pdg_arr[ind] == 1000030070:
            li7 += 1
        elif pdg_arr[ind] == 1000030080:
            li8 += 1
        elif pdg_arr[ind] == 1000030090:
            li9 += 1
        elif pdg_arr[ind] == 1000040070:
            be7 += 1
        elif pdg_arr[ind] == 1000040080:
            be8 += 1
        elif pdg_arr[ind] == 1000040090:
            be9 += 1
        elif pdg_arr[ind] == 1000040100:
            be10 += 1
        elif pdg_arr[ind] == 1000050080:
            b8 += 1
        elif pdg_arr[ind] == 1000050090:
            b9 += 1
        elif pdg_arr[ind] == 1000050100:
            b10 += 1
        elif pdg_arr[ind] == 1000050110:
            b11 += 1
        elif pdg_arr[ind] == 1000060090:
            c9 += 1
        elif pdg_arr[ind] == 1000060100:
            c10 += 1
        elif pdg_arr[ind] == 1000060110:
            c11 += 1
        elif pdg_arr[ind] == 1000060120:
            c12 += 1
        elif (pdg_arr[ind] == -321 or pdg_arr[ind] == 321 or pdg_arr[ind] == 3122 or pdg_arr[ind] == 3112 or
              pdg_arr[ind] == 3222 or pdg_arr[ind] == -311 or pdg_arr[ind] == 311 or pdg_arr[ind] == -2112 or
              pdg_arr[ind] == -2212 or pdg_arr[ind] == -11 or pdg_arr[ind] == 11 or pdg_arr[ind] == 3212 or
              pdg_arr[ind] == 130):
            # Kaon, Sigma, Lambda, ...
            continue
        else:
            print("new PDG: {0:.0f}".format(pdg_arr[ind]))

    return (n, p, piminus, piplus, pi0, gamma, deuteron, triton, he3, alpha, li6, li7, li8, li9, be7, be8, be9, be10,
            b8, b9, b10, b11, c9, c10, c11, c12)


""" information after NC interaction and after deexcitation: """
# Should txt file be saved:
SAVE_TXT = True
# calculate deexcitation channels considering pi0 and gammas:
CALC_PI0_AND_GAMMA = False
# set the path to file gen_NC_onlyC12_250000evts_seed1.root:
file_generator = "/home/astro/blum/juno/atmoNC/data_NC/output_generator/gen_NC_onlyC12_250000evts_seed1.root"
# set output path:
output_path = "/home/astro/blum/juno/atmoNC/data_NC/output_checkoutNC/"
# open root file:
rfile_generator = ROOT.TFile(file_generator)
# get the TTree from the TFile:
rtree_generator = rfile_generator.Get("genEvt")
# get the number of entries, i.e. events, in the ROOT-file:
N_entries_generator = rtree_generator.GetEntries()

""" deexcitation channels with pi0 and gamma : """
if CALC_PI0_AND_GAMMA:
    # preallocate number of event with deexcitation:
    N_evts_deex = 0
    # preallocate number of events with incorrect deexcitation (e.g. number of protons -999, ...):
    N_evts_incorrect = 0
    # C11:
    number_c11_iso = 0
    n_c11_noDeex = 0
    n_c11_c11_Ygamma = 0
    n_c11_c11_Xpi0 = 0
    n_c11_c11_Xpi0_Ygamma = 0
    n_c11_c10_n_Xpi0_Ygamma = 0
    n_c11_c10_n_Ygamma = 0
    n_c11_b10_p_Ygamma = 0
    n_c11_b10_p_Xpi0_Ygamma = 0
    n_c11_b9_n_p_Ygamma = 0
    n_c11_b9_d_Ygamma = 0
    n_c11_b9_d_Xpi0_Ygamma = 0
    n_c11_b9_n_p_Xpi0_Ygamma = 0
    n_c11_be9_2p_Ygamma = 0
    n_c11_be9_2p_Xpi0_Ygamma = 0
    n_c11_be8_p_d_Xpi0_Ygamma = 0
    n_c11_be8_p_d_Ygamma = 0
    n_c11_be8_he3_Ygamma = 0
    n_c11_be8_he3_Xpi0_Ygamma = 0
    n_c11_be8_n_2p_Ygamma = 0
    n_c11_be8_n_2p_Xpi0_Ygamma = 0
    n_c11_be7_alpha_Xpi0_Ygamma = 0
    n_c11_be7_alpha_Ygamma = 0
    n_c11_li6_p_alpha_Ygamma = 0
    n_c11_li6_p_alpha_Xpi0_Ygamma = 0
    n_c11_li5_d_alpha = 0
    n_c11_li5_d_alpha_Xpi0 = 0
    n_c11_li5_n_p_alpha_Xpi0 = 0
    n_c11_li5_n_p_alpha = 0
    n_c11_li5_d_alpha_Ygamma = 0
    n_c11_he3_2alpha = 0
    # C10 deex:
    number_c10_iso = 0
    n_c10_noDeex = 0
    n_c10_c10_Xpi0 = 0
    n_c10_c10_Ygamma = 0
    n_c10_c10_Xpi0_Ygamma = 0
    n_c10_c9_n_Ygamma = 0
    n_c10_c9_n_Xpi0_Ygamma = 0
    n_c10_b9_p_Ygamma = 0
    n_c10_b9_p_Xpi0_Ygamma = 0
    n_c10_b8_n_p_Ygamma = 0
    n_c10_b8_n_p_Xpi0_Ygamma = 0
    n_c10_b8_d_Xpi0_Ygamma = 0
    n_c10_b8_d_Ygamma = 0
    n_c10_be8_2p_Ygamma = 0
    n_c10_be8_2p_Xpi0_Ygamma = 0
    n_c10_be7_p_d_Ygamma = 0
    n_c10_be7_n_2p_Ygamma = 0
    n_c10_be7_p_d_Xpi0_Ygamma = 0
    n_c10_be7_n_2p_Xpi0_Ygamma = 0
    n_c10_be7_he3_Xpi0_Ygamma = 0
    n_c10_be7_he3_Ygamma = 0
    n_c10_li6_2p_d_Xpi0_Ygamma = 0
    n_c10_li6_p_he3_Ygamma = 0
    n_c10_li6_2p_d_Ygamma = 0
    n_c10_li6_n_3p_Ygamma = 0
    n_c10_li6_p_he3_Xpi0_Ygamma = 0
    n_c10_li6_n_3p_Xpi0_Ygamma = 0
    n_c10_li5_n_2p_d = 0
    n_c10_he3_p_d_alpha = 0
    n_c10_b7_2n_p_Xpi0 = 0
    n_c10_be6_n_p_d_Xpi0 = 0
    n_c10_li5_p_2d = 0
    n_c10_be6_p_t = 0
    n_c10_he3_n_2p_alpha = 0
    n_c10_li5_n_p_he3 = 0
    n_c10_b7_2n_p = 0
    n_c10_be6_n_p_d = 0
    n_c10_li5_d_he3 = 0
    n_c10_he3_p_d_alpha_Xpi0 = 0
    n_c10_li5_p_alpha_Xpi0 = 0
    n_c10_li4_n_p_alpha_Xpi0 = 0
    n_c10_rest = 0
    # B11 deexcitation:
    number_b11_iso = 0
    n_b11_noDeex = 0
    n_b11_b11_Ygamma = 0
    n_b11_b11_Xpi0 = 0
    n_b11_b11_Xpi0_Ygamma = 0
    n_b11_b10_n_Ygamma = 0
    n_b11_b10_n_Xpi0_Ygamma = 0
    n_b11_b9_2n_Ygamma = 0
    n_b11_b9_2n_Xpi0_Ygamma = 0
    n_b11_be10_p_Ygamma = 0
    n_b11_be10_p_Xpi0_Ygamma = 0
    n_b11_be9_n_p_Ygamma = 0
    n_b11_be9_d_Ygamma = 0
    n_b11_be9_d_Xpi0_Ygamma = 0
    n_b11_be9_n_p_Xpi0_Ygamma = 0
    n_b11_be8_n_d_Xpi0_Ygamma = 0
    n_b11_be8_n_d_Ygamma = 0
    n_b11_be8_t_Xpi0_Ygamma = 0
    n_b11_be8_t_Ygamma = 0
    n_b11_be8_2n_p_Ygamma = 0
    n_b11_Li7_alpha_Ygamma = 0
    n_b11_Li7_alpha_Xpi0_Ygamma = 0
    n_b11_Li6_n_alpha_Xpi0_Ygamma = 0
    n_b11_Li6_n_alpha_Ygamma = 0
    n_b11_he5_d_alpha = 0
    n_b11_he5_d_alpha_Xpi0 = 0
    n_b11_he5_n_p_alpha = 0
    n_b11_he6_p_alpha = 0
    n_b11_he5_n_p_alpha_Xpi0 = 0
    n_b11_t_2alpha = 0
    n_b11_he6_p_alpha_Xpi0 = 0
    n_b11_t_2alpha_Xpi0 = 0
    # B10 deexcitation:
    number_b10_iso = 0
    n_b10_noDeex = 0
    n_b10_b10_Xpi0 = 0
    n_b10_b10_Ygamma = 0
    n_b10_b10_Xpi0_Ygamma = 0
    n_b10_b9_n_Ygamma = 0
    n_b10_b9_n_Xpi0_Ygamma = 0
    n_b10_be9_p_Xpi0_Ygamma = 0
    n_b10_be9_p_Ygamma = 0
    n_b10_be8_n_p_Ygamma = 0
    n_b10_be8_d_Xpi0_Ygamma = 0
    n_b10_be8_d_Ygamma = 0
    n_b10_be8_n_p_Xpi0_Ygamma = 0
    n_b10_be7_t_Xpi0_Ygamma = 0
    n_b10_be7_t_Ygamma = 0
    n_b10_be7_n_d_Ygamma = 0
    n_b10_li7_he3_Ygamma = 0
    n_b10_li7_he3_Xpi0_Ygamma = 0
    n_b10_li7_p_d_Ygamma = 0
    n_b10_li7_p_d_Xpi0_Ygamma = 0
    n_b10_li6_alpha_Ygamma = 0
    n_b10_li6_alpha_Xpi0_Ygamma = 0
    n_b10_li6_p_t_Xpi0_Ygamma = 0
    n_b10_li6_p_t_Ygamma = 0
    n_b10_li6_n_he3_Ygamma = 0
    n_b10_he5_p_alpha = 0
    n_b10_he5_p_alpha_Xpi0 = 0
    n_b10_li5_n_alpha = 0
    n_b10_li5_n_alpha_Xpi0 = 0
    n_b10_he5_p_alpha_Xpi0_Ygamma = 0
    n_b10_he5_p_alpha_Ygamma = 0
    # Be10 deexcitation:
    number_be10_iso = 0
    n_be10_noDeex = 0
    n_be10_be10_Xpi0 = 0
    n_be10_be10_Ygamma = 0
    n_be10_be10_Xpi0_Ygamma = 0
    n_be10_be9_n_Xpi0_Ygamma = 0
    n_be10_be9_n_Ygamma = 0
    n_be10_be8_2n_Xpi0_Ygamma = 0
    n_be10_be8_2n_Ygamma = 0
    n_be10_li9_p_Ygamma = 0
    n_be10_li9_p_Xpi0_Ygamma = 0
    n_be10_li8_n_p_Ygamma = 0
    n_be10_li8_n_p_Xpi0_Ygamma = 0
    n_be10_li8_d_Xpi0_Ygamma = 0
    n_be10_li8_d_Ygamma = 0
    n_be10_li7_n_d_Ygamma = 0
    n_be10_li7_2n_p_Ygamma = 0
    n_be10_li7_2n_p_Xpi0_Ygamma = 0
    n_be10_li7_t_Ygamma = 0
    n_be10_li7_n_d_Xpi0_Ygamma = 0
    n_be10_li6_n_t_Ygamma = 0
    n_be10_li6_2n_d_Xpi0_Ygamma = 0
    n_be10_li6_3n_p_Ygamma = 0
    n_be10_li6_n_t_Xpi0_Ygamma = 0
    n_be10_li6_2n_d_Ygamma = 0
    n_be10_li6_3n_p_Xpi0_Ygamma = 0
    n_be10_h4_n_p_alpha = 0
    n_be10_he7_p_d = 0
    n_be10_he5_n_p_t_Xpi0 = 0
    n_be10_he5_2n_p_d = 0
    n_be10_he5_n_2d = 0
    n_be10_he7_he3 = 0
    n_be10_h4_n_p_alpha_Xpi0 = 0
    n_be10_he5_n_p_t = 0
    n_be10_alpha_n_d_t = 0
    n_be10_he6_n_p_d = 0
    n_be10_he6_2d = 0
    n_be10_alpha_n_d_t_Xpi0 = 0
    n_be10_rest = 0

    # loop over events:
    for event in range(N_entries_generator):
        # preallocate array, where PDG of final particles are stored:
        pdg_array = np.array([])

        # get the current event from generator tree:
        rtree_generator.GetEntry(event)

        # get channel ID:
        channelID = int(rtree_generator.GetBranch('t_channelID').GetLeaf('t_channelID').GetValue())
        # check channelID:
        if channelID == 2 or channelID == 3:
            # go to next event:
            continue
        # get number of proton, neutron, piminus, piplus from channelID (in channelID no gamma's are stored.
        # K_plus and K_minus are stored in channelID, but low statistics):
        n_p, n_n, n_piminus, n_piplus = NC_background_functions.get_number_of_particles_of_channelid(channelID)

        # get isopdg:
        isopdg_generator = int(rtree_generator.GetBranch('t_isoPdg').GetLeaf('t_isoPdg').GetValue())

        # get Npars:
        Npars_generator = int(rtree_generator.GetBranch('t_Npars').GetLeaf('t_Npars').GetValue())
        # loop over Npars to get pdg:
        for index in range(Npars_generator):
            pdg = int(rtree_generator.GetBranch('t_pdg').GetLeaf('t_pdg').GetValue(index))
            pdg_array = np.append(pdg_array, pdg)

        # get number of particles from pdg_array:
        (N_n, N_p, N_piminus, N_piplus, N_pi0, N_gamma, N_deuteron, N_triton, N_he3, N_alpha, N_li6, N_li7, N_li8, N_li9,
         N_be7, N_be8, N_be9, N_be10, N_b8, N_b9, N_b10, N_b11, N_c9, N_c10, N_c11, N_c12) \
            = get_number_of_particles(pdg_array)

        # calculate the number of particles corresponding only to deexcitation
        # ((after NC interaction and after deexcitation) - (after NC interaction and before deexcitation)):
        N_n = N_n - n_n
        N_p = N_p - n_p
        N_piminus = N_piminus - n_piminus
        N_piplus = N_piplus - n_piplus

        # if N_n + N_p + N_piminus + N_piplus + N_pi0 + N_deuteron + N_triton + N_he3 + N_alpha == 0 and N_gamma != 0:
        #     print("----------")
        #     print(isopdg_generator)
        #     print(N_n)
        #     print(N_p)
        #     print(N_piminus)
        #     print(N_piplus)
        #     print(N_pi0)
        #     print(N_gamma)
        #     print(N_deuteron)
        #     print(N_triton)
        #     print(N_he3)
        #     print(N_alpha)
        #     print(N_li6)
        #     print(N_li7)
        #     print(N_li8)
        #     print(N_li9)
        #     print(N_be7)
        #     print(N_be8)
        #     print(N_be9)
        #     print(N_be10)
        #     print(N_b8)
        #     print(N_b9)
        #     print(N_b10)
        #     print(N_b11)
        #     print(N_c9)
        #     print(N_c10)
        #     print(N_c11)

        # get deexcitation channels:
        if isopdg_generator == 1000060120:
            print("new channel: C12* -> C12 + ..., event = {0:d}".format(event))

        elif isopdg_generator == 1000060110:
            number_c11_iso += 1
            # C11 deexcitation:
            if N_c11 == 1:
                if (N_n == N_p == N_piminus == N_piplus == N_pi0 == N_gamma == N_deuteron == N_triton == N_he3 == N_alpha
                        == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c12 == 0):
                    # no deexcitation:
                    n_c11_noDeex += 1

                elif (N_n == N_p == N_piminus == N_piplus == N_pi0 == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c12 == 0 and N_gamma > 0):
                    # C11* -> C11 + 2gamma:
                    n_c11_c11_Ygamma += 1
                    N_evts_deex += 1

                elif (N_n == N_p == N_piminus == N_piplus == N_gamma == N_deuteron == N_triton == N_he3 == N_alpha ==
                        N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c12 == 0 and N_pi0 > 0):
                    # C11* -> C11 + X*pi0:
                    n_c11_c11_Xpi0 += 1
                    N_evts_deex += 1

                elif (N_n == N_p == N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c12 == 0 and N_pi0 > 0 and N_gamma > 0):
                    # C11* -> C11 + X*pi0 + Y*gamma:
                    n_c11_c11_Xpi0_Ygamma += 1
                    N_evts_deex += 1

                else:
                    if N_p < 0:
                        N_evts_incorrect += 1
                    else:
                        print("new channel: C11* -> C11 +..., event = {0:d}".format(event))

            elif N_c10 == 1:
                if (N_p == N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c11 == N_c12 == 0 and N_n == 1 and N_pi0 > 0 and N_gamma > 0):
                    # C11* -> C10 + n + X*pi0 + Y*gamma:
                    n_c11_c10_n_Xpi0_Ygamma += 1
                    N_evts_deex += 1

                elif (N_p == N_piminus == N_piplus == N_pi0 == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c11 == N_c12 == 0 and N_n == 1 and N_gamma > 0):
                    # C11* -> C10 + n + Y*gamma:
                    n_c11_c10_n_Ygamma += 1
                    N_evts_deex += 1

                else:
                    if N_p < 0:
                        N_evts_incorrect += 1
                    else:
                        print("new channel: C11* -> C10 +..., event = {0:d}".format(event))

            elif N_c9 == 1:
                print("new channel: C11* -> C9 + ..., event = {0:d}".format(event))

            elif N_b11 == 1:
                print("new channel: C11* -> B11 + ..., event = {0:d}".format(event))

            elif N_b10 == 1:
                if (N_n == N_piminus == N_piplus == N_pi0 == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == 1 and N_gamma > 0):
                    # C11* -> B10 + p + Y*gamma:
                    n_c11_b10_p_Ygamma += 1
                    N_evts_deex += 1

                elif (N_n == N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == 1 and N_pi0 > 0 and N_gamma > 0):
                    # C11* -> B10 + p + X*pi0 + Y*gamma:
                    n_c11_b10_p_Xpi0_Ygamma += 1
                    N_evts_deex += 1

                else:
                    if N_p < 0:
                        N_evts_incorrect += 1
                    else:
                        print("new channel: C11* -> B10 +..., event = {0:d}".format(event))

            elif N_b9 == 1:
                if (N_piminus == N_piplus == N_pi0 == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_p == 1 and N_gamma > 0):
                    # C11* -> B9 + n + p + Y*gamma:
                    n_c11_b9_n_p_Ygamma += 1
                    N_evts_deex += 1

                elif (N_n == N_p == N_piminus == N_piplus == N_pi0 == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_deuteron == 1 and N_gamma > 0):
                    # C11* -> B9 + d + Y*gamma:
                    n_c11_b9_d_Ygamma += 1
                    N_evts_deex += 1

                elif (N_n == N_p == N_piminus == N_piplus == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_deuteron == 1 and N_pi0 > 0 and N_gamma > 0):
                    # C11* -> B9 + d + X*pi0 + Y*gamma:
                    n_c11_b9_d_Xpi0_Ygamma += 1
                    N_evts_deex += 1

                elif (N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_p == 1 and N_pi0 > 0 and N_gamma > 0):
                    # C11* -> B9 + n + p + X*pi0 + Y*gamma:
                    n_c11_b9_n_p_Xpi0_Ygamma += 1
                    N_evts_deex += 1

                else:
                    if N_p < 0:
                        N_evts_incorrect += 1
                    else:
                        print("ERR: C11* -> B9 +..., event = {0:d}".format(event))

            elif N_b8 == 1:
                print("new channel: C11* -> B8 +..., event = {0:d}".format(event))

            elif N_be10 == 1:
                print("new channel: C11* -> Be10 +..., event = {0:d}".format(event))

            elif N_be9 == 1:
                if (N_n == N_piminus == N_piplus == N_pi0 == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == 2 and N_gamma > 0):
                    # C11* -> Be9 + 2*p + Y*gamma:
                    n_c11_be9_2p_Ygamma += 1
                    N_evts_deex += 1

                elif (N_n == N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == 2 and N_pi0 > 0 and N_gamma > 0):
                    # C11* -> Be9 + 2*p + X*pi0 + Y*gamma:
                    n_c11_be9_2p_Xpi0_Ygamma += 1
                    N_evts_deex += 1

                else:
                    if N_p < 0:
                        N_evts_incorrect += 1
                    else:
                        print("new channel: C11* -> Be9 +..., event = {0:d}".format(event))

            elif N_be8 == 1:
                if (N_n == N_piminus == N_piplus == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == N_deuteron == 1 and N_pi0 > 0 and N_gamma > 0):
                    # C11* -> Be8 + p + d + X*pi0 + Y*gamma:
                    n_c11_be8_p_d_Xpi0_Ygamma += 1
                    N_evts_deex += 1

                elif (N_n == N_piminus == N_piplus == N_pi0 == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == N_deuteron == 1 and N_gamma > 0):
                    # C11* -> Be8 + p + d + Y*gamma:
                    n_c11_be8_p_d_Ygamma += 1
                    N_evts_deex += 1

                elif (N_p == N_n == N_piminus == N_piplus == N_pi0 == N_deuteron == N_triton == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_he3 == 1 and N_gamma > 0):
                    # C11* -> Be8 + He3 + Y*gamma:
                    n_c11_be8_he3_Ygamma += 1
                    N_evts_deex += 1

                elif (N_p == N_n == N_piminus == N_piplus == N_deuteron == N_triton == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_he3 == 1 and N_pi0 > 0 and N_gamma > 0):
                    # C11* -> Be8 + He3 + X*pi0 Y*gamma:
                    n_c11_be8_he3_Xpi0_Ygamma += 1
                    N_evts_deex += 1

                elif (N_piminus == N_piplus == N_pi0 == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == 1 and N_p == 2 and N_gamma > 0):
                    # C11* -> Be8 + n + 2p + Y*gamma:
                    n_c11_be8_n_2p_Ygamma += 1
                    N_evts_deex += 1

                elif (N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == 1 and N_p == 2 and N_pi0 > 0 and N_gamma > 0):
                    # C11* -> Be8 + n + 2p + X*pi0 + Y*gamma:
                    n_c11_be8_n_2p_Xpi0_Ygamma += 1
                    N_evts_deex += 1

                else:
                    if N_p < 0:
                        N_evts_incorrect += 1
                    else:
                        print("new channel: C11* -> Be8 +..., event = {0:d}".format(event))

            elif N_be7 == 1:
                if (N_n == N_p == N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_alpha == 1 and N_pi0 > 0 and N_gamma > 0):
                    # C11* -> Be7 + alpha + X*pi0 + Y*gamma:
                    n_c11_be7_alpha_Xpi0_Ygamma += 1
                    N_evts_deex += 1

                elif (N_n == N_p == N_piminus == N_piplus == N_pi0 == N_deuteron == N_triton == N_he3 == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_alpha == 1 and N_gamma > 0):
                    # C11* -> Be7 + alpha + Y*gamma:
                    n_c11_be7_alpha_Ygamma += 1
                    N_evts_deex += 1

                else:
                    if N_p < 0:
                        N_evts_incorrect += 1
                    else:
                        print("new channel: C11* -> Be7 +..., event = {0:d}".format(event))

            elif N_li9 == 1:
                print("new channel: C11* -> Li9 +..., event = {0:d}".format(event))

            elif N_li8 == 1:
                print("new channel: C11* -> Li8 +..., event = {0:d}".format(event))

            elif N_li7 == 1:
                print("new channel: C11* -> Li7 +..., event = {0:d}".format(event))

            elif N_li6 == 1:
                if (N_n == N_piminus == N_piplus == N_pi0 == N_deuteron == N_triton == N_he3 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == N_alpha == 1 and N_gamma > 0):
                    # C11* -> Li6 + p + alpha + Y*gamma:
                    n_c11_li6_p_alpha_Ygamma += 1
                    N_evts_deex += 1

                elif (N_n == N_piminus == N_piplus == N_deuteron == N_triton == N_he3 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == N_alpha == 1 and N_pi0 > 0 and N_gamma > 0):
                    # C11* -> Li6 + p + alpha + X*pi0 + Y*gamma:
                    n_c11_li6_p_alpha_Xpi0_Ygamma += 1
                    N_evts_deex += 1

                else:
                    if N_p < 0:
                        N_evts_incorrect += 1
                    else:
                        print("new channel: C11* -> Li6 +..., event = {0:d}".format(event))

            else:
                if (N_n == N_p == N_piminus == N_piplus == N_pi0 == N_gamma == N_triton == N_he3 == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_deuteron == N_alpha == 1):
                    # C11* -> Li5 + d + alpha:
                    n_c11_li5_d_alpha += 1
                    N_evts_deex += 1

                elif (N_n == N_p == N_piminus == N_piplus == N_gamma == N_triton == N_he3 == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_deuteron == N_alpha == 1 and N_pi0 > 0):
                    # C11* -> Li5 + d + alpha + X*pi0:
                    n_c11_li5_d_alpha_Xpi0 += 1
                    N_evts_deex += 1

                elif (N_piminus == N_piplus == N_gamma == N_deuteron == N_triton == N_he3 == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_p == N_alpha == 1 and N_pi0 > 0):
                    # C11* -> Li5 + n + p + alpha + X*pi0:
                    n_c11_li5_n_p_alpha_Xpi0 += 1
                    N_evts_deex += 1

                elif (N_piminus == N_piplus == N_pi0 == N_gamma == N_deuteron == N_triton == N_he3 == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_p == N_alpha == 1):
                    # C11* -> Li5 + n + p + alpha:
                    n_c11_li5_n_p_alpha += 1
                    N_evts_deex += 1

                elif (N_n == N_p == N_piminus == N_piplus == N_pi0 == N_triton == N_he3 == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_deuteron == N_alpha == 1 and N_gamma > 0):
                    # C11* -> Li5 + d + alpha + Y*gamma:
                    n_c11_li5_d_alpha_Ygamma += 1
                    N_evts_deex += 1

                elif (N_n == N_p == N_piminus == N_piplus == N_pi0 == N_gamma == N_deuteron == N_triton == N_he3 == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_alpha == 2):
                    # C11* -> He3 + 2*alpha:
                    n_c11_he3_2alpha += 1
                    N_evts_deex += 1

                else:
                    print("new channel: C11* -> ..., event = {0:d}".format(event))

        elif isopdg_generator == 1000060100:
            number_c10_iso += 1
            # C10 deexcitation:
            if N_c10 == 1:
                if (N_n == N_p == N_piminus == N_piplus == N_pi0 == N_gamma == N_deuteron == N_triton == N_he3 == N_alpha
                        == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c11 == N_c12 == 0):
                    # no deexcitation:
                    n_c10_noDeex += 1

                elif (N_n == N_p == N_piminus == N_piplus == N_gamma == N_deuteron == N_triton == N_he3 == N_alpha
                      == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c11 == N_c12 == 0 and N_pi0 > 0):
                    # C10* -> C10 + X*pi0:
                    n_c10_c10_Xpi0 += 1
                    N_evts_deex += 1

                elif (N_n == N_p == N_piminus == N_piplus == N_pi0 == N_deuteron == N_triton == N_he3 == N_alpha
                      == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c11 == N_c12 == 0 and N_gamma > 0):
                    # C10* -> C10 + Y*gamma:
                    n_c10_c10_Ygamma += 1
                    N_evts_deex += 1

                elif (N_n == N_p == N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c11 == N_c12 == 0 and N_pi0 > 0 and N_gamma > 0):
                    # C10* -> C10 + X*pi0 + Y*gamma:
                    n_c10_c10_Xpi0_Ygamma += 1
                    N_evts_deex += 1

                else:
                    if N_p < 0:
                        N_evts_incorrect += 1
                    else:
                        print("new channel: C10* -> C10 + ..., event = {0:d}".format(event))

            elif N_c9 == 1:
                if (N_p == N_piminus == N_piplus == N_pi0 == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c10 == N_c11 == N_c12 == 0 and N_n == 1 and N_gamma > 0):
                    # C10* -> C9 + n + Y*gamma:
                    n_c10_c9_n_Ygamma += 1
                    N_evts_deex += 1

                elif (N_p == N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c10 == N_c11 == N_c12 == 0 and N_n == 1 and N_pi0 > 0 and N_gamma > 0):
                    # C10* -> C9 + n + X*pi0 + Y*gamma:
                    n_c10_c9_n_Xpi0_Ygamma += 1
                    N_evts_deex += 1

                else:
                    if N_p < 0:
                        N_evts_incorrect += 1
                    else:
                        print("new channel: C10* -> C9 + ..., event = {0:d}".format(event))

            elif N_b11 == 1:
                print("new channel: C10* -> B11 +..., event = {0:d}".format(event))

            elif N_b10 == 1:
                print("new channel: C10* -> B10 +..., event = {0:d}".format(event))

            elif N_b9 == 1:
                if (N_n == N_piminus == N_piplus == N_pi0 == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == 1 and N_gamma > 0):
                    # C10* -> B9 + p + Y*gamma:
                    n_c10_b9_p_Ygamma += 1
                    N_evts_deex += 1

                elif (N_n == N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == 1 and N_pi0 > 0 and N_gamma > 0):
                    # C10* -> B9 + p + X*pi0 + Y*gamma:
                    n_c10_b9_p_Xpi0_Ygamma += 1
                    N_evts_deex += 1

                else:
                    if N_p < 0:
                        N_evts_incorrect += 1
                    else:
                        print("new channel: C10* -> B9 + ..., event = {0:d}".format(event))

            elif N_b8 == 1:
                if (N_piminus == N_piplus == N_pi0 == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_p == 1 and N_gamma > 0):
                    # C10* -> B8 + n + p + Y*gamma:
                    n_c10_b8_n_p_Ygamma += 1
                    N_evts_deex += 1

                elif (N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_p == 1 and N_pi0 > 0 and N_gamma > 0):
                    # C10* -> B8 + n + p + X*pi0 + Y*gamma:
                    n_c10_b8_n_p_Xpi0_Ygamma += 1
                    N_evts_deex += 1

                elif (N_n == N_p == N_piminus == N_piplus == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_deuteron == 1 and N_pi0 > 0 and N_gamma > 0):
                    # C10* -> B8 + d + X*pi0 + Y*gamma:
                    n_c10_b8_d_Xpi0_Ygamma += 1
                    N_evts_deex += 1

                elif (N_n == N_p == N_piminus == N_piplus == N_pi0 == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_deuteron == 1 and N_gamma > 0):
                    # C10* -> B8 + d + Y*gamma:
                    n_c10_b8_d_Ygamma += 1
                    N_evts_deex += 1

                else:
                    if N_p < 0:
                        N_evts_incorrect += 1
                    else:
                        print("new channel: C10* -> B8 + ..., event = {0:d}".format(event))

            elif N_be10 == 1:
                print("new channel: C10* -> Be10 +..., event = {0:d}".format(event))

            elif N_be9 == 1:
                print("new channel: C10* -> Be9 +..., event = {0:d}".format(event))

            elif N_be8 == 1:
                if (N_n == N_piminus == N_piplus == N_pi0 == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == 2 and N_gamma > 0):
                    # C10* -> Be8 + 2p + Y*gamma:
                    n_c10_be8_2p_Ygamma += 1
                    N_evts_deex += 1

                elif (N_n == N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == 2 and N_pi0 > 0 and N_gamma > 0):
                    # C10* -> Be8 + 2p + X*pi0 + Y*gamma:
                    n_c10_be8_2p_Xpi0_Ygamma += 1
                    N_evts_deex += 1

                else:
                    if N_p < 0:
                        N_evts_incorrect += 1
                    else:
                        print("new channel: C10* -> Be8 + ..., event = {0:d}".format(event))

            elif N_be7 == 1:
                if (N_n == N_piminus == N_piplus == N_pi0 == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == N_deuteron == 1 and N_gamma > 0):
                    # C10* -> Be7 + p +d + Y*gamma:
                    n_c10_be7_p_d_Ygamma += 1
                    N_evts_deex += 1

                elif (N_piminus == N_piplus == N_pi0 == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == 1 and N_p == 2 and N_gamma > 0):
                    # C10* -> Be7 + n + 2p + Y*gamma:
                    n_c10_be7_n_2p_Ygamma += 1
                    N_evts_deex += 1

                elif (N_n == N_piminus == N_piplus == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == N_deuteron == 1 and N_pi0 > 0 and N_gamma > 0):
                    # C10* -> Be7 + p + d + X*pi0 + Y*gamma:
                    n_c10_be7_p_d_Xpi0_Ygamma += 1
                    N_evts_deex += 1

                elif (N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == 1 and N_p == 2 and N_pi0 > 0 and N_gamma > 0):
                    # C10* -> Be7 + n + 2p + Y*gamma:
                    n_c10_be7_n_2p_Xpi0_Ygamma += 1
                    N_evts_deex += 1

                elif (N_n == N_p == N_piminus == N_piplus == N_deuteron == N_triton == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_he3 == 1 and N_pi0 > 0 and N_gamma > 0):
                    # C10* -> Be7 + He3 + X*pi0 + Y*gamma:
                    n_c10_be7_he3_Xpi0_Ygamma += 1
                    N_evts_deex += 1

                elif (N_n == N_p == N_piminus == N_piplus == N_pi0 == N_deuteron == N_triton == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_he3 == 1 and N_gamma > 0):
                    # C10* -> Be7 + He3 + Y*gamma:
                    n_c10_be7_he3_Ygamma += 1
                    N_evts_deex += 1

                else:
                    if N_p < 0:
                        N_evts_incorrect += 1
                    else:
                        print("new channel: C10* -> Be7 + ..., event = {0:d}".format(event))

            elif N_li9 == 1:
                print("new channel: C10* -> Li9 +..., event = {0:d}".format(event))

            elif N_li8 == 1:
                print("new channel: C10* -> Li8 +..., event = {0:d}".format(event))

            elif N_li7 == 1:
                print("new channel: C10* -> Li7 +..., event = {0:d}".format(event))

            elif N_li6 == 1:
                if (N_n == N_piminus == N_piplus == N_triton == N_he3 == N_alpha ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == 2 and N_deuteron == 1 and N_pi0 > 0 and
                        N_gamma > 0):
                    # C10* -> Li6 + 2p + d + X*pi0 + Y*gamma:
                    n_c10_li6_2p_d_Xpi0_Ygamma += 1
                    N_evts_deex += 1

                elif (N_n == N_piminus == N_piplus == N_pi0 == N_triton == N_deuteron == N_alpha ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == N_he3 == 1 and N_gamma > 0):
                    # C10* -> Li6 + p + He3 + Y*gamma:
                    n_c10_li6_p_he3_Ygamma += 1
                    N_evts_deex += 1

                elif (N_n == N_piminus == N_piplus == N_pi0 == N_triton == N_he3 == N_alpha ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == 2 and N_deuteron == 1 and N_gamma > 0):
                    # C10* -> Li6 + 2p + d + Y*gamma:
                    n_c10_li6_2p_d_Ygamma += 1
                    N_evts_deex += 1

                elif (N_piminus == N_piplus == N_pi0 == N_deuteron == N_triton == N_he3 == N_alpha ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == 1 and N_p == 3 and N_gamma > 0):
                    # C10* -> Li6 + n + 3p + Y*gamma:
                    n_c10_li6_n_3p_Ygamma += 1
                    N_evts_deex += 1

                elif (N_n == N_piminus == N_piplus == N_deuteron == N_triton == N_alpha ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == N_he3 == 1 and N_pi0 > 0 and N_gamma > 0):
                    # C10* -> Li6 + p + he3 + X*pi0 + Y*gamma:
                    n_c10_li6_p_he3_Xpi0_Ygamma += 1
                    N_evts_deex += 1

                elif (N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == 1 and N_p == 3 and N_pi0 > 0 and N_gamma > 0):
                    # C10* -> Li6 + n + 3p + X*pi0 + Y*gamma:
                    n_c10_li6_n_3p_Xpi0_Ygamma += 1
                    N_evts_deex += 1

                else:
                    if N_p < 0:
                        N_evts_incorrect += 1
                    else:
                        print("new channel: C10* -> Li6 + ..., event = {0:d}".format(event))

            else:
                if (N_piminus == N_piplus == N_pi0 == N_gamma == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == 1 and N_p == 2 and N_deuteron == 1):
                    # C10* -> Li5 + n + 2p + d:
                    n_c10_li5_n_2p_d += 1
                    N_evts_deex += 1

                elif (N_n == N_piminus == N_piplus == N_pi0 == N_gamma == N_triton == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == N_deuteron == N_he3 == 1) or \
                        (N_n == N_piminus == N_piplus == N_pi0 == N_gamma == N_triton == N_he3 == N_li6 ==
                         N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                         N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == N_deuteron == N_alpha == 1):
                    # C10* -> He3 + p + d + alpha:
                    n_c10_he3_p_d_alpha += 1
                    N_evts_deex += 1

                elif (N_piminus == N_piplus == N_gamma == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == 2 and N_p == 1 and N_pi0 > 0):
                    # C10* -> B7 + 2n + p + X*pi0:
                    n_c10_b7_2n_p_Xpi0 += 1
                    N_evts_deex += 1

                elif (N_piminus == N_piplus == N_gamma == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_p == N_deuteron == 1 and N_pi0 > 0):
                    # C10* -> Be6 + n + p + d + X*pi0:
                    n_c10_be6_n_p_d_Xpi0 += 1
                    N_evts_deex += 1

                elif (N_n == N_piminus == N_piplus == N_pi0 == N_gamma == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == 1 and N_deuteron == 2):
                    # C10* -> Li5 + p + 2d:
                    n_c10_li5_p_2d += 1
                    N_evts_deex += 1

                elif (N_n == N_piminus == N_piplus == N_pi0 == N_gamma == N_deuteron == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == 1 and N_triton == 1):
                    # C10* -> Be6 + p + t:
                    n_c10_be6_p_t += 1
                    N_evts_deex += 1

                elif (N_piminus == N_piplus == N_pi0 == N_gamma == N_deuteron == N_triton == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == 1 and N_p == 2 and N_he3 == 1):
                    # C10* -> he3 + n + 2p + alpha:
                    n_c10_he3_n_2p_alpha += 1
                    N_evts_deex += 1

                elif (N_piminus == N_piplus == N_pi0 == N_gamma == N_deuteron == N_triton == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_p == N_he3 == 1):
                    # C10* -> li5 + n + p + he3:
                    n_c10_li5_n_p_he3 += 1
                    N_evts_deex += 1

                elif (N_piminus == N_piplus == N_pi0 == N_gamma == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == 2 and N_p == 1):
                    # C10* -> B7 + 2n + p:
                    n_c10_b7_2n_p += 1
                    N_evts_deex += 1

                elif (N_piminus == N_piplus == N_pi0 == N_gamma == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_p == N_deuteron == 1):
                    # C10* -> Be6 + n + p + d:
                    n_c10_be6_n_p_d += 1
                    N_evts_deex += 1

                elif (N_n == N_p == N_piminus == N_piplus == N_pi0 == N_gamma == N_triton == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_deuteron == N_he3 == 1):
                    # C10* -> Li5 + d + he3:
                    n_c10_li5_d_he3 += 1
                    N_evts_deex += 1

                elif (N_n == N_piminus == N_piplus == N_gamma == N_triton == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == N_deuteron == N_he3 == 1 and N_pi0 > 0):
                    # C10* -> He3 + p + d + alpha + X*pi0:
                    n_c10_he3_p_d_alpha_Xpi0 += 1
                    N_evts_deex += 1

                elif (N_n == N_piminus == N_piplus == N_gamma == N_deuteron == N_triton == N_he3 == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == N_alpha == 1 and N_pi0 > 0):
                    # C10* -> li5 + p + alpha + X*pi0:
                    n_c10_li5_p_alpha_Xpi0 += 1
                    N_evts_deex += 1

                elif (N_piminus == N_piplus == N_gamma == N_deuteron == N_triton == N_he3 == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_p == N_alpha == 1 and N_pi0 > 0):
                    # C10* -> li4 + n + p + alpha + X*pi0:
                    n_c10_li4_n_p_alpha_Xpi0 += 1
                    N_evts_deex += 1

                else:
                    if N_p < 0:
                        N_evts_incorrect += 1

                    else:
                        # more deexcitation channels C10* -> ...
                        # print("new channel: C10* -> ..., event = {0:d}".format(event))
                        n_c10_rest += 1
                        N_evts_deex += 1

        elif isopdg_generator == 1000050110:
            number_b11_iso += 1
            # B11 deexcitation:
            if N_c11 == 1:
                print("new channel: B11* -> C11 + ..., event = {0:d}".format(event))

            elif N_c10 == 1:
                print("new channel: B11* -> C10 + ..., event = {0:d}".format(event))

            elif N_c9 == 1:
                print("new channel: B11* -> C9 + ..., event = {0:d}".format(event))

            elif N_b11 == 1:
                if (N_n == N_p == N_piminus == N_piplus == N_pi0 == N_gamma == N_deuteron == N_triton == N_he3 == N_alpha
                        == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0):
                    # no deexcitation:
                    n_b11_noDeex += 1

                elif (N_n == N_p == N_piminus == N_piplus == N_pi0 == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_gamma > 0):
                    # B11* -> B11 + Y*gamma:
                    n_b11_b11_Ygamma += 1
                    N_evts_deex += 1

                elif (N_n == N_p == N_piminus == N_piplus == N_gamma == N_deuteron == N_triton == N_he3 == N_alpha
                      == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_pi0 > 0):
                    # B11* -> B11 + X*pi0:
                    n_b11_b11_Xpi0 += 1
                    N_evts_deex += 1

                elif (N_n == N_p == N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_pi0 > 0 and N_gamma > 0):
                    # B11* -> B11 + X*pi0 + Y*gamma:
                    n_b11_b11_Xpi0_Ygamma += 1
                    N_evts_deex += 1

                else:
                    if N_p < 0:
                        N_evts_incorrect += 1
                    else:
                        print("new channel: B11* -> B11 + ..., event = {0:d}".format(event))

            elif N_b10 == 1:
                if (N_p == N_piminus == N_piplus == N_pi0 == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == 1 and N_gamma > 0):
                    # B11* -> B10 + n + 2gamma:
                    n_b11_b10_n_Ygamma += 1
                    N_evts_deex += 1

                elif (N_p == N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == 1 and N_pi0 > 0 and N_gamma > 0):
                    # B11* -> B10 + n + X*pi0 + Y*gamma:
                    n_b11_b10_n_Xpi0_Ygamma += 1
                    N_evts_deex += 1

                else:
                    if N_p < 0:
                        N_evts_incorrect += 1
                    else:
                        print("new channel: B11* -> B10 + ..., event = {0:d}".format(event))

            elif N_b9 == 1:
                if (N_p == N_piminus == N_piplus == N_pi0 == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == 2 and N_gamma > 0):
                    # B11* -> B8 + 2n + 2gamma:
                    n_b11_b9_2n_Ygamma += 1
                    N_evts_deex += 1

                elif (N_p == N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == 2 and N_pi0 > 0 and N_gamma > 0):
                    # B11* -> B8 + 2n + X*pi0 + Y*gamma:
                    n_b11_b9_2n_Xpi0_Ygamma += 1
                    N_evts_deex += 1

                else:
                    if N_p < 0:
                        N_evts_incorrect += 1
                    else:
                        print("new channel: B11* -> B9 + ..., event = {0:d}".format(event))

            elif N_b8 == 1:
                print("new channel: B11* -> B8 + ..., event = {0:d}".format(event))

            elif N_be10 == 1:
                if (N_n == N_piminus == N_piplus == N_pi0 == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == 1 and N_gamma > 0):
                    # B11* -> Be10 + p + Y*gamma:
                    n_b11_be10_p_Ygamma += 1
                    N_evts_deex += 1

                elif (N_n == N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == 1 and N_pi0 > 0 and N_gamma > 0):
                    # B11* -> Be10 + p + X*pi0 + Y*gamma:
                    n_b11_be10_p_Xpi0_Ygamma += 1
                    N_evts_deex += 1

                else:
                    if N_p < 0:
                        N_evts_incorrect += 1
                    else:
                        print("new channel: B11* -> Be10 + ..., event = {0:d}".format(event))

            elif N_be9 == 1:
                if (N_piminus == N_piplus == N_pi0 == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_p == 1 and N_gamma > 0):
                    # B11* -> Be9 + n + p + Y*gamma:
                    n_b11_be9_n_p_Ygamma += 1
                    N_evts_deex += 1

                elif (N_n == N_p == N_piminus == N_piplus == N_pi0 == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_deuteron == 1 and N_gamma > 0):
                    # B11* -> Be9 + d + Y*gamma:
                    n_b11_be9_d_Ygamma += 1
                    N_evts_deex += 1

                elif (N_n == N_p == N_piminus == N_piplus == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_deuteron == 1 and N_pi0 > 0 and N_gamma > 0):
                    # B11* -> Be9 + d + X*pi0 + Y*gamma:
                    n_b11_be9_d_Xpi0_Ygamma += 1
                    N_evts_deex += 1

                elif (N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_p == 1 and N_pi0 > 0 and N_gamma > 0):
                    # B11* -> Be9 + n + p + X*pi0 + Y*gamma:
                    n_b11_be9_n_p_Xpi0_Ygamma += 1
                    N_evts_deex += 1

                else:
                    if N_p < 0:
                        N_evts_incorrect += 1
                    else:
                        print("new channel: B11* -> Be9 + ..., event = {0:d}".format(event))

            elif N_be8 == 1:
                if (N_p == N_piminus == N_piplus == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_deuteron == 1 and N_pi0 > 0 and N_gamma > 0):
                    # B11* -> Be8 + n + d + X*pi0 + Y*gamma:
                    n_b11_be8_n_d_Xpi0_Ygamma += 1
                    N_evts_deex += 1

                elif (N_p == N_piminus == N_piplus == N_pi0 == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_deuteron == 1 and N_gamma > 0):
                    # B11* -> Be8 + n + d + Y*gamma:
                    n_b11_be8_n_d_Ygamma += 1
                    N_evts_deex += 1

                elif (N_n == N_p == N_piminus == N_piplus == N_deuteron == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_triton == 1 and N_pi0 > 0 and N_gamma > 0):
                    # B11* -> Be8 + t + X*pi0 + Y*gamma:
                    n_b11_be8_t_Xpi0_Ygamma += 1
                    N_evts_deex += 1

                elif (N_n == N_p == N_piminus == N_piplus == N_pi0 == N_deuteron == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_triton == 1 and N_gamma > 0):
                    # B11* -> Be8 + t + Y*gamma:
                    n_b11_be8_t_Ygamma += 1
                    N_evts_deex += 1

                elif (N_piminus == N_piplus == N_pi0 == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == 2 and N_p == 1 and N_gamma > 0):
                    # B11* -> Be8 + 2n + p + Y*gamma:
                    n_b11_be8_2n_p_Ygamma += 1
                    N_evts_deex += 1

                else:
                    if N_p < 0:
                        N_evts_incorrect += 1
                    else:
                        print("new channel: B11* -> Be8 + ..., event = {0:d}".format(event))

            elif N_be7 == 1:
                print("new channel: B11* -> Be7 + ..., event = {0:d}".format(event))

            elif N_li9 == 1:
                print("new channel: B11* -> Li9 + ..., event = {0:d}".format(event))

            elif N_li8 == 1:
                print("new channel: B11* -> Li8 + ..., event = {0:d}".format(event))

            elif N_li7 == 1:
                if (N_n == N_p == N_piminus == N_piplus == N_pi0 == N_deuteron == N_triton == N_he3 == N_li6 ==
                        N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_alpha == 1 and N_gamma > 0):
                    # B11* -> Li7 + alpha + Y*gamma:
                    n_b11_Li7_alpha_Ygamma += 1
                    N_evts_deex += 1

                elif (N_n == N_p == N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_li6 ==
                        N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_alpha == 1 and N_pi0 > 0 and N_gamma > 0):
                    # B11* -> Li7 + alpha + X*pi0 + Y*gamma:
                    n_b11_Li7_alpha_Xpi0_Ygamma += 1
                    N_evts_deex += 1

                else:
                    if N_p < 0:
                        N_evts_incorrect += 1
                    else:
                        print("new channel: B11* -> Li7 + ..., event = {0:d}".format(event))

            elif N_li6 == 1:
                if (N_p == N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_li7 ==
                        N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_alpha == 1 and N_pi0 > 0 and N_gamma > 0):
                    # B11* -> Li6 + n + alpha + X*pi0 + Y*gamma:
                    n_b11_Li6_n_alpha_Xpi0_Ygamma += 1
                    N_evts_deex += 1

                elif (N_p == N_piminus == N_piplus == N_pi0 == N_deuteron == N_triton == N_he3 == N_li7 ==
                        N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_alpha == 1 and N_gamma > 0):
                    # B11* -> Li6 + n + alpha + Y*gamma:
                    n_b11_Li6_n_alpha_Ygamma += 1
                    N_evts_deex += 1

                else:
                    if N_p < 0:
                        N_evts_incorrect += 1
                    else:
                        print("new channel: B11* -> Li6 + ..., event = {0:d}".format(event))

            else:
                if (N_n == N_p == N_piminus == N_piplus == N_pi0 == N_gamma == N_triton == N_he3 == N_li6 == N_li7 ==
                        N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_deuteron == N_alpha == 1):
                    # B11* -> He5 + d + alpha:
                    n_b11_he5_d_alpha += 1
                    N_evts_deex += 1

                elif (N_n == N_p == N_piminus == N_piplus == N_gamma == N_triton == N_he3 == N_li6 == N_li7 ==
                        N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_deuteron == N_alpha == 1 and N_pi0 > 0):
                    # B11* -> He5 + d + alpha + X*pi0:
                    n_b11_he5_d_alpha_Xpi0 += 1
                    N_evts_deex += 1

                elif (N_piminus == N_piplus == N_pi0 == N_gamma == N_deuteron == N_triton == N_he3 == N_li6 == N_li7 ==
                        N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_p == N_alpha == 1):
                    # B11* -> He5 + n + p + alpha:
                    n_b11_he5_n_p_alpha += 1
                    N_evts_deex += 1

                elif (N_n == N_piminus == N_piplus == N_pi0 == N_gamma == N_deuteron == N_triton == N_he3 == N_li6
                      == N_li7 ==
                        N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == N_alpha == 1):
                    # B11* -> He6 + p + alpha:
                    n_b11_he6_p_alpha += 1
                    N_evts_deex += 1

                elif (N_piminus == N_piplus == N_gamma == N_deuteron == N_triton == N_he3 == N_li6 == N_li7 ==
                        N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_p == N_alpha == 1 and N_pi0 > 0):
                    # B11* -> He5 + n + p + alpha + X*pi0:
                    n_b11_he5_n_p_alpha_Xpi0 += 1
                    N_evts_deex += 1

                elif (N_n == N_p == N_piminus == N_piplus == N_pi0 == N_gamma == N_deuteron == N_he3 == N_li6 == N_li7 ==
                        N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_triton == N_alpha == 1) or \
                        (N_n == N_p == N_piminus == N_piplus == N_pi0 == N_gamma == N_deuteron == N_triton == N_he3
                         == N_li6 == N_li7 ==
                         N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                         N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_alpha == 2):
                    # B11* -> t + 2*alpha:
                    n_b11_t_2alpha += 1
                    N_evts_deex += 1

                elif (N_n == N_p == N_piminus == N_piplus == N_gamma == N_deuteron == N_triton == N_he3 == N_li6 == N_li7 ==
                      N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                      N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_alpha == 2 and N_pi0 > 0):
                    # B11* -> t + 2*alpha + X*pi0:
                    n_b11_t_2alpha_Xpi0 += 1
                    N_evts_deex += 1

                elif (N_n == N_piminus == N_piplus == N_gamma == N_deuteron == N_triton == N_he3 == N_li6 == N_li7 ==
                        N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == N_alpha == 1 and N_pi0 > 0):
                    # B11* -> He6 + p + alpha X*pi0:
                    n_b11_he6_p_alpha_Xpi0 += 1
                    N_evts_deex += 1

                else:
                    if N_p < 0:
                        N_evts_incorrect += 1
                    else:
                        print("new channel: B11* -> ..., event = {0:d}".format(event))

        elif isopdg_generator == 1000050100:
            number_b10_iso += 1
            # B10 deexcitation:
            if N_c11 == 1:
                print("new channel: B10* -> C11 + ..., event = {0:d}".format(event))

            elif N_c10 == 1:
                print("new channel: B10* -> C10 + ..., event = {0:d}".format(event))

            elif N_c9 == 1:
                print("new channel: B10* -> C9 + ..., event = {0:d}".format(event))

            elif N_b11 == 1:
                print("new channel: B10* -> B11 + ..., event = {0:d}".format(event))

            elif N_b10 == 1:
                if (N_n == N_p == N_piminus == N_piplus == N_pi0 == N_gamma == N_deuteron == N_triton == N_he3 == N_alpha
                        == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0):
                    # no deexcitation:
                    n_b10_noDeex += 1

                elif (N_n == N_p == N_piminus == N_piplus == N_gamma == N_deuteron == N_triton == N_he3 == N_alpha
                        == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_pi0 > 0):
                    # B10* -> B10 + X*pi0:
                    n_b10_b10_Xpi0 += 1
                    N_evts_deex += 1

                elif (N_n == N_p == N_piminus == N_piplus == N_pi0 == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_gamma > 0):
                    # B10* -> B10 + Y*gamma:
                    n_b10_b10_Ygamma += 1
                    N_evts_deex += 1

                elif (N_n == N_p == N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_pi0 > 0 and N_gamma > 0):
                    # B10* -> B10 + X*pi0 + Y*gamma:
                    n_b10_b10_Xpi0_Ygamma += 1
                    N_evts_deex += 1

                else:
                    if N_p < 0:
                        N_evts_incorrect += 1
                    else:
                        print("new channel: B10* -> B10 + ..., event = {0:d}".format(event))

            elif N_b9 == 1:
                if (N_p == N_piminus == N_piplus == N_pi0 == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == 1 and N_gamma > 0):
                    # B10* -> B9 + n + Y*gamma:
                    n_b10_b9_n_Ygamma += 1
                    N_evts_deex += 1

                elif (N_p == N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == 1 and N_pi0 > 0 and N_gamma > 0):
                    # B10* -> B9 + n + X*pi0 + Y*gamma:
                    n_b10_b9_n_Xpi0_Ygamma += 1
                    N_evts_deex += 1

                else:
                    if N_p < 0:
                        N_evts_incorrect += 1
                    else:
                        print("new channel: B10* -> B9 + ..., event = {0:d}".format(event))

            elif N_b8 == 1:
                print("new channel: B10* -> B8 + ..., event = {0:d}".format(event))

            elif N_be10 == 1:
                print("new channel: B10* -> Be10 + ..., event = {0:d}".format(event))

            elif N_be9 == 1:
                if (N_n == N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == 1 and N_pi0 > 0 and N_gamma > 0):
                    # B10* -> Be9 + p + X*pi0 + Y*gamma:
                    n_b10_be9_p_Xpi0_Ygamma += 1
                    N_evts_deex += 1

                elif (N_n == N_piminus == N_piplus == N_pi0 == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == 1 and N_gamma > 0):
                    # B10* -> Be9 + p + Y*gamma:
                    n_b10_be9_p_Ygamma += 1
                    N_evts_deex += 1

                else:
                    if N_p < 0:
                        N_evts_incorrect += 1
                    else:
                        print("new channel: B10* -> Be9 + ..., event = {0:d}".format(event))

            elif N_be8 == 1:
                if (N_piminus == N_piplus == N_pi0 == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_p == 1 and N_gamma > 0):
                    # B10* -> Be8 + n + p + Y*gamma:
                    n_b10_be8_n_p_Ygamma += 1
                    N_evts_deex += 1

                elif (N_n == N_p == N_piminus == N_piplus == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_deuteron == 1 and N_pi0 > 0 and N_gamma > 0):
                    # B10* -> Be8 + d + X*pi0 + Y*gamma:
                    n_b10_be8_d_Xpi0_Ygamma += 1
                    N_evts_deex += 1

                elif (N_n == N_p == N_piminus == N_piplus == N_pi0 == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_deuteron == 1 and N_gamma > 0):
                    # B10* -> Be8 + d + X*pi0 + Y*gamma:
                    n_b10_be8_d_Ygamma += 1
                    N_evts_deex += 1

                elif (N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_p == 1 and N_pi0 > 0 and N_gamma > 0):
                    # B10* -> Be8 + n + p + X*pi0 + Y*gamma:
                    n_b10_be8_n_p_Xpi0_Ygamma += 1
                    N_evts_deex += 1

                else:
                    if N_p < 0:
                        N_evts_incorrect += 1
                    else:
                        print("new channel: B10* -> Be8 + ..., event = {0:d}".format(event))

            elif N_be7 == 1:
                if (N_n == N_p == N_piminus == N_piplus == N_deuteron == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_triton == 1 and N_pi0 > 0 and N_gamma > 0):
                    # B10* -> Be7 + t + X*pi0 + Y*gamma:
                    n_b10_be7_t_Xpi0_Ygamma += 1
                    N_evts_deex += 1

                elif (N_n == N_p == N_piminus == N_piplus == N_pi0 == N_deuteron == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_triton == 1 and N_gamma > 0):
                    # B10* -> Be7 + t + Y*gamma:
                    n_b10_be7_t_Ygamma += 1
                    N_evts_deex += 1

                elif (N_p == N_piminus == N_piplus == N_pi0 == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_deuteron == 1 and N_gamma > 0):
                    # B10* -> Be7 + n + d + Y*gamma:
                    n_b10_be7_n_d_Ygamma += 1
                    N_evts_deex += 1

                else:
                    if N_p < 0:
                        N_evts_incorrect += 1
                    else:
                        print("new channel: B10* -> Be7 + ..., event = {0:d}".format(event))

            elif N_li9 == 1:
                print("new channel: B10* -> Li9 + ..., event = {0:d}".format(event))

            elif N_li8 == 1:
                print("new channel: B10* -> Li8 + ..., event = {0:d}".format(event))

            elif N_li7 == 1:
                if (N_n == N_p == N_piminus == N_piplus == N_pi0 == N_deuteron == N_triton == N_alpha == N_li6 ==
                        N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_he3 == 1 and N_gamma > 0):
                    # B10* -> Li7 + he3 + Y*gamma:
                    n_b10_li7_he3_Ygamma += 1
                    N_evts_deex += 1

                elif (N_n == N_p == N_piminus == N_piplus == N_deuteron == N_triton == N_alpha == N_li6 ==
                        N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_he3 == 1 and N_pi0 > 0 and N_gamma > 0):
                    # B10* -> Li7 + he3 + X*pi0 + Y*gamma:
                    n_b10_li7_he3_Xpi0_Ygamma += 1
                    N_evts_deex += 1

                elif (N_n == N_piminus == N_piplus == N_pi0 == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == N_deuteron == 1 and N_gamma > 0):
                    # B10* -> Li7 + p + d + Y*gamma:
                    n_b10_li7_p_d_Ygamma += 1
                    N_evts_deex += 1

                elif (N_n == N_piminus == N_piplus == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == N_deuteron == 1 and N_pi0 > 0 and N_gamma > 0):
                    # B10* -> Li7 + p + d + X*pi0 + Y*gamma:
                    n_b10_li7_p_d_Xpi0_Ygamma += 1
                    N_evts_deex += 1

                else:
                    if N_p < 0:
                        N_evts_incorrect += 1
                    else:
                        print("new channel: B10* -> Li7 + ..., event = {0:d}".format(event))

            elif N_li6 == 1:
                if (N_n == N_p == N_piminus == N_piplus == N_pi0 == N_deuteron == N_triton == N_he3 == N_li7 ==
                        N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_alpha == 1 and N_gamma > 0):
                    # B10* -> Li6 + alpha + Y*gamma:
                    n_b10_li6_alpha_Ygamma += 1
                    N_evts_deex += 1

                elif (N_n == N_p == N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_li7 ==
                        N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_alpha == 1 and N_pi0 > 0 and N_gamma > 0):
                    # B10* -> Li6 + alpha + X*pi0 + Y*gamma:
                    n_b10_li6_alpha_Xpi0_Ygamma += 1
                    N_evts_deex += 1

                elif (N_n == N_piminus == N_piplus == N_deuteron == N_he3 == N_alpha == N_li7 ==
                        N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == N_triton == 1 and N_pi0 > 0 and N_gamma > 0):
                    # B10* -> Li6 + p + t + X*pi0 + Y*gamma:
                    n_b10_li6_p_t_Xpi0_Ygamma += 1
                    N_evts_deex += 1

                elif (N_n == N_piminus == N_piplus == N_pi0 == N_deuteron == N_he3 == N_alpha == N_li7 ==
                        N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == N_triton == 1 and N_gamma > 0):
                    # B10* -> Li6 + p + t + Y*gamma:
                    n_b10_li6_p_t_Ygamma += 1
                    N_evts_deex += 1

                elif (N_p == N_piminus == N_piplus == N_pi0 == N_deuteron == N_triton == N_alpha == N_li7 ==
                        N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_he3 == 1 and N_gamma > 0):
                    # B10* -> Li6 + n + he3 + Y*gamma:
                    n_b10_li6_n_he3_Ygamma += 1
                    N_evts_deex += 1

                else:
                    if N_p < 0:
                        N_evts_incorrect += 1
                    else:
                        print("new channel: B10* -> Li6 + ..., event = {0:d}".format(event))

            else:
                if (N_n == N_piminus == N_piplus == N_pi0 == N_gamma == N_deuteron == N_triton == N_he3 == N_li6 == N_li7 ==
                        N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == N_alpha == 1):
                    # B10* -> He5 + p + alpha:
                    n_b10_he5_p_alpha += 1
                    N_evts_deex += 1

                elif (N_n == N_piminus == N_piplus == N_gamma == N_deuteron == N_triton == N_he3 == N_li6 == N_li7 ==
                        N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == N_alpha == 1 and N_pi0 > 0):
                    # B10* -> He5 + p + alpha + X*pi0:
                    n_b10_he5_p_alpha_Xpi0 += 1
                    N_evts_deex += 1

                elif (N_p == N_piminus == N_piplus == N_pi0 == N_gamma == N_deuteron == N_triton == N_he3 == N_li6
                      == N_li7 ==
                        N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_alpha == 1):
                    # B10* -> Li5 + n + alpha:
                    n_b10_li5_n_alpha += 1
                    N_evts_deex += 1

                elif (N_p == N_piminus == N_piplus == N_gamma == N_deuteron == N_triton == N_he3 == N_li6 == N_li7 ==
                        N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_alpha == 1 and N_pi0 > 0):
                    # B10* -> Li5 + n + alpha + X*pi0:
                    n_b10_li5_n_alpha_Xpi0 += 1
                    N_evts_deex += 1

                elif (N_n == N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_li6 == N_li7 ==
                        N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == N_alpha == 1 and N_pi0 > 0 and N_gamma > 0):
                    # B10* -> He5 + p + alpha + X*pi0 + Y*gamma:
                    n_b10_he5_p_alpha_Xpi0_Ygamma += 1
                    N_evts_deex += 1

                elif (N_n == N_piminus == N_piplus == N_pi0 == N_deuteron == N_triton == N_he3 == N_li6 == N_li7 ==
                        N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == N_alpha == 1 and N_gamma > 0):
                    # B10* -> He5 + p + alpha + Y*gamma:
                    n_b10_he5_p_alpha_Ygamma += 1
                    N_evts_deex += 1

                else:
                    if N_p < 0:
                        N_evts_incorrect += 1
                    else:
                        print("new channel: B10* -> ..., event = {0:d}".format(event))

        elif isopdg_generator == 1000040100:
            number_be10_iso += 1
            # Be10 deexcitation:
            if N_c11 == 1:
                print("new channel: Be10* -> C11 + ..., event = {0:d}".format(event))

            elif N_c10 == 1:
                print("new channel: Be10* -> C10 + ..., event = {0:d}".format(event))

            elif N_c9 == 1:
                print("new channel: Be10* -> C9 + ..., event = {0:d}".format(event))

            elif N_b11 == 1:
                print("new channel: Be10* -> B11 + ..., event = {0:d}".format(event))

            elif N_b10 == 1:
                print("new channel: Be10* -> B10 + ..., event = {0:d}".format(event))

            elif N_b9 == 1:
                print("new channel: Be10* -> B9 + ..., event = {0:d}".format(event))

            elif N_b8 == 1:
                print("new channel: Be10* -> B8 + ..., event = {0:d}".format(event))

            elif N_be10 == 1:
                if (N_n == N_p == N_piminus == N_piplus == N_pi0 == N_gamma == N_deuteron == N_triton == N_he3 == N_alpha
                        == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0):
                    # no deexcitation:
                    n_be10_noDeex += 1

                elif (N_n == N_p == N_piminus == N_piplus == N_gamma == N_deuteron == N_triton == N_he3 == N_alpha
                        == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_pi0 > 0):
                    # Be10* -> Be10 + X*pi0:
                    n_be10_be10_Xpi0 += 1
                    N_evts_deex += 1

                elif (N_n == N_p == N_piminus == N_piplus == N_pi0 == N_deuteron == N_triton == N_he3 == N_alpha
                        == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_gamma > 0):
                    # Be10* -> Be10 + Y*gamma:
                    n_be10_be10_Ygamma += 1
                    N_evts_deex += 1

                elif (N_n == N_p == N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_pi0 > 0 and N_gamma > 0):
                    # Be10* -> Be10 + X*pi0 + Y*gamma:
                    n_be10_be10_Xpi0_Ygamma += 1
                    N_evts_deex += 1

                else:
                    if N_p < 0:
                        N_evts_incorrect += 1
                    else:
                        print("new channel: Be10* -> Be10 + ..., event = {0:d}".format(event))

            elif N_be9 == 1:
                if (N_p == N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == 1 and N_pi0 > 0 and N_gamma > 0):
                    # Be10* -> Be9 + n + X*pi0 + Y*gamma:
                    n_be10_be9_n_Xpi0_Ygamma += 1
                    N_evts_deex += 1

                elif (N_p == N_piminus == N_piplus == N_pi0 == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == 1 and N_gamma > 0):
                    # Be10* -> Be9 + n + Y*gamma:
                    n_be10_be9_n_Ygamma += 1
                    N_evts_deex += 1

                else:
                    if N_p < 0:
                        N_evts_incorrect += 1
                    else:
                        print("new channel: Be10* -> Be9 + ..., event = {0:d}".format(event))

            elif N_be8 == 1:
                if (N_p == N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == 2 and N_pi0 > 0 and N_gamma > 0):
                    # Be10* -> Be8 + 2n + X*pi0 + Y*gamma:
                    n_be10_be8_2n_Xpi0_Ygamma += 1
                    N_evts_deex += 1

                elif (N_p == N_piminus == N_piplus == N_pi0 == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_li9 == N_be7 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == 2 and N_gamma > 0):
                    # Be10* -> Be8 + 2n + Y*gamma:
                    n_be10_be8_2n_Ygamma += 1
                    N_evts_deex += 1

                else:
                    if N_p < 0:
                        N_evts_incorrect += 1
                    else:
                        print("new channel: Be10* -> Be8 + ..., event = {0:d}".format(event))

            elif N_be7 == 1:
                print("new channel: Be10* -> Be7 + ..., event = {0:d}".format(event))

            elif N_li9 == 1:
                if (N_n == N_piminus == N_piplus == N_pi0 == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == 1 and N_gamma > 0):
                    # Be10* -> Li9 + p + Y*gamma:
                    n_be10_li9_p_Ygamma += 1
                    N_evts_deex += 1

                elif (N_n == N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li8 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == 1 and N_pi0 > 0 and N_gamma > 0):
                    # Be10* -> Li9 + p + X*pi0 + Y*gamma:
                    n_be10_li9_p_Xpi0_Ygamma += 1
                    N_evts_deex += 1

                else:
                    if N_p < 0:
                        N_evts_incorrect += 1
                    else:
                        print("new channel: Be10* -> Li9 + ..., event = {0:d}".format(event))

            elif N_li8 == 1:
                if (N_piminus == N_piplus == N_pi0 == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_p == 1 and N_gamma > 0):
                    # Be10* -> Li8 + n + p + Y*gamma:
                    n_be10_li8_n_p_Ygamma += 1
                    N_evts_deex += 1

                elif (N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_p == 1 and N_pi0 > 0 and N_gamma > 0):
                    # Be10* -> Li8 + n + p + X*pi0 + Y*gamma:
                    n_be10_li8_n_p_Xpi0_Ygamma += 1
                    N_evts_deex += 1

                elif (N_n == N_p == N_piminus == N_piplus == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_deuteron == 1 and N_pi0 > 0 and N_gamma > 0):
                    # Be10* -> Li8 + d + X*pi0 + Y*gamma:
                    n_be10_li8_d_Xpi0_Ygamma += 1
                    N_evts_deex += 1

                elif (N_n == N_p == N_piminus == N_piplus == N_pi0 == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li7 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_deuteron == 1 and N_gamma > 0):
                    # Be10* -> Li8 + d + Y*gamma:
                    n_be10_li8_d_Ygamma += 1
                    N_evts_deex += 1

                else:
                    if N_p < 0:
                        N_evts_incorrect += 1
                    else:
                        print("new channel: Be10* -> Li8 + ..., event = {0:d}".format(event))

            elif N_li7 == 1:
                if (N_p == N_piminus == N_piplus == N_pi0 == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_deuteron == 1 and N_gamma > 0):
                    # Be10* -> Li7 + n + d + Y*gamma:
                    n_be10_li7_n_d_Ygamma += 1
                    N_evts_deex += 1

                elif (N_piminus == N_piplus == N_pi0 == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == 2 and N_p == 1 and N_gamma > 0):
                    # Be10* -> Li7 + 2n + p + Y*gamma:
                    n_be10_li7_2n_p_Ygamma += 1
                    N_evts_deex += 1

                elif (N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == 2 and N_p == 1 and N_pi0 > 0 and N_gamma > 0):
                    # Be10* -> Li7 + 2n + p + X*pi0 + Y*gamma:
                    n_be10_li7_2n_p_Xpi0_Ygamma += 1
                    N_evts_deex += 1

                elif (N_n == N_p == N_piminus == N_piplus == N_pi0 == N_deuteron == N_he3 == N_alpha == N_li6 ==
                        N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_triton == 1 and N_gamma > 0):
                    # Be10* -> Li7 + t + Y*gamma:
                    n_be10_li7_t_Ygamma += 1
                    N_evts_deex += 1

                elif (N_p == N_piminus == N_piplus == N_triton == N_he3 == N_alpha == N_li6 ==
                        N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_deuteron == 1 and N_pi0 > 0 and N_gamma > 0):
                    # Be10* -> Li7 + n + d + X*pi0 + Y*gamma:
                    n_be10_li7_n_d_Xpi0_Ygamma += 1
                    N_evts_deex += 1

                else:
                    if N_p < 0:
                        N_evts_incorrect += 1
                    else:
                        print("new channel: Be10* -> Li7 + ..., event = {0:d}".format(event))

            elif N_li6 == 1:
                if (N_p == N_piminus == N_piplus == N_pi0 == N_deuteron == N_he3 == N_alpha == N_li7 ==
                        N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_triton == 1 and N_gamma > 0):
                    # Be10* -> Li6 + n + t + Y*gamma:
                    n_be10_li6_n_t_Ygamma += 1
                    N_evts_deex += 1

                elif (N_p == N_piminus == N_piplus == N_triton == N_he3 == N_alpha == N_li7 ==
                        N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == 2 and N_deuteron == 1 and N_pi0 > 0 and
                        N_gamma > 0):
                    # Be10* -> Li6 + 2n + d + X*pi0 + Y*gamma:
                    n_be10_li6_2n_d_Xpi0_Ygamma += 1
                    N_evts_deex += 1

                elif (N_piminus == N_piplus == N_pi0 == N_deuteron == N_triton == N_he3 == N_alpha == N_li7 ==
                        N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == 3 and N_p == 1 and N_gamma > 0):
                    # Be10* -> Li6 + 3n + p + Y*gamma:
                    n_be10_li6_3n_p_Ygamma += 1
                    N_evts_deex += 1

                elif (N_p == N_piminus == N_piplus == N_deuteron == N_he3 == N_alpha == N_li7 ==
                        N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_triton == 1 and N_pi0 > 0 and N_gamma > 0):
                    # Be10* -> Li6 + n + t + X*pi0 + Y*gamma:
                    n_be10_li6_n_t_Xpi0_Ygamma += 1
                    N_evts_deex += 1

                elif (N_p == N_piminus == N_piplus == N_pi0 == N_triton == N_he3 == N_alpha == N_li7 ==
                        N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == 2 and N_deuteron == 1 and N_gamma > 0):
                    # Be10* -> Li6 + 2n + d + Y*gamma:
                    n_be10_li6_2n_d_Ygamma += 1
                    N_evts_deex += 1

                elif (N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha == N_li7 ==
                        N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == 3 and N_p == 1 and N_pi0 > 0 and N_gamma > 0):
                    # Be10* -> Li6 + 3n + p + X*pi0 + Y*gamma:
                    n_be10_li6_3n_p_Xpi0_Ygamma += 1
                    N_evts_deex += 1

                else:
                    if N_p < 0:
                        N_evts_incorrect += 1
                    else:
                        print("new channel: Be10* -> Li6 + ..., event = {0:d}".format(event))

            else:
                if (N_piminus == N_piplus == N_pi0 == N_gamma == N_deuteron == N_triton == N_he3 == N_li6 == N_li7 ==
                        N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_p == N_alpha == 1):
                    # Be10* -> H4 + n + p + alpha:
                    n_be10_h4_n_p_alpha += 1
                    N_evts_deex += 1

                elif (N_n == N_piminus == N_piplus == N_pi0 == N_gamma == N_triton == N_he3 == N_alpha == N_li6 == N_li7 ==
                        N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == N_deuteron == 1):
                    # Be10* -> He7 + p + d:
                    n_be10_he7_p_d += 1
                    N_evts_deex += 1

                elif (N_piminus == N_piplus == N_gamma == N_deuteron == N_he3 == N_alpha == N_li6 == N_li7 ==
                        N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_p == N_triton == 1 and N_pi0 > 0):
                    # Be10* -> He5 + n + p + t + X*pi0:
                    n_be10_he5_n_p_t_Xpi0 += 1
                    N_evts_deex += 1

                elif (N_piminus == N_piplus == N_pi0 == N_gamma == N_triton == N_he3 == N_alpha == N_li6 == N_li7 ==
                        N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == 2 and N_p == N_deuteron == 1):
                    # Be10* -> He5 + 2n + p + d:
                    n_be10_he5_2n_p_d += 1
                    N_evts_deex += 1

                elif (N_p == N_piminus == N_piplus == N_pi0 == N_gamma == N_triton == N_he3 == N_alpha == N_li6 == N_li7 ==
                        N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == 1 and N_deuteron == 2):
                    # Be10* -> He5 + n + 2d:
                    n_be10_he5_n_2d += 1
                    N_evts_deex += 1

                elif (N_n == N_p == N_piminus == N_piplus == N_pi0 == N_gamma == N_deuteron == N_triton == N_alpha
                      == N_li6 == N_li7 ==
                        N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_he3 == 1):
                    # Be10* -> He7 + He3:
                    n_be10_he7_he3 += 1
                    N_evts_deex += 1

                elif (N_piminus == N_piplus == N_gamma == N_deuteron == N_triton == N_he3 == N_li6 == N_li7 ==
                        N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_p == N_alpha == 1 and N_pi0 > 0):
                    # Be10* -> H4 + n + p + alpha + X*pi0:
                    n_be10_h4_n_p_alpha_Xpi0 += 1
                    N_evts_deex += 1

                elif (N_piminus == N_piplus == N_pi0 == N_gamma == N_deuteron == N_he3 == N_alpha == N_li6 == N_li7 ==
                        N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_p == N_triton == 1):
                    # Be10* -> He5 + n + p + t:
                    n_be10_he5_n_p_t += 1
                    N_evts_deex += 1

                elif (N_p == N_piminus == N_piplus == N_pi0 == N_gamma == N_he3 == N_alpha == N_li6 == N_li7 ==
                        N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_deuteron == N_triton == 1) or \
                        (N_p == N_piminus == N_piplus == N_pi0 == N_gamma == N_triton == N_he3 == N_li6 == N_li7 ==
                         N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                         N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_deuteron == N_alpha == 1):
                    # Be10* -> alpha + n + d + t:
                    n_be10_alpha_n_d_t += 1
                    N_evts_deex += 1

                elif (N_piminus == N_piplus == N_pi0 == N_gamma == N_triton == N_he3 == N_alpha == N_li6 == N_li7 ==
                        N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_p == N_deuteron == 1):
                    # Be10* -> He6 + n + p + d:
                    n_be10_he6_n_p_d += 1
                    N_evts_deex += 1

                elif (N_n == N_p == N_piminus == N_piplus == N_pi0 == N_gamma == N_triton == N_he3 == N_alpha == N_li6
                      == N_li7 ==
                        N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_deuteron == 2):
                    # Be10* -> He6 + 2d:
                    n_be10_he6_2d += 1
                    N_evts_deex += 1

                elif (N_p == N_piminus == N_piplus == N_gamma == N_he3 == N_alpha == N_li6 == N_li7 ==
                        N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                        N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_deuteron == N_triton == 1 and N_pi0 > 0):
                    # Be10* -> alpha + n + d + t + X*pi0:
                    n_be10_alpha_n_d_t_Xpi0 += 1
                    N_evts_deex += 1

                else:
                    if N_p < 0:
                        N_evts_incorrect += 1
                    else:
                        # other deexcitation channels:
                        # print("new channel: Be10* -> ..., event = {0:d}".format(event))
                        n_be10_rest += 1
                        N_evts_deex += 1

                        # print("--------------------")
                        # print(N_n)
                        # print(N_p)
                        # print(N_piminus)
                        # print(N_piplus)
                        # print(N_pi0)
                        # print(N_gamma)
                        # print(N_deuteron)
                        # print(N_triton)
                        # print(N_he3)
                        # print(N_alpha)
                        # print(N_li6)
                        # print(N_li7)
                        # print(N_li8)
                        # print(N_li9)
                        # print(N_be7)
                        # print(N_be8)
                        # print(N_be9)
                        # print(N_be10)
                        # print(N_b8)
                        # print(N_b9)
                        # print(N_b10)
                        # print(N_b11)
                        # print(N_c9)
                        # print(N_c10)
                        # print(N_c11)
                        # print(N_c12)

    print(N_entries_generator)
    print(N_evts_deex)
    print(N_evts_incorrect)
    print("")
    print(number_c11_iso)
    print(n_c11_noDeex)
    print(n_c11_c11_Ygamma)
    print(n_c11_c11_Xpi0)
    print(n_c11_c11_Xpi0_Ygamma)
    print(n_c11_c10_n_Xpi0_Ygamma)
    print(n_c11_c10_n_Ygamma)
    print(n_c11_b10_p_Ygamma)
    print(n_c11_b10_p_Xpi0_Ygamma)
    print(n_c11_b9_n_p_Ygamma)
    print(n_c11_b9_d_Ygamma)
    print(n_c11_b9_d_Xpi0_Ygamma)
    print(n_c11_b9_n_p_Xpi0_Ygamma)
    print(n_c11_be9_2p_Ygamma)
    print(n_c11_be9_2p_Xpi0_Ygamma)
    print(n_c11_be8_p_d_Xpi0_Ygamma)
    print(n_c11_be8_p_d_Ygamma)
    print(n_c11_be8_he3_Ygamma)
    print(n_c11_be8_he3_Xpi0_Ygamma)
    print(n_c11_be8_n_2p_Ygamma)
    print(n_c11_be8_n_2p_Xpi0_Ygamma)
    print(n_c11_be7_alpha_Xpi0_Ygamma)
    print(n_c11_be7_alpha_Ygamma)
    print(n_c11_li6_p_alpha_Ygamma)
    print(n_c11_li6_p_alpha_Xpi0_Ygamma)
    print(n_c11_li5_d_alpha)
    print(n_c11_li5_d_alpha_Xpi0)
    print(n_c11_li5_n_p_alpha_Xpi0)
    print(n_c11_li5_n_p_alpha)
    print(n_c11_li5_d_alpha_Ygamma)
    print(n_c11_he3_2alpha)
    print("")
    print(number_c10_iso)
    print(n_c10_noDeex)
    print(n_c10_c10_Xpi0)
    print(n_c10_c10_Ygamma)
    print(n_c10_c10_Xpi0_Ygamma)
    print(n_c10_c9_n_Ygamma)
    print(n_c10_c9_n_Xpi0_Ygamma)
    print(n_c10_b9_p_Ygamma)
    print(n_c10_b9_p_Xpi0_Ygamma)
    print(n_c10_b8_n_p_Ygamma)
    print(n_c10_b8_n_p_Xpi0_Ygamma)
    print(n_c10_b8_d_Xpi0_Ygamma)
    print(n_c10_b8_d_Ygamma)
    print(n_c10_be8_2p_Ygamma)
    print(n_c10_be8_2p_Xpi0_Ygamma)
    print(n_c10_be7_p_d_Ygamma)
    print(n_c10_be7_n_2p_Ygamma)
    print(n_c10_be7_p_d_Xpi0_Ygamma)
    print(n_c10_be7_n_2p_Xpi0_Ygamma)
    print(n_c10_be7_he3_Xpi0_Ygamma)
    print(n_c10_be7_he3_Ygamma)
    print(n_c10_li6_2p_d_Xpi0_Ygamma)
    print(n_c10_li6_p_he3_Ygamma)
    print(n_c10_li6_2p_d_Ygamma)
    print(n_c10_li6_n_3p_Ygamma)
    print(n_c10_li6_p_he3_Xpi0_Ygamma)
    print(n_c10_li6_n_3p_Xpi0_Ygamma)
    print(n_c10_li5_n_2p_d)
    print(n_c10_he3_p_d_alpha)
    print(n_c10_b7_2n_p_Xpi0)
    print(n_c10_be6_n_p_d_Xpi0)
    print(n_c10_li5_p_2d)
    print(n_c10_be6_p_t)
    print(n_c10_he3_n_2p_alpha)
    print(n_c10_li5_n_p_he3)
    print(n_c10_b7_2n_p)
    print(n_c10_be6_n_p_d)
    print(n_c10_li5_d_he3)
    print(n_c10_he3_p_d_alpha_Xpi0)
    print(n_c10_li5_p_alpha_Xpi0)
    print(n_c10_li4_n_p_alpha_Xpi0)
    print(n_c10_rest)
    print("")
    print(number_b11_iso)
    print(n_b11_noDeex)
    print(n_b11_b11_Ygamma)
    print(n_b11_b11_Xpi0)
    print(n_b11_b11_Xpi0_Ygamma)
    print(n_b11_b10_n_Ygamma)
    print(n_b11_b10_n_Xpi0_Ygamma)
    print(n_b11_b9_2n_Ygamma)
    print(n_b11_b9_2n_Xpi0_Ygamma)
    print(n_b11_be10_p_Ygamma)
    print(n_b11_be10_p_Xpi0_Ygamma)
    print(n_b11_be9_n_p_Ygamma)
    print(n_b11_be9_d_Ygamma)
    print(n_b11_be9_d_Xpi0_Ygamma)
    print(n_b11_be9_n_p_Xpi0_Ygamma)
    print(n_b11_be8_n_d_Xpi0_Ygamma)
    print(n_b11_be8_n_d_Ygamma)
    print(n_b11_be8_t_Xpi0_Ygamma)
    print(n_b11_be8_t_Ygamma)
    print(n_b11_be8_2n_p_Ygamma)
    print(n_b11_Li7_alpha_Ygamma)
    print(n_b11_Li7_alpha_Xpi0_Ygamma)
    print(n_b11_Li6_n_alpha_Xpi0_Ygamma)
    print(n_b11_Li6_n_alpha_Ygamma)
    print(n_b11_he5_d_alpha)
    print(n_b11_he5_d_alpha_Xpi0)
    print(n_b11_he5_n_p_alpha)
    print(n_b11_he6_p_alpha)
    print(n_b11_he5_n_p_alpha_Xpi0)
    print(n_b11_t_2alpha)
    print(n_b11_he6_p_alpha_Xpi0)
    print(n_b11_t_2alpha_Xpi0)
    print("")
    print(number_b10_iso)
    print(n_b10_noDeex)
    print(n_b10_b10_Xpi0)
    print(n_b10_b10_Ygamma)
    print(n_b10_b10_Xpi0_Ygamma)
    print(n_b10_b9_n_Ygamma)
    print(n_b10_b9_n_Xpi0_Ygamma)
    print(n_b10_be9_p_Xpi0_Ygamma)
    print(n_b10_be9_p_Ygamma)
    print(n_b10_be8_n_p_Ygamma)
    print(n_b10_be8_d_Xpi0_Ygamma)
    print(n_b10_be8_d_Ygamma)
    print(n_b10_be8_n_p_Xpi0_Ygamma)
    print(n_b10_be7_t_Xpi0_Ygamma)
    print(n_b10_be7_t_Ygamma)
    print(n_b10_be7_n_d_Ygamma)
    print(n_b10_li7_he3_Ygamma)
    print(n_b10_li7_he3_Xpi0_Ygamma)
    print(n_b10_li7_p_d_Ygamma)
    print(n_b10_li7_p_d_Xpi0_Ygamma)
    print(n_b10_li6_alpha_Ygamma)
    print(n_b10_li6_alpha_Xpi0_Ygamma)
    print(n_b10_li6_p_t_Xpi0_Ygamma)
    print(n_b10_li6_p_t_Ygamma)
    print(n_b10_li6_n_he3_Ygamma)
    print(n_b10_he5_p_alpha)
    print(n_b10_he5_p_alpha_Xpi0)
    print(n_b10_li5_n_alpha)
    print(n_b10_li5_n_alpha_Xpi0)
    print(n_b10_he5_p_alpha_Xpi0_Ygamma)
    print(n_b10_he5_p_alpha_Ygamma)
    print("")
    print(number_be10_iso)
    print(n_be10_noDeex)
    print(n_be10_be10_Xpi0)
    print(n_be10_be10_Ygamma)
    print(n_be10_be10_Xpi0_Ygamma)
    print(n_be10_be9_n_Xpi0_Ygamma)
    print(n_be10_be9_n_Ygamma)
    print(n_be10_be8_2n_Xpi0_Ygamma)
    print(n_be10_be8_2n_Ygamma)
    print(n_be10_li9_p_Ygamma)
    print(n_be10_li9_p_Xpi0_Ygamma)
    print(n_be10_li8_n_p_Ygamma)
    print(n_be10_li8_n_p_Xpi0_Ygamma)
    print(n_be10_li8_d_Xpi0_Ygamma)
    print(n_be10_li8_d_Ygamma)
    print(n_be10_li7_n_d_Ygamma)
    print(n_be10_li7_2n_p_Ygamma)
    print(n_be10_li7_2n_p_Xpi0_Ygamma)
    print(n_be10_li7_t_Ygamma)
    print(n_be10_li7_n_d_Xpi0_Ygamma)
    print(n_be10_li6_n_t_Ygamma)
    print(n_be10_li6_2n_d_Xpi0_Ygamma)
    print(n_be10_li6_3n_p_Ygamma)
    print(n_be10_li6_n_t_Xpi0_Ygamma)
    print(n_be10_li6_2n_d_Ygamma)
    print(n_be10_li6_3n_p_Xpi0_Ygamma)
    print(n_be10_h4_n_p_alpha)
    print(n_be10_he7_p_d)
    print(n_be10_he5_n_p_t_Xpi0)
    print(n_be10_he5_2n_p_d)
    print(n_be10_he5_n_2d)
    print(n_be10_he7_he3)
    print(n_be10_h4_n_p_alpha_Xpi0)
    print(n_be10_he5_n_p_t)
    print(n_be10_alpha_n_d_t)
    print(n_be10_he6_n_p_d)
    print(n_be10_he6_2d)
    print(n_be10_alpha_n_d_t_Xpi0)
    print(n_be10_rest)

""" deexcitation cahnnels with gamma (but without considering pi0): """
# preallocate number of event with deexcitation:
N_evts_deex = 0
# preallocate number of events with incorrect deexcitation (e.g. number of protons -999, ...):
N_evts_incorrect = 0
# C11:
number_c11_iso = 0
n_c11_noDeex = 0
n_c11_c11_Ygamma = 0
n_c11_c10_n_Ygamma = 0
n_c11_b10_p_Ygamma = 0
n_c11_b9_n_p_Ygamma = 0
n_c11_b9_d_Ygamma = 0
n_c11_be9_2p_Ygamma = 0
n_c11_be8_p_d_Ygamma = 0
n_c11_be8_he3_Ygamma = 0
n_c11_be8_n_2p_Ygamma = 0
n_c11_be7_alpha_Ygamma = 0
n_c11_li6_p_alpha_Ygamma = 0
n_c11_li5_d_alpha = 0
n_c11_li5_n_p_alpha = 0
n_c11_li5_d_alpha_Ygamma = 0
n_c11_he3_2alpha = 0
# C10 deex:
number_c10_iso = 0
n_c10_noDeex = 0
n_c10_c10_Ygamma = 0
n_c10_c9_n_Ygamma = 0
n_c10_b9_p_Ygamma = 0
n_c10_b8_n_p_Ygamma = 0
n_c10_b8_d_Ygamma = 0
n_c10_be8_2p_Ygamma = 0
n_c10_be7_p_d_Ygamma = 0
n_c10_be7_n_2p_Ygamma = 0
n_c10_be7_he3_Ygamma = 0
n_c10_li6_p_he3_Ygamma = 0
n_c10_li6_2p_d_Ygamma = 0
n_c10_li6_n_3p_Ygamma = 0
n_c10_li5_n_2p_d = 0
n_c10_he3_p_d_alpha = 0
n_c10_li5_p_2d = 0
n_c10_be6_p_t = 0
n_c10_he3_n_2p_alpha = 0
n_c10_li5_n_p_he3 = 0
n_c10_b7_2n_p = 0
n_c10_be6_n_p_d = 0
n_c10_li5_d_he3 = 0
n_c10_li5_p_alpha = 0
n_c10_li4_n_p_alpha = 0
n_c10_rest = 0
# B11 deexcitation:
number_b11_iso = 0
n_b11_noDeex = 0
n_b11_b11_Ygamma = 0
n_b11_b10_n_Ygamma = 0
n_b11_b9_2n_Ygamma = 0
n_b11_be10_p_Ygamma = 0
n_b11_be9_n_p_Ygamma = 0
n_b11_be9_d_Ygamma = 0
n_b11_be8_n_d_Ygamma = 0
n_b11_be8_t_Ygamma = 0
n_b11_be8_2n_p_Ygamma = 0
n_b11_Li7_alpha_Ygamma = 0
n_b11_Li6_n_alpha_Ygamma = 0
n_b11_he5_d_alpha = 0
n_b11_he5_n_p_alpha = 0
n_b11_he6_p_alpha = 0
n_b11_t_2alpha = 0
# B10 deexcitation:
number_b10_iso = 0
n_b10_noDeex = 0
n_b10_b10_Ygamma = 0
n_b10_b9_n_Ygamma = 0
n_b10_be9_p_Ygamma = 0
n_b10_be8_n_p_Ygamma = 0
n_b10_be8_d_Ygamma = 0
n_b10_be7_t_Ygamma = 0
n_b10_be7_n_d_Ygamma = 0
n_b10_li7_he3_Ygamma = 0
n_b10_li7_p_d_Ygamma = 0
n_b10_li6_alpha_Ygamma = 0
n_b10_li6_p_t_Ygamma = 0
n_b10_li6_n_he3_Ygamma = 0
n_b10_he5_p_alpha = 0
n_b10_li5_n_alpha = 0
n_b10_he5_p_alpha_Ygamma = 0
# Be10 deexcitation:
number_be10_iso = 0
n_be10_noDeex = 0
n_be10_be10_Ygamma = 0
n_be10_be9_n_Ygamma = 0
n_be10_be8_2n_Ygamma = 0
n_be10_li9_p_Ygamma = 0
n_be10_li8_n_p_Ygamma = 0
n_be10_li8_d_Ygamma = 0
n_be10_li7_n_d_Ygamma = 0
n_be10_li7_2n_p_Ygamma = 0
n_be10_li7_t_Ygamma = 0
n_be10_li6_n_t_Ygamma = 0
n_be10_li6_3n_p_Ygamma = 0
n_be10_li6_2n_d_Ygamma = 0
n_be10_h4_n_p_alpha = 0
n_be10_he7_p_d = 0
n_be10_he5_2n_p_d = 0
n_be10_he5_n_2d = 0
n_be10_he7_he3 = 0
n_be10_he5_n_p_t = 0
n_be10_alpha_n_d_t = 0
n_be10_he6_n_p_d = 0
n_be10_he6_2d = 0
n_be10_rest = 0

# loop over events:
for event in range(N_entries_generator):
    # preallocate array, where PDG of final particles are stored:
    pdg_array = np.array([])

    # get the current event from generator tree:
    rtree_generator.GetEntry(event)

    # get channel ID:
    channelID = int(rtree_generator.GetBranch('t_channelID').GetLeaf('t_channelID').GetValue())
    # check channelID:
    if channelID == 2 or channelID == 3:
        # go to next event:
        continue
    # get number of proton, neutron, piminus, piplus from channelID (in channelID no gamma's are stored.
    # K_plus and K_minus are stored in channelID, but low statistics):
    n_p, n_n, n_piminus, n_piplus = NC_background_functions.get_number_of_particles_of_channelid(channelID)

    # get isopdg:
    isopdg_generator = int(rtree_generator.GetBranch('t_isoPdg').GetLeaf('t_isoPdg').GetValue())

    # get Npars:
    Npars_generator = int(rtree_generator.GetBranch('t_Npars').GetLeaf('t_Npars').GetValue())
    # loop over Npars to get pdg:
    for index in range(Npars_generator):
        pdg = int(rtree_generator.GetBranch('t_pdg').GetLeaf('t_pdg').GetValue(index))
        pdg_array = np.append(pdg_array, pdg)

    # get number of particles from pdg_array:
    (N_n, N_p, N_piminus, N_piplus, N_pi0, N_gamma, N_deuteron, N_triton, N_he3, N_alpha, N_li6, N_li7, N_li8, N_li9,
     N_be7, N_be8, N_be9, N_be10, N_b8, N_b9, N_b10, N_b11, N_c9, N_c10, N_c11, N_c12) \
        = get_number_of_particles(pdg_array)

    # calculate the number of particles corresponding only to deexcitation
    # ((after NC interaction and after deexcitation) - (after NC interaction and before deexcitation)):
    N_n = N_n - n_n
    N_p = N_p - n_p
    N_piminus = N_piminus - n_piminus
    N_piplus = N_piplus - n_piplus

    # get deexcitation channels:
    if isopdg_generator == 1000060120:
        print("new channel: C12* -> C12 + ..., event = {0:d}".format(event))

    elif isopdg_generator == 1000060110:
        number_c11_iso += 1
        # C11 deexcitation:
        if N_c11 == 1:
            if (N_n == N_p == N_piminus == N_piplus == N_gamma == N_deuteron == N_triton == N_he3 == N_alpha
                    == N_li6 ==
                    N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c12 == 0):
                # no deexcitation:
                n_c11_noDeex += 1

            elif (N_n == N_p == N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                    N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c12 == 0 and N_gamma > 0):
                # C11* -> C11 + 2gamma:
                n_c11_c11_Ygamma += 1
                N_evts_deex += 1

            else:
                if N_p < 0:
                    N_evts_incorrect += 1
                else:
                    print("new channel: C11* -> C11 +..., event = {0:d}".format(event))

        elif N_c10 == 1:
            if (N_p == N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                    N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c11 == N_c12 == 0 and N_n == 1 and N_gamma > 0):
                # C11* -> C10 + n + Y*gamma:
                n_c11_c10_n_Ygamma += 1
                N_evts_deex += 1

            else:
                if N_p < 0:
                    N_evts_incorrect += 1
                else:
                    print("new channel: C11* -> C10 +..., event = {0:d}".format(event))

        elif N_c9 == 1:
            print("new channel: C11* -> C9 + ..., event = {0:d}".format(event))

        elif N_b11 == 1:
            print("new channel: C11* -> B11 + ..., event = {0:d}".format(event))

        elif N_b10 == 1:
            if (N_n == N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                    N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == 1 and N_gamma > 0):
                # C11* -> B10 + p + Y*gamma:
                n_c11_b10_p_Ygamma += 1
                N_evts_deex += 1

            else:
                if N_p < 0:
                    N_evts_incorrect += 1
                else:
                    print("new channel: C11* -> B10 +..., event = {0:d}".format(event))

        elif N_b9 == 1:
            if (N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                    N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_p == 1 and N_gamma > 0):
                # C11* -> B9 + n + p + Y*gamma:
                n_c11_b9_n_p_Ygamma += 1
                N_evts_deex += 1

            elif (N_n == N_p == N_piminus == N_piplus == N_triton == N_he3 == N_alpha == N_li6 ==
                    N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_deuteron == 1 and N_gamma > 0):
                # C11* -> B9 + d + Y*gamma:
                n_c11_b9_d_Ygamma += 1
                N_evts_deex += 1

            else:
                if N_p < 0:
                    N_evts_incorrect += 1
                else:
                    print("ERR: C11* -> B9 +..., event = {0:d}".format(event))

        elif N_b8 == 1:
            print("new channel: C11* -> B8 +..., event = {0:d}".format(event))

        elif N_be10 == 1:
            print("new channel: C11* -> Be10 +..., event = {0:d}".format(event))

        elif N_be9 == 1:
            if (N_n == N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                    N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == 2 and N_gamma > 0):
                # C11* -> Be9 + 2*p + Y*gamma:
                n_c11_be9_2p_Ygamma += 1
                N_evts_deex += 1

            else:
                if N_p < 0:
                    N_evts_incorrect += 1
                else:
                    print("new channel: C11* -> Be9 +..., event = {0:d}".format(event))

        elif N_be8 == 1:
            if (N_n == N_piminus == N_piplus == N_triton == N_he3 == N_alpha == N_li6 ==
                    N_li7 == N_li8 == N_li9 == N_be7 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == N_deuteron == 1 and N_gamma > 0):
                # C11* -> Be8 + p + d + Y*gamma:
                n_c11_be8_p_d_Ygamma += 1
                N_evts_deex += 1

            elif (N_p == N_n == N_piminus == N_piplus == N_deuteron == N_triton == N_alpha == N_li6 ==
                    N_li7 == N_li8 == N_li9 == N_be7 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_he3 == 1 and N_gamma > 0):
                # C11* -> Be8 + He3 + Y*gamma:
                n_c11_be8_he3_Ygamma += 1
                N_evts_deex += 1

            elif (N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                    N_li7 == N_li8 == N_li9 == N_be7 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == 1 and N_p == 2 and N_gamma > 0):
                # C11* -> Be8 + n + 2p + Y*gamma:
                n_c11_be8_n_2p_Ygamma += 1
                N_evts_deex += 1

            else:
                if N_p < 0:
                    N_evts_incorrect += 1
                else:
                    print("new channel: C11* -> Be8 +..., event = {0:d}".format(event))

        elif N_be7 == 1:
            if (N_n == N_p == N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_li6 ==
                    N_li7 == N_li8 == N_li9 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_alpha == 1 and N_gamma > 0):
                # C11* -> Be7 + alpha + Y*gamma:
                n_c11_be7_alpha_Ygamma += 1
                N_evts_deex += 1

            else:
                if N_p < 0:
                    N_evts_incorrect += 1
                else:
                    print("new channel: C11* -> Be7 +..., event = {0:d}".format(event))

        elif N_li9 == 1:
            print("new channel: C11* -> Li9 +..., event = {0:d}".format(event))

        elif N_li8 == 1:
            print("new channel: C11* -> Li8 +..., event = {0:d}".format(event))

        elif N_li7 == 1:
            print("new channel: C11* -> Li7 +..., event = {0:d}".format(event))

        elif N_li6 == 1:
            if (N_n == N_piminus == N_piplus == N_deuteron == N_triton == N_he3 ==
                    N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == N_alpha == 1 and N_gamma > 0):
                # C11* -> Li6 + p + alpha + Y*gamma:
                n_c11_li6_p_alpha_Ygamma += 1
                N_evts_deex += 1

            else:
                if N_p < 0:
                    N_evts_incorrect += 1
                else:
                    print("new channel: C11* -> Li6 +..., event = {0:d}".format(event))

        else:
            if (N_n == N_p == N_piminus == N_piplus == N_gamma == N_triton == N_he3 == N_li6 ==
                    N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_deuteron == N_alpha == 1):
                # C11* -> Li5 + d + alpha:
                n_c11_li5_d_alpha += 1
                N_evts_deex += 1

            elif (N_piminus == N_piplus == N_gamma == N_deuteron == N_triton == N_he3 == N_li6 ==
                    N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_p == N_alpha == 1):
                # C11* -> Li5 + n + p + alpha:
                n_c11_li5_n_p_alpha += 1
                N_evts_deex += 1

            elif (N_n == N_p == N_piminus == N_piplus == N_triton == N_he3 == N_li6 ==
                    N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_deuteron == N_alpha == 1 and N_gamma > 0):
                # C11* -> Li5 + d + alpha + Y*gamma:
                n_c11_li5_d_alpha_Ygamma += 1
                N_evts_deex += 1

            elif (N_n == N_p == N_piminus == N_piplus == N_gamma == N_deuteron == N_triton == N_he3 == N_li6 ==
                    N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_alpha == 2):
                # C11* -> He3 + 2*alpha:
                n_c11_he3_2alpha += 1
                N_evts_deex += 1

            else:
                print("new channel: C11* -> ..., event = {0:d}".format(event))

    elif isopdg_generator == 1000060100:
        number_c10_iso += 1
        # C10 deexcitation:
        if N_c10 == 1:
            if (N_n == N_p == N_piminus == N_piplus == N_gamma == N_deuteron == N_triton == N_he3 == N_alpha
                    == N_li6 ==
                    N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c11 == N_c12 == 0):
                # no deexcitation:
                n_c10_noDeex += 1

            elif (N_n == N_p == N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha
                  == N_li6 ==
                    N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c11 == N_c12 == 0 and N_gamma > 0):
                # C10* -> C10 + Y*gamma:
                n_c10_c10_Ygamma += 1
                N_evts_deex += 1

            else:
                if N_p < 0:
                    N_evts_incorrect += 1
                else:
                    print("new channel: C10* -> C10 + ..., event = {0:d}".format(event))

        elif N_c9 == 1:
            if (N_p == N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                    N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c10 == N_c11 == N_c12 == 0 and N_n == 1 and N_gamma > 0):
                # C10* -> C9 + n + Y*gamma:
                n_c10_c9_n_Ygamma += 1
                N_evts_deex += 1

            else:
                if N_p < 0:
                    N_evts_incorrect += 1
                else:
                    print("new channel: C10* -> C9 + ..., event = {0:d}".format(event))

        elif N_b11 == 1:
            print("new channel: C10* -> B11 +..., event = {0:d}".format(event))

        elif N_b10 == 1:
            print("new channel: C10* -> B10 +..., event = {0:d}".format(event))

        elif N_b9 == 1:
            if (N_n == N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                    N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == 1 and N_gamma > 0):
                # C10* -> B9 + p + Y*gamma:
                n_c10_b9_p_Ygamma += 1
                N_evts_deex += 1

            else:
                if N_p < 0:
                    N_evts_incorrect += 1
                else:
                    print("new channel: C10* -> B9 + ..., event = {0:d}".format(event))

        elif N_b8 == 1:
            if (N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                    N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_p == 1 and N_gamma > 0):
                # C10* -> B8 + n + p + Y*gamma:
                n_c10_b8_n_p_Ygamma += 1
                N_evts_deex += 1

            elif (N_n == N_p == N_piminus == N_piplus == N_triton == N_he3 == N_alpha == N_li6 ==
                    N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_deuteron == 1 and N_gamma > 0):
                # C10* -> B8 + d + Y*gamma:
                n_c10_b8_d_Ygamma += 1
                N_evts_deex += 1

            else:
                if N_p < 0:
                    N_evts_incorrect += 1
                else:
                    print("new channel: C10* -> B8 + ..., event = {0:d}".format(event))

        elif N_be10 == 1:
            print("new channel: C10* -> Be10 +..., event = {0:d}".format(event))

        elif N_be9 == 1:
            print("new channel: C10* -> Be9 +..., event = {0:d}".format(event))

        elif N_be8 == 1:
            if (N_n == N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                    N_li7 == N_li8 == N_li9 == N_be7 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == 2 and N_gamma > 0):
                # C10* -> Be8 + 2p + Y*gamma:
                n_c10_be8_2p_Ygamma += 1
                N_evts_deex += 1

            else:
                if N_p < 0:
                    N_evts_incorrect += 1
                else:
                    print("new channel: C10* -> Be8 + ..., event = {0:d}".format(event))

        elif N_be7 == 1:
            if (N_n == N_piminus == N_piplus == N_triton == N_he3 == N_alpha == N_li6 ==
                    N_li7 == N_li8 == N_li9 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == N_deuteron == 1 and N_gamma > 0):
                # C10* -> Be7 + p +d + Y*gamma:
                n_c10_be7_p_d_Ygamma += 1
                N_evts_deex += 1

            elif (N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                    N_li7 == N_li8 == N_li9 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == 1 and N_p == 2 and N_gamma > 0):
                # C10* -> Be7 + n + 2p + Y*gamma:
                n_c10_be7_n_2p_Ygamma += 1
                N_evts_deex += 1

            elif (N_n == N_p == N_piminus == N_piplus == N_deuteron == N_triton == N_alpha == N_li6 ==
                    N_li7 == N_li8 == N_li9 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_he3 == 1 and N_gamma > 0):
                # C10* -> Be7 + He3 + Y*gamma:
                n_c10_be7_he3_Ygamma += 1
                N_evts_deex += 1

            else:
                if N_p < 0:
                    N_evts_incorrect += 1
                else:
                    print("new channel: C10* -> Be7 + ..., event = {0:d}".format(event))

        elif N_li9 == 1:
            print("new channel: C10* -> Li9 +..., event = {0:d}".format(event))

        elif N_li8 == 1:
            print("new channel: C10* -> Li8 +..., event = {0:d}".format(event))

        elif N_li7 == 1:
            print("new channel: C10* -> Li7 +..., event = {0:d}".format(event))

        elif N_li6 == 1:
            if (N_n == N_piminus == N_piplus == N_triton == N_deuteron == N_alpha ==
                    N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == N_he3 == 1 and N_gamma > 0):
                # C10* -> Li6 + p + He3 + Y*gamma:
                n_c10_li6_p_he3_Ygamma += 1
                N_evts_deex += 1

            elif (N_n == N_piminus == N_piplus == N_triton == N_he3 == N_alpha ==
                    N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == 2 and N_deuteron == 1 and N_gamma > 0):
                # C10* -> Li6 + 2p + d + Y*gamma:
                n_c10_li6_2p_d_Ygamma += 1
                N_evts_deex += 1

            elif (N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha ==
                    N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == 1 and N_p == 3 and N_gamma > 0):
                # C10* -> Li6 + n + 3p + Y*gamma:
                n_c10_li6_n_3p_Ygamma += 1
                N_evts_deex += 1

            else:
                if N_p < 0:
                    N_evts_incorrect += 1
                else:
                    print("new channel: C10* -> Li6 + ..., event = {0:d}".format(event))

        else:
            if (N_piminus == N_piplus == N_gamma == N_triton == N_he3 == N_alpha == N_li6 ==
                    N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == 1 and N_p == 2 and N_deuteron == 1):
                # C10* -> Li5 + n + 2p + d:
                n_c10_li5_n_2p_d += 1
                N_evts_deex += 1

            elif (N_n == N_piminus == N_piplus == N_gamma == N_triton == N_alpha == N_li6 ==
                    N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == N_deuteron == N_he3 == 1) or \
                    (N_n == N_piminus == N_piplus == N_gamma == N_triton == N_he3 == N_li6 ==
                     N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                     N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == N_deuteron == N_alpha == 1):
                # C10* -> He3 + p + d + alpha:
                n_c10_he3_p_d_alpha += 1
                N_evts_deex += 1

            elif (N_n == N_piminus == N_piplus == N_gamma == N_triton == N_he3 == N_alpha == N_li6 ==
                    N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == 1 and N_deuteron == 2):
                # C10* -> Li5 + p + 2d:
                n_c10_li5_p_2d += 1
                N_evts_deex += 1

            elif (N_n == N_piminus == N_piplus == N_gamma == N_deuteron == N_he3 == N_alpha == N_li6 ==
                    N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == 1 and N_triton == 1):
                # C10* -> Be6 + p + t:
                n_c10_be6_p_t += 1
                N_evts_deex += 1

            elif (N_piminus == N_piplus == N_gamma == N_deuteron == N_triton == N_alpha == N_li6 ==
                    N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == 1 and N_p == 2 and N_he3 == 1):
                # C10* -> he3 + n + 2p + alpha:
                n_c10_he3_n_2p_alpha += 1
                N_evts_deex += 1

            elif (N_piminus == N_piplus == N_gamma == N_deuteron == N_triton == N_alpha == N_li6 ==
                    N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_p == N_he3 == 1):
                # C10* -> li5 + n + p + he3:
                n_c10_li5_n_p_he3 += 1
                N_evts_deex += 1

            elif (N_piminus == N_piplus == N_gamma == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                    N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == 2 and N_p == 1):
                # C10* -> B7 + 2n + p:
                n_c10_b7_2n_p += 1
                N_evts_deex += 1

            elif (N_piminus == N_piplus == N_gamma == N_triton == N_he3 == N_alpha == N_li6 ==
                    N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_p == N_deuteron == 1):
                # C10* -> Be6 + n + p + d:
                n_c10_be6_n_p_d += 1
                N_evts_deex += 1

            elif (N_n == N_p == N_piminus == N_piplus == N_gamma == N_triton == N_alpha == N_li6 ==
                    N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_deuteron == N_he3 == 1):
                # C10* -> Li5 + d + he3:
                n_c10_li5_d_he3 += 1
                N_evts_deex += 1

            elif (N_n == N_piminus == N_piplus == N_gamma == N_deuteron == N_triton == N_he3 == N_li6 ==
                    N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == N_alpha == 1):
                # C10* -> li5 + p + alpha:
                n_c10_li5_p_alpha += 1
                N_evts_deex += 1

            elif (N_piminus == N_piplus == N_gamma == N_deuteron == N_triton == N_he3 == N_li6 ==
                    N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_p == N_alpha == 1):
                # C10* -> li4 + n + p + alpha:
                n_c10_li4_n_p_alpha += 1
                N_evts_deex += 1

            else:
                if N_p < 0:
                    N_evts_incorrect += 1

                else:
                    # more deexcitation channels C10* -> ...
                    # print("new channel: C10* -> ..., event = {0:d}".format(event))
                    n_c10_rest += 1
                    N_evts_deex += 1

    elif isopdg_generator == 1000050110:
        number_b11_iso += 1
        # B11 deexcitation:
        if N_c11 == 1:
            print("new channel: B11* -> C11 + ..., event = {0:d}".format(event))

        elif N_c10 == 1:
            print("new channel: B11* -> C10 + ..., event = {0:d}".format(event))

        elif N_c9 == 1:
            print("new channel: B11* -> C9 + ..., event = {0:d}".format(event))

        elif N_b11 == 1:
            if (N_n == N_p == N_piminus == N_piplus == N_gamma == N_deuteron == N_triton == N_he3 == N_alpha
                    == N_li6 ==
                    N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0):
                # no deexcitation:
                n_b11_noDeex += 1

            elif (N_n == N_p == N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                    N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_gamma > 0):
                # B11* -> B11 + Y*gamma:
                n_b11_b11_Ygamma += 1
                N_evts_deex += 1

            else:
                if N_p < 0:
                    N_evts_incorrect += 1
                else:
                    print("new channel: B11* -> B11 + ..., event = {0:d}".format(event))

        elif N_b10 == 1:
            if (N_p == N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                    N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == 1 and N_gamma > 0):
                # B11* -> B10 + n + 2gamma:
                n_b11_b10_n_Ygamma += 1
                N_evts_deex += 1

            else:
                if N_p < 0:
                    N_evts_incorrect += 1
                else:
                    print("new channel: B11* -> B10 + ..., event = {0:d}".format(event))

        elif N_b9 == 1:
            if (N_p == N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                    N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == 2 and N_gamma > 0):
                # B11* -> B8 + 2n + 2gamma:
                n_b11_b9_2n_Ygamma += 1
                N_evts_deex += 1

            else:
                if N_p < 0:
                    N_evts_incorrect += 1
                else:
                    print("new channel: B11* -> B9 + ..., event = {0:d}".format(event))

        elif N_b8 == 1:
            print("new channel: B11* -> B8 + ..., event = {0:d}".format(event))

        elif N_be10 == 1:
            if (N_n == N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                    N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == 1 and N_gamma > 0):
                # B11* -> Be10 + p + Y*gamma:
                n_b11_be10_p_Ygamma += 1
                N_evts_deex += 1

            else:
                if N_p < 0:
                    N_evts_incorrect += 1
                else:
                    print("new channel: B11* -> Be10 + ..., event = {0:d}".format(event))

        elif N_be9 == 1:
            if (N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                    N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_p == 1 and N_gamma > 0):
                # B11* -> Be9 + n + p + Y*gamma:
                n_b11_be9_n_p_Ygamma += 1
                N_evts_deex += 1

            elif (N_n == N_p == N_piminus == N_piplus == N_triton == N_he3 == N_alpha == N_li6 ==
                    N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_deuteron == 1 and N_gamma > 0):
                # B11* -> Be9 + d + Y*gamma:
                n_b11_be9_d_Ygamma += 1
                N_evts_deex += 1

            else:
                if N_p < 0:
                    N_evts_incorrect += 1
                else:
                    print("new channel: B11* -> Be9 + ..., event = {0:d}".format(event))

        elif N_be8 == 1:
            if (N_p == N_piminus == N_piplus == N_triton == N_he3 == N_alpha == N_li6 ==
                    N_li7 == N_li8 == N_li9 == N_be7 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_deuteron == 1 and N_gamma > 0):
                # B11* -> Be8 + n + d + Y*gamma:
                n_b11_be8_n_d_Ygamma += 1
                N_evts_deex += 1

            elif (N_n == N_p == N_piminus == N_piplus == N_deuteron == N_he3 == N_alpha == N_li6 ==
                    N_li7 == N_li8 == N_li9 == N_be7 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_triton == 1 and N_gamma > 0):
                # B11* -> Be8 + t + Y*gamma:
                n_b11_be8_t_Ygamma += 1
                N_evts_deex += 1

            elif (N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                    N_li7 == N_li8 == N_li9 == N_be7 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == 2 and N_p == 1 and N_gamma > 0):
                # B11* -> Be8 + 2n + p + Y*gamma:
                n_b11_be8_2n_p_Ygamma += 1
                N_evts_deex += 1

            else:
                if N_p < 0:
                    N_evts_incorrect += 1
                else:
                    print("new channel: B11* -> Be8 + ..., event = {0:d}".format(event))

        elif N_be7 == 1:
            print("new channel: B11* -> Be7 + ..., event = {0:d}".format(event))

        elif N_li9 == 1:
            print("new channel: B11* -> Li9 + ..., event = {0:d}".format(event))

        elif N_li8 == 1:
            print("new channel: B11* -> Li8 + ..., event = {0:d}".format(event))

        elif N_li7 == 1:
            if (N_n == N_p == N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_li6 ==
                    N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_alpha == 1 and N_gamma > 0):
                # B11* -> Li7 + alpha + Y*gamma:
                n_b11_Li7_alpha_Ygamma += 1
                N_evts_deex += 1

            else:
                if N_p < 0:
                    N_evts_incorrect += 1
                else:
                    print("new channel: B11* -> Li7 + ..., event = {0:d}".format(event))

        elif N_li6 == 1:
            if (N_p == N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_li7 ==
                    N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_alpha == 1 and N_gamma > 0):
                # B11* -> Li6 + n + alpha + Y*gamma:
                n_b11_Li6_n_alpha_Ygamma += 1
                N_evts_deex += 1

            else:
                if N_p < 0:
                    N_evts_incorrect += 1
                else:
                    print("new channel: B11* -> Li6 + ..., event = {0:d}".format(event))

        else:
            if (N_n == N_p == N_piminus == N_piplus == N_gamma == N_triton == N_he3 == N_li6 == N_li7 ==
                    N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_deuteron == N_alpha == 1):
                # B11* -> He5 + d + alpha:
                n_b11_he5_d_alpha += 1
                N_evts_deex += 1

            elif (N_piminus == N_piplus == N_gamma == N_deuteron == N_triton == N_he3 == N_li6 == N_li7 ==
                    N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_p == N_alpha == 1):
                # B11* -> He5 + n + p + alpha:
                n_b11_he5_n_p_alpha += 1
                N_evts_deex += 1

            elif (N_n == N_piminus == N_piplus == N_gamma == N_deuteron == N_triton == N_he3 == N_li6
                  == N_li7 ==
                    N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == N_alpha == 1):
                # B11* -> He6 + p + alpha:
                n_b11_he6_p_alpha += 1
                N_evts_deex += 1

            elif (N_n == N_p == N_piminus == N_piplus == N_gamma == N_deuteron == N_he3 == N_li6 == N_li7 ==
                    N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_triton == N_alpha == 1) or \
                    (N_n == N_p == N_piminus == N_piplus == N_gamma == N_deuteron == N_triton == N_he3
                     == N_li6 == N_li7 ==
                     N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                     N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_alpha == 2):
                # B11* -> t + 2*alpha:
                n_b11_t_2alpha += 1
                N_evts_deex += 1

            else:
                if N_p < 0:
                    N_evts_incorrect += 1
                else:
                    print("new channel: B11* -> ..., event = {0:d}".format(event))

    elif isopdg_generator == 1000050100:
        number_b10_iso += 1
        # B10 deexcitation:
        if N_c11 == 1:
            print("new channel: B10* -> C11 + ..., event = {0:d}".format(event))

        elif N_c10 == 1:
            print("new channel: B10* -> C10 + ..., event = {0:d}".format(event))

        elif N_c9 == 1:
            print("new channel: B10* -> C9 + ..., event = {0:d}".format(event))

        elif N_b11 == 1:
            print("new channel: B10* -> B11 + ..., event = {0:d}".format(event))

        elif N_b10 == 1:
            if (N_n == N_p == N_piminus == N_piplus == N_gamma == N_deuteron == N_triton == N_he3 == N_alpha
                    == N_li6 ==
                    N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0):
                # no deexcitation:
                n_b10_noDeex += 1

            elif (N_n == N_p == N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                    N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_gamma > 0):
                # B10* -> B10 + Y*gamma:
                n_b10_b10_Ygamma += 1
                N_evts_deex += 1

            else:
                if N_p < 0:
                    N_evts_incorrect += 1
                else:
                    print("new channel: B10* -> B10 + ..., event = {0:d}".format(event))

        elif N_b9 == 1:
            if (N_p == N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                    N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == 1 and N_gamma > 0):
                # B10* -> B9 + n + Y*gamma:
                n_b10_b9_n_Ygamma += 1
                N_evts_deex += 1

            else:
                if N_p < 0:
                    N_evts_incorrect += 1
                else:
                    print("new channel: B10* -> B9 + ..., event = {0:d}".format(event))

        elif N_b8 == 1:
            print("new channel: B10* -> B8 + ..., event = {0:d}".format(event))

        elif N_be10 == 1:
            print("new channel: B10* -> Be10 + ..., event = {0:d}".format(event))

        elif N_be9 == 1:
            if (N_n == N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                    N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == 1 and N_gamma > 0):
                # B10* -> Be9 + p + Y*gamma:
                n_b10_be9_p_Ygamma += 1
                N_evts_deex += 1

            else:
                if N_p < 0:
                    N_evts_incorrect += 1
                else:
                    print("new channel: B10* -> Be9 + ..., event = {0:d}".format(event))

        elif N_be8 == 1:
            if (N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                    N_li7 == N_li8 == N_li9 == N_be7 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_p == 1 and N_gamma > 0):
                # B10* -> Be8 + n + p + Y*gamma:
                n_b10_be8_n_p_Ygamma += 1
                N_evts_deex += 1

            elif (N_n == N_p == N_piminus == N_piplus == N_triton == N_he3 == N_alpha == N_li6 ==
                    N_li7 == N_li8 == N_li9 == N_be7 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_deuteron == 1 and N_gamma > 0):
                # B10* -> Be8 + d + X*pi0 + Y*gamma:
                n_b10_be8_d_Ygamma += 1
                N_evts_deex += 1

            else:
                if N_p < 0:
                    N_evts_incorrect += 1
                else:
                    print("new channel: B10* -> Be8 + ..., event = {0:d}".format(event))

        elif N_be7 == 1:
            if (N_n == N_p == N_piminus == N_piplus == N_deuteron == N_he3 == N_alpha == N_li6 ==
                    N_li7 == N_li8 == N_li9 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_triton == 1 and N_gamma > 0):
                # B10* -> Be7 + t + Y*gamma:
                n_b10_be7_t_Ygamma += 1
                N_evts_deex += 1

            elif (N_p == N_piminus == N_piplus == N_triton == N_he3 == N_alpha == N_li6 ==
                    N_li7 == N_li8 == N_li9 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_deuteron == 1 and N_gamma > 0):
                # B10* -> Be7 + n + d + Y*gamma:
                n_b10_be7_n_d_Ygamma += 1
                N_evts_deex += 1

            else:
                if N_p < 0:
                    N_evts_incorrect += 1
                else:
                    print("new channel: B10* -> Be7 + ..., event = {0:d}".format(event))

        elif N_li9 == 1:
            print("new channel: B10* -> Li9 + ..., event = {0:d}".format(event))

        elif N_li8 == 1:
            print("new channel: B10* -> Li8 + ..., event = {0:d}".format(event))

        elif N_li7 == 1:
            if (N_n == N_p == N_piminus == N_piplus == N_deuteron == N_triton == N_alpha == N_li6 ==
                    N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_he3 == 1 and N_gamma > 0):
                # B10* -> Li7 + he3 + Y*gamma:
                n_b10_li7_he3_Ygamma += 1
                N_evts_deex += 1

            elif (N_n == N_piminus == N_piplus == N_triton == N_he3 == N_alpha == N_li6 ==
                    N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == N_deuteron == 1 and N_gamma > 0):
                # B10* -> Li7 + p + d + Y*gamma:
                n_b10_li7_p_d_Ygamma += 1
                N_evts_deex += 1

            else:
                if N_p < 0:
                    N_evts_incorrect += 1
                else:
                    print("new channel: B10* -> Li7 + ..., event = {0:d}".format(event))

        elif N_li6 == 1:
            if (N_n == N_p == N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_li7 ==
                    N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_alpha == 1 and N_gamma > 0):
                # B10* -> Li6 + alpha + Y*gamma:
                n_b10_li6_alpha_Ygamma += 1
                N_evts_deex += 1

            elif (N_n == N_piminus == N_piplus== N_deuteron == N_he3 == N_alpha == N_li7 ==
                    N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == N_triton == 1 and N_gamma > 0):
                # B10* -> Li6 + p + t + Y*gamma:
                n_b10_li6_p_t_Ygamma += 1
                N_evts_deex += 1

            elif (N_p == N_piminus == N_piplus == N_deuteron == N_triton == N_alpha == N_li7 ==
                    N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_he3 == 1 and N_gamma > 0):
                # B10* -> Li6 + n + he3 + Y*gamma:
                n_b10_li6_n_he3_Ygamma += 1
                N_evts_deex += 1

            else:
                if N_p < 0:
                    N_evts_incorrect += 1
                else:
                    print("new channel: B10* -> Li6 + ..., event = {0:d}".format(event))

        else:
            if (N_n == N_piminus == N_piplus == N_gamma == N_deuteron == N_triton == N_he3 == N_li6 == N_li7 ==
                    N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == N_alpha == 1):
                # B10* -> He5 + p + alpha:
                n_b10_he5_p_alpha += 1
                N_evts_deex += 1

            elif (N_p == N_piminus == N_piplus == N_gamma == N_deuteron == N_triton == N_he3 == N_li6
                  == N_li7 ==
                    N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_alpha == 1):
                # B10* -> Li5 + n + alpha:
                n_b10_li5_n_alpha += 1
                N_evts_deex += 1

            elif (N_n == N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_li6 == N_li7 ==
                    N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == N_alpha == 1 and N_gamma > 0):
                # B10* -> He5 + p + alpha + Y*gamma:
                n_b10_he5_p_alpha_Ygamma += 1
                N_evts_deex += 1

            else:
                if N_p < 0:
                    N_evts_incorrect += 1
                else:
                    print("new channel: B10* -> ..., event = {0:d}".format(event))

    elif isopdg_generator == 1000040100:
        number_be10_iso += 1
        # Be10 deexcitation:
        if N_c11 == 1:
            print("new channel: Be10* -> C11 + ..., event = {0:d}".format(event))

        elif N_c10 == 1:
            print("new channel: Be10* -> C10 + ..., event = {0:d}".format(event))

        elif N_c9 == 1:
            print("new channel: Be10* -> C9 + ..., event = {0:d}".format(event))

        elif N_b11 == 1:
            print("new channel: Be10* -> B11 + ..., event = {0:d}".format(event))

        elif N_b10 == 1:
            print("new channel: Be10* -> B10 + ..., event = {0:d}".format(event))

        elif N_b9 == 1:
            print("new channel: Be10* -> B9 + ..., event = {0:d}".format(event))

        elif N_b8 == 1:
            print("new channel: Be10* -> B8 + ..., event = {0:d}".format(event))

        elif N_be10 == 1:
            if (N_n == N_p == N_piminus == N_piplus == N_gamma == N_deuteron == N_triton == N_he3 == N_alpha
                    == N_li6 ==
                    N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0):
                # no deexcitation:
                n_be10_noDeex += 1

            elif (N_n == N_p == N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha
                    == N_li6 ==
                    N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_gamma > 0):
                # Be10* -> Be10 + Y*gamma:
                n_be10_be10_Ygamma += 1
                N_evts_deex += 1

            else:
                if N_p < 0:
                    N_evts_incorrect += 1
                else:
                    print("new channel: Be10* -> Be10 + ..., event = {0:d}".format(event))

        elif N_be9 == 1:
            if (N_p == N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                    N_li7 == N_li8 == N_li9 == N_be7 == N_be8 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == 1 and N_gamma > 0):
                # Be10* -> Be9 + n + Y*gamma:
                n_be10_be9_n_Ygamma += 1
                N_evts_deex += 1

            else:
                if N_p < 0:
                    N_evts_incorrect += 1
                else:
                    print("new channel: Be10* -> Be9 + ..., event = {0:d}".format(event))

        elif N_be8 == 1:
            if (N_p == N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                    N_li7 == N_li8 == N_li9 == N_be7 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == 2 and N_gamma > 0):
                # Be10* -> Be8 + 2n + Y*gamma:
                n_be10_be8_2n_Ygamma += 1
                N_evts_deex += 1

            else:
                if N_p < 0:
                    N_evts_incorrect += 1
                else:
                    print("new channel: Be10* -> Be8 + ..., event = {0:d}".format(event))

        elif N_be7 == 1:
            print("new channel: Be10* -> Be7 + ..., event = {0:d}".format(event))

        elif N_li9 == 1:
            if (N_n == N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                    N_li7 == N_li8 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == 1 and N_gamma > 0):
                # Be10* -> Li9 + p + Y*gamma:
                n_be10_li9_p_Ygamma += 1
                N_evts_deex += 1

            else:
                if N_p < 0:
                    N_evts_incorrect += 1
                else:
                    print("new channel: Be10* -> Li9 + ..., event = {0:d}".format(event))

        elif N_li8 == 1:
            if (N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                    N_li7 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_p == 1 and N_gamma > 0):
                # Be10* -> Li8 + n + p + Y*gamma:
                n_be10_li8_n_p_Ygamma += 1
                N_evts_deex += 1

            elif (N_n == N_p == N_piminus == N_piplus == N_triton == N_he3 == N_alpha == N_li6 ==
                    N_li7 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_deuteron == 1 and N_gamma > 0):
                # Be10* -> Li8 + d + Y*gamma:
                n_be10_li8_d_Ygamma += 1
                N_evts_deex += 1

            else:
                if N_p < 0:
                    N_evts_incorrect += 1
                else:
                    print("new channel: Be10* -> Li8 + ..., event = {0:d}".format(event))

        elif N_li7 == 1:
            if (N_p == N_piminus == N_piplus == N_triton == N_he3 == N_alpha == N_li6 ==
                    N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_deuteron == 1 and N_gamma > 0):
                # Be10* -> Li7 + n + d + Y*gamma:
                n_be10_li7_n_d_Ygamma += 1
                N_evts_deex += 1

            elif (N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha == N_li6 ==
                    N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == 2 and N_p == 1 and N_gamma > 0):
                # Be10* -> Li7 + 2n + p + Y*gamma:
                n_be10_li7_2n_p_Ygamma += 1
                N_evts_deex += 1

            elif (N_n == N_p == N_piminus == N_piplus == N_deuteron == N_he3 == N_alpha == N_li6 ==
                    N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_triton == 1 and N_gamma > 0):
                # Be10* -> Li7 + t + Y*gamma:
                n_be10_li7_t_Ygamma += 1
                N_evts_deex += 1

            else:
                if N_p < 0:
                    N_evts_incorrect += 1
                else:
                    print("new channel: Be10* -> Li7 + ..., event = {0:d}".format(event))

        elif N_li6 == 1:
            if (N_p == N_piminus == N_piplus == N_deuteron == N_he3 == N_alpha == N_li7 ==
                    N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_triton == 1 and N_gamma > 0):
                # Be10* -> Li6 + n + t + Y*gamma:
                n_be10_li6_n_t_Ygamma += 1
                N_evts_deex += 1

            elif (N_piminus == N_piplus == N_deuteron == N_triton == N_he3 == N_alpha == N_li7 ==
                    N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == 3 and N_p == 1 and N_gamma > 0):
                # Be10* -> Li6 + 3n + p + Y*gamma:
                n_be10_li6_3n_p_Ygamma += 1
                N_evts_deex += 1

            elif (N_p == N_piminus == N_piplus == N_triton == N_he3 == N_alpha == N_li7 ==
                    N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == 2 and N_deuteron == 1 and N_gamma > 0):
                # Be10* -> Li6 + 2n + d + Y*gamma:
                n_be10_li6_2n_d_Ygamma += 1
                N_evts_deex += 1

            else:
                if N_p < 0:
                    N_evts_incorrect += 1
                else:
                    print("new channel: Be10* -> Li6 + ..., event = {0:d}".format(event))

        else:
            if (N_piminus == N_piplus == N_gamma == N_deuteron == N_triton == N_he3 == N_li6 == N_li7 ==
                    N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_p == N_alpha == 1):
                # Be10* -> H4 + n + p + alpha:
                n_be10_h4_n_p_alpha += 1
                N_evts_deex += 1

            elif (N_n == N_piminus == N_piplus == N_gamma == N_triton == N_he3 == N_alpha == N_li6 == N_li7 ==
                    N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_p == N_deuteron == 1):
                # Be10* -> He7 + p + d:
                n_be10_he7_p_d += 1
                N_evts_deex += 1

            elif (N_piminus == N_piplus == N_gamma == N_triton == N_he3 == N_alpha == N_li6 == N_li7 ==
                    N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == 2 and N_p == N_deuteron == 1):
                # Be10* -> He5 + 2n + p + d:
                n_be10_he5_2n_p_d += 1
                N_evts_deex += 1

            elif (N_p == N_piminus == N_piplus == N_gamma == N_triton == N_he3 == N_alpha == N_li6 == N_li7 ==
                    N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == 1 and N_deuteron == 2):
                # Be10* -> He5 + n + 2d:
                n_be10_he5_n_2d += 1
                N_evts_deex += 1

            elif (N_n == N_p == N_piminus == N_piplus == N_gamma == N_deuteron == N_triton == N_alpha
                  == N_li6 == N_li7 ==
                    N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_he3 == 1):
                # Be10* -> He7 + He3:
                n_be10_he7_he3 += 1
                N_evts_deex += 1

            elif (N_piminus == N_piplus == N_gamma == N_deuteron == N_he3 == N_alpha == N_li6 == N_li7 ==
                    N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_p == N_triton == 1):
                # Be10* -> He5 + n + p + t:
                n_be10_he5_n_p_t += 1
                N_evts_deex += 1

            elif (N_p == N_piminus == N_piplus == N_gamma == N_he3 == N_alpha == N_li6 == N_li7 ==
                    N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_deuteron == N_triton == 1) or \
                    (N_p == N_piminus == N_piplus == N_gamma == N_triton == N_he3 == N_li6 == N_li7 ==
                     N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                     N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_deuteron == N_alpha == 1):
                # Be10* -> alpha + n + d + t:
                n_be10_alpha_n_d_t += 1
                N_evts_deex += 1

            elif (N_piminus == N_piplus == N_gamma == N_triton == N_he3 == N_alpha == N_li6 == N_li7 ==
                    N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_n == N_p == N_deuteron == 1):
                # Be10* -> He6 + n + p + d:
                n_be10_he6_n_p_d += 1
                N_evts_deex += 1

            elif (N_n == N_p == N_piminus == N_piplus == N_gamma == N_triton == N_he3 == N_alpha == N_li6
                  == N_li7 ==
                    N_li8 == N_li9 == N_be7 == N_be8 == N_be9 == N_be10 == N_b8 == N_b9 == N_b10 == N_b11 ==
                    N_c9 == N_c10 == N_c11 == N_c12 == 0 and N_deuteron == 2):
                # Be10* -> He6 + 2d:
                n_be10_he6_2d += 1
                N_evts_deex += 1

            else:
                if N_p < 0:
                    N_evts_incorrect += 1
                else:
                    # other deexcitation channels:
                    # print("new channel: Be10* -> ..., event = {0:d}".format(event))
                    n_be10_rest += 1
                    N_evts_deex += 1

                    # print("--------------------")
                    # print(N_n)
                    # print(N_p)
                    # print(N_piminus)
                    # print(N_piplus)
                    # print(N_pi0)
                    # print(N_gamma)
                    # print(N_deuteron)
                    # print(N_triton)
                    # print(N_he3)
                    # print(N_alpha)
                    # print(N_li6)
                    # print(N_li7)
                    # print(N_li8)
                    # print(N_li9)
                    # print(N_be7)
                    # print(N_be8)
                    # print(N_be9)
                    # print(N_be10)
                    # print(N_b8)
                    # print(N_b9)
                    # print(N_b10)
                    # print(N_b11)
                    # print(N_c9)
                    # print(N_c10)
                    # print(N_c11)
                    # print(N_c12)

print(N_entries_generator)
print(N_evts_deex)
print(N_evts_incorrect)
print("")
print(number_c11_iso)
print(n_c11_noDeex)
print(n_c11_c11_Ygamma)
print(n_c11_c10_n_Ygamma)
print(n_c11_b10_p_Ygamma)
print(n_c11_b9_n_p_Ygamma)
print(n_c11_b9_d_Ygamma)
print(n_c11_be9_2p_Ygamma)
print(n_c11_be8_p_d_Ygamma)
print(n_c11_be8_he3_Ygamma)
print(n_c11_be8_n_2p_Ygamma)
print(n_c11_be7_alpha_Ygamma)
print(n_c11_li6_p_alpha_Ygamma)
print(n_c11_li5_d_alpha)
print(n_c11_li5_n_p_alpha)
print(n_c11_li5_d_alpha_Ygamma)
print(n_c11_he3_2alpha)
print("")
print(number_c10_iso)
print(n_c10_noDeex)
print(n_c10_c10_Ygamma)
print(n_c10_c9_n_Ygamma)
print(n_c10_b9_p_Ygamma)
print(n_c10_b8_n_p_Ygamma)
print(n_c10_b8_d_Ygamma)
print(n_c10_be8_2p_Ygamma)
print(n_c10_be7_p_d_Ygamma)
print(n_c10_be7_n_2p_Ygamma)
print(n_c10_be7_he3_Ygamma)
print(n_c10_li6_p_he3_Ygamma)
print(n_c10_li6_2p_d_Ygamma)
print(n_c10_li6_n_3p_Ygamma)
print(n_c10_li5_n_2p_d)
print(n_c10_he3_p_d_alpha)
print(n_c10_li5_p_2d)
print(n_c10_be6_p_t)
print(n_c10_he3_n_2p_alpha)
print(n_c10_li5_n_p_he3)
print(n_c10_b7_2n_p)
print(n_c10_be6_n_p_d)
print(n_c10_li5_d_he3)
print(n_c10_li5_p_alpha)
print(n_c10_li4_n_p_alpha)
print(n_c10_rest)
print("")
print(number_b11_iso)
print(n_b11_noDeex)
print(n_b11_b11_Ygamma)
print(n_b11_b10_n_Ygamma)
print(n_b11_b9_2n_Ygamma)
print(n_b11_be10_p_Ygamma)
print(n_b11_be9_n_p_Ygamma)
print(n_b11_be9_d_Ygamma)
print(n_b11_be8_n_d_Ygamma)
print(n_b11_be8_t_Ygamma)
print(n_b11_be8_2n_p_Ygamma)
print(n_b11_Li7_alpha_Ygamma)
print(n_b11_Li6_n_alpha_Ygamma)
print(n_b11_he5_d_alpha)
print(n_b11_he5_n_p_alpha)
print(n_b11_he6_p_alpha)
print(n_b11_t_2alpha)
print("")
print(number_b10_iso)
print(n_b10_noDeex)
print(n_b10_b10_Ygamma)
print(n_b10_b9_n_Ygamma)
print(n_b10_be9_p_Ygamma)
print(n_b10_be8_n_p_Ygamma)
print(n_b10_be8_d_Ygamma)
print(n_b10_be7_t_Ygamma)
print(n_b10_be7_n_d_Ygamma)
print(n_b10_li7_he3_Ygamma)
print(n_b10_li7_p_d_Ygamma)
print(n_b10_li6_alpha_Ygamma)
print(n_b10_li6_p_t_Ygamma)
print(n_b10_li6_n_he3_Ygamma)
print(n_b10_he5_p_alpha)
print(n_b10_li5_n_alpha)
print(n_b10_he5_p_alpha_Ygamma)
print("")
print(number_be10_iso)
print(n_be10_noDeex)
print(n_be10_be10_Ygamma)
print(n_be10_be9_n_Ygamma)
print(n_be10_be8_2n_Ygamma)
print(n_be10_li9_p_Ygamma)
print(n_be10_li8_n_p_Ygamma)
print(n_be10_li8_d_Ygamma)
print(n_be10_li7_n_d_Ygamma)
print(n_be10_li7_2n_p_Ygamma)
print(n_be10_li7_t_Ygamma)
print(n_be10_li6_n_t_Ygamma)
print(n_be10_li6_3n_p_Ygamma)
print(n_be10_li6_2n_d_Ygamma)
print(n_be10_h4_n_p_alpha)
print(n_be10_he7_p_d)
print(n_be10_he5_2n_p_d)
print(n_be10_he5_n_2d)
print(n_be10_he7_he3)
print(n_be10_he5_n_p_t)
print(n_be10_alpha_n_d_t)
print(n_be10_he6_n_p_d)
print(n_be10_he6_2d)
print(n_be10_rest)


