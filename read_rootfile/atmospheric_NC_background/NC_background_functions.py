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
import glob
import numpy as np
from matplotlib import pyplot as plt

"""
class NCData:
    def __init__(self):  # this method creates the class object.
        # event ID (starts from 0):
        self.eventID = 0
        # PDG ID of the projectile (i.e. which neutrino is interacting):
        self.projectilePDG = 0
        # energy of the incoming neutrino in GeV:
        self.projectileEnergy = 0
        # PDG ID of the target particle (either C12 or proton):
        self.targetPDG = 0
        # Channel ID of the NC interaction, represents which particles are produced via the NC interaction:
        self.NCinteration_ch_ID = 0
        # Channel ID of the deexcitation, represents which particles are produced via the deexication of the produced
        # excited isotope:
        self.deexcitation_ID = 0
        # PDG ID of the isotope after the NC interaction, BUT before the deexcitation:
        self.isotopePDG = 0
        # number of final particles after NC interactions and deexcitation:
        self.Nparticles = 0
        # PDG ID of the final particles. It is an array of length "Nparticles":
        self.finalPDG = 0
        # momentum in x-direction of the final particles. It is an array of length "Nparticles":
        self.finalPx = 0
        # momentum in y-direction of the final particles. It is an array of length "Nparticles":
        self.finalPy = 0
        # momentum in z-direction of the final particles. It is an array of length "Nparticles":
        self.finalPz = 0
"""


def read_sample_detsim_user(rootfile_input, r_cut, e_prompt_min, e_prompt_max, e_delayed_min, e_delayed_max,
                            time_cut_min, time_cut_max, distance_cut, time_resolution, number_entries_input):
    """
    function to read the sample_detsim_user.root file and to get visible energy of the prompt signal of
    the IBD-like signals.

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
    rtree_evt = rfile.Get("evt")
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
        number_events = 0
        print("ERROR: number of events in the Trees are NOT equal!!")

    # check if number_events is equal to number_entries_input (if not, the detector simulation was incorrect!!):
    if number_events != number_entries_input:
        number_events = 0
        print("ERROR: number of events are not equal to {0:d}".format(number_entries_input))
        print("-> Detector Simulation not correct!!")

    # preallocate array of visible energy of prompt signal (energy in MeV) (np.array of float):
    e_vis = np.array([])
    # preallocate array, where event ID of the IBD like events is saved (np.array of float):
    evt_id_ibd = np.array([])

    # preallocate number of events, where 2 prompt signal are added to 1 prompt signal:
    number_2promptadded = 0
    # preallocate number of events, where there are 3 or more prompt signal to 1 corresponding delayed signal:
    number_3prompt_1delayed = 0
    # preallocate number of events, where there are 2 prompt and only 1 corresponding delayed signal:
    number_check1 = 0
    # preallocate number of events, where there are 2 possible IBD signals:
    number_check2 = 0
    # preallocate number of events, where there are more than 3 prompt and more than 3 corresponding delayed signal:
    number_check3 = 0

    # loop over every event, i.e. every event, in the TTree:
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
            evt_id = 0
            print("ERROR: event ID in the Trees are NOT equal!!")

        # check number of initial particles of the Trees:
        if n_par_prmtrkdep == n_par_geninfo:
            n_par = n_par_geninfo
        else:
            n_par = 0
            print("ERROR: number of initial particles in the Trees are NOT equal!!")

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
            """ 1 possible prompt AND 1 possible delayed signal: """
            # only one prompt and one delayed signal. Get first entry of the array:
            index_p = index_prompt[0]
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
            """ more than 1 possible prompt signals, BUT only 1 possible delayed signal: """
            # preallocate value that represents the number of IBD-like signals in this event:
            number_ibd_evts = 0

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

                            # increment number of IBD-lie signals in this event:
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
                    print("ERROR: Number of IBD-like events different to length of 'array_e_vis' (evt_ID = {0:d})"
                          .format(evt_id))

                # append evt_id of the IBD like signal to evt_id_ibd array:
                evt_id_ibd = np.append(evt_id_ibd, evt_id)

                # append Qedep of the prompt signal to the e_vis array:
                e_vis = np.append(e_vis, array_e_vis[0])

            elif number_ibd_evts == 0:
                continue

            elif number_ibd_evts == 2:
                # 2 possible prompt signal and only 1 corresponding delayed signals:
                # calculate the time difference between the prompt signals in ns:
                delta_t_pp = np.absolute(time_exit[int(array_prompt_index[0])] - time_exit[int(array_prompt_index[1])])

                if time_resolution <= delta_t_pp <= time_cut_max:
                    # this means, one prompt signal lies between the other prompt and the delayed signal and can be
                    # separated from the other prompt signal. -> NO IBD-like signal -> go to next event:
                    continue

                elif delta_t_pp < time_resolution:
                    # this means, that the two prompt signals can not be separated clearly -> treat the two prompt
                    # signals as one prompt signal with energy = E_prompt1 + E_prompt2:

                    # append evt_id of the IBD like signal to evt_id_ibd array:
                    evt_id_ibd = np.append(evt_id_ibd, evt_id)

                    # append sum of Qedep of first prompt signal and Qedep of second prompt signal to e_vis array:
                    e_vis = np.append(e_vis, array_e_vis[0] + array_e_vis[1])

                    # check:
                    number_2promptadded = number_2promptadded + 1
                    print(evt_id)

                else:
                    # this means, one prompt signal is after the delayed signal (should be not possible):
                    print("WARNING in len(index_prompt) > 1 and len(index_delayed) == 1 ------ 1 prompt after delayed "
                          "signal")

            else:
                # TODO-me: How to deal with events where there are 3 or more IBD-like events (but only one delayed sig.)
                # to estimate the number of such events, increment variable number_2prompt_1delayed:
                number_3prompt_1delayed = number_3prompt_1delayed + 1
                print("WARNING: more than 1 IBD-like signal in event {2:d}: n_IBD_like = {0:.0f}, prompt = {1:.0f}"
                      .format(number_ibd_evts, check_prompt, evt_id))
                print("-----------> not yet included!")


        elif len(index_prompt) == 1 and len(index_delayed) > 1:
            """ 1 possible prompt signal, BUT more than 1 possible delayed signals"""
            # preallocate value that represents the number of IBD-like signals in this event:
            number_ibd_evts = 0

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
                # this means, that there is no IBD-like signal in this event (either 0 or more than 1 delayed signal)
                continue

            elif len(array_delayed_index) == 1:
                # this means, that there is 1 prompt signal and 1 corresponding delayed signal
                # -> therefore 1 IBD-like signal:

                # append evt_id of the IBD like signal to evt_id_ibd array:
                evt_id_ibd = np.append(evt_id_ibd, evt_id)

                # get the index of the prompt signal in the event:
                which_index = np.where(array_number_ibd == 1)[0]
                # check if there is just one entry equal to 1 in the array_delayed_index:
                if len(which_index) == 1:
                    correct_prompt_index = which_index[0]

                    # append Qedep of this prompt signal to the e_vis array:
                    e_vis = np.append(e_vis, e_qdep[index_prompt[correct_prompt_index]])
                else:
                    print("ERROR in line 528")
                    continue

            elif len(array_delayed_index) == 2:
                # this means, that there are 2 possible IBD-like signals (2 prompt and 2 corresponding delayed signals)
                # there are 2 possibilities:
                # 1. the 'two' delayed signals could be the same. Then: 2 prompt signal and 1 corresponding delayed
                # signal.
                # 2. two different delayed signals. Therefore 2 possible IBD-like signals.
                if array_delayed_index[0] == array_delayed_index[1]:
                    # possibility 1:
                    # to estimate the number of such events, increment variable number_check1:
                    number_check1 = number_check1 + 1
                else:
                    # possibility 2:
                    # to estimate the number of such events, increment variable number_check2:
                    number_check2 = number_check2 + 1


            else:
                # TODO-me: How to deal with events where there are 2 or more IBD-like events (but only one delayed sig.)
                # to estimate the number of such events, increment variable number_check3:
                number_check3 = number_check3 + 1
                print("WARNING: More than 1 possible prompt signal to only 1 possible delayed signal: evt ID = {0:d}"
                      .format(evt_id))
                print("----------> not yet included!!!!")


    if number_2promptadded != 0:
        print("number of events, where 2 prompt signals are added to 1 = {0:d}".format(number_2promptadded))

    if number_3prompt_1delayed != 0:
        print("number of events, with more than 1 prompt and just 1 delayed signal = {0:d}"
              .format(number_3prompt_1delayed))

    if number_check1 != 0:
        print("number of events, with 2 prompt and just 1 delayed signal (prompt>1 AND delayed>1) = {0:d}"
              .format(number_check1))

    if number_check2 != 0:
        print("number of events, with 2 possible delayed signals (prompt>1 AND delayed>1) = {0:d}"
              .format(number_check2))

    if number_check3 != 0:
        print("number of events, with 3 prompt and 3 or less delayed signal (prompt>1 AND delayed>1) = {0:d}"
              .format(number_check3))

    return number_events, evt_id_ibd, e_vis


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
    rfile_output = ROOT.TFile(path_output + "genie_data.root", "recreate")
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
        if nc == 1 and tgt == 1000060120:

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


    return number_events, \
           frac_c12_b11_p, frac_c12_b11_n_piplus, frac_c12_b11_n_piminus_2piplus, frac_c12_b11_p_piminus_piplus, \
           frac_c12_b11_p_2piminus_2piplus, frac_c12_b11_piplus, \
           frac_c12_c11_n, frac_c12_c11_p_piminus, frac_c12_c11_n_piminus_piplus, frac_c12_c11_p_2piminus_piplus, \
           frac_c12_c11_p_3piminus_2piplus, frac_c12_c11_n_2piminus_2piplus, \
           frac_c12_b10_p_n, frac_c12_b10_2p_piminus, frac_c12_b10_p_n_piminus_piplus, frac_c12_b10_2n_piplus, \
           frac_c12_b10_2n_piminus_2piplus, frac_c12_b10_2p_2piminus_piplus, frac_c12_b10_2p_3piminus_2piplus, \
           frac_c12_b10_p_n_2piminus_2piplus, \
           frac_c12_c10_2n, frac_c12_c10_p_n_piminus, frac_c12_c10_p_n_2piminus_piplus, frac_c12_c10_2n_piminus_piplus,\
           frac_c12_c10_2p_2piminus, \
           frac_c12_be10_2p, frac_c12_be10_p_n_piplus, frac_c12_be10_p_n_piminus_2piplus, \
           frac_c12_be10_2p_piminus_piplus, frac_c12_be10_2n_2piplus, frac_c12_be10_p_n_2piminus_3piplus, \
           frac_c12_be10_2p_2piminus_2piplus, frac_c12_be10_2p_3piminus_3piplus, \
           frac_c12_b9_p_2n, frac_c12_b9_p_2n_piminus_piplus, frac_c12_b9_2p_n_3piminus_2piplus, \
           frac_c12_b9_2p_n_piminus, frac_c12_b9_3n_piplus, frac_c12_b9_p_2n_2piminus_2piplus, \
           frac_c12_b9_2p_n_2piminus_piplus, \
           frac_c12_be9_2p_n, frac_c12_be9_p_2n_piplus, frac_c12_be9_3p_piminus, frac_c12_be9_p_2n_piminus_2piplus, \
           frac_c12_be9_2p_n_piminus_piplus, frac_c12_be9_2p_n_3piminus_3piplus, frac_c12_be9_2p_n_2piminus_2piplus, \
           frac_c12_be9_3n_2piplus, frac_c12_be9_3p_2piminus_piplus, \
           frac_c12_be8_2p_2n, frac_c12_be8_3p_n_piminus, frac_c12_be8_p_3n_piplus, \
           frac_c12_be8_2p_2n_2piminus_2piplus, frac_c12_be8_4n_2piplus, frac_c12_be8_2p_2n_piminus_piplus, \
           frac_c12_be8_3p_n_2piminus_piplus, frac_c12_be8_4p_2piminus, \
           frac_c12_c9_p_2n_piminus, frac_c12_c9_3n, frac_c12_c9_2p_n_2piminus, frac_c12_c9_3n_2piminus_2piplus, \
           frac_c12_be7_2p_3n, frac_c12_be7_p_4n_piplus, frac_c12_be7_2p_3n_2piminus_2piplus, \
           frac_c12_be7_3p_2n_piminus, frac_c12_be7_4p_n_2piminus, frac_c12_be7_3p_2n_2piminus_piplus, \
           frac_c12_li6_3p_3n, frac_c12_li6_2p_4n_piplus, frac_c12_li6_5p_n_2piminus, \
           frac_c12_li6_2p_4n_piminus_2piplus, frac_c12_li6_4p_2n_piminus, frac_c12_li6_3p_3n_piminus_piplus, \
           frac_c12_li8_3p_n, frac_c12_li8_4p_piminus, frac_c12_li8_4p_2piminus_piplus, frac_c12_li8_2p_2n_piplus, \
           frac_c12_li8_3p_n_piminus_piplus, \
           frac_c12_li7_2p_3n_piplus, frac_c12_li7_4p_n_piminus, frac_c12_li7_3p_2n, frac_c12_li7_3p_2n_piminus_piplus,\
           frac_c12_li7_4p_n_2piminus_piplus, frac_c12_li7_2p_3n_piminus_2piplus, \
           frac_c12_b8_p_3n, frac_c12_b8_p_3n_piminus_piplus, frac_c12_b8_2p_2n_2piminus_piplus, \
           frac_c12_b8_2p_2n_piminus, frac_c12_b8_4n_piplus, \
           frac_c12_li9_2p_n_piplus, frac_c12_li9_3p, frac_c12_li9_3p_piminus_piplus, \
           frac_c12_li9_2p_n_piminus_2piplus, frac_c12_li9_p_2n_piminus_3piplus, \
           frac_c12_c8_4n, frac_c12_he8_4p, frac_c12_b7_p_4n, frac_c12_he7_4p_n, frac_c12_h7_5p, frac_c12_be6_2p_4n, \
           frac_c12_li5_3p_4n, frac_c12_li4_3p_5n, \
           frac_c12_he6_4p_2n, frac_c12_he5_4p_3n, frac_c12_he4_4p_4n, frac_c12_he3_4p_5n, frac_c12_h6_5p_n, \
           frac_c12_h5_5p_2n, frac_c12_h4_5p_3n, frac_c12_h3_5p_4n, \
           frac_c12_h2_5p_5n, frac_c12_c12,\
           frac_c12_noiso, \
           frac_no_c12, frac_es_p, frac_es_e, frac_es_o16, frac_es_n14, frac_es_s32, frac_c12_missing


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

    return event_id, projectile_pdg, projectile_energy, target_pdg, nc_interaction_ch_id, deexcitation_id, \
           isotope_pdg, n_particles, final_pdg, final_px, final_py, final_pz


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

    return energy_range, event_nu_e_incoming, event_nu_e_bar_incoming, event_nu_mu_incoming, event_nu_mu_bar_incoming, \
           event_nu_tau_incoming, event_nu_tau_bar_incoming, number_nu_e_incoming, number_nu_e_bar_incoming, \
           number_nu_mu_incoming, number_nu_mu_bar_incoming, number_nu_tau_incoming, number_nu_tau_bar_incoming, \
           fraction_nu_e, fraction_nu_e_bar, fraction_nu_mu, fraction_nu_mu_bar, fraction_nu_tau, \
           fraction_nu_tau_bar


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

    # TODO: the event rate and exposure time is NOT included yet!!!

    return energy_range, events_c12, events_proton, events_n14, events_o16, events_electron, events_s32, \
           n_c12, n_proton, n_n14, n_o16, n_electron, n_s32, \
           fraction_c12, fraction_proton, fraction_n14, fraction_o16, fraction_electron, fraction_s32


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


    return frac_c12_b11_p, frac_c12_b11_n_piplus, frac_c12_b11_n_piminus_2piplus, frac_c12_b11_p_piminus_piplus, \
           frac_c12_b11_p_2piminus_2piplus, frac_c12_b11_piplus, frac_c12_b11_other, \
           frac_c12_c11_n, frac_c12_c11_p_piminus, frac_c12_c11_n_piminus_piplus, frac_c12_c11_p_2piminus_piplus, \
           frac_c12_c11_p_3piminus_2piplus, frac_c12_c11_n_2piminus_2piplus, frac_c12_c11_other, \
           frac_c12_b10_p_n, frac_c12_b10_2p_piminus, frac_c12_b10_p_n_piminus_piplus, frac_c12_b10_2n_piplus, \
           frac_c12_b10_2n_piminus_2piplus, frac_c12_b10_2p_2piminus_piplus, frac_c12_b10_2p_3piminus_2piplus, \
           frac_c12_b10_p_n_2piminus_2piplus, frac_c12_b10_other, \
           frac_c12_c10_2n, frac_c12_c10_p_n_piminus, frac_c12_c10_p_n_2piminus_piplus, frac_c12_c10_2n_piminus_piplus,\
           frac_c12_c10_2p_2piminus, frac_c12_c10_other, \
           frac_c12_be10_2p, frac_c12_be10_p_n_piplus, frac_c12_be10_p_n_piminus_2piplus, \
           frac_c12_be10_2p_piminus_piplus, frac_c12_be10_2n_2piplus, frac_c12_be10_p_n_2piminus_3piplus, \
           frac_c12_be10_2p_2piminus_2piplus, frac_c12_be10_2p_3piminus_3piplus, frac_c12_be10_other, \
           frac_c12_b9_p_2n, frac_c12_b9_p_2n_piminus_piplus, frac_c12_b9_2p_n_3piminus_2piplus, \
           frac_c12_b9_2p_n_piminus, frac_c12_b9_3n_piplus, frac_c12_b9_p_2n_2piminus_2piplus, \
           frac_c12_b9_2p_n_2piminus_piplus, frac_c12_b9_other, \
           frac_c12_be9_2p_n, frac_c12_be9_p_2n_piplus, frac_c12_be9_3p_piminus, frac_c12_be9_p_2n_piminus_2piplus, \
           frac_c12_be9_2p_n_piminus_piplus, frac_c12_be9_2p_n_3piminus_3piplus, frac_c12_be9_2p_n_2piminus_2piplus, \
           frac_c12_be9_3n_2piplus, frac_c12_be9_3p_2piminus_piplus, frac_c12_be9_other, \
           frac_c12_be8_2p_2n, frac_c12_be8_3p_n_piminus, frac_c12_be8_p_3n_piplus, \
           frac_c12_be8_2p_2n_2piminus_2piplus, frac_c12_be8_4n_2piplus, frac_c12_be8_2p_2n_piminus_piplus, \
           frac_c12_be8_3p_n_2piminus_piplus, frac_c12_be8_4p_2piminus, frac_c12_be8_other, \
           frac_c12_c9_p_2n_piminus, frac_c12_c9_3n, frac_c12_c9_2p_n_2piminus, frac_c12_c9_3n_2piminus_2piplus, \
           frac_c12_c9_other, \
           frac_c12_be7_2p_3n, frac_c12_be7_p_4n_piplus, frac_c12_be7_2p_3n_2piminus_2piplus, \
           frac_c12_be7_3p_2n_piminus, frac_c12_be7_4p_n_2piminus, frac_c12_be7_3p_2n_2piminus_piplus, \
           frac_c12_be7_other, \
           frac_c12_li6_3p_3n, frac_c12_li6_2p_4n_piplus, frac_c12_li6_5p_n_2piminus, \
           frac_c12_li6_2p_4n_piminus_2piplus, frac_c12_li6_4p_2n_piminus, frac_c12_li6_3p_3n_piminus_piplus, \
           frac_c12_li6_other, \
           frac_c12_li8_3p_n, frac_c12_li8_4p_piminus, frac_c12_li8_4p_2piminus_piplus, frac_c12_li8_2p_2n_piplus, \
           frac_c12_li8_3p_n_piminus_piplus, frac_c12_li8_other, \
           frac_c12_li7_2p_3n_piplus, frac_c12_li7_4p_n_piminus, frac_c12_li7_3p_2n, frac_c12_li7_3p_2n_piminus_piplus,\
           frac_c12_li7_4p_n_2piminus_piplus, frac_c12_li7_2p_3n_piminus_2piplus, frac_c12_li7_other, \
           frac_c12_b8_p_3n, frac_c12_b8_p_3n_piminus_piplus, frac_c12_b8_2p_2n_2piminus_piplus, \
           frac_c12_b8_2p_2n_piminus, frac_c12_b8_4n_piplus, frac_c12_b8_other, \
           frac_c12_li9_2p_n_piplus, frac_c12_li9_3p, frac_c12_li9_3p_piminus_piplus, \
           frac_c12_li9_2p_n_piminus_2piplus, frac_c12_li9_p_2n_piminus_3piplus, frac_c12_li9_other, \
           frac_c12_c8_4n, frac_c12_c8_4n_other, frac_c12_he8_4p, frac_c12_he8_4p_other, \
           frac_c12_b7_p_4n, frac_c12_b7_p_4n_other, frac_c12_he7_4p_n, frac_c12_he7_4p_n_other, \
           frac_c12_h7_5p, frac_c12_h7_5p_other, frac_c12_be6_2p_4n, frac_c12_be6_2p_4n_other, \
           frac_c12_he6_4p_2n, frac_c12_he6_4p_2n_other, frac_c12_h6_5p_n, frac_c12_h6_5p_n_other, \
           frac_c12_mass11u, frac_c12_mass10u, frac_c12_mass9u, frac_c12_mass8u, frac_c12_mass7u, frac_c12_mass6u, \
           frac_c12_mass5orless, \
           frac_c12_c12, frac_c12_noiso, frac_c12_noiso_5p_6n, \
           frac_no_c12, frac_es_p, frac_es_e, frac_es_o16, frac_es_n14, frac_es_s32, frac_c12_faulty


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

    return number_entries, number_target_c12, number_no_c12, number_light_iso, \
           number_c11_notex, number_c11_deex, number_c11_li6_p_alpha, number_c11_be7_alpha, number_c11_b10_p, \
           number_c11_b9_n_p, number_c11_be8_p_d, number_c11_be9_2p, number_c11_b9_d, number_c11_be8_he3, \
           number_c11_c10_n, number_c11_li5_d_alpha, number_c11_li5_n_p_alpha, number_c11_missing, \
           number_b11_notex, number_b11_deex, number_b11_li6_n_alpha, number_b11_b9_2n, number_b11_be8_n_d, \
           number_b11_be9_d, number_b11_be10_p, number_b11_b10_n, number_b11_be9_n_p, number_b11_li7_alpha, \
           number_b11_be8_t, number_b11_he5_d_alpha, number_b11_he6_p_alpha, number_b11_be8_2n_p, number_b11_missing, \
           number_c10_notex, number_c10_deex, number_c10_b9_p, number_c10_be7_p_d, number_c10_li6_p_he3, \
           number_c10_he4_p_d_he3, number_c10_li6_2p_d, number_c10_be8_2p, number_c10_be7_n_2p, number_c10_li6_n_3p, \
           number_c10_be6_n_p_d, number_c10_li5_n_2p_d, number_c10_he3_p_d_alpha, number_c10_li5_d_he3, \
           number_c10_li5_p_2d, number_c10_he3_n_2p_alpha, number_c10_li4_n_p_alpha, number_c10_b8_n_p, \
           number_c10_b8_d, number_c10_be6_p_t, number_c10_he4_n_2p_he3, number_c10_li5_n_p_he3, number_c10_missing, \
           number_b10_notex, number_b10_deex, number_b10_be9_p, number_b10_be8_d, number_b10_b9_n, number_b10_be7_t, \
           number_b10_be8_n_p, number_b10_li7_he3, number_b10_he5_p_alpha, number_b10_li6_alpha, \
           number_b10_li5_n_alpha, number_b10_li7_p_d, number_b10_missing, \
           number_be10_notex, number_be10_deex, number_be10_li6_2n_d, number_be10_li7_2n_p, number_be10_he7_n_2p, \
           number_be10_li6_n_t, number_be10_li6_3n_p, number_be10_he5_n_2d, number_be10_be8_2n, \
           number_be10_he5_n_alpha, number_be10_t_n_d_alpha, number_be10_li8_n_p, number_be10_he4_n_d_t, \
           number_be10_he5_n_p_t, number_be10_li7_n_d, number_be10_h4_n_p_alpha, number_be10_he5_2n_p_d, \
           number_be10_be9_n, number_be10_he6_n_p_d, number_be10_he5_d_t, number_be10_h4_d_alpha, \
           number_be10_he4_2n_p_t, number_be10_missing, \
           number_c9_notex, number_c9_deex, number_c9_be7_2p, number_c9_missing, \
           number_b9_notex, number_b9_deex, number_b9_be8_p, number_b9_missing, \
           number_be9_notex, number_be9_deex, number_be9_li8_p, number_be9_missing, \
           number_li9_notex, number_li9_deex, number_li9_h4_n_alpha, number_li9_he7_d, number_li9_li8_n, \
           number_li9_missing, \
           number_b8_notex, number_b8_deex, number_b8_li6_2p, number_b8_missing, \
           number_li8_notex, number_li8_deex, number_li8_li7_n, number_li8_li6_2n, number_li8_missing, \
           number_be7_notex, number_be7_deex, number_be7_li5_d, number_be7_li6_p, number_be7_missing, \
           number_li7_notex, number_li7_deex, number_li7_li6_n, number_li7_missing


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
