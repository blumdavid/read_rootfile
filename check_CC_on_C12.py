""" script to check the GENIE file gntp.101.gst.root of Julia for possible atmospheric CC background on C12.

    Info about gntp.101.gst.root file from Email of Julia:
    Hierbei ist gntp.101.gst.root das "Original" File welches aus der Simulation kommt: Hier ist folgendes simuliert:
    *   solmax Honda Fluesse fuer nu_e, nu_e_bar, nu_mu, nu_mu_bar
    *   scintillator zusammensetzung: 1000060120[0.87924],1000010010[0.1201],1000080160[0.00034],1000070140[0.00027],
        1000160320[5e-5]
    *   von 0-10 GeV
    *   insgesamt 1 Mio reaktionen

"""
import ROOT
import numpy as np

# path, where GENIE file is saved:
input_path = "/home/astro/blum/juno/atmoNC/data_Julia/"
# file name:
input_file = input_path + "gntp.101.gst.root"

# load gntp.101.gst.root file:
rfile_input = ROOT.TFile(input_file)

# get the "gst"-TTree from the TFile:
rtree_input = rfile_input.Get("gst")

# get the number of events in the 'gst' Tree:
num_events = rtree_input.GetEntries()

# define max. energy, that one event can have (normally 100 MeV is the edge for indirect DM search) in GeV:
energy_max = 0.1

# define number of events that pass the cuts for each neutrino type:
number = 0

# loop over entries of gst tree:
for event in range(num_events):

    # get current event in tree:
    rtree_input.GetEntry(event)

    # total kinetic energy of final state particles in GeV:
    energy = 0.0

    # array, where pdgf is stored:
    pdgf_array = np.array([])
    # array, where pf is stored:
    pf_array = np.array([])

    # set flag for final state particles:
    # no muon in final state:
    flag_no_muon = True
    # no pion in final state
    flag_no_pion = True
    # no proton in final state:
    flag_no_proton = True
    # no lambda, kaon or other exotic particles in final state:
    flag_no_exotic = True
    # neutron in final state:
    flag_neutron = False
    # flag_neutron = True

    # get the target PDG and check if the target is C12:
    tgt = rtree_input.GetBranch('tgt').GetLeaf('tgt').GetValue()
    tgt = int(tgt)

    if tgt == 1000060120:
    # if tgt == 2112:

        # get interaction type and check, if interaction type is charged current (cc==1):
        cc = int(rtree_input.GetBranch('cc').GetLeaf('cc').GetValue())
        nc = int(rtree_input.GetBranch('nc').GetLeaf('nc').GetValue())

        if cc == 1:
        # if nc == 1:
            # event is CC on C12:

            # check the type of the incoming neutrino:
            neu = int(rtree_input.GetBranch('neu').GetLeaf('neu').GetValue())

            if neu < 20:

                # neutrino + C12 -> ...
                # get number of final state particles:
                nf = int(rtree_input.GetBranch('nf').GetLeaf('nf').GetValue())
                # loop over nf_nu_e and add Ef[index] to energy_nu_e:
                for index in range(nf):
                    # momentum of final state particle:
                    pxf = float(rtree_input.GetBranch('pxf').GetLeaf('pxf').GetValue(index))
                    pyf = float(rtree_input.GetBranch('pyf').GetLeaf('pyf').GetValue(index))
                    pzf = float(rtree_input.GetBranch('pzf').GetLeaf('pzf').GetValue(index))

                    # total momentum / kinetic energy of final state particle:
                    pf = np.sqrt(pxf**2 + pyf**2 + pzf**2)

                    energy += pf

                    pf_array = np.append(pf_array, pf)

                    pdgf = int(rtree_input.GetBranch('pdgf').GetLeaf('pdgf').GetValue(index))

                    pdgf_array = np.append(pdgf_array, pdgf)

                    if pdgf == 13 or pdgf == -13:
                        # muon in final state:
                        print("muon")
                        flag_no_muon = False
                    elif pdgf == 11 or pdgf == -11:
                        # electron in final state:
                        print("electron")
                    elif pdgf == 211 or pdgf == -211 or pdgf == 111:
                        # print("pion")
                        # pion in final state:
                        flag_no_pion = False
                    elif (pdgf == 311 or pdgf == -311 or pdgf == 321 or pdgf == -321 or pdgf == 3222 or pdgf == 3122 or
                          pdgf == 3112 or pdgf == 130 or pdgf == 3212 or pdgf == 8285 or pdgf == -2112 or
                          pdgf == -2212 or pdgf == 1000060120 or pdgf == 4122 or pdgf == 4222 or pdgf == 421 or
                          pdgf == 4212 or pdgf == 411 or pdgf == -421 or pdgf == -411 or pdgf == 431):
                        # kaon, lambda or other exotic particles in final state:
                        flag_no_exotic = False
                    # elif pdgf == 2212 and pf > 1.0:
                    #     # proton with Ef > 1 GeV in final state:
                    #     flag_no_proton = False
                    elif pdgf == 2112:
                        # neutron in final state:
                        flag_neutron = True

                    elif pdgf != 2112 and pdgf != 2212 and pdgf != 22 and pdgf != 1000060120:
                        print(pdgf)

        else:
            continue

    else:
        continue

    # check the final state energies:
    if 0 < energy < energy_max and flag_no_muon and flag_no_pion and flag_no_exotic and flag_no_proton and flag_neutron:
        # print(pdgf_array)
        # print(Ef_array)
        number += 1

print(number)
