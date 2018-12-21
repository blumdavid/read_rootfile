""" function to check the variables 'nc' and 'qel' of the 'original' GENIE root-file from Julia.

"""
import ROOT

# set the path of the inputs:
input_path = "/home/astro/blum/juno/atmoNC/data_Julia/"

# file name of the input file:
input_name = input_path + "gntp.101.gst.root"

# load the ROOT file:
rfile_input = ROOT.TFile(input_name)
# get the TTree from the TFile:
rtree_input = rfile_input.Get("gst")

# Info-me: "gst;13" is a copy of meta data of "gst;14", "gst;14" contains correct data and is read

# get the number of entries in the ROOT-file:
# number_entries = rtree_input.GetEntries()
number_entries = 200

""" Read the data from the TTree: """
# loop over every entry, i.e. every event, in the TTree:
for event in range(0, number_entries):

    # get the current event in the TTree:
    rtree_input.GetEntry(event)

    # energy of incoming neutrino in GeV:
    ev = rtree_input.GetBranch('Ev').GetLeaf('Ev').GetValue()

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
    nfpi0 = rtree_input.GetBranch('nfpi0').GetLeaf('nfpi0').GetValue()
    nfpi0 = int(nfpi0)

    # get the value of number of final Kaon_minus:
    nfkm = rtree_input.GetBranch('nfkm').GetLeaf('nfkm').GetValue()
    nfkm = int(nfkm)

    # get the value of number of final Kaon_plus:
    nfkp = rtree_input.GetBranch('nfkp').GetLeaf('nfkp').GetValue()
    nfkp = int(nfkp)

    # get value of number of final K_zero:
    nfk0 = rtree_input.GetBranch('nfk0').GetLeaf('nfk0').GetValue()
    nfk0 = int(nfk0)

    # get value of number of final gamma, electron, positron:
    nfem = rtree_input.GetBranch('nfem').GetLeaf('nfem').GetValue()
    nfem = int(nfem)

    # get value of number of final other hadrons:
    nfother = rtree_input.GetBranch('nfother').GetLeaf('nfother').GetValue()
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
    nipi0 = rtree_input.GetBranch('nipi0').GetLeaf('nipi0').GetValue()
    nipi0 = int(nipi0)

    # get the value of number of primary Kaon_minus:
    nikm = rtree_input.GetBranch('nikm').GetLeaf('nikm').GetValue()
    nikm = int(nikm)

    # get the value of number of primary Kaon_plus:
    nikp = rtree_input.GetBranch('nikp').GetLeaf('nikp').GetValue()
    nikp = int(nikp)

    # get value of number of primary K_zero:
    nik0 = rtree_input.GetBranch('nik0').GetLeaf('nik0').GetValue()
    nik0 = int(nik0)

    # get value of number of primary gamma, electron, positron:
    niem = rtree_input.GetBranch('niem').GetLeaf('niem').GetValue()
    niem = int(niem)

    # get value of number of primary other hadrons:
    niother = rtree_input.GetBranch('niother').GetLeaf('niother').GetValue()
    niother = int(niother)

    # check primary and final particles for target C12:
    if tgt == 1000060120 and nc == 1:

        # B11:
        if nfp == 1 and nfn == 0 and nfpim == 0 and nfpip == 0 and nfkm == 0 and nfkp == 0:
            # interaction channel: nu + C12 -> B11 + proton:
            print("nu + C12 -> nu + B11 + p: nc={0}, qel={1},\n"
                  "nfpi0={2:d}, nfk0={3:d}, nfem={4:d}, nfother={5:d}, "
                  "nip={6:d}, nin={7:d}, nipip={8:d}, nipim={9:d}, nipi0={10:d}, nikp={11:d}, nikm={12:d}, "
                  "nik0={13:d}, niem={14:d}, niother={15:d}, ev={16:.2f}\n"
                  .format(nc, qel, nfpi0, nfk0, nfem, nfother, nip, nin, nipip, nipim, nipi0, nikp, nikm, nik0, niem,
                          niother, ev))

        elif nfp == 0 and nfn == 1 and nfpim == 0 and nfpip == 1 and nfkm == 0 and nfkp == 0:
            # interaction channel: nu + C12 -> B11 + n + pi_plus:
            print("nu + C12 -> nu + B11 + n + pi_plus: nc={0}, qel={1},\n "
                  "nfpi0={2:d}, nfk0={3:d}, nfem={4:d}, nfother={5:d}, "
                  "nip={6:d}, nin={7:d}, nipip={8:d}, nipim={9:d}, nipi0={10:d}, nikp={11:d}, nikm={12:d}, "
                  "nik0={13:d}, niem={14:d}, niother={15:d}, ev={16:.2f}\n"
                  .format(nc, qel, nfpi0, nfk0, nfem, nfother, nip, nin, nipip, nipim, nipi0, nikp, nikm, nik0, niem,
                          niother, ev))
        # C11:
        elif nfp == 0 and nfn == 1 and nfpim == 0 and nfpip == 0 and nfkm == 0 and nfkp == 0:
            # interaction channel: nu + C12 -> C11 + n:
            print("nu + C12 -> nu + C11 + n: nc={0}, qel={1},\n "
                  "nfpi0={2:d}, nfk0={3:d}, nfem={4:d}, nfother={5:d}, "
                  "nip={6:d}, nin={7:d}, nipip={8:d}, nipim={9:d}, nipi0={10:d}, nikp={11:d}, nikm={12:d}, "
                  "nik0={13:d}, niem={14:d}, niother={15:d}, ev={16:.2f}\n"
                  .format(nc, qel, nfpi0, nfk0, nfem, nfother, nip, nin, nipip, nipim, nipi0, nikp, nikm, nik0, niem,
                          niother, ev))

        elif nfp == 1 and nfn == 0 and nfpim == 1 and nfpip == 0 and nfkm == 0 and nfkp == 0:
            # interaction channel: nu + C12 -> C11 + p + pi_minus:
            print("nu + C12 -> nu + C11 + p + pi_minus: nc={0}, qel={1},\n "
                  "nfpi0={2:d}, nfk0={3:d}, nfem={4:d}, nfother={5:d}, "
                  "nip={6:d}, nin={7:d}, nipip={8:d}, nipim={9:d}, nipi0={10:d}, nikp={11:d}, nikm={12:d}, "
                  "nik0={13:d}, niem={14:d}, niother={15:d}, ev={16:.2f}\n"
                  .format(nc, qel, nfpi0, nfk0, nfem, nfother, nip, nin, nipip, nipim, nipi0, nikp, nikm, nik0, niem,
                          niother, ev))

        # B10:
        elif nfp == 1 and nfn == 1 and nfpim == 0 and nfpip == 0 and nfkm == 0 and nfkp == 0:
            # interaction channel: nu + C12 -> B10 + p + n:
            print("nu + C12 -> nu + B10 + p + n: nc={0}, qel={1},\n "
                  "nfpi0={2:d}, nfk0={3:d}, nfem={4:d}, nfother={5:d}, "
                  "nip={6:d}, nin={7:d}, nipip={8:d}, nipim={9:d}, nipi0={10:d}, nikp={11:d}, nikm={12:d}, "
                  "nik0={13:d}, niem={14:d}, niother={15:d}\n".format(nc, qel, nfpi0, nfk0, nfem, nfother, nip, nin,
                                                                      nipip, nipim, nipi0, nikp, nikm, nik0, niem,
                                                                      niother))

        elif nfp == 2 and nfn == 0 and nfpim == 1 and nfpip == 0 and nfkm == 0 and nfkp == 0:
            # interaction channel: nu + C12 -> B10 + 2*p + pi_minus:
            print("nu + C12 -> nu + B10 + 2p + pi_minus: nc={0}, qel={1},\n "
                  "nfpi0={2:d}, nfk0={3:d}, nfem={4:d}, nfother={5:d}, "
                  "nip={6:d}, nin={7:d}, nipip={8:d}, nipim={9:d}, nipi0={10:d}, nikp={11:d}, nikm={12:d}, "
                  "nik0={13:d}, niem={14:d}, niother={15:d}\n".format(nc, qel, nfpi0, nfk0, nfem, nfother, nip, nin,
                                                                      nipip, nipim, nipi0, nikp, nikm, nik0, niem,
                                                                      niother))

        else:
            continue
