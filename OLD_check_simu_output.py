""" script to check the different JUNO simulation parts (detsim, det2elec, elec2calib, calib2rec):

"""
import ROOT
import numpy as np
from matplotlib import pyplot as plt


""" Set file number: """
file_number = 0
""" Set event number: """
event_number = 22

""" get TFile and TTree information """
# input root file from detsim.py: user_atmoNC_0.root:
rootfile_detsim = "/local/scratch1/pipc51/astro/blum/detsim_output_data/user_atmoNC_{0:d}.root".format(file_number)
# input root file from det2elec.py: user_elecsim_atmoNC_0.root:
rootfile_det2elec = \
    "/local/scratch1/pipc51/astro/blum/det2elec_output_data/elecsim_output_TAG/" \
    "user_elecsim_atmoNC_{0:d}.root".format(file_number)
# input root file from draw_waveform.py: elecAnalysis_atmoNC_0.root
rootfile_elecAnalysis = "/local/scratch1/pipc51/astro/blum/det2elec_output_data/elecsim_output_TAG/waveforms/" \
                        "elecAnalysis_atmoNC_{0:d}.root".format(file_number)
# input root file from elec2calib.py: user_calib_atmoNC_0.root
rootfile_elec2calib = "/local/scratch1/pipc51/astro/blum/elec2calib_output_data/user_calib_atmoNC_{0:d}.root"\
    .format(file_number)
# input root file from calib2rec.py: user_rec_atmoNC_0.root
rootfile_calib2rec = "/local/scratch1/pipc51/astro/blum/calib2rec_output_data/user_rec_atmoNC_{0:d}.root"\
    .format(file_number)

# load the ROOT files:
rfile_detsim = ROOT.TFile(rootfile_detsim)
rfile_det2elec = ROOT.TFile(rootfile_det2elec)
rfile_elecAnalysis = ROOT.TFile(rootfile_elecAnalysis)
rfile_elec2calib = ROOT.TFile(rootfile_elec2calib)
rfile_calib2rec = ROOT.TFile(rootfile_calib2rec)

# get TTrees from TFiles:
rtree_detsim_geninfo = rfile_detsim.Get("geninfo")
rtree_detsim_prmtrkdep = rfile_detsim.Get("prmtrkdep")
rtree_det2elec = rfile_det2elec.Get("SIMEVT")
rtree_elecAnalysis = rfile_elecAnalysis.Get("SIMEVT")
rtree_elec2calib = rfile_elec2calib.Get("CALIBEVT")
rtree_calib2rec = rfile_calib2rec.Get("USER_OUTPUT")

# Preallocate interesting values:
# detsim:
nInitParticles_geninfo = 0
InitPDGID_geninfo = np.array([])
InitX_geninfo = np.array([])
InitY_geninfo = np.array([])
InitZ_geninfo = np.array([])
ExitX_geninfo = np.array([])
ExitY_geninfo = np.array([])
ExitZ_geninfo = np.array([])
ExitT_geninfo = np.array([])
edep_prmtrkdep = np.array([])
Qedep_prmtrkdep = np.array([])

# elecAnalysis:
fired_PMT_Num = 0
tdc = []
adc = []

# elec2calib:
Charge = np.array([])
Time = np.array([])

# calib2rec:
xrec = 0
yrec = 0
zrec = 0
Erec = 0


""" first read the "geninfo" Tree"""
# get the current event in geninfo tree:
rtree_detsim_geninfo.GetEntry(event_number)
# get the value of the number of initial particles:
nInitParticles_geninfo = int(rtree_detsim_geninfo.GetBranch('nInitParticles').GetLeaf('nInitParticles').GetValue())

# loop over the number of particles to get information about every particle in the event:
for index in range(nInitParticles_geninfo):

    initpdgid_geninfo = int(rtree_detsim_geninfo.GetBranch('InitPDGID').GetLeaf('InitPDGID').GetValue(index))
    InitPDGID_geninfo = np.append(InitPDGID_geninfo, initpdgid_geninfo)

    initx_geninfo = float(rtree_detsim_geninfo.GetBranch('InitX').GetLeaf('InitX').GetValue(index))
    InitX_geninfo = np.append(InitX_geninfo, initx_geninfo)

    inity_geninfo = float(rtree_detsim_geninfo.GetBranch('InitY').GetLeaf('InitY').GetValue(index))
    InitY_geninfo = np.append(InitY_geninfo, inity_geninfo)

    initz_geninfo = float(rtree_detsim_geninfo.GetBranch('InitZ').GetLeaf('InitZ').GetValue(index))
    InitZ_geninfo = np.append(InitZ_geninfo, initz_geninfo)

    exitx_geninfo = float(rtree_detsim_geninfo.GetBranch('ExitX').GetLeaf('ExitX').GetValue(index))
    ExitX_geninfo = np.append(ExitX_geninfo, exitx_geninfo)

    exity_geninfo = float(rtree_detsim_geninfo.GetBranch('ExitY').GetLeaf('ExitY').GetValue(index))
    ExitY_geninfo = np.append(ExitY_geninfo, exity_geninfo)

    exitz_geninfo = float(rtree_detsim_geninfo.GetBranch('ExitZ').GetLeaf('ExitZ').GetValue(index))
    ExitZ_geninfo = np.append(ExitZ_geninfo, exitz_geninfo)

    exitt_geninfo = float(rtree_detsim_geninfo.GetBranch('ExitT').GetLeaf('ExitT').GetValue(index))
    ExitT_geninfo = np.append(ExitT_geninfo, exitt_geninfo)


""" read prmtrkdep tree: """
# get current event in prmtrkdep tree:
rtree_detsim_prmtrkdep.GetEntry(event_number)

# loop over number of particles:
for index in range(nInitParticles_geninfo):
    edep = float(rtree_detsim_prmtrkdep.GetBranch('edep').GetLeaf('edep').GetValue(index))
    edep_prmtrkdep = np.append(edep_prmtrkdep, edep)

    qedep_prmtrkdep = float(rtree_detsim_prmtrkdep.GetBranch('Qedep').GetLeaf('Qedep').GetValue(index))
    Qedep_prmtrkdep = np.append(Qedep_prmtrkdep, qedep_prmtrkdep)


""" read elecAnalysis file: """
# get current event in elecAnalysis file
rtree_elecAnalysis.GetEntry(event_number)
# get number of fired PMTs:
fired_PMT_Num = int(rtree_elecAnalysis.GetBranch('fired_PMT_Num').GetLeaf('fired_PMT_Num').GetValue())

# loop over the number of fired PMTs:

# for index in range(fired_PMT_Num):
for index in range(5):
    tdc_array = np.array([])
    adc_array = np.array([])

    for index2 in range(index*1250, index*1250+1249):

        tdc_value = float(rtree_elecAnalysis.GetBranch('tdc').GetLeaf('tdc').GetValue(index2))
        tdc_array = np.append(tdc_array, tdc_value)

        adc_value = float(rtree_elecAnalysis.GetBranch('adc').GetLeaf('adc').GetValue(index2))
        adc_array = np.append(adc_array, adc_value)

    tdc = np.append(tdc, tdc_array)
    adc = np.append(adc, adc_array)


""" read elec2calib file: """



""" read calib2rec file: """
# get current event in calib2rec file:
rtree_calib2rec.GetEntry(event_number)
xrec = float(rtree_calib2rec.GetBranch('xrec').GetLeaf('xrec').GetValue())
yrec = float(rtree_calib2rec.GetBranch('yrec').GetLeaf('yrec').GetValue())
zrec = float(rtree_calib2rec.GetBranch('zrec').GetLeaf('zrec').GetValue())
Erec = float(rtree_calib2rec.GetBranch('Erec').GetLeaf('Erec').GetValue())


""" plot waveform: """
h1 = plt.figure(1)
plt.plot(tdc[3750:4999], adc[3750:4999])
plt.xlabel("tdc")
plt.ylabel("adc")
h2 = plt.figure(2)
plt.plot(tdc[2500:3749], adc[2500:3749])
plt.xlabel("tdc")
plt.ylabel("adc")


""" compare detsim with reconstruction: """
print("Event {0:d} in file {1:d}".format(event_number, file_number))
print("-------------------------------------------")
print("PDG ID: {0}".format(InitPDGID_geninfo))
print("-------------------------------------------")
print("Position:")
print("detsim:\tinit:\tx={0},\ny={1},\nz={2}".format(InitX_geninfo, InitY_geninfo, InitZ_geninfo))
print("detsim:\texit:\tx={0},\ny={1},\nz={2}".format(ExitX_geninfo, ExitY_geninfo, ExitZ_geninfo))
print("calib2rec:\tx={0},\ny={1},\nz={2}".format(xrec, yrec, zrec))
print("-------------------------------------------")
print("Energy:")
print("detsim: edep={0},\nQedep={1}".format(edep_prmtrkdep, Qedep_prmtrkdep))
print("calib2rec: Erec={0}".format(Erec))

plt.show()














