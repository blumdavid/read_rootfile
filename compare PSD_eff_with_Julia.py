""" script to check the NC, IBD and fast neutron PSD efficiencies for the energy range between 10 MeV and 30 MeV to
    compare it with the DSNB study results of Julia.
"""
import ROOT
import numpy as np
from matplotlib import pyplot as plt
import NC_background_functions

path_PSD_info = "/home/astro/blum/juno/atmoNC/data_NC/output_detsim_v2/" \
                "DCR_results_16000mm_10MeVto100MeV_1000nsto1ms_mult1_1800keVto2550keV_dist500mm_R17700mm_PSD99/" \
                "test_10to20_20to30_30to40_40to100_final/"

""" IBD efficiency: """
# load simulated IBD spectrum (from 10 to 100 MeV with bin-width 1 MeV):
array_IBD_spectrum = np.loadtxt(path_PSD_info + "IBDspectrum_woPSD_bin1000keV.txt")
# load simulated IBD spectrum after PSD (from 10 to 100 MeV with bin-width 1 MeV):
array_IBD_spectrum_PSD = np.loadtxt(path_PSD_info + "IBDspectrum_wPSD_bin1000keV.txt")

# number of IBD events between 10 and 30 MeV before PSD:
N_IBD = np.sum(array_IBD_spectrum[0:20])
# number of IBD events between 10 and 30 MeV after PSD:
N_IBD_PSD = np.sum(array_IBD_spectrum_PSD[0:20])
# calculate IBD efficiency in percent:
IBD_eff = N_IBD_PSD / N_IBD * 100.0

print("IBD efficiency in precent:")
print(IBD_eff)
print("")

""" NC efficiency of IBD-like events: """
# load simulated NC spectrum (from 10 to 100 MeV with bin-width 0.5 MeV):
array_NC_spectrum = np.loadtxt(path_PSD_info + "NCatmo_onlyC12_woPSD_bin500keV.txt")
# load simulated NC spectrum after PSD (from 10 to 100 MeV with bin-width 0.5 MeV):
array_NC_spectrum_PSD = np.loadtxt(path_PSD_info + "NCatmo_onlyC12_wPSD99_bin500keV.txt")

# number of NC events between 10 and 30 MeV before PSD:
N_NC = np.sum(array_NC_spectrum[0:40])
# number of NC events between 10 and 30 MeV after PSD:
N_NC_PSD = np.sum(array_NC_spectrum_PSD[0:40])
# calculate NC efficiency in percent:
NC_eff = N_NC_PSD / N_NC * 100.0

print("NC efficiency in precent:")
print(NC_eff)
print("")

""" fast neutron efficiency from /home/astro/blum/PhD/work/MeVDM_JUNO/fast_neutrons/DCR/results.txt: """
# number of events between 10 and 20 MeV before PSD:
N_FN_10_20 = 253.0
# number of events between 10 and 20 MeV after PSD:
N_FN_10_20_PSD = 0.0
# number of events between 20 and 30 MeV before PSD:
N_FN_20_30 = 536.0
# number of events between 20 and 30 MeV after PSD:
N_FN_20_30_PSD = 1.0
# mean number of events per 10 MeV, if FN spectrum is flat distribution:
N_FN_mean = 500.5
# number of events between 10 and 20 MeV before PSD after weighting:
N_FN_10_20_w = N_FN_10_20 * N_FN_mean / N_FN_10_20
# number of events between 10 and 20 MeV after PSD after weighting:
N_FN_10_20_PSD_w = N_FN_10_20_PSD * N_FN_mean / N_FN_10_20
# number of events between 20 and 30 MeV before PSD after weighting:
N_FN_20_30_w = N_FN_20_30 * N_FN_mean / N_FN_20_30
# number of events between 20 and 30 MeV after PSD after weighting:
N_FN_20_30_PSD_w = N_FN_20_30_PSD * N_FN_mean / N_FN_20_30
# calculate FN efficiency between 10 and 30 MeV in percent:
FN_eff = (N_FN_10_20_PSD_w + N_FN_20_30_PSD_w) / (N_FN_10_20_w + N_FN_20_30_w) * 100.0

print("FN efficiency in precent:")
print(FN_eff)
print("")

""" Display efficiencies between 10 and 40 MeV for 5 MeV bins: """
# energy range:
energy = [10, 15, 20, 25, 30, 35, 40]

# IBD:
N_IBD_10_15 = np.sum(array_IBD_spectrum[0:5])
N_IBD_15_20 = np.sum(array_IBD_spectrum[5:10])
N_IBD_20_25 = np.sum(array_IBD_spectrum[10:15])
N_IBD_25_30 = np.sum(array_IBD_spectrum[15:20])
N_IBD_30_35 = np.sum(array_IBD_spectrum[20:25])
N_IBD_35_40 = np.sum(array_IBD_spectrum[25:30])
N_IBD_10_15_PSD = np.sum(array_IBD_spectrum_PSD[0:5])
N_IBD_15_20_PSD = np.sum(array_IBD_spectrum_PSD[5:10])
N_IBD_20_25_PSD = np.sum(array_IBD_spectrum_PSD[10:15])
N_IBD_25_30_PSD = np.sum(array_IBD_spectrum_PSD[15:20])
N_IBD_30_35_PSD = np.sum(array_IBD_spectrum_PSD[20:25])
N_IBD_35_40_PSD = np.sum(array_IBD_spectrum_PSD[25:30])

# IBD efficiencies in precent:
array_IBD_eff = [N_IBD_10_15_PSD/N_IBD_10_15*100.0, N_IBD_15_20_PSD/N_IBD_15_20*100.0,
                 N_IBD_20_25_PSD/N_IBD_20_25*100.0, N_IBD_25_30_PSD/N_IBD_25_30*100.0,
                 N_IBD_30_35_PSD/N_IBD_30_35*100.0, N_IBD_35_40_PSD/N_IBD_35_40*100.0,
                 N_IBD_35_40_PSD/N_IBD_35_40*100.0]

# NC:
N_NC_10_15 = np.sum(array_NC_spectrum[0:10])
N_NC_15_20 = np.sum(array_NC_spectrum[10:20])
N_NC_20_25 = np.sum(array_NC_spectrum[20:30])
N_NC_25_30 = np.sum(array_NC_spectrum[30:40])
N_NC_30_35 = np.sum(array_NC_spectrum[40:50])
N_NC_35_40 = np.sum(array_NC_spectrum[50:60])
N_NC_10_15_PSD = np.sum(array_NC_spectrum_PSD[0:10])
N_NC_15_20_PSD = np.sum(array_NC_spectrum_PSD[10:20])
N_NC_20_25_PSD = np.sum(array_NC_spectrum_PSD[20:30])
N_NC_25_30_PSD = np.sum(array_NC_spectrum_PSD[30:40])
N_NC_30_35_PSD = np.sum(array_NC_spectrum_PSD[40:50])
N_NC_35_40_PSD = np.sum(array_NC_spectrum_PSD[50:60])

# NC efficiencies in precent:
array_NC_eff = [N_NC_10_15_PSD/N_NC_10_15*100.0, N_NC_15_20_PSD/N_NC_15_20*100.0,
                N_NC_20_25_PSD/N_NC_20_25*100.0, N_NC_25_30_PSD/N_NC_25_30*100.0,
                N_NC_30_35_PSD/N_NC_30_35*100.0, N_NC_35_40_PSD/N_NC_35_40*100.0,
                N_NC_35_40_PSD/N_NC_35_40*100.0]

# FN:
energy_FN = [10, 20, 30, 40]

# FN efficiencies in precent from /home/astro/blum/PhD/work/MeVDM_JUNO/fast_neutrons/DCR/results.txt:
array_FN_eff = [0.0, 0.16, 0.24, 0.24]

plt.figure(1, figsize=(11, 6))
plt.step(energy, array_IBD_eff, color="r", label="IBD", where='post')
plt.step(energy, array_NC_eff, color="b", label="atmo. NC", where='post')
plt.step(energy_FN, array_FN_eff, color="g", label="fast neutron", where='post')
plt.xlim(xmin=10, xmax=40)
plt.xlabel("Visible energy of prompt signal in MeV")
plt.ylabel("PSD efficiency in % (for bins of 5 MeV)")
plt.title("PSD efficiencies of IBD, atmospheric NC and fast neutron events")
plt.legend()
plt.grid()
plt.show()


