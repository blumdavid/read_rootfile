""" script to calculate the figure of merit of Pulse Shape Discrimination between IBD events and IBD-like
    atmospheric NC events:

    Folder: /home/astro/blum/juno/atmoNC/data_NC/output_detsim_v2/
            DCR_results_16000mm_10MeVto100MeV_500nsto1ms_mult1_2400PEto3400PE_dist500mm_R17700mm_PSD99/

    TTR values of IBD events are taken from TTR_beforePSD_IBDevents_275ns_to_600ns.txt

    TTR values of IBD-like atmospheric NC background events are taken from TTR_IBDlike_NCevents_275ns_to_600ns.txt


"""
import datetime
import numpy as np
from scipy.stats import norm
from matplotlib import pyplot as plt

# get the date and time, when the script was run:
date = datetime.datetime.now()
now = date.strftime("%Y-%m-%d %H:%M")

# set input path, where TTR values are stored:
input_path = "/home/astro/blum/juno/atmoNC/data_NC/output_detsim_v2/" \
             "DCR_results_16000mm_10MeVto100MeV_500nsto1ms_mult1_2400PEto3400PE_dist500mm_R17700mm_PSD99/"

# filename for TTR values of IBD events:
filename_IBD = input_path + "TTR_beforePSD_IBDevents_275ns_to_600ns.txt"
# filename for TTR values of IBD-like events:
filename_IBDlike = input_path + "TTR_IBDlike_NCevents_275ns_to_600ns.txt"

# TTR cut value for NC suppression = 98.98 % and IBD suppression = 11.08 %:
TTR_cut_value = 0.01662

# read files:
array_ttr_IBD = np.loadtxt(filename_IBD)
array_ttr_IBDlike = np.loadtxt(filename_IBDlike)

# get minimum TTR value:
if min(array_ttr_IBD) <= min(array_ttr_IBDlike):
    ttr_min = min(array_ttr_IBD)
else:
    ttr_min = min(array_ttr_IBDlike)
# get maximum TTR value:
if max(array_ttr_IBD) >= max(array_ttr_IBDlike):
    ttr_max = max(array_ttr_IBD)
else:
    ttr_max = max(array_ttr_IBDlike)

# calculate difference ttr_max - ttr_min:
ttr_diff = ttr_max - ttr_min
# divide ttr_diff by 100 to get bin_width that corresponds to 100 bins:
bin_width = ttr_diff / 300

# fit gaussian to array to get mu and sigma:
(mu_IBD, sigma_IBD) = norm.fit(array_ttr_IBD)
(mu_IBDlike, sigma_IBDlike) = norm.fit(array_ttr_IBDlike)

""" build histogram with array_ttr_IBD and array_ttr_IBDlike: """
h1 = plt.figure(1, figsize=(15, 8))
bins_ttr = np.arange(ttr_min, ttr_max + bin_width, bin_width)
# build histogram of array_ttr_IBD:
hist_ttr_IBD, bin_edges, patches = plt.hist(array_ttr_IBD, bins=bins_ttr, density=True, histtype="step", color="r",
                                            label="IBD events (entries = {0:d})".format(len(array_ttr_IBD)))
# build histogram of array_ttr_IBDlike:
hist_ttr_IBDlike, bin_edges, patches = plt.hist(array_ttr_IBDlike, bins=bins_ttr, density=True, histtype="step",
                                                color="b",
                                                label="IBD-like NC events (entries = {0:d})"
                                                .format(len(array_ttr_IBDlike)))

# add 'best fit' line for IBD events:
y_IBD = norm.pdf(bins_ttr, mu_IBD, sigma_IBD)
plt.plot(bins_ttr, y_IBD, 'r--', linewidth=2, label="IBD events: $\\mu$ = {0:.5f}, $\\sigma$ = {1:.5f}"
         .format(mu_IBD, sigma_IBD))
# add 'best fit' line for IBDlike events:
y_IBDlike = norm.pdf(bins_ttr, mu_IBDlike, sigma_IBDlike)
plt.plot(bins_ttr, y_IBDlike, 'b--', linewidth=2, label="IBDlike events: $\\mu$ = {0:.5f}, $\\sigma$ = {1:.5f}"
         .format(mu_IBDlike, sigma_IBDlike))

# calculate the FWHM of the best fit lines:
fwhm_IBD = 2 * np.sqrt(2 * np.log(2)) * sigma_IBD
fwhm_IBDlike = 2 * np.sqrt(2 * np.log(2)) * sigma_IBDlike

# calculate figure of merit (see https://www.sciencedirect.com/science/article/pii/S0168900215004532):
FOM = (mu_IBDlike - mu_IBD) / (fwhm_IBD + fwhm_IBDlike)

plt.xlabel('TTR value')
plt.ylabel('Probability')
plt.title("TTR of prompt signals of IBD events and IBD-like atmospheric NC events\n(figure of merit = {0:.5f})"
          .format(FOM))
plt.grid()
plt.legend()

plt.show()

