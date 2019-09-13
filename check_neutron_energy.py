""" Script to check the simulated neutron files from conversion_nPE_MeV/:

    1.  get total number of events
    2.  get number of events with 10 MeV <= Qedep <= 100 MeV
    3.  store Qedep between 10 MeV and 100 MeV in array and build histogram to check distribution

"""
import numpy as np
from matplotlib import pyplot as plt

# load all files, where neutrons have been simulated:
Qedep_1 = np.loadtxt("/home/astro/blum/juno/atmoNC/data_NC/conversion_nPE_MeV/qedep_n_10MeV.txt")
Qedep_2 = np.loadtxt("/home/astro/blum/juno/atmoNC/data_NC/conversion_nPE_MeV/qedep_n_100MeV.txt")
Qedep_3 = np.loadtxt("/home/astro/blum/juno/atmoNC/data_NC/conversion_nPE_MeV/qedep_n_300MeV.txt")
Qedep_4 = np.loadtxt("/home/astro/blum/juno/atmoNC/data_NC/conversion_nPE_MeV/qedep_n_500MeV.txt")
Qedep_5 = np.loadtxt("/home/astro/blum/juno/atmoNC/data_NC/conversion_nPE_MeV/qedep_n_500MeV_2.txt")
Qedep_6 = np.loadtxt("/home/astro/blum/juno/atmoNC/data_NC/conversion_nPE_MeV/qedep_n_500MeV_3.txt")
Qedep_7 = np.loadtxt("/home/astro/blum/juno/atmoNC/data_NC/conversion_nPE_MeV/qedep_n_500MeV_4.txt")
Qedep_8 = np.loadtxt("/home/astro/blum/juno/atmoNC/data_NC/conversion_nPE_MeV/qedep_n_500MeV_5.txt")
Qedep_9 = np.loadtxt("/home/astro/blum/juno/atmoNC/data_NC/conversion_nPE_MeV/qedep_n_500MeV_6.txt")
Qedep_10 = np.loadtxt("/home/astro/blum/juno/atmoNC/data_NC/conversion_nPE_MeV/qedep_n_500MeV_7.txt")
Qedep_11 = np.loadtxt("/home/astro/blum/juno/atmoNC/data_NC/conversion_nPE_MeV/qedep_n_1GeV.txt")

# total number of events in the files:
number = 0
# number of events with 10 MeV <= Qedep <= 100 MeV:
number_10_to_100 = 0
# array, where Qedep between 10 MeV and 100 MeV are stored:
Qedep_10_to_100 = []

# loop over entries in Qedep_1:
number += len(Qedep_1)
for index in range(len(Qedep_1)):
    if 10.0 <= Qedep_1[index] <= 100.0:
        number_10_to_100 += 1
        Qedep_10_to_100.append(Qedep_1[index])

print("number of analyzed events = {0:d}".format(len(Qedep_1)))
print("max of Qedep = {0:.1f} MeV".format(max(Qedep_1)))
print("min of Qedep = {0:.1f} MeV".format(min(Qedep_1)))

# loop over entries in Qedep_2:
number += len(Qedep_2)
for index in range(len(Qedep_2)):
    if 10.0 <= Qedep_2[index] <= 100.0:
        number_10_to_100 += 1
        Qedep_10_to_100.append(Qedep_2[index])

print("number of analyzed events = {0:d}".format(len(Qedep_2)))
print("max of Qedep = {0:.1f} MeV".format(max(Qedep_2)))
print("min of Qedep = {0:.1f} MeV".format(min(Qedep_2)))

# loop over entries in Qedep_3:
number += len(Qedep_3)
for index in range(len(Qedep_3)):
    if 10.0 <= Qedep_3[index] <= 100.0:
        number_10_to_100 += 1
        Qedep_10_to_100.append(Qedep_3[index])

print("number of analyzed events = {0:d}".format(len(Qedep_3)))
print("max of Qedep = {0:.1f} MeV".format(max(Qedep_3)))
print("min of Qedep = {0:.1f} MeV".format(min(Qedep_3)))

# loop over entries in Qedep_4:
number += len(Qedep_4)
for index in range(len(Qedep_4)):
    if 10.0 <= Qedep_4[index] <= 100.0:
        number_10_to_100 += 1
        Qedep_10_to_100.append(Qedep_4[index])

print("number of analyzed events = {0:d}".format(len(Qedep_4)))
print("max of Qedep = {0:.1f} MeV".format(max(Qedep_4)))
print("min of Qedep = {0:.1f} MeV".format(min(Qedep_4)))

# loop over entries in Qedep_5:
number += len(Qedep_5)
for index in range(len(Qedep_5)):
    if 10.0 <= Qedep_5[index] <= 100.0:
        number_10_to_100 += 1
        Qedep_10_to_100.append(Qedep_5[index])

print("number of analyzed events = {0:d}".format(len(Qedep_5)))
print("max of Qedep = {0:.1f} MeV".format(max(Qedep_5)))
print("min of Qedep = {0:.1f} MeV".format(min(Qedep_5)))

# loop over entries in Qedep_6:
number += len(Qedep_6)
for index in range(len(Qedep_6)):
    if 10.0 <= Qedep_6[index] <= 100.0:
        number_10_to_100 += 1
        Qedep_10_to_100.append(Qedep_6[index])

print("number of analyzed events = {0:d}".format(len(Qedep_6)))
print("max of Qedep = {0:.1f} MeV".format(max(Qedep_6)))
print("min of Qedep = {0:.1f} MeV".format(min(Qedep_6)))

# loop over entries in Qedep_7:
number += len(Qedep_7)
for index in range(len(Qedep_7)):
    if 10.0 <= Qedep_7[index] <= 100.0:
        number_10_to_100 += 1
        Qedep_10_to_100.append(Qedep_7[index])

print("number of analyzed events = {0:d}".format(len(Qedep_7)))
print("max of Qedep = {0:.1f} MeV".format(max(Qedep_7)))
print("min of Qedep = {0:.1f} MeV".format(min(Qedep_7)))

# loop over entries in Qedep_8:
number += len(Qedep_8)
for index in range(len(Qedep_8)):
    if 10.0 <= Qedep_8[index] <= 100.0:
        number_10_to_100 += 1
        Qedep_10_to_100.append(Qedep_8[index])

print("number of analyzed events = {0:d}".format(len(Qedep_8)))
print("max of Qedep = {0:.1f} MeV".format(max(Qedep_8)))
print("min of Qedep = {0:.1f} MeV".format(min(Qedep_8)))

# loop over entries in Qedep_9:
number += len(Qedep_9)
for index in range(len(Qedep_9)):
    if 10.0 <= Qedep_9[index] <= 100.0:
        number_10_to_100 += 1
        Qedep_10_to_100.append(Qedep_9[index])

print("number of analyzed events = {0:d}".format(len(Qedep_9)))
print("max of Qedep = {0:.1f} MeV".format(max(Qedep_9)))
print("min of Qedep = {0:.1f} MeV".format(min(Qedep_9)))

# loop over entries in Qedep_10:
number += len(Qedep_10)
for index in range(len(Qedep_10)):
    if 10.0 <= Qedep_10[index] <= 100.0:
        number_10_to_100 += 1
        Qedep_10_to_100.append(Qedep_10[index])

print("number of analyzed events = {0:d}".format(len(Qedep_10)))
print("max of Qedep = {0:.1f} MeV".format(max(Qedep_10)))
print("min of Qedep = {0:.1f} MeV".format(min(Qedep_10)))

# loop over entries in Qedep_11:
number += len(Qedep_11)
for index in range(len(Qedep_11)):
    if 10.0 <= Qedep_11[index] <= 100.0:
        number_10_to_100 += 1
        Qedep_10_to_100.append(Qedep_11[index])

print("number of analyzed events = {0:d}".format(len(Qedep_11)))
print("max of Qedep = {0:.1f} MeV".format(max(Qedep_11)))
print("min of Qedep = {0:.1f} MeV".format(min(Qedep_11)))

print("---------------")
print("number total = {0:d}".format(number))
print("number between 10 and 100 MeV = {0:d}".format(number_10_to_100))

# build histogram from Qedep_10_to_100:
h1 = plt.figure(1, figsize=(15, 8))
bin_width = 1
bins = np.arange(10, 100 + bin_width, bin_width)
plt.hist(Qedep_10_to_100, bins, label="entries = {0:d}".format(number_10_to_100))
plt.xlim(xmin=10.0, xmax=100.0)
plt.xlabel("quenched deposit energy in MeV")
plt.ylabel("events per bin (bin-width = {0:.1f} MeV)".format(bin_width))
plt.title("energy spectrum of the simulated neutron events \nused for pulse shape analysis of fast neutron background")
plt.legend()
plt.grid()
plt.show()
