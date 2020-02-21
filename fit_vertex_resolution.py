""" script to fit parametrization of vertex resolution to the data points from plot of page 20 of
    'review of vertex reconstruction in 2020' (vertex resolution calculated with dark noise and TTS
    (2.8 ns for Hamamatsu and 18 ns for NNVT):
"""
from matplotlib import pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

# get values of Z_rec - Z_edep resolution in mm for E_edep in MeV:
# deposit energy in MeV:
edep_data = np.array([1.02237, 1.63904, 2.57671, 3.54817, 4.54498, 5.53333, 6.53014, 7.52694])

# Z_rec - Z_edep resolution in mm:
resolution_data = np.array([112.344, 85.6044, 65.8242, 55.3846, 49.8901, 46.7766, 43.6630, 42.9304])


def fit_function(energy, a, b, c):
    """
    fit function (parametrization of vertex resolution) for data points:

    :param energy: Edep in MeV
    :param a: fit parameter a
    :param b: fit parameter b
    :param c: fit parameter c
    :return:
    """
    # parametrization function in mm:
    res = np.sqrt((a / np.sqrt(energy))**2 + b**2 + (c / energy)**2)

    return res


def fit_function_2(energy, a):
    """
    fit function (parametrization of vertex resolution) for data points:

    :param energy: Edep in MeV
    :param a: fit parameter a
    :return:
    """
    # parametrization function in mm:
    res = a / np.sqrt(energy)

    return res


def fit_function_3(energy, a, b):
    """
    fit function (parametrization of vertex resolution) for data points;
    :param energy: Edep in MeV
    :param a: fit parameter a
    :param b: fit parameter b
    :return:
    """
    # parametrization function in mm:
    res = np.sqrt((a / np.sqrt(energy))**2 + b**2)

    return res


# fit data with fit function:
fit_parameter, pcov = curve_fit(fit_function, edep_data, resolution_data, bounds=(0, np.inf))

# fit data with fit function 2:
fit_parameter_2, pcov_2 = curve_fit(fit_function_2, edep_data, resolution_data, bounds=(0, np.inf))

# fit data with fit function 3:
fit_parameter_3, pcov_3 = curve_fit(fit_function_3, edep_data, resolution_data, bounds=(0, np.inf))

# energy array in MeV:
e_array = np.arange(0.5, 10.0, 0.1)

# print fit-parameter:
print("best-fit parameters:")
print("a = {0:.2f}".format(fit_parameter[0]))
print("b = {0:.2f}".format(fit_parameter[1]))
print("c = {0:.2f}".format(fit_parameter[2]))

print("best-fit parameters 2:")
print("a = {0:.2f}".format(fit_parameter_2[0]))

print("best-fit parameters 3:")
print("a = {0:.2f}".format(fit_parameter_3[0]))
print("b = {0:.2f}".format(fit_parameter_3[1]))
print(fit_function_2(8, 110.16))

plt.plot(edep_data, resolution_data, "xk", label="data points")
plt.plot(e_array, fit_function(e_array, *fit_parameter), "-b",
         label="a = {0:.2f},\nb = {1:.2f},\nc = {2:.2f}".format(fit_parameter[0], fit_parameter[1], fit_parameter[2]))
plt.plot(e_array, fit_function_2(e_array, *fit_parameter_2), "-r",
         label="a = {0:.2f}".format(fit_parameter_2[0]))
plt.plot(e_array, fit_function_3(e_array, *fit_parameter_3), "-g",
         label="a = {0:.2f},\nb = {1:.2f}".format(fit_parameter_3[0], fit_parameter_3[1]))
plt.xlabel("Edep in MeV")
plt.ylabel("Z_rec - Z_edep resolution in mm")
plt.title("Vertex resolution")
plt.legend()
plt.grid()
plt.show()



