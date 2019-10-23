""" Script to check for which position inside JUNO detector total reflexion occur.

    If the position of the event is too close to the acrylic sphere, there can be total reflexion between the outer
    side of the acrylic sphere and the water, where the PMT are mounted.

    With this script, you can check different positions (distances to the detector center).

    The geometrical calculation are described in my handwritten notes ("total reflexion", from 27-09-2019)




"""
import numpy as np
from matplotlib import pyplot as plt


def alpha_angle(beta, radius):
    """ function to calculate the angle between the position of the event to the hit-point of the photon on the acrylic
    and the line perpendicular to the acrylic sphere

    :param beta: angle from center of detector between the hit-point on the acrylic and the event position in degree
    :param radius: maximum distance of the event position to the center of the detector in m
    :return: alpha: angle between the position of the event to the hit-point and the perpendicular line in degree
    """
    a_1 = np.sin(beta*np.pi/180.0) * radius

    b_1 = np.sqrt((17.7 * np.sqrt(1 - np.sin(beta*np.pi/180.0)**2) - 16.0)**2 + np.sin(beta*np.pi/180.0)**2 * 17.7**2)

    alpha = np.arcsin(a_1/b_1)

    alpha = alpha * 180.0 / np.pi

    return alpha


# refraction index of LS and acrylic:
n_acrylic = 1.49
# refraction index of water:
n_water = 1.33

# critical angle for total reflexion in degree (if alpha > angle_totalreflexion, there is total reflexion):
angle_totalreflexion = np.arcsin(n_water / n_acrylic) * 180.0 / np.pi

# angle in degree to scan the acrylic sphere:
Beta = np.arange(0, 91, 1)
# maximum distance of the event to the detector center in m:
R_1 = 17.5
R_2 = 17.0
R_3 = 16.0
R_4 = 15.8

# alpha in degree for different R:
Alpha_1 = alpha_angle(Beta, R_1)
Alpha_2 = alpha_angle(Beta, R_2)
Alpha_3 = alpha_angle(Beta, R_3)
Alpha_4 = alpha_angle(Beta, R_4)

plt.figure(1, figsize=(15, 8))
plt.plot(Beta, Alpha_1, "b", label="R = {0:.1f} m".format(R_1))
plt.plot(Beta, Alpha_2, "r", label="R = {0:.1f} m".format(R_2))
plt.plot(Beta, Alpha_3, "g", label="R = {0:.1f} m".format(R_3))
plt.plot(Beta, Alpha_4, color="orange", label="R = {0:.1f} m".format(R_4))
plt.hlines(angle_totalreflexion, min(Beta), max(Beta), "k", label="angle of total reflection = {0:.2f} degrees"
           .format(angle_totalreflexion))
plt.xticks(np.arange(min(Beta), max(Beta)+5, 5))
plt.xlim(xmin=min(Beta), xmax=max(Beta))
plt.xlabel("$\\beta$ in degrees")
plt.ylabel("$\\alpha$ in degrees")
plt.title("Investigation of total reflection depending on the position in JUNO detector")
plt.legend()
plt.grid()
plt.show()
