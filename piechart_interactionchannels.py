# script to display the interaction channels in a pie chart:

import matplotlib.pyplot as plt

# data from file interaction_channels_edit.ods, sheet 4 'Interaction Channels after Generator'

# Make data: I have 6 groups and 17 subgroups
# group_names = ['$^{11}$B', '$^{11}$C', '$^{10}$B', '$^{10}$C', '$^{10}$Be', 'lighter\nisotopes']
group_size = [24.6, 22.7, 16.4, 4.4, 4.5, 27.4]
# subgroup_names = ['+p', '+n+$\\pi^+$', '',
#                   '+n', '+p+$\\pi^-$', '',
#                   '+p+n', '+2p+$\\pi^-$', '+2n+$\\pi^+$', '',
#                   '+2n', '+p+n+$\\pi^-$', '',
#                   '+2p', '+p+n+$\\pi^+$', '',
#                   '']
subgroup_size = [21.0, 2.7, 0.9, 18.8, 3.3, 0.6, 14.0, 0.9, 0.8, 0.7, 3.2, 1.0, 0.2, 3.3, 1.0, 0.2, 27.4]

# Create colors
a, b, c, d, e, f = [plt.cm.Blues, plt.cm.Reds, plt.cm.Greens, plt.cm.Purples, plt.cm.Oranges, plt.cm.Greys]

""" with labels: """
"""
# First Ring (inside)
fig, ax = plt.subplots()
ax.axis('equal')
mypie, _ = ax.pie(group_size, radius=0.8, labels=group_names, labeldistance=0.7,
                  colors=[a(0.7), b(0.7), c(0.7), d(0.7), e(0.7), f(0.7)])
plt.setp(mypie, width=0.7, edgecolor='white', linewidth=2)

# Second Ring (outside)
mypie2, _ = ax.pie(subgroup_size, radius=0.8+0.7, labels=subgroup_names, labeldistance=0.75, rotatelabels=True,
                   colors=[a(0.6), a(0.5), a(0.4), b(0.6), b(0.5), b(0.4), c(0.6), c(0.5), c(0.4), c(0.3), d(0.6),
                           d(0.5), d(0.4), e(0.6), e(0.5), e(0.4), f(0.0)])
plt.setp(mypie2, width=0.7, edgecolor='white')

plt.margins(0, 0)

# show it
plt.show()
"""

""" without labels: """
# First Ring (inside)
fig1, ax1 = plt.subplots(1, figsize=(9, 9))
ax1.axis('equal')
mypie, _ = ax1.pie(group_size, radius=0.8, colors=[a(0.9), b(0.9), c(0.9), d(0.9), e(0.9), f(0.7)])
plt.setp(mypie, width=0.5, edgecolor='white', linewidth=1)

# Second Ring (outside)
mypie2, _ = ax1.pie(subgroup_size, radius=0.8+0.3,
                    colors=[a(0.8), a(0.7), a(0.6), b(0.8), b(0.7), b(0.6), c(0.8), c(0.7), c(0.6), c(0.5), d(0.8),
                            d(0.7), d(0.6), e(0.8), e(0.7), e(0.6), f(0.0)])
plt.setp(mypie2, width=0.3, edgecolor='white', linewidth=0.5)

plt.margins(0, 0)

# show it
plt.show()
