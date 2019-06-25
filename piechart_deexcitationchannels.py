# script to display the deexcitation channels in a pie chart:

import matplotlib.pyplot as plt

# data from deex_channels_info.ods, sheet 1 'NEW!! NC only C12, 250000 evts'

# Make data: I have 4 groups and 7 subgroups
# group_names = ['B11', 'B11*', 'C11', 'C11*', 'B10', 'B10*', 'C10', 'C10*', 'Be10', 'Be10*', 'lighter\nisotopes']
group_size = [16.3, 8.2, 15.1, 7.6, 7.2, 9.2, 1.8, 2.6, 1.8, 2.7, 27.5]
# subgroup_names = ['not excited', 'n+d+Be8', 'n+alpha+Li6', 'n+B10', 'rest',
#                   'not excited', 'p+d+Be8', 'p+alpha+Li6', '2p+Be9', 'rest',
#                   'not excited', 'n+p+Be8', 'p+Be9', 'n+B9', 'rest',
#                   'not excited', 'p+d+Be7', '2p+Be8', 'rest'
#                   'not excited', 'n+d+Li7', '2n+p+Li7', 'rest']
subgroup_size = [16.3, 1.5, 1.5, 1.4, 3.8, 15.1, 1.7, 1.5, 1.2, 3.2, 7.2, 2.2, 2.0, 1.7, 3.3, 1.8, 0.4, 0.3, 1.9,
                 1.8, 0.4, 0.4, 1.9, 27.5]

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
fig, ax = plt.subplots(1, figsize=(9, 9))
ax.axis('equal')
mypie, _ = ax.pie(group_size, radius=0.8, colors=[a(0.9), a(0.8), b(0.9), b(0.8), c(0.9), c(0.8), d(0.9), d(0.8),
                                                  e(0.9), e(0.8), f(0.7)])
plt.setp(mypie, width=0.5, edgecolor='white', linewidth=1)

# Second Ring (outside)
mypie2, _ = ax.pie(subgroup_size, radius=0.8+0.3,
                   colors=[a(0.0), a(0.7), a(0.6), a(0.5), a(0.4), b(0.0), b(0.7), b(0.6), b(0.5), b(0.4),
                           c(0.0), c(0.7), c(0.6), c(0.5), c(0.4), d(0.0), d(0.7), d(0.6), d(0.5),
                           e(0.0), e(0.7), e(0.6), e(0.5), f(0.0)])
plt.setp(mypie2, width=0.3, edgecolor='white', linewidth=0.5)

plt.margins(0, 0)

# plt.title("Rates of atmospheric NC neutrino interaction channels on $^{12}$C")

# show it
plt.show()
