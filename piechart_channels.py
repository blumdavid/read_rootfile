# script to display the interaction channels, deexcitation channels and combined channels in a pie chart:

import matplotlib.pyplot as plt

""" combined channels (NC interaction and deexcitation): """
# data from folder /home/astro/blum/juno/atmoNC/data_NC/output_checkoutNC/ from file combined_channels_250000evts.ods

# Make data: I have 10 groups and 20 subgroups
# group_names = ['$^{11}$B', '$^{11}$C', '$^{10}$B', '$^{10}$C', '$^{10}$Be', 'Be8', 'Be9', 'B9', 'Li6', 'rest']
group_size_1 = [16.4, 15.1, 9.7, 2.1, 2.3, 8.5, 7.7, 6.0, 5.7, 26.5]
# subgroup_names = ['+p', '+n+$\\pi^+$', 'rest',
#                   '+n', '+p+$\\pi^-$', 'rest',
#                   '+p+n', '',
#                   'all',
#                   'all',
#                   '+2n+2p', '+n+p+d', 'rest',
#                   '+n+2p', 'rest',
#                   '+2n+p', 'rest',
#                   '+n+p+alpha', 'rest']
subgroup_size_1 = [14.0, 1.8, 0.6, 12.5, 2.2, 0.4, 8.3, 1.4, 2.1, 2.3, 3.6, 3.3, 1.6, 5.4, 2.3, 4.2, 1.8, 3.0, 2.7,
                   26.5]

# Create colors
a, b, c, d, e, f = [plt.cm.Blues, plt.cm.Reds, plt.cm.Greens, plt.cm.Purples, plt.cm.Oranges, plt.cm.Greys]

""" without labels: """
# First Ring (inside)
fig1, ax1 = plt.subplots(1, figsize=(9, 9))
ax1.axis('equal')
mypie_1, _ = ax1.pie(group_size_1, radius=0.8, colors=[a(0.9), b(0.9), c(0.9), d(0.9), e(0.9), a(0.5), b(0.5), c(0.5),
                                                       d(0.5), f(0.7)])
plt.setp(mypie_1, width=0.5, edgecolor='white', linewidth=1)

# Second Ring (outside)
mypie2_1, _ = ax1.pie(subgroup_size_1, radius=0.8+0.3,
                      colors=[a(0.8), a(0.7), a(0.6), b(0.8), b(0.7), b(0.6), c(0.8), c(0.7), d(0.8), e(0.8), a(0.4),
                              a(0.3), a(0.2), b(0.4), b(0.3), c(0.4), c(0.3), d(0.4), d(0.3), f(0.0)])
plt.setp(mypie2_1, width=0.3, edgecolor='white', linewidth=0.5)

plt.margins(0, 0)

# show it
plt.show()


""" NC interaction channels channels: """
# data from folder /home/astro/blum/juno/atmoNC/ from file interaction_channels_edit.ods, sheet
# "Interaction Channels only 250000 events"

# Make data: I have 6 groups and 17 subgroups
# group_names = ['$^{11}$B', '$^{11}$C', '$^{10}$B', '$^{10}$C', '$^{10}$Be', 'lighter\nisotopes']
group_size_2 = [24.6, 22.7, 16.4, 4.4, 4.5, 27.4]
# subgroup_names = ['+p', '+n+$\\pi^+$', '',
#                   '+n', '+p+$\\pi^-$', '',
#                   '+p+n', '+2p+$\\pi^-$', '+2n+$\\pi^+$', '',
#                   '+2n', '+p+n+$\\pi^-$', '',
#                   '+2p', '+p+n+$\\pi^+$', '',
#                   '']
subgroup_size_2 = [21.0, 2.7, 0.9, 18.8, 3.3, 0.6, 14.0, 0.9, 0.8, 0.7, 3.2, 1.0, 0.2, 3.3, 1.0, 0.2, 27.4]

# Create colors
a, b, c, d, e, f = [plt.cm.Blues, plt.cm.Reds, plt.cm.Greens, plt.cm.Purples, plt.cm.Oranges, plt.cm.Greys]

""" without labels: """
# First Ring (inside)
fig2, ax2 = plt.subplots(1, figsize=(9, 9))
ax2.axis('equal')
mypie_2, _ = ax2.pie(group_size_2, radius=0.8, colors=[a(0.9), b(0.9), c(0.9), d(0.9), e(0.9), f(0.7)])
plt.setp(mypie_2, width=0.5, edgecolor='white', linewidth=1)

# Second Ring (outside)
mypie2_2, _ = ax2.pie(subgroup_size_2, radius=0.8+0.3,
                      colors=[a(0.8), a(0.7), a(0.6), b(0.8), b(0.7), b(0.6), c(0.8), c(0.7), c(0.6), c(0.5), d(0.8),
                              d(0.7), d(0.6), e(0.8), e(0.7), e(0.6), f(0.0)])
plt.setp(mypie2_2, width=0.3, edgecolor='white', linewidth=0.5)

plt.margins(0, 0)

# show it
plt.show()


""" Deexcitation channels: """
# data from folder /home/astro/blum/juno/atmoNC/ from deex_channels_info.ods, sheet
# 'NEW! NC only C12 with GAMMA, NO pi0'

# Make data: I have 11 groups and 20 subgroups
# group_names = ['B11', 'B11*', 'C11', 'C11*', 'B10', 'B10*', 'C10', 'C10*', 'Be10', 'Be10*', 'lighter\nisotopes']
group_size_3 = [13.1, 11.5, 12.0, 10.7, 7.2, 9.2, 1.8, 2.6, 1.8, 2.7, 27.4]
# subgroup_names = ['not excited', 'B11+gamma', 'Be8+n+d+gamma', 'Li6+n+alpha+gamma', 'rest',
#                   'not excited', 'C11 + gamma', 'Be8+p+d+gamma', 'Li6+p+alpha+gamma', 'rest',
#                   'not excited', 'Be8+n+p+gamma', 'Be9+p+gamma', 'B9+n+gamma', 'rest',
#                   'not excited', 'excited',
#                   'not excited', 'excited']
subgroup_size_3 = [13.1, 3.3, 1.5, 1.5, 5.2, 12.0, 3.1, 1.7, 1.5, 4.4, 7.2, 2.2, 2.0, 1.7, 3.3,
                   1.8, 2.6, 1.8, 2.7, 27.4]

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
fig3, ax3 = plt.subplots(1, figsize=(9, 9))
ax3.axis('equal')
mypie_3, _ = ax3.pie(group_size_3, radius=0.8, colors=[a(0.9), a(0.8), b(0.9), b(0.8), c(0.9), c(0.8), d(0.9), d(0.8),
                                                       e(0.9), e(0.8), f(0.7)])
plt.setp(mypie_3, width=0.5, edgecolor='white', linewidth=1)

# Second Ring (outside)
mypie2_3, _ = ax3.pie(subgroup_size_3, radius=0.8+0.3,
                      colors=[a(0.0), a(0.7), a(0.6), a(0.5), a(0.4), b(0.0), b(0.7), b(0.6), b(0.5), b(0.4),
                              c(0.0), c(0.7), c(0.6), c(0.5), c(0.4), d(0.9), d(0.8), e(0.9), e(0.8), f(0.0)])
plt.setp(mypie2_3, width=0.3, edgecolor='white', linewidth=0.5)

plt.margins(0, 0)

# plt.title("Rates of atmospheric NC neutrino interaction channels on $^{12}$C")

# show it
plt.show()
