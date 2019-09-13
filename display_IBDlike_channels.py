""" Script to display the final channels (after deexcitation), the NC interaction channels (before deexcitation) and
    the deexcitation channels of IBD-like events:

    Moreover, also the ratios corresponding to the number of initial neutron and to the number of neutrons after
    deexcitation can be displayed.

    All NC events that mimic an IBD-like signature are stored in folder
    /home/astro/blum/juno/atmoNC/data_NC/output_detsim/.

    The information about the different channels is stored in folder
    /home/astro/blum/juno/atmoNC/data_NC/output_detsim/interaction_channels_IBDlike_signal/ and was analyzed with script
    IBDlike_events_channels.py

    The txt files in /interaction_channels_IBDlike_signal/ are summarized in file channels_IBDlike_signals.ods.

"""
import numpy as np
from matplotlib import pyplot as plt


def label_barh(axx, bars, text_format, is_inside=True, **kwargs):
    """
    Attach a text label to each horizontal bar displaying its y value
    """
    max_y_value = max(bar.get_height() for bar in bars)
    if is_inside:
        distance = max_y_value * 0.05
    else:
        distance = max_y_value * 0.1

    for bar in bars:
        text = text_format.format(bar.get_width())
        if is_inside:
            text_x = bar.get_width() - distance
        else:
            text_x = bar.get_width() + distance
        text_y = bar.get_y() + bar.get_height() / 2

        axx.text(text_x, text_y, text, va='center', **kwargs)


""" insert fractions of final channels of IBD-like events from channels_IBDlike_signals.ods, sheet 'combined channels'.
    (These fractions contain 88.83 % of all channels mimicking an IBD-like event): 
    Fractions in %:
"""
Frac_nu_C11_n_final = 27.10
Frac_nu_B10_n_p_final = 16.43
Frac_nu_Be9_n_2p_final = 11.90
Frac_nu_Be8_n_p_d_final = 8.86
Frac_nu_Li6_n_p_alpha_final = 8.19
Frac_nu_Be8_2n_2p_final = 3.98
Frac_nu_B9_2n_p_final = 3.77
Frac_nu_Li8_n_3p_final = 3.21
Frac_nu_C10_2n_final = 1.51
Frac_nu_B9_n_d_final = 1.51
Frac_nu_He7_n_4p_final = 1.21
Frac_nu_Li7_n_2p_d_final = 1.16

channels_final = ('$^7Li$ + $n$ + $2p$ + $d$', '$^7He$ + $n$ + $4p$', '$^9B$ + $n$ + $d$', '$^{10}C$ + $2n$',
                  '$^8Li$ + $n$ + $3p$', '$^9B$ + $2n$ + $p$', '$^8Be$ + $2n$ + $2p$', '$^6Li$ + $n$ + $p$ + $\\alpha$',
                  '$^8Be$ + $n$ + $p$ + $d$', '$^9Be$ + $n$ + $2p$', '$^{10}B$ + $n$ + $p$', '$^{11}C$ + $n$')

pos_final = np.arange(len(channels_final))

frac_final = np.array([Frac_nu_Li7_n_2p_d_final, Frac_nu_He7_n_4p_final, Frac_nu_B9_n_d_final, Frac_nu_C10_2n_final,
                       Frac_nu_Li8_n_3p_final, Frac_nu_B9_2n_p_final, Frac_nu_Be8_2n_2p_final,
                       Frac_nu_Li6_n_p_alpha_final, Frac_nu_Be8_n_p_d_final, Frac_nu_Be9_n_2p_final,
                       Frac_nu_B10_n_p_final, Frac_nu_C11_n_final])

""" display in bar chart """
fig, ax = plt.subplots()
horizontal_bars = ax.barh(pos_final, frac_final, align='center', alpha=0.9)
label_barh(ax, horizontal_bars, "{:.4}", is_inside=False, fontsize=12)
plt.xlim(xmin=0.0, xmax=100.0)
plt.xticks(fontsize=11)
plt.yticks(pos_final, channels_final, fontsize=12)
plt.xlabel('fraction in %', fontsize=12)
plt.title('Atmospheric NC interactions on $^{12}C$ that mimic an IBD signal\n'
          '(interaction channels: $\\nu_{x}$ + $^{12}C$ $\\rightarrow$ $\\nu_{x}$ + ...)', fontsize=15)
plt.grid(axis='x')

""" insert fractions of IBD-like events corresponding to the number of neutrons after deexcitation from 
    channels_IBDlike_signals.ods, sheet 'combined channels'.
    (These fractions contain 100 % of all channels mimicking an IBD-like event)
    Fractions in %:
"""
Frac_0_neutrons = 1.06
Frac_1_neutron = 85.74
Frac_2_neutrons = 12.13
Frac_3_neutrons = 0.83
Frac_4_or_more_neutrons = 0.25

channels_n = ('more than 3 neutrons', '3 neutrons', '2 neutrons', '1 neutron', '0 neutrons')

pos_n = np.arange(len(channels_n))

frac_n = np.array([Frac_4_or_more_neutrons, Frac_3_neutrons, Frac_2_neutrons, Frac_1_neutron, Frac_0_neutrons])

""" display in bar chart """
fig2, ax2 = plt.subplots()
horizontal_bars2 = ax2.barh(pos_n, frac_n, align='center', alpha=0.9)
label_barh(ax2, horizontal_bars2, "{:.4}", is_inside=False, fontsize=12)
plt.xlim(xmin=0.0, xmax=100.0)
plt.xticks(fontsize=11)
plt.yticks(pos_n, channels_n, fontsize=12)
plt.xlabel('fraction in %', fontsize=12)
plt.title('Number of neutrons in atmospheric NC events that mimic an IBD signal\n'
          '(interaction channels: $\\nu_{x}$ + $^{12}C$ $\\rightarrow$ $\\nu_{x}$ + ...)', fontsize=15)
plt.grid(axis='x')

plt.show()
