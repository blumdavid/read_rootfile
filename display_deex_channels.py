""" script to display the different deexcitation channels from DSNB-NC generation in a bar chart:

    The output root file of the DSNB-NC generator (gen_NC_onlyC12_250000evts_seed1.root) is analyzed with
    get_interaction_and_deex_channel.py and the number of events for different deexcitation channels are saved file
    "deex_channels_info.ods"

    Here, two fractions are calculated for each deexcitation channel:
        - number of events / total number of events on C12
        - number pf events / total number of events, where e.g. C11 is produced


    INFO:   In all events pi0 and Gammas are considered!!!!
            Maybe some event, which directly deexcite into pi0 and/or gamma's are overestimated, because you can
            not specify, if they come from the interaction or deexcitation (not set in t_channelID)

"""
import numpy as np
from matplotlib import pyplot as plt


def label_barh(ax, bars, text_format, is_inside=True, **kwargs):
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

        ax.text(text_x, text_y, text, va='center', **kwargs)


""" insert the fractions of the different deexcitation channels of residual isotopes from 'deex_channels_info.ods', 
    sheet 'NEW! NC only C12 with pi0 and Gamma': """

# C11: fractions in % (number of events / total number of events, where C11 is excited) (fractions below contain 91.6 %
# of all C11 deex channels) (fraction of number of total C11 from 250000 evts = 22.72 %):
frac_C11_C11_gamma = 28.78
frac_C11_p_d_Be8_gamma = 16.03
frac_C11_p_alpha_Li6_gamma = 14.32
frac_C11_2p_Be9_gamma = 11.29
frac_C11_p_B10_gamma = 9.34
frac_C11_d_B9_gamma = 6.77
frac_C11_n_p_B9_gamma = 5.10
frac_C11_sum = frac_C11_C11_gamma + frac_C11_p_d_Be8_gamma + frac_C11_p_alpha_Li6_gamma + frac_C11_2p_Be9_gamma + \
               frac_C11_p_B10_gamma + frac_C11_d_B9_gamma + frac_C11_n_p_B9_gamma

ch_C11 = ('$^9B$ + $n$ + $p$ + $\\gamma$', '$^9B$ + $d$ + $\\gamma$', '$^{10}B$ + $p$ + $\\gamma$',
          '$^9Be$ + $2p$ + $\\gamma$', '$^6Li$ + $p$ + $\\alpha$ + $\\gamma$', '$^8Be$ + $p$ + $d$ + $\\gamma$',
          '$^{11}C$ + $\\gamma$')

pos_C11 = np.arange(len(ch_C11))

frac_C11 = np.array([frac_C11_n_p_B9_gamma, frac_C11_d_B9_gamma, frac_C11_p_B10_gamma, frac_C11_2p_Be9_gamma,
                     frac_C11_p_alpha_Li6_gamma, frac_C11_p_d_Be8_gamma, frac_C11_C11_gamma])


# B11: fractions in % (number of events / total number of events, where B11 is excited) (fractions below contain 86.0 %
# of all B11 deex channels) (fraction of number of total B11 from 250000 evts = 24.55 %):
frac_B11_B11_gamma = 28.96
frac_B11_n_d_Be8_gamma = 13.38
frac_B11_n_alpha_Li6_gamma = 12.80
frac_B11_n_B10_gamma = 12.32
frac_B11_n_p_Be9_gamma = 9.94
frac_B11_d_Be9_gamma = 8.63

frac_B11_sum = frac_B11_B11_gamma + frac_B11_d_Be9_gamma + frac_B11_n_p_Be9_gamma + frac_B11_n_B10_gamma + \
               frac_B11_n_alpha_Li6_gamma + frac_B11_n_d_Be8_gamma

ch_B11 = ('$^9Be$ + $d$ + $\\gamma$', '$^9Be$ + $n$ + $p$ + $\\gamma$', '$^{10}B$ + $n$ + $\\gamma$',
          '$^6Li$ + $n$ + $\\alpha$ + $\\gamma$', '$^8Be$ + $n$ + $d$ + $\\gamma$', '$^{11}B$ + $\\gamma$')

pos_B11 = np.arange(len(ch_B11))

frac_B11 = np.array([frac_B11_d_Be9_gamma, frac_B11_n_p_Be9_gamma, frac_B11_n_B10_gamma, frac_B11_n_alpha_Li6_gamma,
                     frac_B11_n_d_Be8_gamma, frac_B11_B11_gamma])


# B10: fractions in % (number of events / total number of events, where B10 is excited) (fractions below contain 92.7 %
# of all B10 deex channels) (fraction of number of total B10 from 250000 evts = 16.41 %):
frac_B10_n_p_Be8_gamma = 24.09
frac_B10_p_Be9_gamma = 22.08
frac_B10_n_B9_gamma = 18.52
frac_B10_d_Be8_gamma = 7.58
frac_B10_He3_Li7_gamma = 7.00
frac_B10_p_alpha_He5_gamma = 6.93
frac_B10_t_Be7_gamma = 6.56

frac_B10_sum = frac_B10_t_Be7_gamma + frac_B10_He3_Li7_gamma + frac_B10_p_alpha_He5_gamma + frac_B10_d_Be8_gamma + \
               frac_B10_n_B9_gamma + frac_B10_p_Be9_gamma + frac_B10_n_p_Be8_gamma

ch_B10 = ('$^7Be$ + $t$ + $\\gamma$', '$^5He$ + $p$ + $\\alpha$ + $\\gamma$', '$^7Li$ + $^3He$ + $\\gamma$',
          '$^8Be$ + $d$ + $\\gamma$', '$^9B$ + $n$ + $\\gamma$', '$^9Be$ + $p$ + $\\gamma$',
          '$^8Be$ + $n$ + $p$ + $\\gamma$')

pos_B10 = np.arange(len(ch_B10))

frac_B10 = np.array([frac_B10_t_Be7_gamma, frac_B10_p_alpha_He5_gamma, frac_B10_He3_Li7_gamma, frac_B10_d_Be8_gamma,
                     frac_B10_n_B9_gamma, frac_B10_p_Be9_gamma, frac_B10_n_p_Be8_gamma])

# C10: fraction of number of total C10 from 250000 evts = 4.46 %
# Be10: fraction of number of total C10 from 250000 evts = 4.47 %
# lighter isotopes: fraction of number of total isotope from 250000 evts < 2.7 %


""" display in bar chart """
# C11:
fig_c11, ax_c11 = plt.subplots(figsize=(11, 6))
horizontal_bars_c11 = ax_c11.barh(pos_C11, frac_C11, align='center', alpha=0.9)
label_barh(ax_c11, horizontal_bars_c11, "{:.4}", is_inside=False, fontsize=12)
plt.xticks(fontsize=11)
plt.yticks(pos_C11, ch_C11, fontsize=12)
plt.xlabel('fraction in %', fontsize=12)
plt.title('{0:.1f}'.format(frac_C11_sum) + ' % of all de-excitation channels of residual $^{11}C^*$ nuclei',
          fontsize=15)
plt.grid(axis='x')

# B11
fig_b11, ax_b11 = plt.subplots(figsize=(11, 6))
horizontal_bars_b11 = ax_b11.barh(pos_B11, frac_B11, align='center', alpha=0.9)
label_barh(ax_b11, horizontal_bars_b11, "{:.4}", is_inside=False, fontsize=12)
plt.xticks(fontsize=11)
plt.yticks(pos_B11, ch_B11, fontsize=12)
plt.xlabel('fraction in %', fontsize=12)
plt.title('{0:.1f}'.format(frac_B11_sum) + ' % of all de-excitation channels of residual $^{11}B^*$ nuclei',
          fontsize=15)
plt.grid(axis='x')

# B10
fig_b10, ax_b10 = plt.subplots(figsize=(11, 6))
horizontal_bars_b10 = ax_b10.barh(pos_B10, frac_B10, align='center', alpha=0.9)
label_barh(ax_b10, horizontal_bars_b10, "{:.4}", is_inside=False, fontsize=12)
plt.xticks(fontsize=11)
plt.yticks(pos_B10, ch_B10, fontsize=12)
plt.xlabel('fraction in %', fontsize=12)
plt.title('{0:.1f}'.format(frac_B10_sum) + ' % of all de-excitation channels of residual $^{10}B^*$ nuclei',
          fontsize=15)
plt.grid(axis='x')


plt.show()
