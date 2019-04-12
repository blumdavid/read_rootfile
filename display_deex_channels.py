""" script to display the different deexcitation channels from DSNB-NC generation in a bar chart:

    The output root file of the DSNB-NC generator (gen_NC_onlyC12_250000evts_seed1.root) is analyzed with
    checkout_NCgen.py and the number of events for different deexcitation channels are saved in txt file
    (deex_channels_NC_onlyC12_from_generator_250000evts.txt).

    This txt file is edited in "deex_channels_info.ods" (also the files
    deex_channels_NC_onlyC12_from_generator_250000evts_old.txt
    "deex_channels_NC_onlyC12_from_generator_11330evts_old.txt" and
    "deex_channels_QEL_NC_from_generator_11330evts_old.txt" are
    edited in this file).
    Here, two fractions are calculated for each deexcitation channel:
        - number of events / total number of events on C12
        - number pf events / total number of events, where e.g. C11 is produced


    INFO:   all events can have additional gammas. Gammas are not specified in deex_id, but saved in hepevt file
            (therefore gammas will be simulated in detector simulation)


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
    sheet 'NEW!! NC only C12, 250000 evts': """

# C11: fractions in % (number of events / total number of events, where C11 is excited) (fractions below contain 88.2 %
# of all C11 deex channels) (fraction of number of total C11 from 250000 evts = 22.72 %):
frac_C11_p_d_Be8 = 22.51
frac_C11_p_alpha_Li6 = 20.09
frac_C11_2p_Be9 = 15.87
frac_C11_p_B10 = 13.12
frac_C11_d_B9 = 9.50
frac_C11_n_p_B9 = 7.14
frac_C11_sum = frac_C11_p_d_Be8 + frac_C11_p_alpha_Li6 + frac_C11_2p_Be9 + frac_C11_p_B10 + frac_C11_d_B9 + \
               frac_C11_n_p_B9

ch_C11 = ('$n$ + $p$ + $^9B$', '$d$ + $^9B$', '$p$ + $^{10}B$', '$2p$ + $^9Be$', '$p$ + $\\alpha$ + $^6Li$',
          '$p$ + $d$ + $^8Be$')

pos_C11 = np.arange(len(ch_C11))

frac_C11 = np.array([frac_C11_n_p_B9, frac_C11_d_B9, frac_C11_p_B10, frac_C11_2p_Be9, frac_C11_p_alpha_Li6,
                     frac_C11_p_d_Be8])


# B11: fractions in % (number of events / total number of events, where B11 is excited) (fractions below contain 87.4 %
# of all B11 deex channels) (fraction of number of total B11 from 250000 evts = 24.43 %):
frac_B11_n_d_Be8 = 18.83
frac_B11_n_alpha_Li6 = 18.34
frac_B11_n_B10 = 17.34
frac_B11_n_p_Be9 = 13.99
frac_B11_d_Be9 = 12.14
frac_B11_p_Be10 = 6.72

frac_B11_sum = frac_B11_p_Be10 + frac_B11_d_Be9 + frac_B11_n_p_Be9 + frac_B11_n_B10 + frac_B11_n_alpha_Li6 + \
               frac_B11_n_d_Be8

ch_B11 = ('$p$ + $^{10}Be$', '$d$ + $^9Be$', '$n$ + $p$ + $^9Be$', '$n$ + $^{10}B$', '$n$ + $\\alpha$ + $^6Li$',
          '$n$ + $d$ + $^8Be$')

pos_B11 = np.arange(len(ch_B11))

frac_B11 = np.array([frac_B11_p_Be10, frac_B11_d_Be9, frac_B11_n_p_Be9, frac_B11_n_B10, frac_B11_n_alpha_Li6,
                     frac_B11_n_d_Be8])


# B10: fractions in % (number of events / total number of events, where B10 is excited) (fractions below contain 86.5 %
# of all B10 deex channels) (fraction of number of total B10 from 250000 evts = 16.38 %):
frac_B10_n_p_Be8 = 24.1621
frac_B10_p_Be9 = 22.16
frac_B10_n_B9 = 18.57
frac_B10_d_Be8 = 7.60
frac_B10_He3_Li7 = 7.03
frac_B10_p_alpha_He5 = 6.97

frac_B10_sum = frac_B10_He3_Li7 + frac_B10_p_alpha_He5 + frac_B10_d_Be8 + frac_B10_n_B9 + frac_B10_p_Be9 + \
               frac_B10_n_p_Be8

ch_B10 = ('$p$ + $\\alpha$ + $^5He$', '$^3He$ + $^7Li$', '$d$ + $^8Be$', '$n$ + $^9B$', '$p$ + $^9Be$',
          '$n$ + $p$ + $^8Be$')

pos_B10 = np.arange(len(ch_B10))

frac_B10 = np.array([frac_B10_p_alpha_He5, frac_B10_He3_Li7, frac_B10_d_Be8, frac_B10_n_B9, frac_B10_p_Be9,
                     frac_B10_n_p_Be8])



# C10: fraction of number of total C10 from 250000 evts = 4.46 %
# Be10: fraction of number of total C10 from 250000 evts = 4.47 %
# lighter isotopes: fraction of number of total isotope from 250000 evts < 2.7 %


""" display in bar chart """
# C11:
fig_c11, ax_c11 = plt.subplots()
horizontal_bars_c11 = ax_c11.barh(pos_C11, frac_C11, align='center', alpha=0.9)
label_barh(ax_c11, horizontal_bars_c11, "{:.4}", is_inside=False, fontsize=12)
plt.xticks(fontsize=11)
plt.yticks(pos_C11, ch_C11, fontsize=12)
plt.xlabel('fraction in %', fontsize=12)
plt.title('{0:.0f}'.format(frac_C11_sum) + ' % of all de-excitation channels of residual $^{11}C^*$ nuclei',
          fontsize=15)
plt.grid(axis='x')

# B11
fig_b11, ax_b11 = plt.subplots()
horizontal_bars_b11 = ax_b11.barh(pos_B11, frac_B11, align='center', alpha=0.9)
label_barh(ax_b11, horizontal_bars_b11, "{:.4}", is_inside=False, fontsize=12)
plt.xticks(fontsize=11)
plt.yticks(pos_B11, ch_B11, fontsize=12)
plt.xlabel('fraction in %', fontsize=12)
plt.title('{0:.0f}'.format(frac_B11_sum) + ' % of all de-excitation channels of residual $^{11}B^*$ nuclei',
          fontsize=15)
plt.grid(axis='x')

# B10
fig_b10, ax_b10 = plt.subplots()
horizontal_bars_b10 = ax_b10.barh(pos_B10, frac_B10, align='center', alpha=0.9)
label_barh(ax_b10, horizontal_bars_b10, "{:.4}", is_inside=False, fontsize=12)
plt.xticks(fontsize=11)
plt.yticks(pos_B10, ch_B10, fontsize=12)
plt.xlabel('fraction in %', fontsize=12)
plt.title('{0:.0f}'.format(frac_B10_sum) + ' % of all de-excitation channels of residual $^{10}B^*$ nuclei',
          fontsize=15)
plt.grid(axis='x')


plt.show()
