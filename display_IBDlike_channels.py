""" Script to display the final channels (after deexcitation), the NC interaction channels (before deexcitation) and
    the deexcitation channels of IBD-like events:

    Moreover, also the ratios corresponding to the number of initial neutron and to the number of neutrons after
    deexcitation can be displayed.

    All NC events that mimic an IBD-like signature are stored in folder
    /home/astro/blum/juno/atmoNC/data_NC/output_detsim_v2/.

    The information about the different channels is stored in folder
    /home/astro/blum/juno/atmoNC/data_NC/output_detsim_v2/DCR_results_16000mm_10MeVto100MeV_1000nsto1ms_mult1_
    1800keVto2550keV_dist500mm_R17700mm_PSD99/
    and was analyzed with script IBDlike_events_channels_v2.py

    The txt files are summarized in file channels_IBDlike_signals.ods.

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


""" DCR_results_16000mm_10MeVto100MeV_1000nsto1ms_mult1_1800keVto2550keV_dist500mm_R17700mm: """

# insert fractions of final channels of IBD-like events from channels_IBDlike_signals_v2.ods, sheet
# 'combined channels' (DCR_results_16000mm_10MeVto100MeV_1000nsto1ms_mult1_1800keVto2550keV_dist500mm_R17700mm).
# (These fractions contain 88.67 % of all channels mimicking an IBD-like event):
# Fractions in %:

Frac_nu_C11_n_final_1 = 23.71
Frac_nu_B10_n_p_final_1 = 16.92
Frac_nu_Be9_n_2p_final_1 = 12.50
Frac_nu_Be8_n_p_d_final_1 = 9.21
Frac_nu_Li6_n_p_alpha_final_1 = 8.61
Frac_nu_Be8_2n_2p_final_1 = 4.46
Frac_nu_B9_2n_p_final_1 = 4.29
Frac_nu_Li8_n_3p_final_1 = 3.57
Frac_nu_C10_2n_final_1 = 1.48
Frac_nu_B9_n_d_final_1 = 1.26
Frac_nu_He7_n_4p_final_1 = 1.36
Frac_nu_Li7_n_2p_d_final_1 = 1.31

channels_final_1 = ('$^7Li$ + $n$ + $2p$ + $d$', '$^7He$ + $n$ + $4p$', '$^9B$ + $n$ + $d$', '$^{10}C$ + $2n$',
                    '$^8Li$ + $n$ + $3p$', '$^9B$ + $2n$ + $p$', '$^8Be$ + $2n$ + $2p$', '$^6Li$ + $n$ + $p$ + '
                                                                                         '$\\alpha$',
                    '$^8Be$ + $n$ + $p$ + $d$', '$^9Be$ + $n$ + $2p$', '$^{10}B$ + $n$ + $p$', '$^{11}C$ + $n$')

pos_final_1 = np.arange(len(channels_final_1))

frac_final_1 = np.array([Frac_nu_Li7_n_2p_d_final_1, Frac_nu_He7_n_4p_final_1, Frac_nu_B9_n_d_final_1,
                         Frac_nu_C10_2n_final_1, Frac_nu_Li8_n_3p_final_1, Frac_nu_B9_2n_p_final_1,
                         Frac_nu_Be8_2n_2p_final_1, Frac_nu_Li6_n_p_alpha_final_1, Frac_nu_Be8_n_p_d_final_1,
                         Frac_nu_Be9_n_2p_final_1, Frac_nu_B10_n_p_final_1, Frac_nu_C11_n_final_1])

""" display in bar chart """
fig_1, ax_1 = plt.subplots(figsize=(11, 6))
horizontal_bars_1 = ax_1.barh(pos_final_1, frac_final_1, align='center', alpha=0.9)
label_barh(ax_1, horizontal_bars_1, "{:.4}", is_inside=False, fontsize=12)
plt.xlim(xmin=0.0, xmax=25.0)
plt.xticks(fontsize=11)
plt.yticks(pos_final_1, channels_final_1, fontsize=12)
plt.xlabel('fraction in %', fontsize=12)
plt.title('Atmospheric NC interactions on $^{12}C$ that mimic an IBD signal\n'
          '(interaction channels: $\\nu_{x}$ + $^{12}C$ $\\rightarrow$ $\\nu_{x}$ + ...)', fontsize=15)
plt.grid(axis='x')

# insert fractions of IBD-like events corresponding to the number of neutrons after deexcitation from
# channels_IBDlike_signals_v2.ods, sheet 'combined channels'.
# DCR_results_16000mm_10MeVto100MeV_1000nsto1ms_mult1_1800keVto2550keV_dist500mm_R17700mm/
# (These fractions contain 100 % of all channels mimicking an IBD-like event)
# Fractions in %:

Frac_0_neutrons_1 = 1.10
Frac_1_neutron_1 = 84.61
Frac_2_neutrons_1 = 13.09
Frac_3_neutrons_1 = 0.91
Frac_4_or_more_neutrons_1 = 0.29

channels_n_1 = ('more than 3 neutrons', '3 neutrons', '2 neutrons', '1 neutron', '0 neutrons')

pos_n_1 = np.arange(len(channels_n_1))

frac_n_1 = np.array([Frac_4_or_more_neutrons_1, Frac_3_neutrons_1, Frac_2_neutrons_1, Frac_1_neutron_1,
                     Frac_0_neutrons_1])

""" display in bar chart """
fig2_1, ax2_1 = plt.subplots(figsize=(11, 6))
horizontal_bars2_1 = ax2_1.barh(pos_n_1, frac_n_1, align='center', alpha=0.9)
label_barh(ax2_1, horizontal_bars2_1, "{:.4}", is_inside=False, fontsize=12)
plt.xlim(xmin=0.0, xmax=100.0)
plt.xticks(fontsize=11)
plt.yticks(pos_n_1, channels_n_1, fontsize=12)
plt.xlabel('fraction in %', fontsize=12)
plt.title('Number of neutrons in atmospheric NC events that mimic an IBD signal\n'
          '(interaction channels: $\\nu_{x}$ + $^{12}C$ $\\rightarrow$ $\\nu_{x}$ + ...)', fontsize=15)
plt.grid(axis='x')

# TODO-me: These value must be updated!!

""" PSD is considered:
    DCR_results_16000mm_10MeVto100MeV_500nsto1ms_mult1_2400PEto3400PE_dist500mm_R17700mm_PSD99: """

# insert fractions of final channels of IBD-like events from channels_IBDlike_signals_v2.ods, sheet
# 'combined channels with PSD'
# (DCR_results_16000mm_10MeVto100MeV_500nsto1ms_mult1_2400PEto3400PE_dist500mm_R17700mm_PSD99).
# (These fractions contain 88.61 % of all channels mimicking an IBD-like event):
# Fractions in %:

Frac_nu_C11_n_final_PSD = 11.39
Frac_nu_B10_n_p_final_PSD = 10.13
Frac_nu_Be9_n_2p_final_PSD = 11.39
Frac_nu_Be8_n_p_d_final_PSD = 6.33
Frac_nu_Li6_n_p_alpha_final_PSD = 3.80
Frac_nu_Be8_2n_2p_final_PSD = 8.86
Frac_nu_B9_2n_p_final_PSD = 1.27
Frac_nu_Li8_n_3p_final_PSD = 8.86
Frac_nu_C10_2n_final_PSD = 1.27
Frac_nu_B9_n_d_final_PSD = 2.53
Frac_nu_He7_n_4p_final_PSD = 0.00
Frac_nu_Li7_n_2p_d_final_PSD = 13.92
Frac_nu_Be7_2n_p_d_final_PSD = 8.86

channels_final_PSD = ('$^7Be$ + $2n$ + $p$ + $d$', '$^7Li$ + $n$ + $2p$ + $d$', '$^7He$ + $n$ + $4p$',
                      '$^9B$ + $n$ + $d$', '$^{10}C$ + $2n$',
                      '$^8Li$ + $n$ + $3p$', '$^9B$ + $2n$ + $p$', '$^8Be$ + $2n$ + $2p$', '$^6Li$ + $n$ + $p$ + '
                                                                                         '$\\alpha$',
                      '$^8Be$ + $n$ + $p$ + $d$', '$^9Be$ + $n$ + $2p$', '$^{10}B$ + $n$ + $p$', '$^{11}C$ + $n$')

pos_final_PSD = np.arange(len(channels_final_PSD))

frac_final_PSD = np.array([Frac_nu_Be7_2n_p_d_final_PSD, Frac_nu_Li7_n_2p_d_final_PSD, Frac_nu_He7_n_4p_final_PSD,
                           Frac_nu_B9_n_d_final_PSD,
                           Frac_nu_C10_2n_final_PSD, Frac_nu_Li8_n_3p_final_PSD, Frac_nu_B9_2n_p_final_PSD,
                           Frac_nu_Be8_2n_2p_final_PSD, Frac_nu_Li6_n_p_alpha_final_PSD, Frac_nu_Be8_n_p_d_final_PSD,
                           Frac_nu_Be9_n_2p_final_PSD, Frac_nu_B10_n_p_final_PSD, Frac_nu_C11_n_final_PSD])

""" display in bar chart """
fig_PSD, ax_PSD = plt.subplots(figsize=(11, 6))
horizontal_bars_PSD = ax_PSD.barh(pos_final_PSD, frac_final_PSD, align='center', alpha=0.9)
label_barh(ax_PSD, horizontal_bars_PSD, "{:.4}", is_inside=False, fontsize=12)
plt.xlim(xmin=0.0, xmax=100.0)
plt.xticks(fontsize=11)
plt.yticks(pos_final_PSD, channels_final_PSD, fontsize=12)
plt.xlabel('fraction in %', fontsize=12)
plt.title('Atmospheric NC interactions on $^{12}C$ that mimic an IBD signal and pass PSD cut\n'
          '(interaction channels: $\\nu_{x}$ + $^{12}C$ $\\rightarrow$ $\\nu_{x}$ + ...)', fontsize=15)
plt.grid(axis='x')

# insert fractions of IBD-like events corresponding to the number of neutrons after deexcitation from
# channels_IBDlike_signals_v2.ods, sheet 'combined channels with PSD'.
# DCR_results_16000mm_10MeVto100MeV_500nsto1ms_mult1_2400PEto3400PE_dist500mm_R17700mm_PSD99/
# (These fractions contain 100 % of all channels mimicking an IBD-like event)
# Fractions in %:

Frac_0_neutrons_PSD = 0.0
Frac_1_neutron_PSD = 70.89
Frac_2_neutrons_PSD = 22.78
Frac_3_neutrons_PSD = 6.33
Frac_4_or_more_neutrons_PSD = 0.0

channels_n_PSD = ('more than 3 neutrons', '3 neutrons', '2 neutrons', '1 neutron', '0 neutrons')

pos_n_PSD = np.arange(len(channels_n_PSD))

frac_n_PSD = np.array([Frac_4_or_more_neutrons_PSD, Frac_3_neutrons_PSD, Frac_2_neutrons_PSD, Frac_1_neutron_PSD,
                       Frac_0_neutrons_PSD])

""" display in bar chart """
fig2_PSD, ax2_PSD = plt.subplots(figsize=(11, 6))
horizontal_bars2_PSD = ax2_PSD.barh(pos_n_PSD, frac_n_PSD, align='center', alpha=0.9)
label_barh(ax2_PSD, horizontal_bars2_PSD, "{:.4}", is_inside=False, fontsize=12)
plt.xlim(xmin=0.0, xmax=100.0)
plt.xticks(fontsize=11)
plt.yticks(pos_n_PSD, channels_n_PSD, fontsize=12)
plt.xlabel('fraction in %', fontsize=12)
plt.title('Number of neutrons in atmospheric NC events that mimic an IBD signal\n'
          '(interaction channels: $\\nu_{x}$ + $^{12}C$ $\\rightarrow$ $\\nu_{x}$ + ...)', fontsize=15)
plt.grid(axis='x')

plt.show()
