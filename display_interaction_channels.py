""" script to display the different interaction channels from Julia's GENIE simulation in a bar chart:

    The original GENIE root file of Julia is read with "read_GENIE_file.py" and the fractions of the interaction
    channels for the different isotopes are saved in file "interaction_channels_qel_NC_201617evts.txt".

    This txt file is edited in "interaction_channels_edit.ods" (also the file "interaction_channels_NC_318097evts.txt"
    is edited in this libre office file). The calculated fraction are defined as
    number_of_events/number_of_total_events, where the total events are all NC and QEL events (channels with targets
    C12, proton, electron, N14, O16, S32). In "interaction_channels_edit.ods", these fractions are converted to
    fractions, which are defined as number_of_events/number_of_events_with_target_C12 (only channels with C12 as target
    are taken into account).
    These are the only channel that can mimic an IBD event (the other channels are just elastic scattering on the target
    and therefore there is no neutron).

    In the second sheet of the libre office file ('sorted information'), the interaction channels are sorted with
    descending fraction.

    INFO:   all channels can also include Pion_0, Kaon_0, Kaon_plus, Kaon_minus (or heavier hadrons), but this is very
            unlikely, because events like this should not be quasi-elastic and therefore not stored.


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


""" insert the fractions of the different interaction channels from 'interaction_channels_edit.ods', 
    sheet 'Sorted_Information'. These fractions contain 89.7480 % of all interaction channels 
    (fractions of nu + C12 -> ... in %): """
Frac_nu_B11_p = 29.0523
Frac_nu_C11_n = 24.9883
Frac_nu_B10_p_n = 18.2382
Frac_nu_Be10_2p = 4.1941
Frac_nu_C10_2n = 4.0271
Frac_nu_Be9_2p_n = 1.1153
Frac_nu_B9_p_2n = 1.0579
Frac_nu_Be8_2p_2n = 0.9999
Frac_nu_Li6_3p_3n = 0.8657
Frac_nu_Li7_3p_2n = 0.8510
Frac_nu_Be7_2p_3n = 0.8258
Frac_nu_Li9_3p = 0.7461
Frac_nu_Li8_3p_n = 0.7191
Frac_nu_C9_3n = 0.6975
Frac_nu_B8_p_3n = 0.6892
Frac_nu_He5_4p_3n = 0.6805
Frac_other = 10.2520


channels = ('other rare channels', '$^5He$ + $4p$ + $3n$', '$^8B$ + $p$ + $3n$', '$^9C$ + $3n$', '$^8Li$ + $3p$ + $n$',
            '$^9Li$ + $3p$',
            '$^7Be$ + $2p$ + $3n$', '$^7Li$ + $3p$ + $2n$', '$^6Li$ + $3p$ + $3n$', '$^8Be$ + $2p$ + $2n$',
            '$^9B$ + $p$ + $2n$', '$^9Be$ + $2p$ + $n$', '$^{10}C$ + $2n$', '$^{10}Be$ + $2p$', '$^{10}B$ + $p$ + $n$',
            '$^{11}C$ + $n$', '$^{11}B$ + $p$')

pos = np.arange(len(channels))

fractions = np.array([Frac_other, Frac_nu_He5_4p_3n, Frac_nu_B8_p_3n, Frac_nu_C9_3n, Frac_nu_Li8_3p_n, Frac_nu_Li9_3p,
                      Frac_nu_Be7_2p_3n, Frac_nu_Li7_3p_2n, Frac_nu_Li6_3p_3n, Frac_nu_Be8_2p_2n, Frac_nu_B9_p_2n,
                      Frac_nu_Be9_2p_n, Frac_nu_C10_2n, Frac_nu_Be10_2p, Frac_nu_B10_p_n, Frac_nu_C11_n, Frac_nu_B11_p])


""" display in bar chart """
fig, ax = plt.subplots()
horizontal_bars = ax.barh(pos, fractions, align='center', alpha=0.9)
label_barh(ax, horizontal_bars, "{:.4}", is_inside=False, fontsize=12)
plt.xticks(fontsize=11)
plt.yticks(pos, channels, fontsize=12)
plt.xlabel('fraction in %', fontsize=12)
plt.title('Atmospheric NC QEL neutrino interactions on $^{12}C$\n'
          '(interaction channels: $\\nu_{x}$ + $^{12}C$ $\\rightarrow$ $\\nu_{x}$ + ...)', fontsize=15)
plt.grid(axis='x')


""" insert the fractions of the different interaction channels from 'interaction_channels_edit.ods', 
    sheet 'Sorted_Information 2'. These fractions contain 74.80581% of all interaction channels 
    (fractions of nu + C12 -> ... in %): """
F_nu_B11_p = 20.9711
F_nu_C11_n = 18.7591
F_nu_B10_p_n = 13.9626
F_nu_C11_p_piminus = 3.3507
F_nu_Be10_2p = 3.3229
F_nu_C10_2n = 3.2678
F_nu_B11_n_piplus = 2.6327
F_nu_B9_p_2n = 2.0312
F_nu_Be9_2p_n = 1.6938
F_nu_noiso = 1.6814
F_nu_Be8_2p_2n = 1.1800
F_nu_C10_p_n_piminus = 1.0105
F_nu_Be10_p_n_piplus = 0.9422

channels_2 = ('$^{10}Be$ + $p$ + $n$ + $\\pi^+$', '$^{10}C$ + $p$ + $n$ + $\\pi^-$', '$^8Be$ + $2p$ + $2n$',
              'no isotope (only $p$, $n$, $\\pi$, ...)', '$^9Be$ + $2p$ + $n$', '$^9B$ + $p$ + $2n$',
              '$^{11}B$ + $n$ + $\\pi^+$', '$^{10}C$ + $2n$', '$^{10}Be$ + $2p$',
              '$^{11}C$ + $p$ + $\\pi^-$', '$^{10}B$ + $p$ + $n$', '$^{11}C$ + $n$', '$^{11}B$ + $p$')

pos_2 = np.arange(len(channels_2))

fractions_2 = np.array([F_nu_Be10_p_n_piplus, F_nu_C10_p_n_piminus, F_nu_Be8_2p_2n, F_nu_noiso, F_nu_Be9_2p_n,
                        F_nu_B9_p_2n, F_nu_B11_n_piplus, F_nu_C10_2n, F_nu_Be10_2p, F_nu_C11_p_piminus, F_nu_B10_p_n,
                        F_nu_C11_n, F_nu_B11_p])

""" display in bar chart """
fig2, ax2 = plt.subplots()
horizontal_bars_2 = ax2.barh(pos_2, fractions_2, align='center', alpha=0.9)
label_barh(ax2, horizontal_bars_2, "{:.4}", is_inside=False, fontsize=12)
plt.xticks(fontsize=11)
plt.yticks(pos_2, channels_2, fontsize=12)
plt.xlabel('fraction in %', fontsize=12)
plt.title('75% of all atmospheric NC neutrino interactions on $^{12}C$\n'
          '(interaction channels: $\\nu_{x}$ + $^{12}C$ $\\rightarrow$ $\\nu_{x}$ + ...)', fontsize=15)
plt.grid(axis='x')


""" insert the fractions of the different interaction channels from 'interaction_channels_edit.ods', 
    sheet 'Interaction Channels after Generator'. These fraction are calculated with checkout_NCgen.py and saved in 
    file NC_onlyC12_interaction_channels_250000evts.txt.
    
    IMPORTANT: fractions of isotopes are only included, when there exists the corresponding TALYS simulation of the 
    isotope!    """
# total fractions of each isotope (in %), if total fraction > 2 % (sum of these fraction 77.65 %):
frac_nu_B11 = 24.55
frac_nu_C11 = 22.72
frac_nu_B10 = 16.41
frac_nu_Be10 = 4.50
frac_nu_C10 = 4.43
frac_nu_B9 = 2.68
frac_nu_Be9 = 2.36

channels_3 = ('$^{9}Be$ + ...', '$^{9}B$ + ...', '$^{10}C$ + ...', '$^{10}Be$ + ...', '$^{10}B$ + ...',
              '$^{11}C$ + ...', '$^{11}B$ + ...')

pos_3 = np.arange(len(channels_3))

fractions_3 = np.array([frac_nu_Be9, frac_nu_B9, frac_nu_C10, frac_nu_Be10, frac_nu_B10, frac_nu_C11, frac_nu_B11])

""" display in bar chart """
fig3, ax3 = plt.subplots()
horizontal_bars_3 = ax3.barh(pos_3, fractions_3, align='center', alpha=0.9)
label_barh(ax3, horizontal_bars_3, "{:.4}", is_inside=False, fontsize=12)
plt.xticks(fontsize=11)
plt.yticks(pos_3, channels_3, fontsize=12)
plt.xlabel('fraction in %', fontsize=12)
plt.title('78% of all atmospheric NC neutrino interactions on $^{12}C$\n'
          '(interaction channels: $\\nu_{x}$ + $^{12}C$ $\\rightarrow$ $\\nu_{x}$ + ...)', fontsize=15)
plt.grid(axis='x')

# leading fractions in % (sum of these fractions 70.02 %):
frac_nu_B11_p = 21.01
frac_nu_C11_n = 18.79
frac_nu_B10_p_n = 14.03
frac_nu_Be10_2p = 3.30
frac_nu_C11_p_piminus = 3.25
frac_nu_C10_2n = 3.24
frac_nu_B11_n_piplus = 2.67
frac_nu_B9_p_2n = 2.02
frac_nu_Be9_2p_n = 1.71

channels_4 = ('$^9Be$ + $2p$ + $n$', '$^9B$ + $p$ + $2n$', '$^{11}B$ + $n$ + $\\pi^+$', '$^{10}C$ + $2n$',
              '$^{11}C$ + $p$ + $\\pi^-$', '$^{10}Be$ + $2p$', '$^{10}B$ + $p$ + $n$', '$^{11}C$ + $n$',
              '$^{11}B$ + $p$')

pos_4 = np.arange(len(channels_4))

fractions_4 = np.array([frac_nu_Be9_2p_n, frac_nu_B9_p_2n, frac_nu_B11_n_piplus, frac_nu_C10_2n, frac_nu_C11_p_piminus,
                        frac_nu_Be10_2p, frac_nu_B10_p_n, frac_nu_C11_n, frac_nu_B11_p])

""" display in bar chart """
fig4, ax4 = plt.subplots()
horizontal_bars_4 = ax4.barh(pos_4, fractions_4, align='center', alpha=0.9)
label_barh(ax4, horizontal_bars_4, "{:.4}", is_inside=False, fontsize=12)
plt.xticks(fontsize=11)
plt.yticks(pos_4, channels_4, fontsize=12)
plt.xlabel('fraction in %', fontsize=12)
plt.title('70% of all atmospheric NC neutrino interactions on $^{12}C$\n'
          '(interaction channels: $\\nu_{x}$ + $^{12}C$ $\\rightarrow$ $\\nu_{x}$ + ...)', fontsize=15)
plt.grid(axis='x')


plt.show()
