import ROOT
import numpy as np
from matplotlib import pyplot as plt
import NC_background_functions

"""
efficiency_volume_cut = 100.019
efficiency_prompt_energy_cut = 98.873
efficiency_time_cut = 100.004
efficiency_delayed_energy_cut = 100.087
efficiency_neutron_multiplicity_cut = 100.0
efficiency_distance_cut = 100.0

error_efficiency_volume_cut = 0.368
error_efficiency_prompt_energy_cut = 2.006
error_efficiency_time_cut = 0.449
error_efficiency_delayed_energy_cut = 0.452
error_efficiency_neutron_multiplicity_cut = 0.454
error_efficiency_distance_cut = 0.454

efficiency_muon_veto = 100.00
error_efficiency_muon_veto = 0.00

# consider the cut efficiencies for volume, prompt energy, time, delayed energy, multiplicity, distance cut and muon
# veto cut in percent:
cut_efficiency = efficiency_volume_cut/100.0 * efficiency_prompt_energy_cut/100.0 * efficiency_time_cut/100.0 * \
                 efficiency_delayed_energy_cut/100.0 * efficiency_neutron_multiplicity_cut/100.0 \
                 * efficiency_distance_cut/100.0 * efficiency_muon_veto/100.0 * 100.0
# calculate statistical error of cut_efficiency with Gaussian error propagation in percent:
error_cut_efficiency = (efficiency_prompt_energy_cut/100.0 * efficiency_time_cut/100.0 *
                        efficiency_delayed_energy_cut/100.0 * efficiency_neutron_multiplicity_cut/100.0 *
                        efficiency_distance_cut/100.0 * efficiency_muon_veto/100.0 * error_efficiency_volume_cut +
                        efficiency_volume_cut/100.0 * efficiency_time_cut/100.0 * efficiency_delayed_energy_cut/100.0 *
                        efficiency_neutron_multiplicity_cut/100.0 * efficiency_distance_cut/100.0 *
                        efficiency_muon_veto/100.0 * error_efficiency_prompt_energy_cut +
                        efficiency_volume_cut/100.0 * efficiency_prompt_energy_cut/100.0 *
                        efficiency_delayed_energy_cut/100.0 * efficiency_neutron_multiplicity_cut/100.0 *
                        efficiency_distance_cut/100.0 * efficiency_muon_veto/100.0 * error_efficiency_time_cut +
                        efficiency_volume_cut/100.0 * efficiency_prompt_energy_cut/100.0 *
                        efficiency_time_cut/100.0 * efficiency_neutron_multiplicity_cut/100.0 *
                        efficiency_distance_cut/100.0 * efficiency_muon_veto/100.0 *
                        error_efficiency_delayed_energy_cut +
                        efficiency_volume_cut/100.0 * efficiency_prompt_energy_cut/100.0 *
                        efficiency_time_cut/100.0 * efficiency_delayed_energy_cut/100.0 * efficiency_distance_cut/100.0
                        * efficiency_muon_veto/100.0 * error_efficiency_neutron_multiplicity_cut +
                        efficiency_volume_cut/100.0 * efficiency_prompt_energy_cut/100.0 *
                        efficiency_time_cut/100.0 * efficiency_delayed_energy_cut/100.0 *
                        efficiency_neutron_multiplicity_cut/100.0 * efficiency_muon_veto/100.0 *
                        error_efficiency_distance_cut +
                        efficiency_volume_cut/100.0 * efficiency_prompt_energy_cut/100.0 *
                        efficiency_time_cut/100.0 * efficiency_delayed_energy_cut/100.0 *
                        efficiency_neutron_multiplicity_cut/100.0 * efficiency_distance_cut/100.0 *
                        error_efficiency_muon_veto)

print("cut_efficiency = {0:.3f} %".format(cut_efficiency))
print("stat. error of cut efficiency = {0:.3f} %".format(error_cut_efficiency))
"""








