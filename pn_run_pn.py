# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 11:28:17 2020

@author: Matteo
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

from pn_classes import *
from pn_functions import *
from pn_generate_pn import get_pn

def run_pn():
# =============================================================================
# # Import pore network
# =============================================================================
    pn = import_pickle('input/pn')
    soil = import_pickle('input/soil')
    test = import_pickle('input/test')
    macro = import_pickle('input/macro')
    fluid = import_pickle('input/fluid')

    [d, f] = pn.get_frequency_distribution()
    pn.plot_diameters_psd(d, f)
    saturation = pn.get_pn_saturation()

# =============================================================================
# # Import experiments
# =============================================================================
    exp = np.loadtxt('experimental data/kaolin.txt', skiprows=1)
    exp_s = exp[:, 0]
    exp_e = exp[:, 1]
    exp_ew = exp[:, 2]
    exp_sr = exp[:, 3]
    exp_w = exp[:, 4]


# =============================================================================
# # Initialize plot window
# =============================================================================
    # Main plot
    figure = MainFigure()
    figure.init_macro_figures(macro, exp)
    figure.init_fig1(saturation)
    figure.init_fig2(d, f)

# =============================================================================
# # Iteration
# =============================================================================
    iter_plot = 500
    delta_for_plot = 0
    delta_for_psd = 0
    for suction in test.suction_steps:
        print('suction '+str(suction)+'kPa')
        delta_for_plot += suction-test.suction_previous_step
        delta_for_psd += suction-test.suction_previous_step

# Desaturation law
        washburn_limit = 2 * fluid.washburn_constant / (suction)*10**3  # um

# Strain
        strain_sat = ((soil.lambda_water / macro.void_ratio[macro.idx]) *
                      (suction-test.suction_previous_step) / suction)
        strain_desat = (soil.kappa_air / macro.void_ratio[macro.idx]) * suction
        strain = np.array([strain_sat[0], strain_desat[0]])

# Update diameters
        dim = np.array([pn.x_dim, pn.y_dim, pn.z_dim])
        alf1 = cycle(dim, pn.saturation_ext, pn.diameters, washburn_limit,
                     strain, cord_number='26')
        [pn.diameters, pn.saturation_ext] = alf1

# Update macroscopic variables
        test.suction_previous_step = suction
        saturation = pn.get_pn_saturation()
        macro.update_macro(pn.diameters, pn.saturation, suction)
        ddist, freq = pn.get_frequency_distribution()

# Save psd
        # if (delta_for_psd >= 1000 or
        #     suction == test.suction_steps[-1] or
        #     suction == test.suction_steps[0]):
        #     delta_for_psd = 0
        #     save_psd(ddist, frequency, suction)

# Update figure
        if delta_for_plot >= iter_plot or suction == test.suction_steps[-1]:
            figure.upadate_main_figure(suction, macro, saturation, ddist, freq)
            delta_for_plot = 0
            plt.pause(0.1)

# Early exit
        if np.sum(saturation) == 0:  # Early exit
            print('Early exit: Pore network desaturated at suction:'
                  '{} kPa'.format(suction))
            break

    return pn, macro
