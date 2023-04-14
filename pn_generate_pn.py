# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 17:33:37 2020

@author: Matteo
"""

import numpy as np
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__),
                             os.path.pardir, 'micropy'))
from micropy import pore_distribution as pdist

from pn_classes import *
from pn_functions import save_pickle

def get_pn(x_dim=10, y_dim=10, z_dim=10):
# =============================================================================
# # Import experimental data
# =============================================================================
# Import PSD input
    distribution_input_parameters = {'intervals': 50,
                                     'dmax': 161, 'dmin': 0.0022}
    distribution = pdist.DataElaboration(distribution_input_parameters)
    distribution.get_cpd_from_file('experimental data/cpd_input_w47.txt', True)
    alf = distribution.psd.get_psd_from_cpd(distribution.cpd.d,
                                            distribution.cpd.e)
    distribution.psd.d, distribution.psd.e = alf

# =============================================================================
# # Initialize PORE NETWORK
# =============================================================================
# Initialize soil, fluid and test
    soil = Soil(Gs =2.605, lamda_water=-0.168173896, kappa_air = -6.92541E-05)
    fluid = Fluid(contact_angle=0, surface_tension= 0.072)
    test = Test(s_max=10000, s_min=94, s_delta=10)

# Initialize pn
    pn = PN(distribution, x_dim_pn=x_dim, y_dim_pn=y_dim, z_dim_pn=z_dim)
    [d, f] = pn.get_frequency_distribution()
    pn.plot_diameters_psd(d, f)

# Boundary conditions: HYDRAULIC
    pn.saturation_ext[0, :, :] = 0  # water=1, air=0

# Initialize macro values
    saturation = pn.get_pn_saturation()
    macro = Macro(soil.Gs, pn.diameters, np.max(pn.cpd_input.e),
                  test.suction_steps)
    macro.update_macro(pn.diameters, pn.saturation, test.suction_steps[0])

# =============================================================================
# # Save pore network
# =============================================================================
    save_pickle(pn, 'input/pn')
    save_pickle(soil, 'input/soil')
    save_pickle(test, 'input/test')
    save_pickle(macro, 'input/macro')
    save_pickle(fluid, 'input/fluid')
