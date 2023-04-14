# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 17:56:04 2020

@author: Matteo
"""

from pn_generate_pn import *
from pn_run_pn import *

dim = 50
get_pn(dim, dim, dim)
run_pn()

# se kappa_air metti -6 invece che -4 lo swelling funziona meglio. Controlla se e' lo stesso per Sr
# TODO: create module
