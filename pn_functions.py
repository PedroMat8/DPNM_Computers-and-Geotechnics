# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 14:46:04 2020

@author: Matteo
"""
import numpy as np
from numba import jit
import pickle

# =============================================================================
# # PN loop
# =============================================================================
@jit(nopython=True)
def cycle(dim, saturation_ext, diameters, washburn_limit, strain, cord_number):

    diameters += strain[0]*np.multiply(diameters,saturation_ext[1:-1, 1:-1, 1:-1])
    control = True
    while control:
        # TODO: do it in a recursive manner: diameter-->mask-->create a list with the coordinate number
        control = False
        d = np.multiply(diameters,saturation_ext[1:-1, 1:-1, 1:-1])
        index = np.where(d > washburn_limit)
        for i in range(len(index[0])):
            x = index[0][i]
            y = index[1][i]
            z = index[2][i]
            counter = coordination_count(saturation_ext, x, y, z, cord_number)

            if counter < 26:
                saturation_ext[x+1,y+1,z+1] = 0
                # diameters[x, y, z] -= delta[1]*diameters[x, y, z]
                diameters[x, y, z] -= strain[1]*diameters[x, y, z]
                control = True
    # print('in cycle', diameters[10,10,10])
    return diameters, saturation_ext

@jit(nopython=True)
def coordination_count(saturation_ext, x, y, z, number):

    if number == '26':
        return (saturation_ext[x, y, z] +
                saturation_ext[x, y+1, z] + saturation_ext[x, y+2, z] +
                saturation_ext[x+1, y, z] + saturation_ext[x+1, y+2, z] +
                saturation_ext[x+2, y, z] + saturation_ext[x+2, y+1, z] +
                saturation_ext[x+2, y+2, z] +

                saturation_ext[x, y, z+1] +
                saturation_ext[x, y+1, z+1] + saturation_ext[x, y+2, z+1] +

                saturation_ext[x+1, y, z+1] + saturation_ext[x+1, y+2, z+1] +

                saturation_ext[x+2, y, z+1] +
                saturation_ext[x+2, y+1, z+1] + saturation_ext[x+2, y+2, z+1] +

                saturation_ext[x, y, z+2] +
                saturation_ext[x, y+1, z+2] + saturation_ext[x, y+2, z+2] +
                saturation_ext[x+1, y, z+2] + saturation_ext[x+1, y+2, z+2] +

                saturation_ext[x+2, y, z+2] +
                saturation_ext[x+2, y+1, z+2] + saturation_ext[x+2, y+2, z+2] +

                saturation_ext[x+1, y+1, z] + saturation_ext[x+1, y+1, z+2])

    elif number == '6':
            return (saturation_ext[x, y+1, z+1] +
                    saturation_ext[x+2, y+1, z+1] +
                    saturation_ext[x+1, y, z+1] +
                    saturation_ext[x+1, y+2, z+1] +
                    saturation_ext[x+1, y+1, z] +
                    saturation_ext[x+1, y+1, z+2])
    else:
        print('coordination number is neither 6 or 26')

# =============================================================================
# # Import/Export
# =============================================================================
def save_psd(ddist, frequency, suction):
        matrix = np.column_stack((ddist, frequency))
        header=('diameters_psd\tfrequency_psd')
        namefile = 'output/'+str(suction)+'.txt'
        np.savetxt(namefile, matrix,
                header=header, delimiter='\t', fmt='%s')

def save_pickle(whattosave, filename):
    with open(filename, 'wb') as file:
        pickle.dump(whattosave, file, pickle.HIGHEST_PROTOCOL)
    return

def import_pickle(whattoopen):
    with open(whattoopen, 'rb') as f:
        return pickle.load(f)