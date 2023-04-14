# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 10:20:08 2020

@author: Matteo
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), os.path.pardir, 'micropy'))
from micropy import pore_distribution as pdist


class Soil:
    def __init__(self, Gs, lamda_water, kappa_air):
        self.Gs = Gs
        self.lambda_water = lamda_water  # Cam Clay lambda slurry
        self.kappa_air = kappa_air  # Cam Clay kappa dry powder

class Fluid:
    def __init__(self, contact_angle, surface_tension):
        self.contact_angle = contact_angle  # degrees
        self.surface_tension = surface_tension  # N/m
        self.washburn_constant = (self.surface_tension *
                                  np.cos(np.radians(self.contact_angle)))

class Test:
    def __init__(self, s_max, s_min, s_delta):
        self.suction_max = s_max  # [kPa]
        # 94 if w47, 230 if w40
        self.suction_min = s_min  # kPa
        self.suction_delta = s_delta  # kPa
        self.suction_steps = np.arange(
            self.suction_min, self.suction_max, self.suction_delta)
        self.suction_previous_step = self.suction_steps[0]

class Macro():
        def __init__(self, Gs, diameters, cpd_e_max, suction_steps):
            self.A = 999 # um^2 --> does not matter

            e_init = cpd_e_max
            number_particles = (np.size(diameters)*np.mean(diameters)/e_init)

            self.volume_solids = self.A * number_particles
            self.mass_solids = self.volume_solids * Gs #  10^-3 pg

            array_dim = np.size(suction_steps)
            self.void_ratio = np.zeros(array_dim)
            self.porosity = np.zeros(array_dim)
            self.saturated_pores = np.zeros(array_dim)
            self.degree_saturation = np.zeros(array_dim)
            self.water_content = np.zeros(array_dim)
            self.suction_plot = suction_steps
            self.idx = 0

        def update_macro(self, diameters, saturation, suction):
# TODO: controlla che ogni volta che lancio la funzione saturation sia updated
            volume_voids = self.A*(np.sum(diameters))  # um^3
            volume_water = self.A * (np.sum(diameters * saturation))  # um^3
            mass_water = volume_water*0.997  # 10^-3 pg

            self.idx = np.where(self.suction_plot == suction)
            self.saturated_pores[self.idx] =(np.sum(saturation)/
                                             np.size(diameters)) * 100
            self.degree_saturation[self.idx] = volume_water / volume_voids * 100
            self.void_ratio[self.idx] = (volume_voids / self.volume_solids)
            self.porosity[self.idx] = volume_voids /(
                                        volume_voids+self.volume_solids)
            self.water_content[self.idx] = mass_water / self.mass_solids
            return self

class PN():
    def __init__(self, distribution, x_dim_pn, y_dim_pn, z_dim_pn):
        """ x--> downward
            y --> leftward
            z --> towards you """
        self.x_dim = x_dim_pn
        self.y_dim = y_dim_pn
        self.z_dim = z_dim_pn

        self.dist_input = {'intervals': distribution.inputs.intervals,
                           'dmax': distribution.inputs.dmax,
                           'dmin': distribution.inputs.dmin}
        self.psd_input = distribution.psd
        self.cpd_input = distribution.cpd

        self.diameters = np.zeros(shape=(x_dim_pn, y_dim_pn, z_dim_pn))
        self.saturation_ext = np.full(
                (x_dim_pn+2, y_dim_pn+2, z_dim_pn+2), 1)  # water=1, air=0

        self.generate_3D_pore_network()


    def get_pn_saturation(self):
        self.saturation = self.saturation_ext[
                1:self.x_dim+1, 1:self.y_dim+1, 1:self.z_dim+1]
        return self.saturation

    def generate_3D_pore_network(self):
        d_starting = self.cpd_input.d
        e_starting = self.cpd_input.e
        d_min = np.min(d_starting)
        d_max = np.max(d_starting)
        x=self.x_dim
        y=self.y_dim
        z=self.z_dim
        intervals = x*y*z
        [d_new, e_new] = pdist.DataElaboration.interpolate_e(
                d_min, d_max, d_starting, e_starting, intervals)
        [d_psd, e_psd] = pdist.DataElaboration.PSD.get_psd_from_cpd(d_new, e_new)
        f = e_psd/d_psd
        frequency = f/np.sum(f)
        foo = np.random.choice(d_psd, size=intervals, p=frequency)
        self.diameters = foo.reshape((x, y, z))

    def get_frequency_distribution(self):
        '''Generate frequency distribution of the pore network'''
        # Check pn diameters distribution
        unique_diameters, diameters_counts = np.unique(self.diameters,
                                                       return_counts=True)
        cum = np.cumsum(diameters_counts)

        dist = pdist.DataElaboration(self.dist_input)
        dist.get_cpd_from_array(unique_diameters, cum)
        [dist.psd.d, dist.psd.e] = dist.psd.get_psd_from_cpd(
            dist.cpd.d, dist.cpd.e)
# TODO: Swtich off the following two rows if you want to plot number frequency
        dist.psd.e *= dist.psd.d
        dist.psd.e /= np.sum(dist.psd.e)
        return dist.psd.d, dist.psd.e

    def plot_diameters_psd(self, dist_d, dist_f):
        frequency_psd_input = self.psd_input.e/np.sum(self.psd_input.e)
        plt.figure(1)
        plt.semilogx(self.psd_input.d, frequency_psd_input)
        plt.semilogx(dist_d, dist_f)
        plt.legend(('Input PSD', 'Pore Network PSD'))
        plt.xlabel('Diameter [um]')
        plt.ylabel('Frequency Distribution [-]')
        return self

class MainFigure():
    def __init__(self):
        self.fig, ([self.ax1, self.ax2, self.ax3],
                   [self.ax4, self.ax5, self.ax6]) = plt.subplots(2, 3)

        self.ax1.title.set_text('Pore saturation')
        self.ax2.set_xlabel('pore diameter [um]')
        self.ax2.set_ylabel('pore frequency [%]')
        self.ax3.set_xlabel('suction [kPa]')
        self.ax3.set_ylabel('saturated pores ratio [%]')
        self.ax4.set_xlabel('suction [kPa]')
        self.ax4.set_ylabel('water content [-]')
        self.ax5.set_xlabel('suction [kPa]')
        self.ax5.set_ylabel('void ratio [-]')
        self.ax6.set_xlabel('suction [kPa]')
        self.ax6.set_ylabel('degree of saturation [%]')

    def init_fig1(self, saturation):
        cmap_saturation = plt.cm.Blues
        norm_saturation = colors.Normalize(vmin=0., vmax=1.)
        self.im_saturation = self.ax1.imshow(
                            saturation[:, :, int(saturation.shape[2]/2)],
                            cmap=cmap_saturation,
                            norm=norm_saturation)
        cbar_saturation = self.fig.colorbar(self.im_saturation, ax=self.ax1)
        cbar_saturation.set_ticks([0, 1])
        cbar_saturation.set_ticklabels(['Dry', 'Saturated'])
        cbar_saturation.set_label = ('Pore saturation')
        return self

    def init_fig2(self, ddist, frequency):
        self.ax2.semilogx(ddist, frequency, label='0 kPa')
        self.legend_size= 6
        self.ax2.legend( prop={'size': self.legend_size})
        return self

    def update_fig1(self, saturation):
        self.im_saturation.set_data(
                                saturation[:, :, int(saturation.shape[2]/2)])
        return self

    def update_fig2(self, ddist, frequency, suction):
        self.ax2.semilogx(ddist, frequency, label='{} kPa'.format(suction))
        self.ax2.legend(prop={'size': self.legend_size})

    def init_macro_figures(self, macro, exp):
        exp_s = exp[:, 0]
        exp_e = exp[:, 1]
        exp_ew = exp[:, 2]
        exp_sr = exp[:, 3]
        exp_w = exp[:, 4]

        self.line42, = self.ax4.semilogx(exp_s, exp_w, 'go', label = 'exp')
        self.line52, = self.ax5.semilogx(exp_s, exp_e, 'go', label = 'exp')
        self.line62, = self.ax6.semilogx(exp_s, exp_sr*100, 'go', label = 'exp')

        self.line31, = self.ax3.semilogx(macro.suction_plot[0],
                             macro.saturated_pores[0], 'r', label = 'DPNM')
        self.line41, = self.ax4.semilogx(macro.suction_plot[0],
                             macro.water_content[0], 'r', label = 'DPNM')
        self.line51, = self.ax5.semilogx(macro.suction_plot[0],
                               macro.void_ratio[0], 'r', label = 'DPNM')
        self.line61, = self.ax6.semilogx(macro.suction_plot[0],
                              macro.degree_saturation[0], 'r',
                              label = 'DPNM')

        axes = [self.ax3,self.ax4, self.ax5, self.ax6]
        for ax in axes:
            ax.relim()
            ax.autoscale_view()
            ax.legend()

        return self

    def update_macro_figures(self, macro, idx):
        self.line31.set_xdata(macro.suction_plot[0:idx])
        self.line31.set_ydata(macro.saturated_pores[0:idx])
        self.line41.set_xdata(macro.suction_plot[0:idx])
        self.line41.set_ydata(macro.water_content[0:idx])
        self.line51.set_xdata(macro.suction_plot[0:idx])
        self.line51.set_ydata(macro.void_ratio[0:idx])
        self.line61.set_xdata(macro.suction_plot[0:idx])
        self.line61.set_ydata(macro.degree_saturation[0:idx])

        axes = [self.ax3, self.ax4, self.ax5, self.ax6]
        for ax in axes:
            ax.relim()
            ax.autoscale_view()
        return self

    def upadate_main_figure(self, suction, macro, saturation,
                            ddist, frequency):

        idx = macro.idx[0][-1]
        w = round(macro.water_content[idx], 2)
        sat_por = round(macro.saturated_pores[idx], 2)
        sr = round(macro.degree_saturation[idx])
        subtitle = [suction, w, sat_por, sr]
        subtitle_string = 's {a[0]} kPa; w {a[1]}\nsaturated pores {a[2]}%; Sr {a[3]}%'
        self.fig.suptitle(subtitle_string.format(a=subtitle))

        self.update_fig1(saturation)
        self.update_fig2(ddist, frequency, suction)
        self.update_macro_figures(macro, idx)
