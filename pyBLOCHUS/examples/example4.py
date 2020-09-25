"""Create example4 plot.

This script demonstrates how to calculate lookup-tables for pre-polarization
switch-off ramps that can subsequently be used with MRSMatlab to forward
model MRS sounding curves.
Basically this script allows to create the data that is shown in Fig. 4
of the paper:
    Hiller, T., Dlugosch, R. and Müller-Petke, M., "Utilizing pre-
    polarization to enhance SNMR signals - effect of imperfect
    switch-off", Geophysical Journal International Vol. 222(2), p.815-826, 2020

    general settings:
        - magnetic field B0 is set to 48 µT
        - T1, T2 - relaxation is ignored
    lookup table parameters:
        - pre-polarization field Bp ranges between 4.8µT and 480mT
        - initial orientation between Bp and B0 varies between 0° and 180°
          (angle between B0 and Bp)

See also: BLOCHUS
Author: Thomas Hiller
email: thomas.hiller[at]leibniz-liag.de
License: GNU GPLv3 (at end)
"""
import time
import numpy as np
import matplotlib.pyplot as plt

from pyBLOCHUS import BlochusPrepol as prepol

SAVEDATA = False
SAVEPNG = False


def test(save_png=False, save_data=False):
    """Perform the lookup table calculation."""
    # the different ramps
    my_ramp = ['exp', 'linexp', 'halfcos', 'lin']
    my_ramp_time = np.array([1e-4, 1e-3, 4e-3])

    # IMPORTANT NOTE:
    # This very coarse (and basically useless) fac_num=10 x theta_num=10 grid
    # takes for the 4x3 ramp combinations (shapes x times) roughly 75min! The
    # resolution in Fig.4 of the paper is fac_num=51 x theta_num=360 ... so
    # you don't want to start this calculation on your normal machine.
    # NOTE 2:
    # The Matlab version of this script only takes 30min

    # local lookup table settings (careful with the discretization)
    fac_num = 10
    theta_num = 10
    # pre-polarization factor [B0]
    factor = np.logspace(-1, 4, fac_num)
    # angle between Bp and B0 [deg]
    theta = np.linspace(0, 180, theta_num)
    # initialize prepol Bloch class
    bloch = prepol(B0=4.8e-5)

    start_time = time.time()
    # create output dictionary
    data = [dict() for x in range(len(my_ramp)*len(my_ramp_time))]
    count = 0
    for i_1 in np.arange(len(my_ramp)):
        for i_2 in np.arange(len(my_ramp_time)):
            ramp = my_ramp[i_1]
            t_ramp = my_ramp_time[i_2]

            # set ramp parameter
            bloch.Ramp.shape = ramp
            bloch.Ramp.t_ramp = t_ramp
            # set basic t_sim
            bloch.t_sim = t_ramp

            # output variables
            adiab_quality = np.zeros((len(theta), len(factor)))
            mag_final = np.zeros((len(theta), len(factor), 3))
            for t_1 in np.arange(len(theta)):
                for f_1 in np.arange(len(factor)):

                    bloch.Ramp.factor = factor[f_1]

                    # different ramps need different settings
                    if ramp == 'exp':
                        bloch.Ramp.switch_factor = 1
                        bloch.Ramp.t_slope = t_ramp/10
                    elif ramp == 'linexp':
                        bloch.Ramp.switch_factor = factor[f_1]/10
                        bloch.Ramp.t_slope = t_ramp/2
                    else:
                        bloch.Ramp.switch_factor = 1
                        bloch.Ramp.t_slope = t_ramp

                    # orientation
                    bloch.Ramp.theta = theta[t_1]

                    # 'Ramp.orient' gets updated automatically
                    # pre-polarization field into the direction of
                    # 'Ramp.orient' and Ramp.factor x B0 stronger than B0
                    b_p = bloch.Ramp.orient*bloch.Ramp.factor*bloch.B0
                    # Earth's magnetic field (as vector)
                    b_e = bloch.B0*bloch.zunit
                    # initial magnetization points into the sum of both
                    minit = b_p + b_e
                    # for simplicity equilibrium magnetization is B0
                    # (obtaining the real magnetization values is simply a
                    # multiplication of the magnetic fields with the "Curie"
                    # factor)
                    bloch.M0 = bloch.B0*bloch.zunit
                    # solve
                    bloch.solve(m_init=minit, use_numba=True)

                    # command line output
                    string = 'running: Ramp = {:7s} '.format(ramp)
                    string += '({} / {}) | '.format(i_1+1, len(my_ramp))
                    string += 'Tramp = {:3.1f}ms '.format(t_ramp*1e3)
                    string += '({} / {}) | '.format(i_2+1, len(my_ramp_time))
                    string += 'theta = {:5.1f}° '.format(theta[t_1])
                    string += '({:2d} / {:2d}) | '.format(t_1+1, len(theta))
                    string += 'factor = {:5.1e} '.format(factor[f_1])
                    string += '({:2d} / {:2d})'.format(f_1+1, len(factor))
                    print(string)

                    adiab_quality[t_1, f_1] = bloch.Ramp.adiab_qual
                    mag_final[t_1, f_1, :] = bloch.m[:, -1]
            # save the data for the current ramp - time - set
            dat = {'t_ramp': t_ramp, 'ramp': ramp, 'factor': factor,
                   'theta': theta, 'adiab_quality': adiab_quality,
                   'mag_final': mag_final}
            data[count] = dat
            count = count + 1

    end_time = time.time()
    print(end_time-start_time)
    # create the output figure
    fig = plt.figure(figsize=(11, 11))
    count = 0
    for i_1 in np.arange(len(my_ramp)):
        for i_2 in np.arange(len(my_ramp_time)):
            count = count + 1
            # grab the data from the dict
            adiab_quality = data[count-1]['adiab_quality']
            # create x and y vector
            # NOTE: due to the peculiarity of 'pcolor' the center of the
            # patches does not conform to the actual values used in the
            # simulation
            xvec = np.linspace(theta.min(), theta.max(), theta_num+1)
            yvec = np.logspace(np.log10(factor.min()),
                               np.log10(factor.max()), fac_num+1)
            # 2d vectors for pcolor
            x2d, y2d = np.meshgrid(xvec, yvec)
            ax1 = fig.add_subplot(4, 3, count)
            # plot data
            ax1.pcolor(x2d, y2d, adiab_quality.T, clim=(-1, 1))
            # adjust axes
            ax1.set_yscale('log')
            low, bot, wid, hei = ax1.get_position().bounds
            ax1.set_position([low, bot, wid*0.9, hei*0.8])
            ax1.set_xlim(theta.min(), theta.max())
            ax1.set_xticks([0, 45, 90, 135, 180])
            ax1.set_ylim(factor.min(), factor.max())
            ax1.set_xlabel(r'angle $\theta$ ($\angle B_0 B_p$) [deg]')
            ax1.set_ylabel(r'$B_p$ [$B_0$]')
            ax1.set_title('{} | {:.1f} ms'.format(my_ramp[i_1],
                                                  my_ramp_time[i_2]*1e3))

    # save data
    if save_data:
        np.save('example4_lookup-data', data, allow_pickle=True)
    # save figure
    if save_png:
        fig.savefig('example4.png', dpi=300)


# actually call the main function
# ask the user to proceed
QSTR = 'This script runs approx. 75min on a (Core i5-6200) Laptop. '
QSTR += 'Continue? [y/n]: '
answer = input(QSTR)
if answer in {"y", "ye", "yes"}:
    test(SAVEPNG, SAVEDATA)
else:
    print("Canceled. Bye.")

# License:
# GNU GPLv3
#
# blochUS
# Copyright (C) 2019 Thomas Hiller
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
