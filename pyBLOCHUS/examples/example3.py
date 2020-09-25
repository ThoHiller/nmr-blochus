"""Create example3 plots.

This script demonstrates how p-polarization switch-off ramps are simulated
with BLOCHUS.
Both examples were also used as benchmark cases in the paper:
Hiller, T., Dlugosch, R. and Müller-Petke, M., "Utilizing pre-
polarization to enhance SNMR signals - effect of imperfect
switch-off", Geophysical Journal International Vol. 222(2), p.815-826, 2020

Two p-polarization switch-off benchmarks are shown:

1.) from Melton et al., 1995, J Mag Res A, Vol. 117, p.164-170
    general settings:
     - magnetic field B0 is set to 50 µT
     - pre-polarization field is a factor 100 larger
     - uses a linear switch-off ramp with different switch-off rate (slopes)

2.) Conradi et al., 2017, J Mag Res, Vol. 281, p.241-245
    general settings:
        - magnetic field B0 is set to 50 µT
        - pre-polarization field is a factor 50 larger
        - uses a linexp switch-off ramp with different initial angles theta
          (angle between B0 and Bp) and different switch-over fields Bstar

See also: BLOCHUS
Author: Thomas Hiller
email: thomas.hiller[at]leibniz-liag.de
License: GNU GPLv3 (at end)
"""
import numpy as np
# we need several matplotlib libs for plotting
import matplotlib as mp
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
# although pylint is complaining, that these classes are not used (which
# technically is true), they have to be imported to make the 3d plots
from mpl_toolkits.mplot3d import axes3d, Axes3D  # pylint: disable=W0611

from pyBLOCHUS import BlochusPrepol as prepol

# activate different tests
# 1 - Melton et al. 1995 benchmark - NOTE: takes ~1min
# 2 - Conradi et al. 2017 benchmark - NOTE: takes ~7min
usetest = [True, True]
SAVEPNG = False


def test_one(save_png=False):
    """Perform the first test: 'Melton et al. 1995' benchmark."""
    # initialize prepol Bloch class
    bloch = prepol()
    # the ramp shape is linear -> 'lin'
    bloch.Ramp.shape = 'lin'
    # initial orientation is towards y-axis (phi=90
    bloch.Ramp.phi = 90
    # and in the x-y-plane (theta=90)
    bloch.Ramp.theta = 90

    # 'Ramp.orient' gets updated automatically
    # pre-polarization field into the direction of 'Ramp.orient'
    # and Ramp.factor x B0 stronger than B0
    b_p = bloch.Ramp.orient*bloch.Ramp.factor*bloch.B0
    # Earth's magnetic field (as vector)
    b_e = bloch.B0*bloch.zunit
    # initial magnetization points into the sum of both
    m_init = b_p + b_e
    # for simplicity equilibrium magnetization is B0 (obtaining the real
    # magnetization values is simply a multiplication of the magnetic fields
    # with the "Curie" factor)
    bloch.M0 = bloch.B0*bloch.zunit

    # predefined relaxation rates (k/gamma) from Melton et al. 1995
    relax_rates = np.array([1/16, 1/8, 1/4, 1/2, 1.0, 2.0, 4.0, 8.0, 16.0])
    # initialize output variables
    m_final = np.zeros((3, len(relax_rates)), dtype=float)
    phi_final = np.zeros((len(relax_rates),), dtype=float)
    theta_final = np.zeros((len(relax_rates),), dtype=float)
    # loop over relaxation rates
    print('start calculation of relaxation rates ...')
    for i in np.arange(0, len(relax_rates)):
        # current relaxation rate
        rate = relax_rates[i]
        # calculate switch-off ramp time
        gamma = bloch.Ramp.factor/rate
        bloch.Ramp.t_ramp = gamma/(bloch.B0*bloch.gamma)
        # update t_slope
        bloch.Ramp.t_slope = bloch.Ramp.t_ramp
        # update simulation time t_sim
        bloch.t_sim = bloch.Ramp.t_ramp
        # solve
        bloch.solve(m_init, use_numba=True)

        # calculate final angles describing the orientation of M
        # see p. 167 of the paper
        m_end = bloch.m[:, -1]
        m_end_n = np.linalg.norm(m_end)
        theta = 90 - np.rad2deg(np.arccos(m_end[2]/m_end_n))
        phi = 90 - np.rad2deg(np.arctan2(m_end[1], m_end[0]))

        # gather output data
        m_final[:, i] = m_end/m_end_n
        phi_final[i] = phi
        theta_final[i] = theta
        # inform the user about the progress
        print('{}/{} done'.format(i+1, len(relax_rates)))

    # create output figure
    fig = plt.figure(figsize=(9, 5))
    ax1 = fig.add_subplot(121, projection='3d')
    ax2 = fig.add_subplot(122)

    # draw only a part of the Bloch sphere
    bloch.plot_bloch_sphere(ax1, (30, 30), 1, (0, 90), (0, 90))

    # create an index vector for discrete colors out of the jet colormap
    index = range(len(relax_rates))
    jet = plt.get_cmap('jet')
    c_norm = colors.Normalize(vmin=0, vmax=index[-1])
    scalar_map = cmx.ScalarMappable(norm=c_norm, cmap=jet)

    # loop over all relaxation rates and plot the result with the
    # corresponding color
    for i in np.arange(0, len(relax_rates)):
        color_val = scalar_map.to_rgba(index[i])
        ax1.quiver(0, 0, 0, m_final[0, i], m_final[1, i], m_final[2, i],
                   color=color_val)
        ax2.plot(relax_rates[i], theta_final[i], marker='o', fillstyle='none',
                 color=color_val)
        ax2.plot(relax_rates[i], phi_final[i], marker='s', fillstyle='none',
                 color=color_val)

    # plot the trace of all final positions on the Bloch sphere surface
    m_x = np.append(np.append(np.array([0.0]), m_final[0, :]), np.array([0.0]))
    m_y = np.append(np.append(np.array([0.0]), m_final[1, :]), np.array([1.0]))
    m_z = np.append(np.append(np.array([1.0]), m_final[2, :]), np.array([0.0]))
    ax1.plot3D(m_x, m_y, m_z, 'k--')
    # adjust axis position
    ax1.set_position([0.0, -0.11, 0.35*1.5, 0.77*1.4])

    # plot the fit from the Melton et al. paper (eq. 7)
    x_1 = relax_rates[3:len(relax_rates)]
    y_1 = 50.8/np.sqrt(x_1)
    ax2.plot(x_1, y_1, 'k--')
    ax2.annotate(r'$50.8^{\circ}\sqrt{k/\Gamma}$', xy=(2.8, 30),
                 xycoords='data', xytext=(0.3, 20), textcoords='data',
                 arrowprops=dict(arrowstyle="->", connectionstyle="arc3"))

    # add a legend
    lgd_labels = [r'$\theta_f$', r'$\phi_f$']
    ax2.legend(lgd_labels, loc='upper right')

    # set axis properties
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_yticks([10, 100])
    ax2.get_yaxis().set_minor_formatter(mp.ticker.NullFormatter())
    ax2.set_xlim(0.04, 30)
    ax2.set_ylim(10, 100)
    ax2.grid(color='lightgray', which='both', axis='both')
    ax2.set_xlabel(r'relaxation rate k/$\Gamma$')
    ax2.set_ylabel('final angles [deg]')

    # save figure
    if save_png:
        fig.savefig('example3a.png', dpi=300)


def test_two(save_png=False):
    """Perform the second test: 'Conradi et al. 2017' benchmark."""
    # choose which parameters should be used
    conradi_fig = '4b'
    if conradi_fig == '4a':
        conradi_factor = 320
        conradi_switch = np.array([40.0, 5.0, 2.4])
        conradi_t_slope = 0.006
        lgd_labels = ['2000', '250', '120']
    else:
        conradi_factor = 50
        conradi_switch = np.array([5.0, 1.2, 0.02])
        conradi_t_slope = 0.01
        lgd_labels = ['250', '6', '1']

    # initialize prepol Bloch class
    bloch = prepol()
    # the ramp shape is linear & exponential -> 'linexp'
    bloch.Ramp.shape = 'linexp'
    # the factor is 50
    bloch.Ramp.factor = conradi_factor
    # initial orientation is towards y-axis (phi=90)
    bloch.Ramp.phi = 90

    bloch.Ramp.t_slope = conradi_t_slope
    bloch.Ramp.t_ramp = bloch.Ramp.t_slope * 2
    bloch.t_sim = bloch.Ramp.t_ramp

    bloch.M0 = bloch.B0*bloch.zunit

    theta = np.linspace(0, 200, 101, endpoint=True)
    adiab_quality = np.zeros((len(theta), 3))

    for i in np.arange(0, len(conradi_switch)):
        bloch.Ramp.switch_factor = conradi_switch[i]
        for j in np.arange(0, len(theta)):
            # running theta value
            bloch.Ramp.theta = theta[j]

            # 'Ramp.orient' gets updated automatically
            # pre-polarization field into the direction of 'Ramp.orient'
            # and Ramp.factor x B0 stronger than B0
            b_p = bloch.Ramp.orient*bloch.Ramp.factor*bloch.B0
            # Earth's magnetic field (as vector)
            b_e = bloch.B0*bloch.zunit
            # initial magnetization points into the sum of both
            m_init = b_p + b_e
            # solve
            bloch.solve(m_init, use_numba=True)

            adiab_quality[j, i] = bloch.Ramp.adiab_qual
            # inform the user about the progress
            print('B* = {} [µT] | Theta = {} [deg]'.format(lgd_labels[i],
                                                           theta[j]))

    # create output figure
    fig = plt.figure(figsize=(7, 5))
    ax1 = fig.add_subplot(111)
    col = ['r', 'g', 'b']
    for i in np.arange(0, len(conradi_switch)):
        ax1.plot(theta, adiab_quality[:, i], color=col[i])

    ax1.set_xlim(-20, 210)
    ax1.set_xticks([0, 45, 90, 135, 180])
    ax1.set_ylim(-1.05, 1.05)
    ax1.set_yticks([-1, -0.5, 0, 0.5, 1])
    ax1.grid(color='lightgray', axis='both')
    ax1.set_xlabel(r'angle $\theta$ ($\angle B_0 B_p$) [deg]')
    ax1.set_ylabel(r'adiabatic quality $p$')
    ax1.legend(lgd_labels, title=r'$B^\ast$ [µT]', loc='lower left')

    # save figure
    if save_png:
        fig.savefig('example3b.png', dpi=300)


# actually call the main function
# ask the user to proceed
QSTR = 'This script runs approx. 8min on a (Core i5-6200) Laptop. '
QSTR += 'Continue? [y/n]: '
answer = input(QSTR)
if answer in {"y", "ye", "yes"}:
    if usetest[0]:
        test_one(SAVEPNG)
    if usetest[1]:
        test_two(SAVEPNG)
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
