"""Create example1 plots.

This script demonstrates how basic pulses can be simulated with BLOCHUS
everything is done with a Hydrogen proton (gyromagnetic ratio > 0)

general settings:
-> magnetic field B0 is set to 50 µT
  -> fL(H)  = -2128.9 Hz
-> pulse length is set to 5 ms
-> different pulses are used (pi/2, pi)
    -> effect of different pulse axis
    -> effect of frequency offsets
    -> T1 and T2 relaxation times are 1 ms and 0.5 ms, respectively

See also: blochUS
Author: Thomas Hiller
email: thomas.hiller[at]leibniz-liag.de
License: GNU GPLv3 (at end)
"""
import numpy as np
import matplotlib.pyplot as plt
# although pylint is complaining, that these classes are not used (which
# technically is true), they have to be imported to make the 3d plots
from mpl_toolkits.mplot3d import axes3d, Axes3D  # pylint: disable=W0611

from pyBLOCHUS import BlochusPulse as pulse

# activate different tests
# 1 - pi/2 - pulse
#     this will create a figure corresponding to "example2a_ref"
# 2 - pi - pulse
#     this will create a figure corresponding to "example2b_ref"
# 3 - pi/2 - pulses with different pulse axes
#     this will create a figure corresponding to "example2c_ref"
# 4 - pi/2 pulses with different frequency offsets
#     this will create a figure corresponding to "example2d_ref"
usetest = [True, True, True, True]
SAVEPNG = False


def test_one(save_png=False):
    """Perform the first test demonstrating a pi/2 pulse."""
    # create output figure
    fig = plt.figure(figsize=(7, 7))
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222, projection='3d')
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224, projection='3d')
    # initialize pulse Bloch class
    bloch = pulse(T1=0.001, T2=0.0005, t_sim=0.005)
    # solve
    bloch.solve(m_init=bloch.zunit)
    # create lab-frame plots
    bloch.plot_magnetization(ax1, version='lab')
    bloch.plot_magnetization3d(ax2, version='lab')
    # create rot-frame plots
    bloch.plot_magnetization(ax3, version='rot')
    bloch.plot_magnetization3d(ax4, version='rot')
    # adjust axes position
    low, bot, wid, hei = ax1.get_position().bounds
    ax1.set_position([low, bot*1.1, wid, hei-0.05])
    low, bot, wid, hei = ax3.get_position().bounds
    ax3.set_position([low, bot*1.1, wid, hei-0.05])
    # adjust titles
    ax1.set_title(r'$^1$H | B0 = {:.1f}$\mu$T'.format(bloch.B0*1e6))
    ax2.set_title(r'$\gamma>0$')
    ax3.set_title(r'$^1$H | B1 ~ {:.2f}$\mu$T'.format(bloch.Pulse.factor *
                                                      bloch.B0*1e6))
    ax4.set_title(r'$(\pi/2)_{+x}$ - pulse')
    # adjust legends
    ax1.legend(prop={"size": 8}, loc='lower left')
    ax3.legend(prop={"size": 8}, loc='lower right')
    # save figure
    if save_png:
        fig.savefig('example2a.png', dpi=300)


def test_two(save_png=False):
    """Perform the second test demonstrating a pi pulse."""
    # create output figure
    fig = plt.figure(figsize=(7, 7))
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222, projection='3d')
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224, projection='3d')
    # initialize pulse Bloch class
    bloch = pulse(T1=0.001, T2=0.0005, t_sim=0.005)
    # set pi pulse
    bloch.Pulse.pulse_type = 'pi'
    # solve
    bloch.solve(m_init=bloch.zunit)
    # create lab-frame plots
    bloch.plot_magnetization(ax1, version='lab')
    bloch.plot_magnetization3d(ax2, version='lab')
    # create rot-frame plots
    bloch.plot_magnetization(ax3, version='rot')
    bloch.plot_magnetization3d(ax4, version='rot')
    # adjust axes position
    low, bot, wid, hei = ax1.get_position().bounds
    ax1.set_position([low, bot*1.1, wid, hei-0.05])
    low, bot, wid, hei = ax3.get_position().bounds
    ax3.set_position([low, bot*1.1, wid, hei-0.05])
    # adjust titles
    ax1.set_title(r'$^1$H | B0 = {:.1f}$\mu$T'.format(bloch.B0*1e6))
    ax2.set_title(r'$\gamma>0$')
    ax3.set_title(r'$^1$H | B1 ~ {:.2f}$\mu$T'.format(bloch.Pulse.factor *
                                                      bloch.B0*1e6))
    ax4.set_title(r'$(\pi)_{+x}$ - pulse')
    # adjust legends
    ax1.legend(prop={"size": 8}, loc='lower left')
    ax3.legend(prop={"size": 8}, loc='lower right')
    # save figure
    if save_png:
        fig.savefig('example2b.png', dpi=300)


def test_three(save_png=False):
    """Perform the third test demonstrating different pulse axes."""
    # create output figure
    fig = plt.figure(figsize=(7, 7))
    ax1 = fig.add_subplot(221, projection='3d')
    ax2 = fig.add_subplot(222, projection='3d')
    ax3 = fig.add_subplot(223, projection='3d')
    ax4 = fig.add_subplot(224, projection='3d')
    # initialize pulse Bloch class
    bloch = pulse(T1=0.001, T2=0.0005, t_sim=0.005)
    # solve #1
    bloch.solve(m_init=bloch.zunit)
    # create rot-frame plot
    bloch.plot_magnetization3d(ax1, version='rot')
    # solve #2
    bloch.Pulse.axis = '+y'
    bloch.solve(m_init=bloch.zunit)
    # create rot-frame plot
    bloch.plot_magnetization3d(ax2, version='rot')
    # solve #3
    bloch.Pulse.axis = '-x'
    bloch.solve(m_init=bloch.zunit)
    # create rot-frame plot
    bloch.plot_magnetization3d(ax3, version='rot')
    # solve #4
    bloch.Pulse.axis = '-y'
    bloch.solve(m_init=bloch.zunit)
    # create rot-frame plot
    bloch.plot_magnetization3d(ax4, version='rot')

    # adjust titles
    ax1.set_title(r'$(\pi/2)_{+x}$ - pulse')
    ax2.set_title(r'$(\pi/2)_{+y}$ - pulse')
    ax3.set_title(r'$(\pi/2)_{-x}$ - pulse')
    ax4.set_title(r'$(\pi/2)_{-y}$ - pulse')
    # save figure
    if save_png:
        fig.savefig('example2c.png', dpi=300)


def test_four(save_png=False):
    """Perform the fourth test demonstrating different off-resonances.

    here we basically re-plot Fig. 10.28 (p. 255) from
    Levitt, 2008, "spin dynamics" 2nd ed
    positive ratio means the pulse frequency is "higher" as absolute value
    (i.e. instead of -2000 Hz it is -2100 Hz)
    because of the negative Larmor frequency of H-protons the rotation axis
    dips downwards for negative ratios
    """
    # initialize pulse Bloch class
    bloch = pulse(T1=0.001, T2=0.0005, t_sim=0.005)
    # pulse axis
    bloch.Pulse.axis = '+y'
    # df/omega_nut ratio
    ratio = np.arange(-4., 5., 1.)
    # nutation frequency for a pi/2 pulse with 5 ms is 50 Hz
    omega_nut = (np.pi/2.0)/bloch.Pulse.t_pulse/2.0/np.pi  # [Hz]

    # create output figure
    fig = plt.figure(figsize=(7, 7))

    for i in np.arange(len(ratio)):
        # frequency offsets [Hz]
        freq_offset = ratio[i]*omega_nut
        # update frequency offset parameter
        bloch.Pulse.fmod.v_range[1] = freq_offset
        # solve
        bloch.solve(m_init=bloch.zunit)
        # create rot-frame plot
        ax1 = fig.add_subplot(331+i, projection='3d')
        bloch.plot_magnetization3d(ax1, version='rot')
        ax1.set_title(r'$\Omega_0 / \omega_{{{}}}$ = {}'
                      .format('nut', ratio[i]))

    # save figure
    if save_png:
        fig.savefig('example2d.png', dpi=300)


# actually call the main function
if usetest[0]:
    test_one(SAVEPNG)
if usetest[1]:
    test_two(SAVEPNG)
if usetest[2]:
    test_three(SAVEPNG)
if usetest[3]:
    test_four(SAVEPNG)

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
