"""Create example1 plots.

This script demonstrates the basic use of the blochUS ode-solver. Two
different protons (H & He) are used to show the effect of the sign of the
gyromagnetic ratio gamma
1.) Hydrogen proton (gyromagnetic ratio > 0)
2.) Helium proton (gyromagnetic ratio < 0)

general settings:
-> magnetic field B0 is set to 50 µT
  -> fL(H)  = -2128.9 Hz
  -> fL(He) =  1621.7 Hz
-> simulation time is 5 ms
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

from pyBLOCHUS import BlochusBasic as basic

# activate different tests
# 1 - using different protons with different signs for gamma
#     this will create two figures corresponding to "example1a_ref" and
#       "example1b_ref"
# 2 - comparison of relaxation time effect
#     this will create one figure corresponding to "example1c_ref"
usetest = [True, True]
SAVEPNG = False


def test_one(save_png=False):
    """Perform the first test demonstrating the effect of the sign of gamma."""
    # create output figure
    fig = plt.figure(figsize=(7, 7))
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222, projection='3d')
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224, projection='3d')
    # initialize basic Bloch class
    bloch = basic(T1=0.001, T2=0.0005, t_sim=0.005)
    # initial magnetization
    m_init = np.array([1.0, 0.0, 0.0])
    # solve
    bloch.solve(m_init)
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
    ax3.set_title(r'$^1$H | fL = {:.1f}Hz'.format(bloch.larmor_f))
    # adjust legends
    ax1.legend(prop={"size": 8}, loc='lower right')
    ax3.legend(prop={"size": 8}, loc='lower right')
    # save figure
    if save_png:
        fig.savefig('example1a.png', dpi=300)

    # second example with He proton
    # create output figure
    fig = plt.figure(figsize=(7, 7))
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222, projection='3d')
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224, projection='3d')
    # initialize basic Bloch class
    bloch = basic(nucleus='3He', T1=0.001, T2=0.0005, t_sim=0.005)
    # initial magnetization
    m_init = np.array([1.0, 0.0, 0.0])
    # solve
    bloch.solve(m_init)
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
    ax1.set_title(r'$^3$He | B0 = {:.1f}$\mu$T'.format(bloch.B0*1e6))
    ax2.set_title(r'$\gamma<0$')
    ax3.set_title(r'$^3$He | fL = {:.1f}Hz'.format(bloch.larmor_f))
    ax1.legend(prop={"size": 8}, loc='lower right')
    ax3.legend(prop={"size": 8}, loc='lower right')
    # save figure
    if save_png:
        fig.savefig('example1b.png', dpi=300)


def test_two(save_png=False):
    """Perform the second test demonstrating the effect T1 and T2."""
    # create output figure
    fig = plt.figure(figsize=(6, 9))
    ax1 = fig.add_subplot(321)
    ax2 = fig.add_subplot(322, projection='3d')
    ax3 = fig.add_subplot(323)
    ax4 = fig.add_subplot(324, projection='3d')
    ax5 = fig.add_subplot(325)
    ax6 = fig.add_subplot(326, projection='3d')
    # initialize basic Bloch class
    bloch = basic(T1=0.001, T2=0.001, t_sim=0.005)
    # initial magnetization
    m_init = np.array([1.0, 0.0, 0.0])
    # solve
    bloch.solve(m_init)
    # create lab-frame plots
    bloch.plot_magnetization(ax1, version='lab')
    bloch.plot_magnetization3d(ax2, version='lab')

    bloch.T2 = 2*bloch.T1
    # solve
    bloch.solve(m_init)
    # create lab-frame plots
    bloch.plot_magnetization(ax3, version='lab')
    bloch.plot_magnetization3d(ax4, version='lab')

    bloch.T2 = bloch.T1/2
    # solve
    bloch.solve(m_init)
    # create lab-frame plots
    bloch.plot_magnetization(ax5, version='lab')
    bloch.plot_magnetization3d(ax6, version='lab')

    # adjust axes position
    low, bot, wid, hei = ax1.get_position().bounds
    ax1.set_position([low, bot*1.1, wid, hei*0.9])
    low, bot, wid, hei = ax3.get_position().bounds
    ax3.set_position([low, bot*1.1, wid, hei*0.9])
    low, bot, wid, hei = ax5.get_position().bounds
    ax5.set_position([low, bot*1.1, wid, hei*0.9])
    # adjust titles
    ax1.set_title(r'$^1$H | B0 = {:.1f}$\mu$T'.format(bloch.B0*1e6))
    ax3.set_title(r'$^1$H | B0 = {:.1f}$\mu$T'.format(bloch.B0*1e6))
    ax5.set_title(r'$^1$H | B0 = {:.1f}$\mu$T'.format(bloch.B0*1e6))
    ax2.set_title(r'T1 = T2')
    ax4.set_title(r'T2 = 2T1')
    ax6.set_title(r'T1 = 2T2')
    # adjust legends
    ax1.legend(prop={"size": 8}, loc='lower right')
    ax3.legend(prop={"size": 8}, loc='lower right')
    ax5.legend(prop={"size": 8}, loc='lower right')
    # save figure
    if save_png:
        fig.savefig('example1c.png', dpi=300)


# actually call the main function
if usetest[0]:
    test_one(SAVEPNG)
if usetest[1]:
    test_two(SAVEPNG)

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
