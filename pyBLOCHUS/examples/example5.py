"""Create example5 plot.

This script demonstrates how to calculate lookup-tables for an adiabatic
half passage excitation pulse that can subsequently be used with MRSMatlab
to forward model MRS sounding curves.
Basically this script allows to create the data that is shown in Fig. 4b
of the paper:
    Grunewald, E., Grombacher, D. and Walsh, D., "Adiabatic pulses enhance
    surface nuclear magnetic resonance measurement and survey speed for
    groundwater investigations", Geophysics Vol. 81(4), WB85-WB96, 2016
    DOI: 10.1190/GEO2015-0527.1

NOTE: In this simulation the "My" component has a switched sign compared
to Grunewald et al., 2016 due to the implementation of the reference phase
and therewith accounting for the sign of the Larmor frequency (see e.g.
Levitt, 2002). In practice, this has no effect on the result.

    general settings:
        - magnetic field B0 is set to 50 µT
        - T1, T2 - relaxation is ignored

    lookup table parameters for SWEEP1:
        - B1 ranges between 0.01µT and 100µT
        - frequency sweeps linearly from 2100 Hz to 2300 Hz, so from -200
          to 0 Hz off-resonance
        - amplitude tuning with quality factor Q=10

See also: BLOCHUS
Author: Thomas Hiller
email: thomas.hiller[at]leibniz-liag.de
License: GNU GPLv3 (at end)
"""

import numpy as np
import matplotlib.pyplot as plt

from pyBLOCHUS import BlochusPulse as pulse

SAVEDATA = False
SAVEPNG = True


def test(save_png=False, save_data=False):
    """Perform the lookup table calculation."""
    # LOOKUP TABLE SETTINGS
    # global lookup table settings (change as you like):
    # B1 excitation amplitude [T]
    # in the Grunewald et al., 2016 paper, there are in total 5000 points which
    # would increase the calculation time to more than 10h!
    pulse_num = 1000
    b_pulse = np.logspace(-8, -4, pulse_num)
    # pulse length [s]
    t_pulse = 0.08

    # initialize pulse Bloch class
    bloch = pulse(t_sim=t_pulse)
    # set Larmor freq.
    bloch.larmor_f = -2300

    # pulse type
    bloch.Pulse.pulse_type = 'AHP'
    bloch.Pulse.axis = '+y'
    # frequency modulation
    bloch.Pulse.fmod.mod_type = 'lin'
    bloch.Pulse.fmod.v_range[0] = -200
    bloch.Pulse.fmod.v_range[1] = 0
    # current modulation
    bloch.Pulse.imod.mod_type = 'const'
    bloch.Pulse.imod.qual_fac = 10

    # create output dictionary
    mag_final = np.zeros((pulse_num, 3))
    # loop over all B1 values
    for i_1 in reversed(np.arange(pulse_num)):
        # because the B1 values are given in [T] they need to by divided by
        # B0 because the pulse factor is given in units of [B0]
        bloch.Pulse.factor = b_pulse[i_1]/bloch.B0
        # solve
        bloch.solve(m_init=bloch.zunit, use_numba=True)

        # command line output
        string = '{:d} / {:d}'.format(i_1+1, pulse_num)
        print(string)
        # save final magnetization orientation in rotating frame of reference
        mag_final[i_1, :] = bloch.mrot[:, -1]

    # create the output figure
    fig = plt.figure(figsize=(7, 5))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    # plot data
    ax1.plot(mag_final[:, 0], b_pulse*1e6)
    ax2.plot(mag_final[:, 1], b_pulse*1e6)
    # adjust axes
    ax1.set_yscale('log')
    ax2.set_yscale('log')
    ax1.set_xlim(-1, 1)
    ax2.set_xlim(-1, 1)
    ax1.set_ylim(b_pulse.min()*1e6, b_pulse.max()*1e6)
    ax2.set_ylim(b_pulse.min()*1e6, b_pulse.max()*1e6)
    ax1.set_xlabel(r'$M_x/M_0$')
    ax2.set_xlabel(r'$M_y/M_0$')
    ax1.set_ylabel(r'$B_1$ [$\mu T$]')
    ax1.grid(color='lightgray')
    ax2.grid(color='lightgray')

    # save data
    if save_data:
        data = {'B1': b_pulse, 'M': mag_final}
        np.save('example5_lookup-data', data, allow_pickle=True)
    # save figure
    if save_png:
        fig.savefig('example5.png', dpi=300)


# actually call the main function
# ask the user to proceed
QSTR = 'This script runs approx. 2h on a (Core i5-6200) Laptop. '
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
