function data = BLOCHUS_loadDefaults
%BLOCHUS_loadDefaults loads default GUI data values
%
% Syntax:
%       BLOCHUS_loadDefaults
%
% Inputs:
%       none
%
% Outputs:
%       out - default data structure
%
% Example:
%       out = BLOCHUS_loadDefaults
%
% Other m-files required:
%       none
%
% Subfunctions:
%       getInitData
%
% MAT-files required:
%       none
%
% See also BLOCHUS
% Author: Thomas Hiller
% email: thomas.hiller[at]leibniz-liag.de
% License: GNU GPLv3 (at end)

%------------- BEGIN CODE --------------

% aux data
data.info.ToolTips = 1;
data.info.Timer = 0;

% get init data
% Note to user: use the "getInitData" function to adjust default
% settings and parameter range
init = getInitData;
data.init = init;

% basic data
data.basic.nucleus = init.nucleus;
data.basic.gamma = init.gamma;
data.basic.type = init.type;
data.basic.M0 = init.M0;
data.basic.Minit = init.Minit;
data.basic.B0 = init.B0(1);
data.basic.Omega0 = init.Omega0(1);
data.basic.T1relax = init.T1relax(1);
data.basic.T2relax = init.T2relax(1);
data.basic.Tsim = init.Tsim(1);

% pre-polarization data
data.prepol.Ramp = init.PrePolRamp;
data.prepol.RDS = init.PrePolRDS;
data.prepol.Factor = init.PrePolFactor(1);
data.prepol.Theta = init.PrePolTheta(1);
data.prepol.Phi = init.PrePolPhi(1);
data.prepol.SwitchFactor = init.PrePolSwitchFactor(1);
data.prepol.Tramp = init.PrePolTramp(1);
data.prepol.Tslope = init.PrePolTslope(1);

% pulse data
data.pulse.Type = init.PulseType;
data.pulse.Axis = init.PulseAxis;
data.pulse.Polarization = init.PulsePolarization;
data.pulse.RDP = init.PulseRDP;
data.pulse.Ttau = init.PulseTtau(1);
data.pulse.B1Factor = init.PulseB1Factor(1);
data.pulse.DFmode = init.PulseDFmode;
data.pulse.DFstart = init.PulseDFstart(1);
data.pulse.DFend = init.PulseDFend(1);
data.pulse.DFA = init.PulseDFA(1);
data.pulse.DFB = init.PulseDFB(1);
data.pulse.Imode = init.PulseImode;
data.pulse.Istart = init.PulseIstart(1);
data.pulse.Iend = init.PulseIend(1);
data.pulse.IA = init.PulseIA(1);
data.pulse.IB = init.PulseIB(1);
data.pulse.Q = init.PulseQ(1);
data.pulse.Qdf = init.PulseQdf(1);
data.pulse.Twait = init.PulseTwait(1);
data.pulse.MIDINP = init.PulseMIDINP;
data.pulse.MIDIsf = init.PulseMIDIsf;

end

% define init values and range
function init = getInitData
% --- BASIC ---
% proton [string]
init.nucleus = '1H';
% gyromagnetic ratio [rad/s/T]
init.gamma = getGyroRatio(init.nucleus);
% simulation type [string]
init.type = 'std';
% equilibrium magnetization [A/m]
init.M0 = [0 0 1];
% initial magnetization [A/m]
init.Minit = [1 0 0];
% Earth's magnetic field [T] default: 48킫 min: 1fT max: 1T
init.B0 = [48/1e6 1e-12 1];
% corresponding Larmor frequencies [Hz]
init.Omega0 = [getOmega0(init.gamma,init.B0(1))/2/pi -1e9 1e9];
% T1 relaxation time [ms] default: 100ms min: 1탎 max: 1000s
init.T1relax = [100 1e-3 1e6];
% T2 relaxation time [ms] default: 50ms min: 1탎 max: 1000s
init.T2relax = [50 1e-3 1e6];
% simulation time [ms] default: 50ms min: 1탎 max: 1000s
init.Tsim = [50 1e-3 1e6];

% --- PRE-POLARIZATION ---
% switch-off ramp shape [string]
init.PrePolRamp = 'exp';
% relaxation during switch-off [0/1]
init.PrePolRDS = 0;
% pre-polarization B-field amplitude in units of [B0]
init.PrePolFactor = [100 1e-6 1e6];
% polar angle of pre-polarization field [deg]
init.PrePolTheta = [90 1e-3 360];
% azimuthal angle of pre-polarization field [deg]
init.PrePolPhi = [0 0 360];
% pre-polarization B-field switch amplitude in units of [B0]
init.PrePolSwitchFactor = [1 1e-6 1e6];
% switch-off ramp time [ms] default: 1ms min: 1탎 max: 1000s
init.PrePolTramp = [1 1e-6 1e6];
% switch-off ramp slope time [ms] default: 0.1ms min: 1탎 max: 1000s
init.PrePolTslope = [0.1 1e-6 1e6];

% --- PULSE ---
% pulse type [string]
init.PulseType = 'pi_half';
% pulse axis [string]
init.PulseAxis = '+x';
% pulse polarization [string]
init.PulsePolarization = 'circular';
% relaxation during pulse [0/1]
init.PulseRDP = 0;
% pulse length [ms] default: 20ms min: 1탎 max: 1000s
init.PulseTtau = [20 1e-6 1e6];
% pulse amplitude in units of [B0] default: 0.0061163 min: 1e-12 max: 1e6
init.PulseB1Factor = [0.0061163 1e-12 1e6];
% frequency modulation function [string]
init.PulseDFmode = 'const';
% start of frequency sweep [Hz] default: 0Hz min: -1MHz max: 1MHz
init.PulseDFstart = [0 -1e6 1e6];
% end of frequency sweep [Hz] default: 0Hz min: -1MHz max: 1MHz
init.PulseDFend = [0 -1e6 1e6];
% slope parameter A
init.PulseDFA = [1 -1e6 1e6];
% slope parameter B
init.PulseDFB = [0 -1e6 1e6];
% current modulation function [string]
init.PulseImode = 'const';
% start of current sweep [A] default: 1A min: 0A max: 1A
init.PulseIstart = [1 0 1];
% end of current sweep [A] default: 1A min: 0A max: 1A
init.PulseIend = [1 0 1];
% slope parameter A
init.PulseIA = [1 -1e6 1e6];
% slope parameter B
init.PulseIB = [0 -1e6 1e6];
% quality factor Q
init.PulseQ = [0 0 100];
% quality factor off-resonance [Hz]
init.PulseQdf = [0 0 1e6];
% wait time between pre-polarization switch-off and pulse [ms]
init.PulseTwait = [2 0 1e6];
% number of periods (only for MIDI)
init.PulseMIDINP = 1;
% MIDI sampling frequency (fixed)
init.PulseMIDIsf = 50000;

end

%------------- END OF CODE --------------

%% License:
% GNU GPLv3
%
% BLOCHUS
% Copyright (C) 2019 Thomas Hiller
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
