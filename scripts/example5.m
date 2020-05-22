% script: example5.m
% This script demonstrates how to calculate lookup-tables for an adiabatic
% half passage excitation pulse that can subsequently be used with MRSMatlab
% to forward model MRS sounding curves.
% Basically this script allows to create the data that is shown in Fig. 4b
% of the paper:
%   Grunewald, E., Grombacher, D. and Walsh, D., "Adiabatic pulses enhance
%   surface nuclear magnetic resonance measurement and survey speed for
%   groundwater investigations", Geophysics Vol. 81(4), WB85-WB96, 2016
%   DOI: 10.1190/GEO2015-0527.1
%
% NOTE: In this simulation the "My" component has a switched sign compared
% to Grunewald et al., 2016 due to the implementation of the reference phase
% and therewith accounting for the sign of the Larmor frequency (see e.g.
% Levitt, 2002). In practice, this has no effect on the result.
%
%   general settings:
%       - magnetic field B0 is set to 50 µT
%       - T1, T2 - relaxation is ignored
%
%   lookup table parameters for SWEEP1:
%       - B1 ranges between 0.01µT and 100µT
%       - frequency sweeps linearly from 2100 Hz to 2300 Hz, so from -200
%         to 0 Hz
%       - amplitude tuning with quality factor Q=10
%
% See also: BLOCHUS
% Author: Thomas Hiller
% email: thomas.hiller[at]leibniz-liag.de
% License: GNU GPLv3 (at end)

%------------- BEGIN CODE --------------

% clear data
clear variables; clc;

% because this script takes approx. 75min, ask the user before continuing
answer = questdlg({'This script needs about 75 min. to run.',...
    'Do you want to continue?'},'Dessert Menu','Yes','No','Yes');
% handle response
switch answer
    case 'Yes'
        % nothing to do
    case {'No',''}
        return
end

% output switches
save_data = true;
plot_data = true;

% --- LOOKUP TABLE SETTINGS -----------------------------------------------
% global lookup table settings (change as you like):
% B1 excitation amplitude [T]
% in the Grunewald et al., 2016 paper, there are in total 5000 points which
% would increase the calculation time to more than 6h!
N = 1000;
B1 = logspace(-8,-4,N);
% pulse length [s]
Ttau = 0.08;
% -------------------------------------------------------------------------

% basic settings
zunit = [0 0 1]'; % z unit vector
% ODE solver error tolerance
tol = 1e-9;
% ODE solver options
options = odeset('RelTol',tol,'AbsTol',[tol tol tol]);

% hydrogen proton
nucleus = '1H';

% parameter needed for the ODE solver
% simulation type [string]
odeparam.type = 'pulse';
% equilibrium magnetization [A/m]
odeparam.M0 = zunit;
% relaxation time T1 [s]
odeparam.T1 = 2;
% relaxation time T2 [s]
odeparam.T2 = 1;
% gyromagnetic ratio [rad/s/T]
odeparam.gamma = getGyroRatio(nucleus);
% primary (Earth's) magnetic field B0 [T]
odeparam.B0 = getB0(odeparam.gamma,2300*2*pi);

% pulse settings for BLOCHUS
% dummy time vector with sampling freq 50 kHz
tt = 0:1/50000:Ttau; % [s]
PulseType = 'AHP';

% pulse type [string]
fmod.PulseType = PulseType;
% shape of the modulation function [string]
fmod.shape = 'lin';
% time vector [s]
fmod.t = tt;
% pulse start time [s]
fmod.t0 = tt(1);
% pulse end time [s]
fmod.t1 = tt(end);
% start frequency [Hz]
fmod.v0 = -200;
% end frequency [Hz]
fmod.v1 = 0;
% modulation parameter (MIDI)
fmod.A = 1;
% modulation parameter (MIDI)
fmod.B = 0;

% pulse type [string]
Imod.PulseType = PulseType;
% shape of the modulation function [string]
Imod.shape = 'lin';
% check if quality factor tuning is activated [0/1]
Imod.useQ = 1;
% quality factor
Imod.Q = 10;
% quality factor off-resonance frequency [Hz]
Imod.Qdf = 0;
% time vector [s]
Imod.t = tt;
% pulse start time [s]
Imod.t0 = tt(1);
% pulse end time [s]
Imod.t1 = tt(end);
% start current [A]
Imod.v0 = 1;
% end current [A]
Imod.v1 = 1;
% modulation parameter (MIDI)
Imod.A = 1;
% modulation parameter (MIDI)
Imod.B = 0;

% output data
Mfinal = zeros(numel(B1),6);
for i = numel(B1):-1:1
    
    % relaxation during pulse [0/1]
    odeparam.RDP = 0;
    % pulse length [s]
    odeparam.Ttau = tt(end);
    
    % excitation pulse parameter
    % pulse type [string]
    pulseparam.PulseType = PulseType;
    % gyromagnetic ratio [rad/s/T]
    pulseparam.gamma = odeparam.gamma;
    % Larmor frequency [rad/s]
    pulseparam.omega0 = getOmega0(odeparam.gamma,odeparam.B0);
    % pulse amplitude [T]
    pulseparam.Amp = B1(i);
    % pulse frequency modulation [struct]    
    pulseparam.fmod = fmod;
    % pulse current modulation [struct]
    pulseparam.Imod = Imod;
    % auxiliary pulse phase [rad]
    pulseparam.phi = 0;
    % pulse axis [string]
    pulseparam.PulseAxis = '+y';
    % pulse polarization [string]
    pulseparam.PulsePolarization = 'circular';
    
    % save pulse parameter settings
    odeparam.pulseparam = pulseparam;
    
    % initial magnetization pointing towards B0 (+z) [A/m]
    Minit = zunit;  
    
    % ODE solver call
    % OUTPUT: time T and magnetization M
    [T,M] = ode45(@(t,m) fcn_BLOCHUS_ode(t,m,odeparam),[0 odeparam.Ttau],Minit,options);
    
    % in order to get M in the rotating frame of reference, get the pulse
    % phase theta for all simulation time steps
    pulseparam.fmod.t = T;
    pulseparam.Imod.t = T;
    pulseparam.t = T;
    [Bpulse,~,~,theta]= getPulseTimeSeries(pulseparam);
    
    % get M(T) in rotating frame of reference 
    Mrot = getMrot(M,theta);
    
    % final magnetization in lab and rot fame
    Mfinal(i,1:6) = [M(end,:) Mrot(end,:)];
    
    % show status
    disp([num2str(i),' / ',num2str(numel(B1)),' done']);    
end

% save lookup table
if save_data
    data.odeparam = odeparam;
    data.B1 = B1;
    data.Mfinal = Mfinal;
    save('example5_lookup_table.mat','data');
end

% plot data
if plot_data
    figure;
    subplot(121);
    semilogy(Mfinal(:,4),B1.*1e6);
    xlim([-1 1]);
    xlabel('Mx');
    ylabel('B1 [µT]');
    grid on;
    
    subplot(122);
    semilogy(Mfinal(:,5),B1.*1e6);
    xlim([-1 1]);
    xlabel('My');
    grid on;
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
