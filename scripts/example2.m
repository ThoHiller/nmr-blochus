% script: example2.m
% this script demonstrates how basic pulses are used with BLOCHUS
% everything is done with a Hydrogen proton (gyromagnetic ratio > 0 )
%
% settings:
% -> magnetic field B0 is set to 50 µT
%   -> fL(H)  = -2128.9 Hz
% -> pulse length is set to 5 ms
%   -> different pulses are used (pi/2, pi)
%   -> effect of different pulse axis
%   -> effect of frequency offsets
% -> T1 and T2 relaxation times are 1 ms and 0.5 ms, respectively
%
% See also: BLOCHUS
% Author: Thomas Hiller
% email: thomas.hiller[at]leibniz-liag.de
% License: GNU GPLv3 (at end)

%------------- BEGIN CODE --------------

% clear data
clear variables; clc;

% activate different tests
% 1 - pi/2 - pulse
%     this will create a figure corresponding to "example2a_ref"
% 2 - pi - pulse
%     this will create a figure corresponding to "example2b_ref"
% 3 - pi/2 - pulses with different pulse axes
%     this will create a figure corresponding to "example2c_ref"
% 4 - pi/2 pulses with different frequency offsets
%     this will create a figure corresponding to "example2d_ref"
usetest = [1 1 1 1];

%% basic settings
zunit = [0 0 1]'; % z unit vector
% ODE solver error tolerance
tol = 1e-9;
% ODE solver options
options = odeset('RelTol',tol,'AbsTol',[tol tol tol]);

%% standard parameter
nucleus = '1H'; % hydrogen proton
% parameter needed for the ODE
odeparam.type = 'pulse';
odeparam.M0 = zunit; % equilibrium magnetization [A/m]
odeparam.B0 = 5e-5; % primary magnetic field B0 [T]
odeparam.T1 = 0.0010; % relaxation time T1 [s]
odeparam.T2 = 0.0005; % relaxation time T2 [s]
odeparam.gamma = getGyroRatio(nucleus); % gyromagnetic ratio [rad/s/T]

%% pulse modulation standard values
fmod.PulseType = 'pi_half';
fmod.shape = 'const';
fmod.t0 = 0;
fmod.t1 = 0.005;
fmod.t = 0;
fmod.v0 = 0;
fmod.v1 = 0;
fmod.A = 1;
fmod.B = 0;

Imod.shape = 'const';
Imod.useQ = 0;
Imod.Q = 0;
Imod.Qdf = 0;
Imod.t0 = 0;
Imod.t1 = 0.005;
Imod.t = 0;
Imod.v0 = 1;
Imod.v1 = 1;
Imod.A = 1;
Imod.B = 0;

%% 1.) pi/2 pulse
if usetest(1)
    % excitation pulse parameter
    odeparam.RDP = 0;
    odeparam.Ttau = 0.005; % [s]
    pulseparam.PulseType = fmod.PulseType;
    pulseparam.Amp = odeparam.B0*abs((pi/2/(odeparam.gamma*odeparam.B0*odeparam.Ttau)));
    pulseparam.gamma = odeparam.gamma;
    pulseparam.omega0 = getOmega0(odeparam.gamma,odeparam.B0);
    pulseparam.fmod = fmod;
    pulseparam.Imod = Imod;
    pulseparam.phi = 0;
    pulseparam.PulseAxis = '+x';
    pulseparam.PulsePolarization = 'circular';
    
    odeparam.pulseparam = pulseparam;
    
    % initial magnetization pointing towards B0 (+z)
    Minit = zunit; % [A/m]
    % simulation time
    Tsim = odeparam.Ttau; % [s]
    
    % ODE solver call
    % OUTPUT: time T and magnetization M
    [T,M] = ode45(@(t,m) fcn_BLOCHUS_ode(t,m,odeparam),[0 Tsim],Minit,options);
    
    % rotate M(T) vector into rotating reference frame
    Mrot = getMrot(M,getOmega0(odeparam.gamma,odeparam.B0).*T);
    
    f1 = figure; set(f1,'Name','example2a_ref');
    ax1 = subplot(2,2,1); hold(ax1,'on');
    plot(T.*1000,M(:,1),'r','Parent',ax1);
    plot(T.*1000,M(:,2),'g','Parent',ax1);
    plot(T.*1000,M(:,3),'b','Parent',ax1);
    lgh = legend('Mx','My','Mz','Location','SouthEast');
    set(lgh,'FontSize',10);
    set(ax1,'XLim',[0 5],'XTick',0:1:5,'XTickLabel',{'0','1','2','3','4','5 ms'});
    set(ax1,'YLim',[-1 1],'YTick',-1:0.25:1,'YTickLabel',{'-1','','-0.5','','0','','0.5','','1'});
    set(get(ax1,'Title'),'String','lab frame | ^{1}H | B0 = 50µT');
    
    ax2 = subplot(2,2,3); hold(ax2,'on');
    plot(T.*1000,Mrot(:,1),'r','Parent',ax2);
    plot(T.*1000,Mrot(:,2),'g','Parent',ax2);
    plot(T.*1000,Mrot(:,3),'b','Parent',ax2);
    set(ax2,'XLim',[0 5],'XTick',0:1:5,'XTickLabel',{'0','1','2','3','4','5 ms'});
    set(ax2,'YLim',[-1 1],'YTick',-1:0.25:1,'YTickLabel',{'-1','','-0.5','','0','','0.5','','1'});
    set(get(ax2,'Title'),'String','rot frame | ^{1}H | B1 ~ 1.17µT');
    
    ax3 = subplot(2,2,2); hold(ax3,'on');
    plot3(M(:,1),M(:,2),M(:,3),'k');
    plotBSphere(30,30,ax3);
    view([134 30]); hold off;
    set(ax3,'Color','w','XColor','none','YColor','none',...
        'ZColor','none');
    axis equal tight
    set(get(ax3,'Title'),'String','\gamma > 0');
    
    ax4 = subplot(2,2,4); hold(ax4,'on');
    plot3(Mrot(:,1),Mrot(:,2),Mrot(:,3),'k');
    plotBSphere(30,30,ax4);
    view([134 30]); hold off;
    set(ax4,'Color','w','XColor','none','YColor','none',...
        'ZColor','none');
    axis equal tight;
    set(get(ax4,'Title'),'String','(\pi/2)_{+x} - pulse');
end

%% 2.) pi pulse
if usetest(2)
    fmod.PulseType = 'pi';
    % excitation pulse parameter
    odeparam.RDP = 0;
    odeparam.Ttau = 0.005; % [s]
    pulseparam.PulseType = fmod.PulseType;
    pulseparam.Amp = odeparam.B0*abs((pi/(odeparam.gamma*odeparam.B0*odeparam.Ttau)));
    pulseparam.gamma = odeparam.gamma;
    pulseparam.omega0 = getOmega0(odeparam.gamma,odeparam.B0);
    pulseparam.fmod = fmod;
    pulseparam.Imod = Imod;
    pulseparam.phi = 0;
    pulseparam.PulseAxis = '+x';
    pulseparam.PulsePolarization = 'circular';
    
    odeparam.pulseparam = pulseparam;
    
    % initial magnetization pointing towards B0 (+z)
    Minit = zunit; % [A/m]
    % simulation time
    Tsim = odeparam.Ttau; % [s]
    
    % ODE solver call
    % OUTPUT: time T and magnetization M
    [T,M] = ode45(@(t,m) fcn_BLOCHUS_ode(t,m,odeparam),[0 Tsim],Minit,options);
    
    % rotate M(T) vector into rotating reference frame
    Mrot = getMrot(M,getOmega0(odeparam.gamma,odeparam.B0).*T);
    
    f2 = figure; set(f2,'Name','example2b_ref');
    ax1 = subplot(2,2,1); hold(ax1,'on');
    plot(T.*1000,M(:,1),'r','Parent',ax1);
    plot(T.*1000,M(:,2),'g','Parent',ax1);
    plot(T.*1000,M(:,3),'b','Parent',ax1);
    lgh = legend('Mx','My','Mz','Location','SouthEast');
    set(lgh,'FontSize',10);
    set(ax1,'XLim',[0 5],'XTick',0:1:5,'XTickLabel',{'0','1','2','3','4','5 ms'});
    set(ax1,'YLim',[-1 1],'YTick',-1:0.25:1,'YTickLabel',{'-1','','-0.5','','0','','0.5','','1'});
    set(get(ax1,'Title'),'String','lab frame | ^{1}H | B0 = 50µT');
    
    ax2 = subplot(2,2,3); hold(ax2,'on');
    plot(T.*1000,Mrot(:,1),'r','Parent',ax2);
    plot(T.*1000,Mrot(:,2),'g','Parent',ax2);
    plot(T.*1000,Mrot(:,3),'b','Parent',ax2);
    set(ax2,'XLim',[0 5],'XTick',0:1:5,'XTickLabel',{'0','1','2','3','4','5 ms'});
    set(ax2,'YLim',[-1 1],'YTick',-1:0.25:1,'YTickLabel',{'-1','','-0.5','','0','','0.5','','1'});
    set(get(ax2,'Title'),'String','rot frame | ^{1}H | B1 ~ 2.35µT');
    
    ax3 = subplot(2,2,2); hold(ax3,'on');
    plot3(M(:,1),M(:,2),M(:,3),'k');
    plotBSphere(30,30,ax3);
    view([134 30]); hold off;
    set(ax3,'Color','w','XColor','none','YColor','none',...
        'ZColor','none');
    axis equal tight
    set(get(ax3,'Title'),'String','\gamma > 0');
    
    ax4 = subplot(2,2,4); hold(ax4,'on');
    plot3(Mrot(:,1),Mrot(:,2),Mrot(:,3),'k');
    plotBSphere(30,30,ax4);
    view([134 30]); hold off;
    set(ax4,'Color','w','XColor','none','YColor','none',...
        'ZColor','none');
    axis equal tight;
    set(get(ax4,'Title'),'String','(\pi)_{+x} - pulse');
end

%% 3.) pi/2 pulses with different pulse axes
if usetest(3)
    fmod.PulseType = 'pi_half';
    % excitation pulse parameter
    odeparam.RDP = 0;
    odeparam.Ttau = 0.005; % [s]
    pulseparam.PulseType = fmod.PulseType;
    pulseparam.Amp = odeparam.B0*abs((pi/2/(odeparam.gamma*odeparam.B0*odeparam.Ttau)));
    pulseparam.gamma = odeparam.gamma;
    pulseparam.omega0 = getOmega0(odeparam.gamma,odeparam.B0);
    pulseparam.fmod = fmod;
    pulseparam.Imod = Imod;
    pulseparam.phi = 0;
    pulseparam.PulseAxis = '+x';
    pulseparam.PulsePolarization = 'circular';
    
    odeparam.pulseparam = pulseparam;
    
    % initial magnetization pointing towards B0 (+z)
    Minit = zunit; % [A/m]
    % simulation time
    Tsim = odeparam.Ttau; % [s]
    
    % ODE solver call
    % OUTPUT: time T and magnetization M
    [T,M] = ode45(@(t,m) fcn_BLOCHUS_ode(t,m,odeparam),[0 Tsim],Minit,options);
    
    % rotate M(T) vector into rotating reference frame
    Mrot = getMrot(M,getOmega0(odeparam.gamma,odeparam.B0).*T);
    
    f3 = figure; set(f3,'Name','example2c_ref');
    ax1 = subplot(2,2,1); hold(ax1,'on');
    plot3(Mrot(:,1),Mrot(:,2),Mrot(:,3),'k');
    plotBSphere(30,30,ax1);
    view([134 30]); hold off;
    set(ax1,'Color','w','XColor','none','YColor','none','ZColor','none');
    axis equal tight;
    set(get(ax1,'Title'),'String','(\pi/2)_{+x} - pulse');
    
    % +y pulse
    pulseparam.PulseAxis = '+y';
    odeparam.pulseparam = pulseparam;
    [T,M] = ode45(@(t,m) fcn_BLOCHUS_ode(t,m,odeparam),[0 Tsim],Minit,options);
    Mrot = getMrot(M,getOmega0(odeparam.gamma,odeparam.B0).*T);
    
    ax2 = subplot(2,2,2); hold(ax2,'on');
    plot3(Mrot(:,1),Mrot(:,2),Mrot(:,3),'k');
    plotBSphere(30,30,ax2);
    view([134 30]); hold off;
    set(ax2,'Color','w','XColor','none','YColor','none','ZColor','none');
    axis equal tight;
    set(get(ax2,'Title'),'String','(\pi/2)_{+y} - pulse');
    
    % -x pulse
    pulseparam.PulseAxis = '-x';
    odeparam.pulseparam = pulseparam;
    [T,M] = ode45(@(t,m) fcn_BLOCHUS_ode(t,m,odeparam),[0 Tsim],Minit,options);
    Mrot = getMrot(M,getOmega0(odeparam.gamma,odeparam.B0).*T);
    
    ax3 = subplot(2,2,3); hold(ax3,'on');
    plot3(Mrot(:,1),Mrot(:,2),Mrot(:,3),'k');
    plotBSphere(30,30,ax3);
    view([134 30]); hold off;
    set(ax3,'Color','w','XColor','none','YColor','none','ZColor','none');
    axis equal tight;
    set(get(ax3,'Title'),'String','(\pi/2)_{-x} - pulse');
    
    % -y pulse
    pulseparam.PulseAxis = '-y';
    odeparam.pulseparam = pulseparam;
    [T,M] = ode45(@(t,m) fcn_BLOCHUS_ode(t,m,odeparam),[0 Tsim],Minit,options);
    Mrot = getMrot(M,getOmega0(odeparam.gamma,odeparam.B0).*T);
    
    ax4 = subplot(2,2,4); hold(ax4,'on');
    plot3(Mrot(:,1),Mrot(:,2),Mrot(:,3),'k');
    plotBSphere(30,30,ax4);
    view([134 30]); hold off;
    set(ax4,'Color','w','XColor','none','YColor','none','ZColor','none');
    axis equal tight;
    set(get(ax4,'Title'),'String','(\pi/2)_{-y} - pulse');
end

%% 4.) pi/2 pulses with different frequency offsets
% here we basically re-plot Fig. 10.28 (p. 255) from
% Levitt, 2008, "spin dynamics" 2nd ed
% positive ratio means the pulse frequency is "higher" as absolute value
% (i.e. instead of -2000 Hz it is -2100 Hz)
% because of the negative Larmor frequency of H-protons the rotation axis
% dips downwards for negative ratios
if usetest(4)    
    % df/omega_nut ratio
    ratio = -4:1:4;
    
    fmod.PulseType = 'pi_half';
    % excitation pulse parameter
    odeparam.RDP = 0;
    odeparam.Ttau = 0.005; % [s]
    pulseparam.PulseType = fmod.PulseType;
    pulseparam.Amp = odeparam.B0*abs((pi/2/(odeparam.gamma*odeparam.B0*odeparam.Ttau)));
    pulseparam.gamma = odeparam.gamma;
    pulseparam.omega0 = getOmega0(odeparam.gamma,odeparam.B0);
    pulseparam.fmod = fmod;
    pulseparam.Imod = Imod;
    pulseparam.phi = 0;
    pulseparam.PulseAxis = '+y';
    pulseparam.PulsePolarization = 'circular';
    
    % nutation frequency for a pi/2 pulse with 5 ms is 50 Hz
    omega_nut = (pi/2)/odeparam.Ttau/2/pi; % [Hz]
    
    f4 = figure; set(f4,'Name','example2d_ref');
    for i = 1:numel(ratio)
        
        df = ratio(i)*omega_nut; % frequency offsets [Hz]
        pulseparam.fmod.v0 = df;
        pulseparam.fmod.v1 = df;
        
        odeparam.pulseparam = pulseparam;
        
        % initial magnetization pointing towards B0 (+z)
        Minit = zunit; % [A/m]
        % simulation time
        Tsim = odeparam.Ttau; % [s]
        
        % ODE solver call
        % OUTPUT: time T and magnetization M
        [T,M] = ode45(@(t,m) fcn_BLOCHUS_ode(t,m,odeparam),[0 Tsim],Minit,options);
        
        % in order to get M in the rotating frame get the pulse phase theta
        % for all time steps
        pulseparam.fmod.t = T;
        pulseparam.Imod.t = T;
        pulseparam.t = T;
        [Bpulse,~,~,theta]= getPulseTimeSeries(pulseparam);
        
        % rotate M(T) vector into rotating reference frame
        Mrot = getMrot(M,theta);
        
        ax(i) = subplot(3,3,i); hold(ax(i),'on'); %#ok<*SAGROW>
        plot3(Mrot(:,1),Mrot(:,2),Mrot(:,3),'k');
        plotBSphere(30,30,ax(i));
        view([134 30]); hold off;
        set(ax(i),'Color','w','XColor','none','YColor','none','ZColor','none');
        axis equal;
        set(get(ax(i),'Title'),'String',['\Omega_0 / \omega_{nut} = ',num2str(ratio(i))]);
    end
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
