% script: example1.m
% This script demonstrates the basic use of the BLOCHUS ode-solver. Two
% different protons (H & He) are used to show the effect of the sign of the
% gyromagnetic ratio \gamma
% 1.) Hydrogen proton (gyromagnetic ratio > 0 )
% 2.) Helium proton (gyromagnetic ratio < 0 )
%
% general settings:
% -> magnetic field B0 is set to 50 µT
%   -> fL(H)  = -2128.9 Hz
%   -> fL(He) =  1621.7 Hz
% -> simulation time is 5 ms
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
% 1 - using different protons with different signs in \gamma
%     this will create two figures corresponding to "example1a_ref" and "example1b_ref"
% 2 - comparison of relaxation time effect
%     this will create one figure corresponding to "example1c_ref"
usetest = [1 1];

%% basic settings
% get default values
% defaults = BLOCHUS_loadDefaults('basic');
zunit = [0 0 1]'; % z unit vector
% ODE solver error tolerance
tol = 1e-9;
% ODE solver options
options = odeset('RelTol',tol,'AbsTol',[tol tol tol]);

%% 1.) Hydrogen proton vs Helium proton
if usetest(1)
    nucleus = '1H'; % hydrogen proton
    
    % parameter needed for the ODE
    odeparam.type = 'std';
    odeparam.M0 = zunit; % equilibrium magnetization [A/m]
    odeparam.B0 = 5e-5; % primary magnetic field B0 [T]
    odeparam.T1 = 0.0010; % relaxation time T1 [s]
    odeparam.T2 = 0.0005; % relaxation time T2 [s]
    odeparam.gamma = getGyroRatio(nucleus); % gyromagnetic ratio [rad/s/T]
    
    % initial magnetization pointing towards +x
    Minit = [1 0 0]; % [A/m]
    % simulation time
    Tsim = 0.005; % [s]
    
    % ODE solver call
    % OUTPUT: time T and magnetization M
    [T,M] = ode45(@(t,m) fcn_BLOCHUS_ode(t,m,odeparam),[0 Tsim],Minit,options);
    
    % rotate M(T) vector into rotating reference frame
    Mrot = getMrot(M,getOmega0(odeparam.gamma,odeparam.B0).*T);
    
    f1 = figure; set(f1,'Name','example1a_ref');
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
    set(get(ax2,'Title'),'String','rot frame | ^{1}H');
    
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
    
    % --- Helium ---
    nucleus = '3He'; % helium proton
    odeparam.gamma = getGyroRatio(nucleus); % gyromagnetic ratio [rad/s/T]
    
    % ODE solver call
    % OUTPUT: time T and magnetization M
    [T,M] = ode45(@(t,m) fcn_BLOCHUS_ode(t,m,odeparam),[0 Tsim],Minit,options);
    
    % rotate M(T) vector into rotating reference frame
    Mrot = getMrot(M,getOmega0(odeparam.gamma,odeparam.B0).*T);
    
    f2 = figure; set(f2,'Name','example1b_ref');
    ax5 = subplot(2,2,1); hold(ax5,'on');
    plot(T.*1000,M(:,1),'r','Parent',ax5);
    plot(T.*1000,M(:,2),'g','Parent',ax5);
    plot(T.*1000,M(:,3),'b','Parent',ax5);
    lgh = legend('Mx','My','Mz','Location','SouthEast');
    set(lgh,'FontSize',10);
    set(ax5,'XLim',[0 5],'XTick',0:1:5,'XTickLabel',{'0','1','2','3','4','5 ms'});
    set(ax5,'YLim',[-1 1],'YTick',-1:0.25:1,'YTickLabel',{'-1','','-0.5','','0','','0.5','','1'});
    set(get(ax5,'Title'),'String','lab frame | ^{3}He | B0 = 50µT');
    
    ax6 = subplot(2,2,3); hold(ax6,'on');
    plot(T.*1000,Mrot(:,1),'r','Parent',ax6);
    plot(T.*1000,Mrot(:,2),'g','Parent',ax6);
    plot(T.*1000,Mrot(:,3),'b','Parent',ax6);
    set(ax6,'XLim',[0 5],'XTick',0:1:5,'XTickLabel',{'0','1','2','3','4','5 ms'});
    set(ax6,'YLim',[-1 1],'YTick',-1:0.25:1,'YTickLabel',{'-1','','-0.5','','0','','0.5','','1'});
    set(get(ax6,'Title'),'String','rot frame | ^{3}He');
    
    ax7 = subplot(2,2,2); hold(ax7,'on');
    plot3(M(:,1),M(:,2),M(:,3),'k');
    plotBSphere(30,30,ax7);
    view([134 30]); hold off;
    set(ax7,'Color','w','XColor','none','YColor','none',...
        'ZColor','none');
    axis equal tight;
    set(get(ax7,'Title'),'String','\gamma < 0');
    
    ax8 = subplot(2,2,4); hold(ax8,'on');
    plot3(Mrot(:,1),Mrot(:,2),Mrot(:,3),'k');
    plotBSphere(30,30,ax8);
    view([134 30]); hold off;
    set(ax8,'Color','w','XColor','none','YColor','none',...
        'ZColor','none');
    axis equal tight;
end

%% 2.) relaxation time effect
% here we basically re-plot the Figs. 11.31 - 11.33 (p. 286/287)
% from Levitt, 2008, "spin dynamics" 2nd ed, to show the effect of varying
% T1/T2 ratios
if usetest(2)
    nucleus = '1H'; % hydrogen proton
    
    % --- T2 = T1
    % parameter needed for the ODE
    odeparam.type = 'std';
    odeparam.M0 = zunit; % equilibrium magnetization [A/m]
    odeparam.B0 = 5e-5; % primary magnetic field B0 [T]
    odeparam.T1 = 0.0010; % relaxation time T1 [s]
    odeparam.T2 = 0.0010; % relaxation time T2 [s]
    odeparam.gamma = getGyroRatio(nucleus); % gyromagnetic ratio [rad/s/T]
    
    % initial magnetization pointing towards +x
    Minit = [1 0 0]; % [A/m]
    % simulation time
    Tsim = 0.005; % [s]
    
    % ODE solver call
    % OUTPUT: time T and magnetization M
    [T,M] = ode45(@(t,m) fcn_BLOCHUS_ode(t,m,odeparam),[0 Tsim],Minit,options);
    
    % rotate M(T) vector into rotating reference frame
    Mrot = getMrot(M,getOmega0(odeparam.gamma,odeparam.B0).*T);
    
    f3 = figure; set(f3,'Name','example1c_ref');
    ax1 = subplot(3,2,1); hold(ax1,'on');
    plot(T.*1000,M(:,1),'r','Parent',ax1);
    plot(T.*1000,M(:,2),'g','Parent',ax1);
    plot(T.*1000,M(:,3),'b','Parent',ax1);
    set(ax1,'XLim',[0 5],'XTick',0:1:5,'XTickLabel',{'0','1','2','3','4','5 ms'});
    set(ax1,'YLim',[-1 1],'YTick',-1:0.25:1,'YTickLabel',{'-1','','-0.5','','0','','0.5','','1'});
    set(get(ax1,'Title'),'String','lab frame | ^{1}H | B0 = 50µT');
    
    ax2 = subplot(3,2,2); hold(ax2,'on');
    plot3(M(:,1),M(:,2),M(:,3),'k');
    plotBSphere(30,30,ax2);
    view([134 30]); hold off;
    set(ax2,'Color','w','XColor','none','YColor','none',...
        'ZColor','none');
    axis equal tight
    set(get(ax2,'Title'),'String','T1 = T2');
    
    % --- T2 = 2*T1
    odeparam.T2 = 0.0020; % relaxation time T2 [s]

    % ODE solver call
    % OUTPUT: time T and magnetization M
    [T,M] = ode45(@(t,m) fcn_BLOCHUS_ode(t,m,odeparam),[0 Tsim],Minit,options);
    
    % rotate M(T) vector into rotating reference frame
    Mrot = getMrot(M,getOmega0(odeparam.gamma,odeparam.B0).*T);
    
    ax3 = subplot(3,2,3); hold(ax3,'on');
    plot(T.*1000,M(:,1),'r','Parent',ax3);
    plot(T.*1000,M(:,2),'g','Parent',ax3);
    plot(T.*1000,M(:,3),'b','Parent',ax3);
    lgh = legend('Mx','My','Mz','Location','SouthEast');
    set(lgh,'FontSize',10);
    set(ax3,'XLim',[0 5],'XTick',0:1:5,'XTickLabel',{'0','1','2','3','4','5 ms'});
    set(ax3,'YLim',[-1 1],'YTick',-1:0.25:1,'YTickLabel',{'-1','','-0.5','','0','','0.5','','1'});
    set(get(ax3,'Title'),'String','lab frame | ^{1}H | B0 = 50µT');
    
    ax4 = subplot(3,2,4); hold(ax4,'on');
    plot3(M(:,1),M(:,2),M(:,3),'k');
    plotBSphere(30,30,ax4);
    view([134 30]); hold off;
    set(ax4,'Color','w','XColor','none','YColor','none',...
        'ZColor','none');
    axis equal tight
    set(get(ax4,'Title'),'String','T2 = 2T1');
    
    % --- T1 = 2*T2
    odeparam.T2 = 0.0005; % relaxation time T2 [s]

    % ODE solver call
    % OUTPUT: time T and magnetization M
    [T,M] = ode45(@(t,m) fcn_BLOCHUS_ode(t,m,odeparam),[0 Tsim],Minit,options);
    
    % rotate M(T) vector into rotating reference frame
    Mrot = getMrot(M,getOmega0(odeparam.gamma,odeparam.B0).*T);
    
    ax5 = subplot(3,2,5); hold(ax5,'on');
    plot(T.*1000,M(:,1),'r','Parent',ax5);
    plot(T.*1000,M(:,2),'g','Parent',ax5);
    plot(T.*1000,M(:,3),'b','Parent',ax5);
    set(ax5,'XLim',[0 5],'XTick',0:1:5,'XTickLabel',{'0','1','2','3','4','5 ms'});
    set(ax5,'YLim',[-1 1],'YTick',-1:0.25:1,'YTickLabel',{'-1','','-0.5','','0','','0.5','','1'});
    set(get(ax5,'Title'),'String','lab frame | ^{1}H | B0 = 50µT');
    
    ax6 = subplot(3,2,6); hold(ax6,'on');
    plot3(M(:,1),M(:,2),M(:,3),'k');
    plotBSphere(30,30,ax6);
    view([134 30]); hold off;
    set(ax6,'Color','w','XColor','none','YColor','none',...
        'ZColor','none');
    axis equal tight
    set(get(ax6,'Title'),'String','T1 = 2T2');
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
