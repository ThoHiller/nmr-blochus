% script: example3.m
% This script demonstrates how pre-polarization switch-off ramps are simulated
% with BLOCHUS.
% Both examples were also used as benchmark cases in the paper:
%	Hiller, T., Dlugosch, R. and Müller-Petke, M., "Utilizing pre-
%	polarization to enhance SNMR signals - effect of imperfect
%	switch-off", Geophysical Journal International Vol. 222(2), p.815-826, 2020
%
% Two pre-polarization switch-off benchmarks are shown:
%
% 1.) from Melton et al., 1995, J Mag Res A, Vol. 117, p.164-170
%
%     general settings:
%     - magnetic field B0 is set to 50 µT
%     - pre-polarization field is a factor 100 larger
%     - uses a linear switch-off ramp with different switch-off rate (slopes)
%
% 2.) Conradi et al., 2017, J Mag Res, Vol. 281, p.241-245
%
%     general settings:
%     - magnetic field B0 is set to 50 µT
%     - pre-polarization field is a factor 50 larger
%     - uses a linexp switch-off ramp with different initial angles theta
%     (angle between B0 and Bp) and different switch-over fields Bstar
%
% See also: BLOCHUS
% Author: Thomas Hiller
% email: thomas.hiller[at]leibniz-liag.de
% License: GNU GPLv3 (at end)

%------------- BEGIN CODE --------------

% clear data
clear variables; clc;

% activate different tests
% 1 - Melton et al. 1995 benchmark
% 2 - Conradi et al. 2017 benchmark

% because this script takes a few minutes ask the user before continuing
answer = questdlg({'This script needs about 5 min. to run.',...
    'Do you want to continue?'},'Dessert Menu','Yes','No','Yes');
% handle response
switch answer
    case 'Yes'
        usetest = [1 1];
    case {'No',''}
        usetest = [0 0];
end

if sum(usetest) > 0
    % basic settings
    zunit = [0 0 1]'; % z unit vector
    % ODE solver error tolerance
    tol = 1e-9;
    % ODE solver options
    options = odeset('RelTol',tol,'AbsTol',[tol tol tol]);
end
%% 1.) Melton et al.(1995)
if usetest(1)
    % hydrogen proton
    nucleus = '1H';
    % pre-polarization factor [B0]
    PrePolFactor = 100;
    
    % parameter needed for the ODE solver
    % simulation type [string]
    odeparam.type = 'prepol';
    % equilibrium magnetization [A/m]
    odeparam.M0 = zunit;
    % primary (Earth's) magnetic field B0 [T]
    odeparam.B0 = 5e-5;
    % relaxation time T1 [s]
    odeparam.T1 = 2;
    % relaxation time T2 [s]
    odeparam.T2 = 1;
    % gyromagnetic ratio [rad/s/T]
    odeparam.gamma = getGyroRatio(nucleus);
    
    % pre-polarization switch-off ramp parameter
    % relaxation during switch-off [0/1]
    odeparam.RDS = 0;
    % switch-off ramp type [string]
    rampparam.ramp = 'lin';
    % amplitude of pre-polarization field
    rampparam.Bmax = odeparam.B0*PrePolFactor;
    % switch over magnetization (for linexp case)
    rampparam.Bstar = odeparam.B0;
    
    % initial orientation is toards y-axis
    orient = getRotationMatrixFromAngleandAxis(pi/2,[0 1 0])*zunit;
    orient = getRotationMatrixFromAngleandAxis(pi/2,[0 0 1])*orient;
    % initial magnetization in the direction of B0+Bp
    Morient = orient*rampparam.Bmax + odeparam.B0*zunit;
    % Minit does not get normalized in this case
    Minit = Morient;
    % because Minit is not normalized to 1, M0 needs to be adjusted
    odeparam.M0 = odeparam.B0*zunit;
    
    % orientation of the pre-polarization pulse axis (y-axis)
    odeparam.orient = orient;
    
    % predefined relaxation rates (k/GAMMA)
    swrates = [1/16 1/8 1/4 1/2 1 2 4 8 16];
    % initialize output variables
    M_final = zeros(numel(swrates),3) ;
    Phi_final = zeros(numel(swrates),1) ;
    Theta_final = zeros(numel(swrates),1) ;
    for r = 1:numel(swrates)
        
        % current relaxation rate
        rate = swrates(r);
        % calculate switch-off ramp time
        GAMMA = PrePolFactor/rate;
        Tramp = GAMMA/(odeparam.B0*odeparam.gamma);
        
        % switch-off ramp time [s]
        rampparam.Tramp = Tramp;
        % switch over ramp time [s] (for linexp case)
        rampparam.Tslope = Tramp;
        
        % add the ramp parameter to the ode parameter struct
        odeparam.rampparam = rampparam;
        
        % ODE solver call
        % OUTPUT: time T and magnetization M
        [T,M] = ode45(@(t,m) fcn_BLOCHUS_ode(t,m,odeparam),[0 Tramp],Minit,options);
        
        % calculate final angles describing the orientation of M
        % see p. 167 of the paper
        theta = 90 - acosd(M(end,3)/norm(M(end,:)));
        phi = 90 - atan2d(M(end,2),M(end,1));
        
        % output data
        M_final(r,:) = M(end,:)./odeparam.B0;
        Phi_final(r) = phi;
        Theta_final(r) = theta;
    end
    
    f1 = figure; set(f1,'Name','example3a_ref');
    % color the different switch-off rates
    col = jet(numel(swrates));
    ax1 = subplot(1,2,1);
    hold(ax1,'on');
    % plot a colored arrow showing the endpoint of the magnetization
    % this corresponds to Fig. 3 in Melton et al.(1995)
    for r = 1:numel(swrates)
        % normalize the arrow to 1
        M_final(r,:) = M_final(r,:)./norm(M_final(r,:));
        quiver3(0,0,0,M_final(r,1),M_final(r,2),M_final(r,3),0,...
            'Color',col(r,:),'LineWidth',2);
    end
    % plot the trace on the sphere surface
    plot3([0; M_final(:,1); 0],[0; M_final(:,2); 1],[1; M_final(:,3); 0],...
        'k--','LineWidth',2);
    % plot the Bloch sphere
    bsh = plotBSphere(30,30,ax1);
    text(-0.1,0.45,0.8,'30°','Parent',ax1);
    text(-0.1,0.8,0.4,'60°','Parent',ax1);
    set(ax1,'Color','w','XColor','none','YColor','none','ZColor','none');
    axis equal tight
    view([135 45]);
    set(ax1,'XLim',[0 1.2],'YLim',[0 1.2],'ZLim',[0 1.2]);
    
    % final angles for all switch-off rates
    % this corresponds to Fig. 4 in Melton et al.(1995)
    ax2 = subplot(1,2,2);
    hold(ax2,'on');
    for r = 1:numel(swrates)
        plot(swrates(r),Phi_final(r),'s','Color',col(r,:),'MarkerSize',8);
        plot(swrates(r),Theta_final(r),'o','Color',col(r,:),'MarkerSize',8);
    end
    % the fit from the paper eq.7
    xI = 0.5:16; yI = 50.8./sqrt(xI);
    plot(xI,yI,'k--');
    x = [0.675 0.775];
    y = [0.4 0.5];
    ah = annotation('textarrow',x,y,'String','$$50.8^\circ \sqrt(k/\Gamma)$$');
    set(ah,'Interpreter','latex','TextBackgroundColor','w','TextLineWidth',1);
    set(ax2,'XScale','log','YScale','log');
    set(ax2,'XLim',[0.04 30],'XTick',[0.1 1 10],'XTickLabel',{'0.1','1','10'})
    set(ax2,'YLim',[10 100])
    grid on; box on; axis square;
    legend('\phi_f','\theta_f');
    xlabel('relaxation rate k/\Gamma');ylabel('final angles [°]');
    
end

%% 2.) Conradi et al.(2017)
if usetest(2)
    
    % Fig. 4a of Conradi et al.(2017)
    %     conradi_Bmax = 16/1e3; % [T]
    %     conradi_Bstar = [2000 250 120]./1e6; % [T]
    %     conradi_Tslope = 6/1e3; % [s]
    %     conradi_lgdstr = {'2000','250','120'};
    
    % Fig. 4b of Conradi et al.(2017)
    conradi_Bmax = 2.5/1e3; % [T]
    conradi_Bstar = [250 60 1]./1e6; % [T]
    conradi_Tslope = 10/1e3; % [s]
    conradi_lgdstr = {'250','6','1'};
    
    % hydrogen proton
    nucleus = '1H';
    % pre-polarization factor [B0]
    PrePolFactor = 50;
    
    % parameter needed for the ODE solver
    % simulation type [string]
    odeparam.type = 'prepol';
    % equilibrium magnetization [A/m]
    odeparam.M0 = zunit;
    % primary (Earth's) magnetic field B0 [T]
    odeparam.B0 = 5e-5;
    % relaxation time T1 [s]
    odeparam.T1 = 2;
    % relaxation time T2 [s]
    odeparam.T2 = 1;
    % gyromagnetic ratio [rad/s/T]
    odeparam.gamma = getGyroRatio(nucleus);
    
    % pre-polarization switch-off ramp parameter
    % relaxation during switch-off [0/1]
    odeparam.RDS = 0;
    % switch-off ramp type [string]
    rampparam.ramp = 'linexp';
    % amplitude of pre-polarization field [T]
    rampparam.Bmax = conradi_Bmax;
    
    % switch-off ramp time [s]
    rampparam.Tramp = conradi_Tslope*2;
    % switch over ramp time [s] (for linexp case)
    rampparam.Tslope = conradi_Tslope;
    
    % because Minit is not normalized to 1, M0 needs to be adjusted
    odeparam.M0 = odeparam.B0*zunit;
    
    % predefined angles theta [deg]
    theta = linspace(0,200,101);
    % predefined switch over magnetization [T]
    Bstar = conradi_Bstar;
    % initialize output variable
    p = zeros(numel(theta),numel(Bstar));
    for b = 1:numel(Bstar)
        for r = 1:numel(theta)
            
            % initial orientation is towards y-axis
            orient = getRotationMatrixFromAngleandAxis(deg2rad(theta(r)),[0 1 0])*zunit;
            % initial magnetization in the direction of B0+Bp
            Morient = orient*rampparam.Bmax + odeparam.B0*zunit;
            % Minit does not get normalized in this case
            Minit = Morient;
            % orientation of the pre-polarization pulse axis (y-axis)
            odeparam.orient = orient;
            
            % switch over magnetization (for linexp case) [T]
            rampparam.Bstar = Bstar(b);
            
            % add the ramp parameter to the ode parameter struct
            odeparam.rampparam = rampparam;
            
            % ODE solver call
            % OUTPUT: time T and magnetization M
            [T,M] = ode45(@(t,m) fcn_BLOCHUS_ode(t,m,odeparam),[0 rampparam.Tramp],Minit,options);
            
            % normalized magnetization vector at end of switch-off ramp
            MMn = M(end,:)./norm(M(end,:));
            % "adiabatic quality" p of the switch-off ramp
            % -> orientation of M with respect to z_unit
            p(r,b) = dot(MMn(end,:),zunit)./norm(zunit);
            
            clc;
            disp(['Bstar: ',sprintf('%d',Bstar(b)*1e6),'µT',...
                ' Theta: ',sprintf('%d',theta(r)),'°']);
        end
    end
    
    %% this corresponds to Fig. 4b in Conradi et al.(2017)
    f2 = figure; set(f2,'Name','example3b_ref');
    col = [1 0 0;0 1 0;0 0 1];
    ax1 = axes;
    hold(ax1,'on');
    for b = 1:numel(Bstar)
        plot(theta,p(:,b),'Color',col(b,:));
    end
    set(ax1,'XLim',[-20 210],'XTick',0:45:180);
    set(ax1,'YLim',[-1.05 1.05]);
    grid on; box on;
    lh = legend(conradi_lgdstr,'Location','SouthWest');
    set(get(lh,'Title'),'String','B^\ast [µT]')
    xlabel('angle \theta (\angle B_0 B_p) [°]');
    ylabel('adiabatic quality p');
    
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
