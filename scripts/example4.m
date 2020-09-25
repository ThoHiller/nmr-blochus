% script: example4.m
% This script demonstrates how to calculate lookup-tables for pre-polarization
% switch-off ramps that can subsequently be used with MRSMatlab to forward
% model MRS sounding curves.
% Basically this script allows to create the data that is shown in Fig. 4
% of the paper:
%   Hiller, T., Dlugosch, R. and Müller-Petke, M., "Utilizing pre-
%   polarization to enhance SNMR signals - effect of imperfect
%   switch-off", Geophysical Journal International Vol. 222(2), p.815-826, 2020
%
%   general settings:
%       - magnetic field B0 is set to 48 µT
%       - T1, T2 - relaxation is ignored
%
%   lookup table parameters:
%       - pre-polarization field Bp ranges between 4.8µT and 480mT
%       - initial orientation between Bp and B0 varies between 0° and 180°
%
% See also: BLOCHUS
% Author: Thomas Hiller
% email: thomas.hiller[at]leibniz-liag.de
% License: GNU GPLv3 (at end)

%------------- BEGIN CODE --------------

% clear data
clear variables; clc;

% because this script takes approx. 30min, ask the user before continuing
answer = questdlg({'This script needs about 30 min. to run.',...
    'Do you want to continue?'},'Continue?','Yes','No','Yes');
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
% all four major ramp shapes
lookup_ramps = {'exp','linexp','halfcos','lin'};
% only three distinct switch-off times [s]
lookup_Tramps = [0.1 1 4]/1e3;

% IMPORTANT NOTE:
% This very coarse (and basically useless) Np=10 x Nt=10 grid takes for the
% 4x3 ramp combinations (shapes x times) roughly 30min! The resolution in
% Fig.4 of the paper is Np=51 x Nt=360 ... so you don't want to start this
% calculation on your normal machine.
% local lookup table settings (careful with the discretization)
Np = 10;
Nt = 10;
% pre-polarization factor [B0]
PPfac = logspace(-1,4,Np);
% angle between Bp and B0 [deg]
theta = linspace(0,180,Nt);
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
odeparam.type = 'prepol';
% equilibrium magnetization [A/m]
odeparam.M0 = zunit;
% primary (Earth's) magnetic field B0 [T]
odeparam.B0 = 4.8e-5;
% relaxation time T1 [s]
odeparam.T1 = 2;
% relaxation time T2 [s]
odeparam.T2 = 1;
% gyromagnetic ratio [rad/s/T]
odeparam.gamma = getGyroRatio(nucleus);

% pre-polarization switch-off ramp parameter
% relaxation during switch-off [0/1]
odeparam.RDS = 0;

% because Minit is not normalized to 1, M0 needs to be adjusted
odeparam.M0 = odeparam.B0*zunit;

data = cell(numel(lookup_ramps),numel(lookup_Tramps));
for i1 = 1:numel(lookup_ramps)
    for i2 = 1:numel(lookup_Tramps)
        % current switch-off ramp and time
        ramp = lookup_ramps{i1};
        Tramp = lookup_Tramps(i2);
        
        infostr1 = ['running: Ramp=',ramp,' | Tramp=',...
            sprintf('%3.1f',Tramp*1e3),' ms'];
        disp(infostr1);
        
        % create the output data matrix
        % adiabatic quality p
        p = zeros(numel(theta),numel(PPfac));
        % final magnetization vector
        Mfinal = zeros(numel(theta),numel(PPfac),3);
        % loop over the local parameters
        for tt = 1:numel(theta)
            for pp = 1:numel(PPfac)
                
                % switch-off ramp type [string]
                rampparam.ramp = ramp;
                % gyromagnetic ratio [rad/s/T]
                rampparam.gamma = odeparam.gamma;
                % primary (Earth's) magnetic field B0 [T]
                rampparam.B0 = odeparam.B0;
                % amplitude of pre-polarization field [T]
                rampparam.Bmax = PPfac(pp)*odeparam.B0;
                % switch-off ramp time [s]
                rampparam.Tramp = Tramp;
                
                % different ramps need different settings
                switch ramp
                    case 'exp'
                        SwitchFactor = 1;                        
                        % switch over ramp time [s]
                        rampparam.Tslope = Tramp/10;
                    case 'linexp'
                        SwitchFactor = PPfac(pp)/10;
                        rampparam.Tslope = Tramp/2;
                    case {'halfcos','lin'}
                        SwitchFactor = 1;
                        rampparam.Tslope = Tramp;
                end
                % switch-over amplitude for the "linexp" ramp (factor*B0) [T]
                rampparam.Bstar = SwitchFactor*odeparam.B0;

                % initial orientation is towards y-axis
                orient = getRotationMatrixFromAngleandAxis(deg2rad(theta(tt)),[0 1 0])*zunit;
                % initial magnetization in the direction of B0+Bp
                Morient = orient*rampparam.Bmax + odeparam.B0*zunit;
                % Minit does not get normalized in this case
                Minit = Morient;
                % orientation of the pre-polarization pulse axis (y-axis)
                odeparam.orient = orient;

                % add the ramp parameter to the ode parameter struct
                odeparam.rampparam = rampparam;
                
                % ODE solver call
                % OUTPUT: time T and magnetization M
                [T,M] = ode45(@(t,m) fcn_BLOCHUS_ode(t,m,odeparam),[0 rampparam.Tramp],Minit,options);
                
                % normalized magnetization vector at end of switch-off ramp
                MMn = M(end,:)./norm(M(end,:));
                % "adiabatic quality" p of the switch-off ramp
                % -> orientation of M with respect to z_unit
                p(tt,pp) = dot(MMn(end,:),zunit)./norm(zunit);
                % final orientation of M
                Mfinal(tt,pp,1:3) = M(end,:);
                
            end
            infostr2 = [' | theta=',num2str(theta(tt)),'°'];
            disp([infostr1,infostr2]);
        end
        disp([infostr1,' ... finished.']);
        
        % save the local lookup table for this particular combination of
        % ramp shape and ramp time
        data{i1,i2}.Tramp = Tramp;
        data{i1,i2}.ramp = ramp;
        data{i1,i2}.PPfac = PPfac;
        data{i1,i2}.theta = theta;
        data{i1,i2}.odeparam = odeparam;
        data{i1,i2}.M = Mfinal;
        data{i1,i2}.p = p;        
    end
end

% save lookup table
if save_data
    save('example4_lookup_table.mat','data');
end

% plot data
if plot_data
    [xx,yy] = meshgrid(linspace(0,180,Nt+1),logspace(-1,4,Np+1));
    figure;
    count = 0;
    for i1 = 1:numel(lookup_ramps)
        for i2 = 1:numel(lookup_Tramps)
            count = count + 1;
            p = data{i1,i2}.p;
            subplot(numel(lookup_ramps),numel(lookup_Tramps),count);
            surf(xx,yy,ones(size(xx)),p'); view(2);
            shading flat;
            set(gca,'XLim',[0 180],'XTick',linspace(0,180,5),'YScale','log',...
                'YLim',[0.1 1e4],'YTick',logspace(-1,4,6),'FontSize',10);
            xlabel('\theta [deg]');
            ylabel('Bp [B0]');
            title([lookup_ramps{i1},' | ',sprintf('%3.1f',lookup_Tramps(i2)*1e3),'ms']);
        end
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
