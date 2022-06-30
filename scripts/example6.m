% script: example6.m
% This script demonstrates how to calculate a lookup-table for an arbitrarily
% shaped pre-polarization switch-off ramp with BLOCHUS.
%
% The switch-off ramp is imported as 'time' and 'current' data and used
% inside BLOCHUS via interpolation
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
% three custom ramps
lookup_name = {'ramp_2ms','ramp_3ms','ramp_3ms_raw'};

% IMPORTANT NOTE:
% This very coarse (and basically useless) Np=10 x Nt=10 grid takes for the
% 3 ramps roughly 30min!
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

data = cell(1,numel(lookup_name));
for i1 = 1:numel(lookup_name)
    % load ramp file
    dataR = load([lookup_name{i1},'.mat']);
    dataR = dataR.data;
    Tramp = dataR.t(end);
    
    % create the output data matrix
    % adiabatic quality p
    p = zeros(numel(theta),numel(PPfac));
    % final magnetization vector
    Mfinal = zeros(numel(theta),numel(PPfac),3);
    % loop over the local parameters
    for tt = 1:numel(theta)
        for pp = 1:numel(PPfac)
            
            % switch-off ramp type [string]
            rampparam.ramp = 'custom';
            % gyromagnetic ratio [rad/s/T]
            rampparam.gamma = odeparam.gamma;
            % primary (Earth's) magnetic field B0 [T]
            rampparam.B0 = odeparam.B0;
            % amplitude of pre-polarization field [T]
            rampparam.Bmax = PPfac(pp)*odeparam.B0;
            % switch-off ramp time [s]
            rampparam.Tramp = Tramp;
            % ramp interpolation data
            rampparam.interp.t = dataR.t;
            rampparam.interp.I = dataR.I;
            
            % switch-off slope time is NOT NEEDED [s]
            rampparam.Tslope = Tramp;
            % switch-over amplitude is NOT NEEDED [T]
            rampparam.Bstar = odeparam.B0;
            
            % initial orientation is towards x-axis
            orient = getRotationMatrixFromAngleandAxis(deg2rad(theta(tt)),[0 1 0])*zunit;
            % initial magnetization in the direction of B0+Bp
            Morient = orient*rampparam.Bmax + odeparam.B0*zunit;
            % Minit does not get normalized in this case
            Minit = Morient;
            % orientation of the pre-polarization pulse axis (x-axis)
            odeparam.orient = orient;
            
            % add the ramp parameter to the ode parameter struct
            odeparam.rampparam = rampparam;
            
            % ODE solver call
            % OUTPUT: time T and magnetization M
            [T,M] = ode45(@(t,m) fcn_BLOCHUS_ode(t,m,odeparam),rampparam.interp.t,Minit,options);
            
            % normalized magnetization vector at end of switch-off ramp
            MMn = M(end,:)./norm(M(end,:));
            % "adiabatic quality" p of the switch-off ramp
            % -> orientation of M with respect to z_unit
            p(tt,pp) = dot(MMn(end,:),zunit)./norm(zunit);
            % final orientation of M
            Mfinal(tt,pp,1:3) = M(end,:);
            
        end
        infostr1 = lookup_name{i1};
        infostr2 = ['theta=',num2str(theta(tt)),'°'];
        disp([infostr1,' | ',infostr2]);
    end
    % save the local lookup table for this particular combination of
    % ramp shape and ramp time
    data{i1}.ramp = lookup_name{i1};
    data{i1}.PPfac = PPfac;
    data{i1}.theta = theta;
    data{i1}.odeparam = odeparam;
    data{i1}.M = Mfinal;
    data{i1}.p = p;
end

% save lookup table
if save_data
    save('example6_lookup_table.mat','data');
end

% plot data
if plot_data
    [xx,yy] = meshgrid(linspace(0,180,Nt+1),logspace(-1,4,Np+1));
    figure;
    for i1 = 1:numel(lookup_name)
        p = data{i1}.p;
        subplot(1,numel(lookup_name),i1);
        surf(xx,yy,ones(size(xx)),p'); view(2);
        shading flat;
        set(gca,'XLim',[0 180],'XTick',linspace(0,180,5),'YScale','log',...
            'YLim',[0.1 1e4],'YTick',logspace(-1,4,6),'FontSize',12);
        xlabel('\theta [deg]');
        ylabel('Bp [B0]');
        tstr = lookup_name{i1};
        tstr = strrep(tstr,'_',' ');
        title(tstr);
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
