function onPushRun(src,~)
%onPushRun starts the calculation
%
% Syntax:
%       onPushRun(src)
%
% Inputs:
%       src - handle of the calling object
%
% Outputs:
%       none
%
% Example:
%       onPushRun(src)
%
% Other m-files required:
%       fcn_BLOCHUS_ode
%       getFFT
%       getMrot
%       getOmega0
%       getPulsePhase
%       getPulseTimeSeries
%       plotResults
%
% Subfunctions:
%       none
%
% MAT-files required:
%       none
%
% See also: BLOCHUS
% Author: Thomas Hiller
% email: thomas.hiller[at]leibniz-liag.de
% License: GNU GPLv3 (at end)

%------------- BEGIN CODE --------------

% get GUI handle
fig = ancestor(src,'figure','toplevel');

if ~isempty(fig) && strcmp(get(fig,'Tag'),'BLOCHUS')
    % get GUI data
    gui = getappdata(fig,'gui');
    data = getappdata(fig,'data');
    
    % change the pushbutton color to indicate that a calculation is running
    set(src,'BackGroundColor','r');
    % reset the color of the animate button
    set(gui.push_handles.Animate,'BackGroundColor',[.94 .94 .94]);
    
    % basic parameter that are always needed
    % z unit vector
    zunit = [0 0 1]';
    % initial magnetization [A/m]
    Minit = data.basic.Minit(:);
    % total simulation time Tsim [s]
    Tsim = data.basic.Tsim/1e3; 
    
    % parameter needed for the ODE solver
    % simulation type [string]
    odeparam.type = data.basic.type;
    % equilibrium magnetization [A/m]
    odeparam.M0 = data.basic.M0(:);
    % primary (Earth's) magnetic field B0 [T]
    odeparam.B0 = data.basic.B0;
    % longitudinal relaxation time T1 [s]
    odeparam.T1 = data.basic.T1relax/1e3;
    % transversal relaxation time T2 [s]
    odeparam.T2 = data.basic.T2relax/1e3;
    % gyromagnetic ratio [rad/s/T]
    odeparam.gamma = data.basic.gamma; 
    
    % ODE solver error tolerance
    tol = 1e-9;
    % ODE solver options
    options = odeset('RelTol',tol,'AbsTol',[tol tol tol]);
    
    % time the calculation
    t0 = tic;
    switch data.basic.type
        case 'std'
            % update Info field
            set(gui.text_handles.Status,'String','Calculation of relaxation ...');
            pause(0.01);
            
            % normalized initial magnetization
            Minit = Minit./norm(Minit);
            % update the corresponding GUI fields
            data.basic.Minit = Minit;
            set(gui.edit_handles.Minitx,'String',sprintf('%4.3f',Minit(1)));
            set(gui.edit_handles.Minity,'String',sprintf('%4.3f',Minit(2)));
            set(gui.edit_handles.Minitz,'String',sprintf('%4.3f',Minit(3)));
            
            % ODE solver call
            % OUTPUT: time T and magnetization M
            [TT,MM] = ode45(@(t,m) fcn_BLOCHUS_ode(t,m,odeparam),[0 Tsim],Minit,options);
            
            % save data
            data.results.basic.T = TT;
            data.results.basic.M = MM;
            
            % rotate M into rotating frame of reference 
            Mrot = getMrot(MM,getOmega0(odeparam.gamma,odeparam.B0).*TT);
            data.results.basic.Mrot = Mrot;
            
            % get FFTof M in the laboratory frame of reference
            [Xm,fx] = getFFT(TT,MM(:,1:2));
            data.results.basic.Mspec.fx = fx;
            data.results.basic.Mspec.X = Xm;
            
            % update Info field
            set(gui.text_handles.Status,'String','Calculation of relaxation ... finished.');
            pause(0.01);
            
        case 'prepol'
            % update Info field
            set(gui.text_handles.Status,'String','Calculation of pre-polarization switch-off ...');
            pause(0.01);
            
            % initial magnetization in the direction of B0+Bp
            Morient = data.results.prepol.orient*data.results.prepol.Bmax + odeparam.B0*zunit;
            % Minit does not get normalized in this case
            Minit = Morient;%./norm(Morient);
            % because Minit is not normalized to 1, M0 needs to be adjusted
            odeparam.M0 = odeparam.B0*zunit;
            % update the corresponding GUI fields
            set(gui.edit_handles.Minitx,'String',sprintf('%4.3f',Minit(1)/norm(Minit)));
            set(gui.edit_handles.Minity,'String',sprintf('%4.3f',Minit(2)/norm(Minit)));
            set(gui.edit_handles.Minitz,'String',sprintf('%4.3f',Minit(3)/norm(Minit)));
            
            % orientation of the pre-polarization pulse axis
            odeparam.orient = data.results.prepol.orient;
            
            % pre-polarization switch-off ramp parameter
            % relaxation during switch-off [0/1]
            odeparam.RDS = data.prepol.RDS;
            % switch-off ramp type [string]
            rampparam.ramp = data.prepol.Ramp;            
            % amplitude of pre-polarization field
            rampparam.Bmax = data.results.prepol.Bmax;
            % switch over magnetization (for linexp case)
            rampparam.Bstar = data.results.prepol.Bstar;
            % switch-off ramp time [s]
            rampparam.Tramp = data.prepol.Tramp/1e3;
            % switch over ramp time [s] (for linexp case)
            rampparam.Tslope = data.prepol.Tslope/1e3;
            
            % add the ramp parameter to the ode parameter struct
            odeparam.rampparam = rampparam;
            
            % ODE solver call
            % OUTPUT: time T and magnetization M
            [TT,MM] = ode45(@(t,m) fcn_BLOCHUS_ode(t,m,odeparam),[0 Tsim],Minit,options);
            
            % normalized magnetization vector at end of switch-off ramp
            % because the simulation can be longer than the ramp, find the
            % last point of the ramp to calculate the adiabatic quality p
            indt = find(TT<=rampparam.Tramp,1,'last');
            MMn = MM(indt,:)./norm(MM(indt,:));
            % "adiabatic quality" p of the switch-off ramp
            % -> orientation of M with respect to z_unit          
            p = dot(MMn,zunit)./norm(zunit);
            
            % save data
            data.results.basic.T = TT;
            data.results.basic.M = MM;
            data.results.prepol.p = p;
            
            % update Info field
            set(gui.text_handles.Status,'String','Calculation of pre-polarization switch-off ... finished.');
            pause(0.01);
            
            % activate the lab-frame panels to show the results
            set(gui.panels.Plot.Mag,'Selection',1);
            set(gui.panels.Plot.Sphere,'Selection',1);
            
        case 'pulse'
            % update Info field
            set(gui.text_handles.Status,'String','Calculation of excitation pulse ...');
            pause(0.01);
            
            % normalized initial magnetization
            Minit = Minit./norm(Minit);
            % update the corresponding GUI fields
            data.basic.Minit = Minit;
            set(gui.edit_handles.Minitx,'String',sprintf('%4.3f',Minit(1)));
            set(gui.edit_handles.Minity,'String',sprintf('%4.3f',Minit(2)));
            set(gui.edit_handles.Minitz,'String',sprintf('%4.3f',Minit(3)));
            
            % excitation pulse parameter
            % relaxation during pulse [0/1]
            odeparam.RDP = data.pulse.RDP;
            % pulse length [s]
            odeparam.Ttau = data.pulse.Ttau/1e3;
            % pulse type [string]
            pulseparam.PulseType = data.pulse.Type;
            % gyromagnetic ratio [rad/s/T]
            pulseparam.gamma = odeparam.gamma;
            % Larmor frequency [rad]
            pulseparam.omega0 = getOmega0(odeparam.gamma,odeparam.B0);
            % pulse amplitude [B0]
            pulseparam.Amp = odeparam.B0*data.pulse.B1Factor;
            % pulse frequency modulation [struct]
            pulseparam.fmod = data.results.pulse.fmod;
            % pulse current modulation [struct]
            pulseparam.Imod = data.results.pulse.Imod;
            % auxiliary pulse phase [rad]
            pulseparam.phi = 0;
            % pulse axis [string]
            pulseparam.PulseAxis = data.pulse.Axis;
            % pulse polarization [string]
            pulseparam.PulsePolarization = data.pulse.Polarization;
            % if the discrete MIDI pulses are used add the corresponding
            % data
            if isfield(data,'pulse_MIDI')
                % MIDI pulse data [struct]
                pulseparam.MIDI = data.pulse_MIDI;
            end
            
            % add the pulse parameter to the ode parameter struct
            odeparam.pulseparam = pulseparam;
            
            % ODE solver call
            % OUTPUT: time T and magnetization M
            switch data.pulse.Type
                case {'MIDI_OR','MIDI_AP'}
                    % if Tsim > Ttau extend the discrete time steps
                    % this is needed because for the discrete pulses
                    % specific time steps are fed into the solver
                    Tsim = unique([data.pulse_MIDI.t;(0:1/50e3:Tsim)']);
                    [TT,MM] = ode45(@(t,m) fcn_BLOCHUS_ode(t,m,odeparam),Tsim,Minit,options);
                otherwise
                    [TT,MM] = ode45(@(t,m) fcn_BLOCHUS_ode(t,m,odeparam),[0 Tsim],Minit,options);
            end
            
            % save data
            data.results.basic.T = TT;
            data.results.basic.M = MM;
            
            % in order to transform M in the rotating frame of reference
            % get the pulse phase theta for all simulated time steps
            pulseparam.fmod.t = TT;
            pulseparam.Imod.t = TT;
            pulseparam.t = TT;
            % get also the pulse amplitudes at these time steps
            [Bpulse,~,~,theta] = getPulseTimeSeries(pulseparam);
            
            % if the simulation is longer then the pulse (additional
            % relaxation afterwards), correct for the pulse phase for
            % points AFTER the pulse by simply calculating the phase at
            % the end of the pulse
            if data.basic.Tsim>data.pulse.Ttau
                % correction phase
                dphi = zeros(size(theta));
                % pulse length
                Ttau = odeparam.Ttau;
                % get pulse phase at the end of the pulse
                dphival = getPulsePhase(Ttau,0,pulseparam,1);
                % and add it to the correction phase for all time steps after
                % the pulse
                indT = TT>Ttau;
                dphi(indT) = dphival;
                % set theta for all time steps after the pulse simply to
                % omega0*t
                theta(indT) = getOmega0(odeparam.gamma,odeparam.B0).*TT(indT);
                % rotate M into rotating frame of reference
                Mrot = getMrot(MM,theta,dphi);
                % get FFT of M in the laboratory frame of reference
                [Xm,fmx] = getFFT(TT,MM(:,1:2));
                % get FFT of the pulse
                [Xb,fbx] = getFFT(TT(TT<=Ttau),Bpulse(TT<=Ttau,1:2));
                
                % update Info field
                set(gui.text_handles.Status,'String','Calculation of excitation pulse & relaxation ... finished.');
                pause(0.01);
            else
                % rotate M into rotating frame of reference
                Mrot = getMrot(MM,theta);
                % get FFT of M in the laboratory frame of reference
                [Xm,fmx] = getFFT(TT,MM(:,1:2));
                % get FFT of the pulse
                [Xb,fbx] = getFFT(TT,Bpulse(:,1:2));
                
                % update Info field
                set(gui.text_handles.Status,'String','Calculation of excitation pulse ... finished.');
                pause(0.01);
            end
            
            % save data
            data.results.basic.Mrot = Mrot;
            data.results.basic.Mspec.fx = fmx;
            data.results.basic.Mspec.X = Xm;
            data.results.pulse.Bspec.fx = fbx;
            data.results.pulse.Bspec.X = Xb;
            
        case 'prepolpulse'
            
            % this is basically the combination of all of the above
            % simulation types
			%  --- IMPORTANT NOTE: ----------------------------------------
			% The pre-polarization switch-off only "lives" in the laboratory
			% frame of reference. However, for convenience reasons I plot
			% the lab-frame pre-polarization data also in the rotating
			% frame of reference! This is u-n-p-h-y-s-i-c-a-l and pure "eye candy"
			%  -------------------------------------------------------------
            
            % these times we will later on, so init them already here
            % total simulation time [s]
            Tsim = data.basic.Tsim/1e3;
            % switch-off ramp time [s]
            Tramp = data.prepol.Tramp/1e3;
            % wait time between switch-off ramp and pulse [s]
            Twait = data.pulse.Twait/1e3;
            % pulse length [s]
            Ttau = data.pulse.Ttau/1e3;
            
            % check if the simulation time is long enough to process all
            % different stages, if not fix it
            % NOTE: of course the simulation time can be longer (relaxation)
            % after the pulse
            if Tsim-(Tramp+Twait+Ttau) < 0
                Tsim = Tramp+Twait+Ttau;
                % update the corresponding GUI data and field
                data.basic.Tsim = Tsim*1e3; % [ms]
                set(gui.edit_handles.Tsim,'String',num2str(data.basic.Tsim));
                % no relaxation after excitation pulse
                Trelax = 0;
            else
                % relaxation after excitation pulse [s]
                Trelax = Tsim-(Tramp+Twait+Ttau);
            end
            
            % -------------------------------------------------------------
            %% --- 1.) pre-polarization switch-off ------------------------
            % -------------------------------------------------------------
            % update Info field
            set(gui.text_handles.Status,'String','Calculation of pre-polarization switch-off ...');
            pause(0.01);
            
            % initial magnetization in the direction of B0+Bp
            Morient = data.results.prepol.orient*data.results.prepol.Bmax + odeparam.B0*zunit;
             % Minit does not get normalized in this case
            Minit = Morient;%./norm(Morient);
            % because Minit is not normalized to 1, M0 needs to be adjusted
            odeparam.M0 = odeparam.B0*zunit;
            % update the corresponding GUI fields
            set(gui.edit_handles.Minitx,'String',sprintf('%4.3f',Minit(1)/norm(Minit)));
            set(gui.edit_handles.Minity,'String',sprintf('%4.3f',Minit(2)/norm(Minit)));
            set(gui.edit_handles.Minitz,'String',sprintf('%4.3f',Minit(3)/norm(Minit)));

            % update the simulation type
            odeparam.type = 'prepol';
            
            % orientation of the pre-polarization pulse axis
            odeparam.orient = data.results.prepol.orient;
            
            % pre-polarization switch-off ramp parameter
            % relaxation during switch-off [0/1]
            odeparam.RDS = data.prepol.RDS;
            % switch-off ramp type [string]
            rampparam.ramp = data.prepol.Ramp;            
            % amplitude of pre-polarization field
            rampparam.Bmax = data.results.prepol.Bmax;
            % switch over magnetization (for linexp case)
            rampparam.Bstar = data.results.prepol.Bstar;
            % switch-off ramp time [s]
            rampparam.Tramp = Tramp;
            % switch over ramp time [s] (for linexp case)
            rampparam.Tslope = data.prepol.Tslope/1e3;
            
            % add the ramp parameter to the ode parameter struct
            odeparam.rampparam = rampparam;
            
            % ODE solver call
            % OUTPUT: time T and magnetization M
            [TT,MM] = ode45(@(t,m) fcn_BLOCHUS_ode(t,m,odeparam),[0 Tramp],Minit,options);
            
            % data of the first stage
            TT1 = TT;
            MM1 = MM;

            % normalized magnetization vector at end of switch-off ramp
            MMn = MM(end,:)./norm(MM(end,:));
            % "adiabatic quality" p of the switch-off ramp
            % -> orientation of M with respect to z_unit
            data.results.prepol.p = dot(MMn(end,:),zunit)./norm(zunit);
            
            % -------------------------------------------------------------
            %% --- 2.) wait time (if any) ---------------------------------
            % -------------------------------------------------------------
            if Twait > 0
                % update Info field
                set(gui.text_handles.Status,'String','Calculation of wait time ...');
                pause(0.01);
                % initial magnetization is the end of the switch-off ramp
                Minit2 = MM1(end,:)';
                % update the simulation type
                odeparam.type = 'std';
                
                % ODE solver call
                % OUTPUT: time T and magnetization M
                [TT,MM] = ode45(@(t,m) fcn_BLOCHUS_ode(t,m,odeparam),[0 Twait],Minit2,options);
                
                % rotate M into rotating frame of reference 
                MM2rot = getMrot(MM,getOmega0(odeparam.gamma,odeparam.B0).*TT);
                % get the additional phase at the end of the wait time
                phiWait = getOmega0(odeparam.gamma,odeparam.B0).*TT(end);
                
                % because the simulation was done in "local" time
                % coordinates, shift the time vector to the end of the
                % switch-off ramp
                TT2 = TT + Tramp;
                % save data
                MM2 = MM;
                % initial magnetization for the pulse is the end of the
                % wait time period
                Minit3 = MM2(end,:)';
            else
                % if wait time is 0, use dummy values
                TT2 = [];
                MM2 = [];
                MM2rot = [];
                phiWait = 0;
                % initial magnetization for the pulse is the end of the
                % switch-off ramp
                Minit3 = MM1(end,:)';
            end
            
            % -------------------------------------------------------------
            %% --- 3.) excitation pulse -----------------------------------
            % -------------------------------------------------------------
            
            % update the simulation type
            odeparam.type = 'pulse';

            % excitation pulse parameter
            % relaxation during pulse [0/1]
            odeparam.RDP = data.pulse.RDP;
            % pulse length [s]
            odeparam.Ttau = Ttau;
            
            % pulse type [string]
            pulseparam.PulseType = data.pulse.Type;
            % gyromagnetic ratio [rad/s/T]
            pulseparam.gamma = odeparam.gamma;
            % Larmor frequency [rad]
            pulseparam.omega0 = getOmega0(odeparam.gamma,odeparam.B0);
            % pulse amplitude [B0]
            pulseparam.Amp = odeparam.B0*data.pulse.B1Factor;
            % pulse frequency modulation [struct]
            pulseparam.fmod = data.results.pulse.fmod;
            % pulse current modulation [struct]
            pulseparam.Imod = data.results.pulse.Imod;
            % auxiliary pulse phase [rad]
            pulseparam.phi = 0;
            % pulse axis [string]
            pulseparam.PulseAxis = data.pulse.Axis;
            % pulse polarization [string]
            pulseparam.PulsePolarization = data.pulse.Polarization;
            % if the discrete MIDI pulses are used add the corresponding
            % data
            if isfield(data,'pulse_MIDI')
                % MIDI pulse data [struct]
                pulseparam.MIDI = data.pulse_MIDI;
            end

            % add the pulse parameter to the ode parameter struct
            odeparam.pulseparam = pulseparam;
            
            % update Info field
            set(gui.text_handles.Status,'String','Calculation of excitation pulse ...');
            pause(0.01);
            
            % if there is a time span after the pulse
            if abs(Trelax) > eps
                % the "local" simulation time is prolonged
                Tsim1 = Ttau+Trelax;
            else
                % otherwise the "local" simulation time is simply the pulse
                % length
                Tsim1 = Ttau;
            end
            
            % ODE solver call
            % OUTPUT: time T and magnetization M
            switch data.pulse.Type
                case {'MIDI_OR','MIDI_AP'}
                    % if Tsim1 > Ttau extend the discrete time steps
                    % this is needed because for the discrete pulses
                    % specific time steps are fed into the solver
                    Tsim1 = unique([data.pulse_MIDI.t;(0:1/50e3:Tsim1)']);
                    [TT,MM] = ode45(@(t,m) fcn_BLOCHUS_ode(t,m,odeparam),Tsim1,Minit3,options);
                otherwise
                    [TT,MM] = ode45(@(t,m) fcn_BLOCHUS_ode(t,m,odeparam),[0 Tsim1],Minit3,options);
            end
            
            % in order to transform M in the rotating frame of reference
            % get the pulse phase theta for all simulated time steps
            pulseparam.fmod.t = TT;
            pulseparam.Imod.t = TT;
            pulseparam.t = TT;
            % get also the pulse amplitudes at these time steps
            [Bpulse,~,~,theta] = getPulseTimeSeries(pulseparam);
            
            % if the simulation is longer then the pulse (additional
            % relaxation afterwards), correct for the pulse phase for
            % points AFTER the pulse by simply calculating the phase at
            % the end of the pulse
            if abs(Trelax) > eps%max(Tsim1)>Ttau
                % correction phase
                dphi = zeros(size(theta));
                % get pulse phase at the end of the pulse
                dphival = getPulsePhase(Ttau,0,pulseparam,1);
                % and add it to the correction phase for all time steps after
                % the pulse
                indT = TT>Ttau;
                dphi(indT) = dphival;
                % set theta for all time steps after the pulse simply to
                % omega0*t
                theta(indT) = getOmega0(odeparam.gamma,odeparam.B0).*TT(indT);
                % rotate M into rotating frame of reference
                % and account for the phase from the wait time before the
                % pulse
                MM3rot = getMrot(MM,theta,dphi+phiWait);
                % get FFT of the pulse
                [Xb,fbx] = getFFT(TT(TT<=Ttau),Bpulse(TT<=Ttau,1:2));
                
                % update Info field
                set(gui.text_handles.Status,'String',['Calculation of excitation pulse',...
                    ' & relaxation ... finished.']);
                pause(0.01);
            else
                % rotate M into rotating frame of reference
                % and account for the phase from the wait time before the
                % pulse
                MM3rot = getMrot(MM,theta,phiWait);
                % get FFT of the pulse
                [Xb,fbx] = getFFT(TT,Bpulse(:,1:2));
                
                % update Info field
                set(gui.text_handles.Status,'String','Calculation of excitation pulse ... finished.');
                pause(0.01);
            end
            
            % -------------------------------------------------------------
            %% --- 4.) final data merging ---------------------------------
            % -------------------------------------------------------------
            % because the simulation was done in "local" time
            % coordinates, shift the time vector to the end of the
            % wait time after the switch-off            
            TT3 = TT + Tramp + Twait;
            % save data
            MM3 = MM;
            
            % combine data from all stages
            TT = [TT1;TT2;TT3];
            MM = [MM1;MM2;MM3];
            % because there is no M in the rotating frame of reference
            % during the switch-off ramp, use the lab-frame data instead
			%  --- IMPORTANT NOTE: ----------------------------------------
			% The pre-polarization switch-off only "lives" in the laboratory
			% frame of reference. However, for convenience reasons I plot
			% the lab-frame pre-polarization data also in the rotating
			% frame of reference! This is u-n-p-h-y-s-i-c-a-l and pure "eye candy"
			% MM1 is lab-frame data!
			%  -------------------------------------------------------------
            MMrot = [MM1;MM2rot;MM3rot];
            
            % remove possible duplicate points
            [TT,ix] = unique(TT);
            MM = MM(ix,:);
            MMrot = MMrot(ix,:);
            
            % save combined data
            data.results.basic.T = TT;
            data.results.basic.M = MM;
            data.results.basic.Mrot = MMrot;
            
            % get FFT for the combined M
            [Xm,fmx] = getFFT(TT,MM(:,1:2));
            data.results.basic.Mspec.fx = fmx;
            data.results.basic.Mspec.X = Xm;
            % save pulse FFT
            data.results.pulse.Bspec.fx = fbx;
            data.results.pulse.Bspec.X = Xb;

    end
    % time the calculation
    data.info.Timer = toc(t0);
    
    % update GUI data
    setappdata(fig,'data',data);
    % plot results
    plotResults(fig);
    % update status bar
    updateStatusInformation(fig);
    % activate animation button
    set(gui.push_handles.Animate,'Enable','on');
    
    % change the pushbutton color back to green
    set(src,'BackGroundColor','g');
    
else
    warndlg({'onPushRun:','There is no figure with the BLOCHUS Tag open.'},...
        'BLOCHUS error');
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
