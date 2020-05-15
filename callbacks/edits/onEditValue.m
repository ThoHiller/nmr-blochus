function onEditValue(src,~)
%onEditValue updates all edit field values, checks for wrong inputs and
%restores a default value if necessary
%
% Syntax:
%       onEditValue(src)
%
% Inputs:
%       src - handle of the calling object
%
% Outputs:
%       none
%
% Example:
%       onEditValue(src)
%
% Other m-files required:
%       getB0
%       getOmega0
%       getPulseParameters
%       getRampParameters
%       plotBpulse
%       plotRamp
%
% Subfunctions:
%       createDataString
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
    
    % get the value of the field
    val = str2double(get(src,'String'));
    % get the tag of the edit field
    tag = get(src,'Tag');
    % get the user data of the field
    ud = get(src,'UserData');
    
    if isstruct(ud) % for new fields
        % get the default values [default min max]
        defaults = ud.defaults;
        
        % check if the value is numeric
        % if not reset to defaults stored in user data
        if isnan(val)
            set(src,'String',num2str(defaults(1)));
            val = str2double(get(src,'String'));
        end
        % check if the value is out of bounds
        % if yes reset to default
        if val < defaults(2) || val > defaults(3)
            set(src,'String',num2str(defaults(1)));
            val = str2double(get(src,'String')); %#ok<*NASGU>
        end
        
    else % old style (Needs to be removed after the refactoring)
        % check if the value is numeric
        % if not reset to defaults stored in user data
        if isnan(val)
            set(src,'String',num2str(ud(3)));
            val = str2double(get(src,'String'));
        end
        % check if the value is out of bounds
        % if yes reset to default
        if val < ud(1) || val > ud(2)
            set(src,'String',num2str(ud(3)));
            val = str2double(get(src,'String')); %#ok<*NASGU>
        end
        disp(['there is still an old style edit field:',tag])
    end
    
    % get the data field to update from the field tag
    out = createDataString(tag);
    % update the corresponding data field
    updstr = [out.updstr,'=val;'];
    eval(updstr);
    
    % update the data inside the GUI
    setappdata(fig,'data',data);
    
    % depending on the particular edit field, further actions are
    % necessary
    switch tag
        % -----------------------------------------------------------------
        % --- BASICS ------------------------------------------------------
        % -----------------------------------------------------------------
        case 'basic_B0'
            % new Larmor frequency due to changed B0
            data.basic.Omega0 = getOmega0(data.basic.gamma,data.basic.B0)/2/pi;
            % update the GUI fields
            set(gui.edit_handles.Omega0,'String',sprintf('%7.2f',data.basic.Omega0));
            setappdata(fig,'data',data);
            % if "Pulse" is activated, update the pulse as if the
            % pulse length (or number of periods) would have been changed
            % manually
            if get(gui.check_handles.Pulse,'Value') == 1
                onEditValue(gui.edit_handles.PulseTtau);
                data = getappdata(fig,'data');
            end
            % if "PrePol" is activated and the active ramp is "melton1995"
            % update the switch-off ramp settings due to the changed Larmor
            % frequency
            if get(gui.check_handles.PrePol,'Value') == 1 && ...
                    get(gui.popup_handles.PrePolRamp,'Value') == 5
                setappdata(fig,'data',data);
                onPopupPrePolRamp(gui.popup_handles.PrePolRamp);
                data = getappdata(fig,'data');
            end
            if get(gui.check_handles.PrePol,'Value') == 1
                setappdata(fig,'data',data);
                getRampParameters(fig);
                data = getappdata(fig,'data');
            end
            
        case 'basic_Omega0'
            % fix sign if it is wrong
            if (data.basic.gamma > 0 && data.basic.Omega0 > 0) || ...
                    (data.basic.gamma < 0 && data.basic.Omega0 < 0)
                data.basic.Omega0 = -data.basic.Omega0;
            end
            % new B0 due to changed Larmor frequency
            data.basic.B0 = getB0(data.basic.gamma,data.basic.Omega0*2*pi);
            % update the GUI fields
            set(gui.edit_handles.Omega0,'String',sprintf('%7.2f',data.basic.Omega0));
            set(gui.edit_handles.B0,'String',sprintf('%5.4e',data.basic.B0));
            setappdata(fig,'data',data);
            % if "Pulse" is activated, update the pulse as if the
            % pulse length (or number of periods) would have been changed
            % manually
            if get(gui.check_handles.Pulse,'Value') == 1
                onEditValue(gui.edit_handles.PulseTtau);
                data = getappdata(fig,'data');
            end
            % if "PrePol" is activated and the active ramp is "melton1995"
            % update the switch-off ramp settings due to the changed Larmor
            % frequency
            if get(gui.check_handles.PrePol,'Value') == 1 && ...
                    get(gui.popup_handles.PrePolRamp,'Value') == 5
                setappdata(fig,'data',data);
                onPopupPrePolRamp(gui.popup_handles.PrePolRamp);
                data = getappdata(fig,'data');
            end
            if get(gui.check_handles.PrePol,'Value') == 1
                setappdata(fig,'data',data);
                getRampParameters(fig);
                data = getappdata(fig,'data');
            end
            
        case 'basic_Minitx'
            % x-component of init magnetization
            data.basic.Minit(1) = val;
            
        case 'basic_Minity'
            % y-component of init magnetization
            data.basic.Minit(2) = val;
            
        case 'basic_Minitz'
            % z-component of init magnetization
            data.basic.Minit(3) = val;
            
            % -------------------------------------------------------------
            % --- PRE-POLARIZATION ----------------------------------------
            % -------------------------------------------------------------
        case {'prepol_Factor','prepol_SwitchFactor'}
            % check if the ramp type is "melton1995", if yes update the
            % switch-off time to maintain the same switch-off rate
            if get(gui.popup_handles.PrePolRamp,'Value') == 5
                % k / rate
                GAMMA = data.prepol.Factor/data.prepol.SwitchFactor;
                % corresponding ramp time
                Tramp = GAMMA/(data.basic.B0*data.basic.gamma)*1e3; % [ms]
                data.prepol.Tramp = Tramp;
                data.prepol.Tslope = Tramp;
                switch data.basic.type
                    case 'prepol'
                        data.basic.Tsim = Tramp;
                        set(gui.edit_handles.Tsim,'String',num2str(data.basic.Tsim));
                    otherwise
                        % nothing to do
                end
                setappdata(fig,'data',data);
                set(gui.edit_handles.PrePolTramp,'String',num2str(data.prepol.Tramp));
                set(gui.edit_handles.PrePolTslope,'String',num2str(data.prepol.Tslope));
            end
            % update switch-off ramp settings
            getRampParameters(fig);
            data = getappdata(fig,'data');
            plotRamp(fig);
            
        case {'prepol_Theta','prepol_Phi','prepol_Tslope'}
            % update switch-off ramp settings
            getRampParameters(fig);
            data = getappdata(fig,'data');
            plotRamp(fig);
            
        case 'prepol_Tramp'
            % for "linear" and "half cosine" ramp set the Tslope value,
            % which is not need here equal to Tramp
            if strcmp(data.prepol.Ramp,'lin') || strcmp(data.prepol.Ramp,'halfcos')
                data.prepol.Tslope = data.prepol.Tramp;
                set(gui.edit_handles.PrePolTslope,'String',num2str(data.prepol.Tslope));
            end
            % in case of pure "PrePol" update the total simulation time
            switch data.basic.type
                case 'prepol'
                    data.basic.Tsim = data.prepol.Tramp;
                    set(gui.edit_handles.Tsim,'String',num2str(data.basic.Tsim));
                otherwise
                    % nothing to do
            end
            setappdata(fig,'data',data);
            % update switch-off ramp settings
            getRampParameters(fig);
            data = getappdata(fig,'data');
            plotRamp(fig);
            
            % -------------------------------------------------------------
            % --- PULSE ---------------------------------------------------
            % -------------------------------------------------------------
        case 'pulse_B1Factor'
            % if the pulse amplitude changes some settings have to be
            % adjusted, depending on the chosen pulse type
            switch data.pulse.Type
                case {'pi_half','pi'}
                    switch data.pulse.Type
                        case 'pi_half'
                            % in case of "pi_half" the pulse length gets updated
                            data.pulse.Ttau = abs(1e3*(pi/2/(data.basic.gamma*data.basic.B0*data.pulse.B1Factor)));
                        case 'pi'
                            % in case of "pi" the pulse length gets updated
                            data.pulse.Ttau = abs(1e3*(pi/(data.basic.gamma*data.basic.B0*data.pulse.B1Factor)));
                    end
                    set(gui.edit_handles.PulseTtau,'String',num2str(data.pulse.Ttau));
                    % in case of pure "Pulse" update the total simulation time
                    switch data.basic.type
                        case 'pulse'
                            data.basic.Tsim = data.pulse.Ttau;
                            set(gui.edit_handles.Tsim,'String',num2str(data.basic.Tsim));
                        otherwise
                            % nothing to do
                    end
                    setappdata(fig,'data',data);
                    % update pulse settings
                    getPulseParameters(fig);
                    data = getappdata(fig,'data');
                    plotBpulse(fig);
                    
                otherwise
                    % update pulse settings
                    getPulseParameters(fig);
                    data = getappdata(fig,'data');
                    plotBpulse(fig);
            end
            
        case 'pulse_Ttau'
            % if the pulse length changes some settings have to be
            % adjusted, depending on the chosen pulse type
            switch data.pulse.Type
                case {'pi_half','pi'}
                    switch data.pulse.Type
                        case 'pi_half'
                            % in case of "pi_half" the pulse amplitude gets updated
                            data.pulse.B1Factor = abs(1e3*(pi/2/(data.basic.gamma*data.basic.B0*data.pulse.Ttau)));
                        case 'pi'
                            % in case of "pi" the pulse amplitude gets updated
                            data.pulse.B1Factor = abs(1e3*(pi/(data.basic.gamma*data.basic.B0*data.pulse.Ttau)));
                    end
                    set(gui.edit_handles.PulseB1Factor,'String',num2str(data.pulse.B1Factor));
                    % in case of pure "Pulse" update the total simulation time
                    switch data.basic.type
                        case 'pulse'
                            data.basic.Tsim = data.pulse.Ttau;
                            set(gui.edit_handles.Tsim,'String',num2str(data.basic.Tsim));
                        otherwise
                            % nothing to do
                    end
                    setappdata(fig,'data',data);
                    % update pulse settings
                    getPulseParameters(fig);
                    data = getappdata(fig,'data');
                    plotBpulse(fig);
                    
                case 'MIDI_OR'
                    % in case of "MIDI_OR" the PulseTtau-field holds the
                    % number of periods
                    data.pulse.MIDINP = val;
                    % update pulse length
                    data.pulse.Ttau = abs(data.pulse.MIDINP/data.basic.Omega0)*1e3;
                    % update pulse amplitude to pi/2 value
                    data.pulse.B1Factor = abs(1e3*(pi/2/(data.basic.gamma*data.basic.B0*data.pulse.Ttau)));
                    set(gui.edit_handles.PulseB1Factor,'String',num2str(data.pulse.B1Factor));
                    % in case of pure "Pulse" update the total simulation time
                    switch data.basic.type
                        case 'pulse'
                            data.basic.Tsim = data.pulse.Ttau;
                            set(gui.edit_handles.Tsim,'String',num2str(data.basic.Tsim));
                        otherwise
                            % nothing to do
                    end
                    setappdata(fig,'data',data);
                    % update pulse settings
                    getPulseParameters(fig);
                    data = getappdata(fig,'data');
                    plotBpulse(fig);
                    
                case 'MIDI_AP'
                    % in case of "MIDI_AP" the PulseTtau-field holds the
                    % number of periods
                    data.pulse.MIDINP = val;
                    % update pulse length assuming a constant frequency
                    data.pulse.Ttau = abs(data.pulse.MIDINP/data.basic.Omega0)*1e3;
                    % now calculate the actual pulse settings
                    setappdata(fig,'data',data);
                    getPulseParameters(fig);
                    data = getappdata(fig,'data');
                    % adjust the pulse length to the actual value
                    data.pulse.Ttau = data.pulse_MIDI.t(end)*1e3;
                    % in case of pure "Pulse" update the total simulation time
                    switch data.basic.type
                        case 'pulse'
                            data.basic.Tsim = data.pulse.Ttau;
                            set(gui.edit_handles.Tsim,'String',num2str(data.basic.Tsim));
                        otherwise
                            % nothing to do
                    end
                    setappdata(fig,'data',data);
                    % update pulse plot
                    plotBpulse(fig);
                    
                otherwise
                    % for all other pulse types ("free" & "AHP")
                    % in case of pure "Pulse" update the total simulation time
                    switch data.basic.type
                        case 'pulse'
                            data.basic.Tsim = data.pulse.Ttau;
                            set(gui.edit_handles.Tsim,'String',num2str(data.basic.Tsim));
                        otherwise
                            % nothing to do
                    end
                    % update GUI data
                    setappdata(fig,'data',data);
                    % update pulse settings
                    getPulseParameters(fig);
                    data = getappdata(fig,'data');
                    plotBpulse(fig);
            end
            
        otherwise
            switch data.basic.type
                % update all changes related to the pulse modulation
                case {'pulse','prepolpulse'}
                    setappdata(fig,'data',data);
                    getPulseParameters(fig);
                    data = getappdata(fig,'data');
                    plotBpulse(fig);
                    
                otherwise
                    % nothing to do for all other fields e.g. T1, T2, Tsim
            end
    end
    
    % update GUI data
    setappdata(fig,'data',data);
else
    warndlg({'onEditValue:','There is no figure with the BLOCHUS Tag open.'},...
        'BLOCHUS error');
end

end

%% helper function to create the update string
function out = createDataString(tag)
% find the underscore
ind = strfind(tag,'_');
% the panel name is before the underscore
out.panel = tag(1:ind(1)-1);
% the field name afterwards
out.field = tag(ind(1)+1:end);
% replace the underscore with a dot
tag(ind) = '.';
% create the update string
out.updstr = ['data.',tag];

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
