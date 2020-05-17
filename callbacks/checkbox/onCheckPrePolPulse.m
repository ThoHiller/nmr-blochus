function onCheckPrePolPulse(src,~)
%onCheckPrePolPulse activates / deactivates all control elements needed for
%either the pre-polarization switch-off ramp or a B1-pulse or both
%
% Syntax:
%       onCheckPrePolPulse(src)
%
% Inputs:
%       src - handle of the calling object
%
% Outputs:
%       none
%
% Example:
%       onCheckPrePolPulse(src)
%
% Other m-files required:
%       onPopupPrePolRamp
%       onPopupPulseType
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
    gui  = getappdata(fig,'gui');
    data = getappdata(fig,'data');
    
    % check status of the check boxes
    isPrePol = get(gui.check_handles.PrePol,'Value');
    isPulse = get(gui.check_handles.Pulse,'Value');
    
    if isPrePol == 0 && isPulse == 0 % just standard relaxation
        data.basic.type = 'std';
        
        % menu settings
        set(gui.menu.view_figures_ramp,'Enable','off');
        set(gui.menu.view_figures_pulse,'Enable','off');
        
        % std settings
        data.basic.Minit = [1 0 0];
        set(gui.edit_handles.Minitx,'String',num2str(data.basic.Minit(1)));
        set(gui.edit_handles.Minity,'String',num2str(data.basic.Minit(2)));
        set(gui.edit_handles.Minitz,'String',num2str(data.basic.Minit(3)));
        set(gui.edit_handles.Minitx,'Enable','on');
        set(gui.edit_handles.Minity,'Enable','on');
        set(gui.edit_handles.Minitz,'Enable','on');
        setappdata(fig,'data',data);
        
        % clear magnetization FFT axis
        clearSingleAxis(gui.axes_handles.MagFFT);
        
        % deactivate PrePol settings
        set(gui.check_handles.PrePolRDS,'Enable','off');
        set(gui.popup_handles.PrePolRamp,'Enable','off');
        set(gui.edit_handles.PrePolFactor,'Enable','off');
        set(gui.edit_handles.PrePolTheta,'Enable','off');
        set(gui.edit_handles.PrePolPhi,'Enable','off');
        set(gui.edit_handles.PrePolSwitchFactor,'Enable','off');
        set(gui.edit_handles.PrePolTramp,'Enable','off');
        set(gui.edit_handles.PrePolTslope,'Enable','off');
        
        % clear PrePol axes
        clearSingleAxis(gui.axes_handles.Bpre);
        clearSingleAxis(gui.axes_handles.alpha);
        clearSingleAxis(gui.axes_handles.dadt);
        clearSingleAxis(gui.axes_handles.wda);
        
        % deactivate Pulse settings
        set(gui.check_handles.PulseRDP,'Enable','off');
        set(gui.popup_handles.PulseType,'Enable','off');
        set(gui.popup_handles.PulseAxis,'Enable','off');
        set(gui.popup_handles.PulsePolarization,'Enable','off');        
        set(gui.edit_handles.PulseB1Factor,'Enable','off');
        set(gui.edit_handles.PulseTtau,'Enable','off');
        set(gui.edit_handles.PulseTwait,'Enable','off');
        set(gui.popup_handles.PulseDFmode,'Enable','off');
        set(gui.edit_handles.PulseDFstart,'Enable','off');
        set(gui.edit_handles.PulseDFend,'Enable','off');
        set(gui.edit_handles.PulseDFA,'Enable','off');
        set(gui.edit_handles.PulseDFB,'Enable','off');
        set(gui.popup_handles.PulseImode,'Enable','off');
        set(gui.edit_handles.PulseIstart,'Enable','off');
        set(gui.edit_handles.PulseIend,'Enable','off');
        set(gui.edit_handles.PulseIA,'Enable','off');
        set(gui.edit_handles.PulseIB,'Enable','off');
        set(gui.check_handles.PulseQ,'Enable','off');
        set(gui.edit_handles.PulseQ,'Enable','off');
        set(gui.edit_handles.PulseQdf,'Enable','off');
        
        % clear Pulse axes
        clearSingleAxis(gui.axes_handles.PulseB);
        clearSingleAxis(gui.axes_handles.PulseSetupF);
        clearSingleAxis(gui.axes_handles.PulseSetupI);
        clearSingleAxis(gui.axes_handles.PulseFFT);
        
    elseif isPrePol == 1 && isPulse == 0 % only pre-polarization switch-off
        data.basic.type = 'prepol';
        
        % menu settings
        set(gui.menu.view_figures_ramp,'Enable','on');
        set(gui.menu.view_figures_pulse,'Enable','off');
        
        % std settings
        data.basic.Tsim = data.prepol.Tramp;
        set(gui.edit_handles.Tsim,'String',num2str(data.basic.Tsim));
        data.basic.Minit = [1 0 0];
        set(gui.edit_handles.Minitx,'String',num2str(data.basic.Minit(1)));
        set(gui.edit_handles.Minity,'String',num2str(data.basic.Minit(2)));
        set(gui.edit_handles.Minitz,'String',num2str(data.basic.Minit(3)));
        set(gui.edit_handles.Minitx,'Enable','off');
        set(gui.edit_handles.Minity,'Enable','off');
        set(gui.edit_handles.Minitz,'Enable','off');
        
        % clear magnetization FFT axis
        clearSingleAxis(gui.axes_handles.MagFFT);
        
        % activate PrePol settings
        set(gui.check_handles.PrePolRDS,'Enable','on');
        set(gui.popup_handles.PrePolRamp,'Enable','on');
        set(gui.edit_handles.PrePolFactor,'Enable','on');
        set(gui.edit_handles.PrePolTheta,'Enable','on');
        set(gui.edit_handles.PrePolPhi,'Enable','on');
        set(gui.edit_handles.PrePolSwitchFactor,'Enable','on');
        set(gui.edit_handles.PrePolTramp,'Enable','on');
        set(gui.edit_handles.PrePolTslope,'Enable','on');
        setappdata(fig,'data',data);
        onPopupPrePolRamp(gui.popup_handles.PrePolRamp);
        data = getappdata(fig,'data');
        
        % deactivate Pulse settings
        set(gui.check_handles.PulseRDP,'Enable','off');
        set(gui.popup_handles.PulseType,'Enable','off');
        set(gui.popup_handles.PulseAxis,'Enable','off');
        set(gui.popup_handles.PulsePolarization,'Enable','off');
        set(gui.edit_handles.PulseB1Factor,'Enable','off');
        set(gui.edit_handles.PulseTtau,'Enable','off');
        set(gui.edit_handles.PulseTwait,'Enable','off');
        set(gui.popup_handles.PulseDFmode,'Enable','off');
        set(gui.edit_handles.PulseDFstart,'Enable','off');
        set(gui.edit_handles.PulseDFend,'Enable','off');
        set(gui.edit_handles.PulseDFA,'Enable','off');
        set(gui.edit_handles.PulseDFB,'Enable','off');
        set(gui.popup_handles.PulseImode,'Enable','off');
        set(gui.edit_handles.PulseIstart,'Enable','off');
        set(gui.edit_handles.PulseIend,'Enable','off');
        set(gui.edit_handles.PulseIA,'Enable','off');
        set(gui.edit_handles.PulseIB,'Enable','off');
        set(gui.check_handles.PulseQ,'Enable','off');
        set(gui.edit_handles.PulseQ,'Enable','off');
        set(gui.edit_handles.PulseQdf,'Enable','off');
        
        % clear Pulse axes
        clearSingleAxis(gui.axes_handles.PulseB);
        clearSingleAxis(gui.axes_handles.PulseSetupF);
        clearSingleAxis(gui.axes_handles.PulseSetupI);
        clearSingleAxis(gui.axes_handles.PulseFFT);
        
    elseif isPrePol == 0 && isPulse == 1 % only Pulse
        data.basic.type = 'pulse';
        
        % menu settings
        set(gui.menu.view_figures_ramp,'Enable','off');
        set(gui.menu.view_figures_pulse,'Enable','on');
        
        % std settings
        data.basic.Tsim = data.pulse.Ttau;
        set(gui.edit_handles.Tsim,'String',num2str(data.basic.Tsim));
        data.basic.Minit = [0 0 1];
        set(gui.edit_handles.Minitx,'String',num2str(data.basic.Minit(1)));
        set(gui.edit_handles.Minity,'String',num2str(data.basic.Minit(2)));
        set(gui.edit_handles.Minitz,'String',num2str(data.basic.Minit(3)));
        set(gui.edit_handles.Minitx,'Enable','on');
        set(gui.edit_handles.Minity,'Enable','on');
        set(gui.edit_handles.Minitz,'Enable','on');
        setappdata(fig,'data',data);
        
        % clear magnetization FFT axis
        clearSingleAxis(gui.axes_handles.MagFFT);
        
        % deactivate PrePol settings
        set(gui.check_handles.PrePolRDS,'Enable','off');
        set(gui.popup_handles.PrePolRamp,'Enable','off');
        set(gui.edit_handles.PrePolFactor,'Enable','off');
        set(gui.edit_handles.PrePolTheta,'Enable','off');
        set(gui.edit_handles.PrePolPhi,'Enable','off');
        set(gui.edit_handles.PrePolSwitchFactor,'Enable','off');
        set(gui.edit_handles.PrePolTramp,'Enable','off');
        set(gui.edit_handles.PrePolTslope,'Enable','off');
        
        % clear PrePol axes
        clearSingleAxis(gui.axes_handles.Bpre);
        clearSingleAxis(gui.axes_handles.alpha);
        clearSingleAxis(gui.axes_handles.dadt);
        clearSingleAxis(gui.axes_handles.wda);
        
        % activate Pulse settings
        set(gui.check_handles.PulseRDP,'Enable','on');
        set(gui.popup_handles.PulseType,'Enable','on');
        set(gui.popup_handles.PulseAxis,'Enable','on');
        set(gui.popup_handles.PulsePolarization,'Enable','on');        
        set(gui.edit_handles.PulseB1Factor,'Enable','on');
        set(gui.edit_handles.PulseTtau,'Enable','on');
        set(gui.edit_handles.PulseDFstart,'Enable','on');
        set(gui.edit_handles.PulseTwait,'Enable','off');
        onPopupPulseType(gui.popup_handles.PulseType);
        data = getappdata(fig,'data');
        
    elseif isPrePol == 1 && isPulse == 1 % pre-polarization + pulse
        data.basic.type = 'prepolpulse';
        
        % menu settings
        set(gui.menu.view_figures_ramp,'Enable','on');
        set(gui.menu.view_figures_pulse,'Enable','on');
        
        % std settings
        data.basic.Minit = [1 0 0];
        set(gui.edit_handles.Minitx,'String',num2str(data.basic.Minit(1)));
        set(gui.edit_handles.Minity,'String',num2str(data.basic.Minit(2)));
        set(gui.edit_handles.Minitz,'String',num2str(data.basic.Minit(3)));
        set(gui.edit_handles.Minitx,'Enable','off');
        set(gui.edit_handles.Minity,'Enable','off');
        set(gui.edit_handles.Minitz,'Enable','off');
        setappdata(fig,'data',data);
        
        % clear magnetization FFT axis
        clearSingleAxis(gui.axes_handles.MagFFT);
        clearSingleAxis(gui.axes_handles.MagR);
        clearSingleAxis(gui.axes_handles.SphereR);
        
        % activate PrePol settings
        set(gui.check_handles.PrePolRDS,'Enable','on');
        set(gui.popup_handles.PrePolRamp,'Enable','on');
        set(gui.edit_handles.PrePolFactor,'Enable','on');
        set(gui.edit_handles.PrePolTheta,'Enable','on');
        set(gui.edit_handles.PrePolPhi,'Enable','on');
        set(gui.edit_handles.PrePolSwitchFactor,'Enable','on');
        set(gui.edit_handles.PrePolTramp,'Enable','on');
        set(gui.edit_handles.PrePolTslope,'Enable','on');
        setappdata(fig,'data',data);
        onPopupPrePolRamp(gui.popup_handles.PrePolRamp);
        data = getappdata(fig,'data');
        
        % activate Pulse settings
        set(gui.check_handles.PulseRDP,'Enable','on');
        set(gui.popup_handles.PulseType,'Enable','on');
        set(gui.popup_handles.PulseAxis,'Enable','on');
        set(gui.popup_handles.PulsePolarization,'Enable','on');        
        set(gui.edit_handles.PulseB1Factor,'Enable','on');
        set(gui.edit_handles.PulseTtau,'Enable','on');
        set(gui.edit_handles.PulseDFstart,'Enable','on');
        set(gui.edit_handles.PulseTwait,'Enable','on');
        setappdata(fig,'data',data);
        onPopupPulseType(gui.popup_handles.PulseType);
        data = getappdata(fig,'data');
    end
    
    % because the settings changed, deactivate the "Animate" button
    set(gui.push_handles.Animate,'Enable','off');
    % update all data inside the GUI
    setappdata(fig,'data',data);
    % update status bar
    updateStatusInformation(fig);
else
    warndlg({'onCheckPrePolPulse:','There is no figure with the BLOCHUS Tag open.'},...
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
