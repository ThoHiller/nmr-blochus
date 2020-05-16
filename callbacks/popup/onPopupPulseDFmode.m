function onPopupPulseDFmode(src,~)
%onPopupPulseDFmode selects the frequency modulation for an adiabatic pulse
%
% Syntax:
%       onPopupPulseDFmode(src)
%
% Inputs:
%       src - handle of the calling object
%
% Outputs:
%       none
%
% Example:
%       onPopupPulseDFmode(src)
%
% Other m-files required:
%       getPulseParameters;
%       plotPulse;
%
% Subfunctions:
%       none
%
% MAT-files required:
%       none
%
% See also BLOCHUS
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
    
    % get the popup menu entry
    val = get(src,'Value');
    
    % set the corresponding frequency modulation
    switch val
        case 1 % linear
            data.pulse.DFmode = 'lin';
            set(gui.edit_handles.PulseDFA,'Enable','off');
            set(gui.edit_handles.PulseDFB,'Enable','off');
            
        case 2 % tanh MIDI
            data.pulse.DFmode = 'tanhMIDI';
            data.pulse.DFstart = -200;
            data.pulse.DFA = 0.7;
            set(gui.edit_handles.PulseDFstart,'String',num2str(data.pulse.DFstart));
            set(gui.edit_handles.PulseDFA,'String',num2str(data.pulse.DFA));
            set(gui.edit_handles.PulseDFA,'Enable','on');
            set(gui.edit_handles.PulseDFB,'Enable','on');
            
        case 3 % tanh GMR
            data.pulse.DFmode = 'tanhGMR';
            data.pulse.DFstart = -300;
            set(gui.edit_handles.PulseDFstart,'String',num2str(data.pulse.DFstart));
            set(gui.edit_handles.PulseDFA,'Enable','off');
            set(gui.edit_handles.PulseDFB,'Enable','off');
            
        case 4 % exp
            data.pulse.DFmode = 'exp';
            data.pulse.DFstart = -300;
            data.pulse.DFA = 10;
            set(gui.edit_handles.PulseDFstart,'String',num2str(data.pulse.DFstart));
            set(gui.edit_handles.PulseDFA,'String',num2str(data.pulse.DFA));
            set(gui.edit_handles.PulseDFA,'Enable','on');
            set(gui.edit_handles.PulseDFB,'Enable','off');
    end
    
    % update GUI data
    setappdata(fig,'data',data);
    % update pulse settings
    getPulseParameters(gui.figh);
    % get GUI data
    data = getappdata(fig,'data');
    % update GUI data
    setappdata(fig,'data',data);
    % plot pulse
    plotPulse(gui.figh);
    
else
    warndlg({'onPopupPulseDFmode:','There is no figure with the BLOCHUS Tag open.'},...
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
