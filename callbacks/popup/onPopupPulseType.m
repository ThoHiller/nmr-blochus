function onPopupPulseType(src,~)
%onPopupPulseType selects the pulse type
%
% Syntax:
%       onPopupPulseType(src)
%
% Inputs:
%       src - handle of the calling object
%
% Outputs:
%       none
%
% Example:
%       onPopupPulseType(src)
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
    
    % set the corresponding pulse type
    switch val
        case 1 % pi/2
            data.pulse.Type = 'pi_half';
            set(gui.text_handles.PulseTtau,'String',['Pulse length ',char(hex2dec('3C4')),' [ms]']);
            set(gui.edit_handles.PulseTtau,'String',num2str(data.pulse.Ttau));
            % adjust B1 amplitude
            data.pulse.B1Factor = abs(1e3*(pi/2/(data.basic.gamma*data.basic.B0*data.pulse.Ttau)));
            set(gui.edit_handles.PulseB1Factor,'String',num2str(data.pulse.B1Factor));
            % set offset frequency
            data.pulse.DFstart = 0;
            data.pulse.DFend = 0;
            set(gui.text_handles.PulseDF,'String','Offset [Hz]');
            set(gui.edit_handles.PulseDFstart,'String',num2str(data.pulse.DFstart),'Enable','on');
            set(gui.edit_handles.PulseDFend,'String',num2str(data.pulse.DFend),'Enable','off');
            % disable all other frequency stuff
            set(gui.popup_handles.PulseDFmode,'Enable','off');
            set(gui.edit_handles.PulseDFA,'Enable','off');
            set(gui.edit_handles.PulseDFB,'Enable','off');
            % set current
            data.pulse.Istart = 1;
            data.pulse.Iend = 1;
            set(gui.edit_handles.PulseIstart,'String',num2str(data.pulse.Istart),'Enable','off');
            set(gui.edit_handles.PulseIend,'String',num2str(data.pulse.Iend),'Enable','off');
            set(gui.popup_handles.PulseImode,'Enable','off');
            set(gui.edit_handles.PulseIA,'Enable','off');
            set(gui.edit_handles.PulseIB,'Enable','off');
            set(gui.check_handles.PulseQ,'Enable','off');
            set(gui.edit_handles.PulseQ,'Enable','off');
            set(gui.edit_handles.PulseQdf,'Enable','off');
            
        case 2 % pi
            data.pulse.Type = 'pi';
            set(gui.text_handles.PulseTtau,'String',['Pulse length ',char(hex2dec('3C4')),' [ms]']);
            set(gui.edit_handles.PulseTtau,'String',num2str(data.pulse.Ttau));
            % adjust B1 amplitude
            data.pulse.B1Factor = abs(1e3*(pi/(data.basic.gamma*data.basic.B0*data.pulse.Ttau)));
            set(gui.edit_handles.PulseB1Factor,'String',num2str(data.pulse.B1Factor));
            % set offset frequency
            data.pulse.DFstart = 0;
            data.pulse.DFend = 0;
            set(gui.text_handles.PulseDF,'String','Offset [Hz]');
            set(gui.edit_handles.PulseDFstart,'String',num2str(data.pulse.DFstart),'Enable','on');
            set(gui.edit_handles.PulseDFend,'String',num2str(data.pulse.DFend),'Enable','off');
            % disable all other frequency stuff
            set(gui.popup_handles.PulseDFmode,'Enable','off');
            set(gui.edit_handles.PulseDFA,'Enable','off');
            set(gui.edit_handles.PulseDFB,'Enable','off');
            % set current
            data.pulse.Istart = 1;
            data.pulse.Iend = 1;
            set(gui.edit_handles.PulseIstart,'String',num2str(data.pulse.Istart),'Enable','off');
            set(gui.edit_handles.PulseIend,'String',num2str(data.pulse.Iend),'Enable','off');
            set(gui.popup_handles.PulseImode,'Enable','off');
            set(gui.edit_handles.PulseIA,'Enable','off');
            set(gui.edit_handles.PulseIB,'Enable','off');
            set(gui.check_handles.PulseQ,'Enable','off');
            set(gui.edit_handles.PulseQ,'Enable','off');
            set(gui.edit_handles.PulseQdf,'Enable','off');
            
        case 3 % free
            data.pulse.Type = 'free';
            set(gui.text_handles.PulseTtau,'String',['Pulse length ',char(hex2dec('3C4')),' [ms]']);
            set(gui.edit_handles.PulseTtau,'String',num2str(data.pulse.Ttau));
            % set offset frequency
            data.pulse.DFstart = 0;
            data.pulse.DFend = 0;
            set(gui.text_handles.PulseDF,'String','Offset [Hz]');
            set(gui.edit_handles.PulseDFstart,'String',num2str(data.pulse.DFstart),'Enable','on');
            set(gui.edit_handles.PulseDFend,'String',num2str(data.pulse.DFend),'Enable','off');
            % disable all other frequency stuff
            set(gui.popup_handles.PulseDFmode,'Enable','off');
            set(gui.edit_handles.PulseDFA,'Enable','off');
            set(gui.edit_handles.PulseDFB,'Enable','off');
            % set current
            data.pulse.Istart = 1;
            data.pulse.Iend = 1;
            set(gui.edit_handles.PulseIstart,'String',num2str(data.pulse.Istart),'Enable','off');
            set(gui.edit_handles.PulseIend,'String',num2str(data.pulse.Iend),'Enable','off');
            set(gui.popup_handles.PulseImode,'Enable','off');
            set(gui.edit_handles.PulseIA,'Enable','off');
            set(gui.edit_handles.PulseIB,'Enable','off');
            set(gui.check_handles.PulseQ,'Enable','off');
            set(gui.edit_handles.PulseQ,'Enable','off');
            set(gui.edit_handles.PulseQdf,'Enable','off');            
            
        case 4 % AHP
            data.pulse.Type = 'AHP';
            set(gui.text_handles.PulseTtau,'String',['Pulse length ',char(hex2dec('3C4')),' [ms]']);
            set(gui.edit_handles.PulseTtau,'String',num2str(data.pulse.Ttau));
            
            set(gui.text_handles.PulseDF,'String','start [Hz] end [Hz] A B');
            set(gui.popup_handles.PulseDFmode,'Enable','on');
            set(gui.edit_handles.PulseDFstart,'Enable','on');
            set(gui.edit_handles.PulseDFend,'Enable','on');
            set(gui.edit_handles.PulseDFA,'Enable','off');
            set(gui.edit_handles.PulseDFB,'Enable','off');
            
            if get(gui.check_handles.PulseQ,'Value') == 0
                data.pulse.Istart = 0;
            else
                data.pulse.Istart = 1;
            end
            data.pulse.Iend = 1;
            set(gui.popup_handles.PulseImode,'Enable','on');
            set(gui.edit_handles.PulseIstart,'Enable','on','String',num2str(data.pulse.Istart));
            set(gui.edit_handles.PulseIend,'Enable','on','String',num2str(data.pulse.Iend));
            set(gui.edit_handles.PulseIA,'Enable','off');
            set(gui.edit_handles.PulseIB,'Enable','off');
            set(gui.check_handles.PulseQ,'Enable','on');
            
            setappdata(fig,'data',data);
            onCheckPulseQ(gui.check_handles.PulseQ);
            data = getappdata(fig,'data');
            setappdata(fig,'data',data);
            onPopupPulseDFmode(gui.popup_handles.PulseDFmode);
            data = getappdata(fig,'data');
            setappdata(fig,'data',data);
            onPopupPulseImode(gui.popup_handles.PulseImode);
            data = getappdata(fig,'data');
  
        case 5 % MIDI OR
            data.pulse.Type = 'MIDI_OR';
            % calculate number of periods based on current pulse time and fL
            data.pulse.MIDINP = ceil(abs(data.pulse.Ttau/1000 * data.basic.Omega0));
            set(gui.edit_handles.PulseTtau,'String',num2str(data.pulse.MIDINP));
            set(gui.text_handles.PulseTtau,'String','Number of periods');         
            set(gui.text_handles.PulseDF,'String','Offset [Hz]');
            set(gui.popup_handles.PulseDFmode,'Enable','off');
            data.pulse.DFstart = 0;
            data.pulse.DFend = 0;
            set(gui.edit_handles.PulseDFstart,'Enable','off','String',num2str(data.pulse.DFstart));
            set(gui.edit_handles.PulseDFend,'Enable','off','String',num2str(data.pulse.DFend));
            set(gui.edit_handles.PulseDFA,'Enable','off');
            set(gui.edit_handles.PulseDFB,'Enable','off');
            
            set(gui.popup_handles.PulseImode,'Enable','off');
            set(gui.edit_handles.PulseIstart,'Enable','on');
            data.pulse.Istart = 0.6;
            data.pulse.Iend = 0.6;
            set(gui.edit_handles.PulseIstart,'String',num2str(data.pulse.Istart));
            set(gui.edit_handles.PulseIend,'Enable','off');
            set(gui.edit_handles.PulseIA,'Enable','off');
            set(gui.edit_handles.PulseIB,'Enable','off');
            set(gui.check_handles.PulseQ,'Enable','off');
            set(gui.edit_handles.PulseQ,'Enable','off');
            set(gui.edit_handles.PulseQdf,'Enable','off');
            
            % update the GUI as if the number of periods would have been
            % changed manually
            setappdata(fig,'data',data);
            onEditValue(gui.edit_handles.PulseTtau);
            data = getappdata(fig,'data');
            
        case 6 % MIDI AP
            data.pulse.Type = 'MIDI_AP';
            % calculate number of periods based on current pulse time and fL
            data.pulse.MIDINP = ceil(abs(data.pulse.Ttau/1000 * data.basic.Omega0));
            set(gui.edit_handles.PulseTtau,'String',num2str(data.pulse.MIDINP));
            set(gui.text_handles.PulseTtau,'String','Number of periods');         
            set(gui.text_handles.PulseDF,'String','start [Hz] end [Hz] A B');
            data.pulse.DFmode = 'tanhMIDI';
            set(gui.popup_handles.PulseDFmode,'Enable','on','Value',2);
            data.pulse.DFstart = -200;
            data.pulse.DFend = 0;
            set(gui.edit_handles.PulseDFstart,'Enable','on','String',num2str(data.pulse.DFstart));
            set(gui.edit_handles.PulseDFend,'Enable','on','String',num2str(data.pulse.DFend));
            data.pulse.DFA = 0.7;
            data.pulse.DFB = 0;
            set(gui.edit_handles.PulseDFA,'Enable','on','String',num2str(data.pulse.DFA));
            set(gui.edit_handles.PulseDFB,'Enable','on','String',num2str(data.pulse.DFB));
            
            data.pulse.Imode = 'tanhMIDI';
            set(gui.popup_handles.PulseImode,'Enable','on','Value',2);
            data.pulse.Istart = 0.15;
            data.pulse.Iend = 0.6;
            set(gui.edit_handles.PulseIstart,'Enable','on','String',num2str(data.pulse.Istart));
            set(gui.edit_handles.PulseIend,'Enable','on','String',num2str(data.pulse.Iend));
            data.pulse.IA = 0.5;
            data.pulse.IB = 1;
            set(gui.edit_handles.PulseIA,'Enable','on','String',num2str(data.pulse.IA));
            set(gui.edit_handles.PulseIB,'Enable','on','String',num2str(data.pulse.IB));
            set(gui.check_handles.PulseQ,'Enable','off');
            set(gui.edit_handles.PulseQ,'Enable','off');
            set(gui.edit_handles.PulseQdf,'Enable','off');
            
            % update the GUI as if the number of periods would have been
            % changed manually
            setappdata(fig,'data',data);
            onEditValue(gui.edit_handles.PulseTtau);
            data = getappdata(fig,'data');
    end
    
    % update GUI data
    setappdata(fig,'data',data);
    % update pulse settings
    getPulseParameters(fig);
    % get GUI data
    data = getappdata(fig,'data');
    % update GUI data
    setappdata(fig,'data',data);
    % plot pulse
    plotPulse(fig);
    % update status bar
    updateStatusInformation(fig);
else
    warndlg({'onPopupPulseType:','There is no figure with the BLOCHUS Tag open.'},...
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
