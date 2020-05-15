function onMenuImport(src,~)
%onMenuImport handles the extra menu entries
%
% Syntax:
%       onMenuImport(src)
%
% Inputs:
%       src - handle of the calling object
%
% Outputs:
%       none
%
% Example:
%       onMenuImport(src)
%
% Other m-files required:
%       switchToolTips
%       updateStatusInformation
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
    
    % after the import these values should be strings
    Sessionpath = -1;
    Sessionfile = -1;
    % 'pathstr' hold s the name of the chosen data path
    [pathstr,~,~] = fileparts(pwd);
    % get the file name
    [Sessionfile,Sessionpath] = uigetfile(pathstr,...
        'Choose BLOCHUS session file');
    
    % only continue if user didn't cancel
    if sum(Sessionpath) > 0
        % check if it is a valid session file
        tmp = load(fullfile(Sessionpath,Sessionfile),'savedata');
        if isfield(tmp,'savedata') && isfield(tmp.savedata,'data') && ...
                isfield(tmp.savedata,'isPulse') && isfield(tmp.savedata,'isPrePol')
            savedata = tmp.savedata;
            
            % copy data
            data.info = savedata.data.info;
            data.init = savedata.data.init;
            data.basic = savedata.data.basic;
            data.prepol = savedata.data.prepol;
            data.pulse = savedata.data.pulse;
            if isfield(savedata.data,'results')
                data.results = savedata.data.results;
            end
            % update GUI data
            setappdata(fig,'data',data);
            
            % update all edit-fields
            set(gui.edit_handles.B0,'String',num2str(data.basic.B0));
            set(gui.edit_handles.Omega0,'String',num2str(data.basic.Omega0));
            set(gui.edit_handles.T1relax,'String',num2str(data.basic.T1relax));
            set(gui.edit_handles.T2relax,'String',num2str(data.basic.T2relax));
            set(gui.edit_handles.Tsim,'String',num2str(data.basic.Tsim));
            set(gui.edit_handles.M0x,'String',num2str(data.basic.M0(1)));
            set(gui.edit_handles.M0y,'String',num2str(data.basic.M0(2)));
            set(gui.edit_handles.M0z,'String',num2str(data.basic.M0(3)));
            set(gui.edit_handles.Minitx,'String',num2str(data.basic.Minit(1)));
            set(gui.edit_handles.Minity,'String',num2str(data.basic.Minit(2)));
            set(gui.edit_handles.Minitz,'String',num2str(data.basic.Minit(3)));
            
            set(gui.edit_handles.PrePolFactor,'String',num2str(data.prepol.Factor));
            set(gui.edit_handles.PrePolTheta,'String',num2str(data.prepol.Theta));
            set(gui.edit_handles.PrePolPhi,'String',num2str(data.prepol.Phi));
            set(gui.edit_handles.PrePolSwitchFactor,'String',num2str(data.prepol.SwitchFactor));
            set(gui.edit_handles.PrePolTramp,'String',num2str(data.prepol.Tramp));
            set(gui.edit_handles.PrePolTslope,'String',num2str(data.prepol.Tslope));
            
            set(gui.edit_handles.PulseB1Factor,'String',num2str(data.pulse.B1Factor));
            set(gui.edit_handles.PulseTtau,'String',num2str(data.pulse.Ttau));
            set(gui.edit_handles.PulseDFstart,'String',num2str(data.pulse.DFstart));
            set(gui.edit_handles.PulseDFend,'String',num2str(data.pulse.DFend));
            set(gui.edit_handles.PulseDFA,'String',num2str(data.pulse.DFA));
            set(gui.edit_handles.PulseDFB,'String',num2str(data.pulse.DFB));
            set(gui.edit_handles.PulseIstart,'String',num2str(data.pulse.Istart));
            set(gui.edit_handles.PulseIend,'String',num2str(data.pulse.Iend));
            set(gui.edit_handles.PulseIA,'String',num2str(data.pulse.IA));
            set(gui.edit_handles.PulseIB,'String',num2str(data.pulse.IB));
            set(gui.edit_handles.PulseQ,'String',num2str(data.pulse.Q));
            set(gui.edit_handles.PulseQdf,'String',num2str(data.pulse.Qdf));
            set(gui.edit_handles.PulseTwait,'String',num2str(data.pulse.Twait));
            
            % update nucleus popup
            set(gui.popup_handles.Nuc,'Value',savedata.Nucleus);
            
            % update PrePol uicontrols
            set(gui.check_handles.PrePolRDS,'Value',data.prepol.RDS);
            set(gui.popup_handles.PrePolRamp,'Value',savedata.Ramp);
            set(gui.check_handles.PrePol,'Value',savedata.isPrePol);
            
            % update Pulse uicontrols
            set(gui.check_handles.PulseRDP,'Value',data.pulse.RDP);
            set(gui.popup_handles.PulseType,'Value',savedata.PulseType);
            set(gui.popup_handles.PulseAxis,'Value',savedata.PulseAxis);
            set(gui.popup_handles.PulsePolarization,'Value',savedata.PulsePolarization);
            set(gui.popup_handles.PulseDFmode,'Value',savedata.PulseDFmode);
            set(gui.popup_handles.PulseImode,'Value',savedata.PulseImode);
            set(gui.check_handles.PulseQ,'Value',savedata.PulseQ);
            onCheckPulseQ(gui.check_handles.PulseQ);
            set(gui.check_handles.Pulse,'Value',savedata.isPulse);
            
            % activate the panels (if any)
            onCheckPrePolPulse(gui.check_handles.PrePol);
            
            % plot results (if any)
            if isfield(savedata.data,'results')
                plotResults(fig);
            end
            
        else
            helpdlg({'onMenuImport:';...
                'This seems to be not a valid BLOCHUS session file'},...
                'No session data found');
        end        
    end
    
else
    warndlg({'onMenuImport:','There is no figure with the BLOCHUS Tag open.'},...
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
