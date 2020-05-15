function onPopupPrePolRamp(src,~)
%onPopupPrePolRamp selects the pre-polarization switch-off ramp
%
% Syntax:
%       onPopupPrePolRamp(src)
%
% Inputs:
%       src - handle of the calling object
%
% Outputs:
%       none
%
% Example:
%       onPopupPrePolRamp(src)
%
% Other m-files required:
%       getRampParameters;
%       plotRamp;
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
    
    % set the corresponding switch-off ramp
    switch val
        case 1 % exponential
            data.prepol.Ramp = 'exp';
            data.prepol.Tslope = data.prepol.Tramp/10;
            data.prepol.Factor = 100;
            data.prepol.SwitchFactor = 1;
            
            set(gui.edit_handles.PrePolFactor,'Enable','on');
            set(gui.edit_handles.PrePolTheta,'Enable','on');
            set(gui.edit_handles.PrePolSwitchFactor,'Enable','off');
            set(gui.edit_handles.PrePolTramp,'Enable','on');
            set(gui.edit_handles.PrePolTslope,'Enable','on');
            
            set(gui.edit_handles.PrePolFactor,'String',num2str(data.prepol.Factor));
            set(gui.edit_handles.PrePolTheta,'String',num2str(data.prepol.Theta));
            set(gui.edit_handles.PrePolSwitchFactor,'String',num2str(data.prepol.SwitchFactor));
            set(gui.edit_handles.PrePolTramp,'String',num2str(data.prepol.Tramp));
            set(gui.edit_handles.PrePolTslope,'String',num2str(data.prepol.Tslope));
            
            set(gui.text_handles.PrePolSwitchFactor,'String','Switch factor [B0]');
            
        case 2 % linear & exponential (generic)
            data.prepol.Ramp = 'linexp';
            data.prepol.Tslope = data.prepol.Tramp/2;
            data.prepol.Factor = 100;
            data.prepol.SwitchFactor = 10;
            
            set(gui.edit_handles.PrePolFactor,'Enable','on');
            set(gui.edit_handles.PrePolTheta,'Enable','on');
            set(gui.edit_handles.PrePolSwitchFactor,'Enable','on');
            set(gui.edit_handles.PrePolTramp,'Enable','on');
            set(gui.edit_handles.PrePolTslope,'Enable','on');
            
            set(gui.edit_handles.PrePolFactor,'String',num2str(data.prepol.Factor));
            set(gui.edit_handles.PrePolTheta,'String',num2str(data.prepol.Theta));
            set(gui.edit_handles.PrePolSwitchFactor,'String',num2str(data.prepol.SwitchFactor));
            set(gui.edit_handles.PrePolTramp,'String',num2str(data.prepol.Tramp));
            set(gui.edit_handles.PrePolTslope,'String',num2str(data.prepol.Tslope));
            
            set(gui.text_handles.PrePolSwitchFactor,'String','Switch factor [B0]');
            
        case 3 % half cosine (generic)
            data.prepol.Ramp = 'halfcos';
            data.prepol.Tslope = data.prepol.Tramp;
            data.prepol.Factor = 100;
            data.prepol.SwitchFactor = 1;
            
            set(gui.edit_handles.PrePolFactor,'Enable','on');
            set(gui.edit_handles.PrePolTheta,'Enable','on');
            set(gui.edit_handles.PrePolSwitchFactor,'Enable','off');
            set(gui.edit_handles.PrePolTramp,'Enable','on');
            set(gui.edit_handles.PrePolTslope,'Enable','off');
            
            set(gui.edit_handles.PrePolFactor,'String',num2str(data.prepol.Factor));
            set(gui.edit_handles.PrePolTheta,'String',num2str(data.prepol.Theta));
            set(gui.edit_handles.PrePolSwitchFactor,'String',num2str(data.prepol.SwitchFactor));
            set(gui.edit_handles.PrePolTramp,'String',num2str(data.prepol.Tramp));
            set(gui.edit_handles.PrePolTslope,'String',num2str(data.prepol.Tslope));
            
            set(gui.text_handles.PrePolSwitchFactor,'String','Switch factor [B0]');    
            
        case 4 % linear (generic)
            data.prepol.Ramp = 'lin';
            data.prepol.Tslope = data.prepol.Tramp;
            data.prepol.Factor = 100;
            data.prepol.SwitchFactor = 1;
            
            set(gui.edit_handles.PrePolFactor,'Enable','on');
            set(gui.edit_handles.PrePolTheta,'Enable','on');
            set(gui.edit_handles.PrePolSwitchFactor,'Enable','off');
            set(gui.edit_handles.PrePolTramp,'Enable','on');
            set(gui.edit_handles.PrePolTslope,'Enable','off');
            
            set(gui.edit_handles.PrePolFactor,'String',num2str(data.prepol.Factor));
            set(gui.edit_handles.PrePolTheta,'String',num2str(data.prepol.Theta));
            set(gui.edit_handles.PrePolSwitchFactor,'String',num2str(data.prepol.SwitchFactor));
            set(gui.edit_handles.PrePolTramp,'String',num2str(data.prepol.Tramp));
            set(gui.edit_handles.PrePolTslope,'String',num2str(data.prepol.Tslope));
            
            set(gui.text_handles.PrePolSwitchFactor,'String','Switch factor [B0]');
            
        case 5 % linear (Melton 1995)
            data.prepol.Ramp = 'lin';
            
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
            
            set(gui.edit_handles.PrePolFactor,'Enable','on');
            set(gui.edit_handles.PrePolTheta,'Enable','off');
            set(gui.edit_handles.PrePolSwitchFactor,'Enable','on');
            set(gui.edit_handles.PrePolTramp,'Enable','off');
            set(gui.edit_handles.PrePolTslope,'Enable','off');
            
            set(gui.edit_handles.PrePolFactor,'String',num2str(data.prepol.Factor));
            set(gui.edit_handles.PrePolTheta,'String',num2str(data.prepol.Theta));
            set(gui.edit_handles.PrePolSwitchFactor,'String',num2str(data.prepol.SwitchFactor));
            set(gui.edit_handles.PrePolTramp,'String',num2str(data.prepol.Tramp));
            set(gui.edit_handles.PrePolTslope,'String',num2str(data.prepol.Tslope));
            
            set(gui.text_handles.PrePolSwitchFactor,'String','cutoff rate');
            
        case 6 % half cosine (MIDI)
            data.prepol.Ramp = 'halfcos';
            
            Tramp = 1; % [ms]
            data.prepol.Tramp = Tramp;
            data.prepol.Tslope = Tramp;            
            data.prepol.Factor = 100;
            data.prepol.SwitchFactor = 1;
            
            switch data.basic.type
                case 'prepol'
                    data.basic.Tsim = Tramp;
                    set(gui.edit_handles.Tsim,'String',num2str(data.basic.Tsim));
                otherwise
                    % nothing to do
            end
            
            set(gui.edit_handles.PrePolFactor,'Enable','on');
            set(gui.edit_handles.PrePolTheta,'Enable','on');
            set(gui.edit_handles.PrePolSwitchFactor,'Enable','off');
            set(gui.edit_handles.PrePolTramp,'Enable','off');
            set(gui.edit_handles.PrePolTslope,'Enable','off');
            
            set(gui.edit_handles.PrePolFactor,'String',num2str(data.prepol.Factor));
            set(gui.edit_handles.PrePolTheta,'String',num2str(data.prepol.Theta));
            set(gui.edit_handles.PrePolSwitchFactor,'String',num2str(data.prepol.SwitchFactor));
            set(gui.edit_handles.PrePolTramp,'String',num2str(data.prepol.Tramp));
            set(gui.edit_handles.PrePolTslope,'String',num2str(data.prepol.Tslope));
            
            set(gui.text_handles.PrePolSwitchFactor,'String','Switch factor [B0]');
    end
    
    % update GUI data
    setappdata(fig,'data',data);    
    % get ramp parameters
    getRampParameters(fig);
    % update ramp plot
    plotRamp(fig);
    
else
    warndlg({'onPopupPrePolRamp:','There is no figure with the BLOCHUS Tag open.'},...
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
