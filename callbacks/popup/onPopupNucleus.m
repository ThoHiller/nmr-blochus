function onPopupNucleus(src,~)
%onPopupNucleus selects the nucleus (proton) to use for the simulation
%
% Syntax:
%       onPopupNucleus(src)
%
% Inputs:
%       src - handle of the calling object
%
% Outputs:
%       none
%
% Example:
%       onPopupNucleus(src)
%
% Other m-files required:
%       getGyroRatio
%       getOmega0
%       onEditValue
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
    gui = getappdata(fig,'gui');
    data = getappdata(fig,'data');
    
    % get the popup menu entry
    val = get(src,'Value');
    
    % set the corresponding nucleus
    switch val
        case 1
            data.basic.nucleus = '1H';
        case 2
            data.basic.nucleus = '2H';
        case 3
            data.basic.nucleus = '3He';
        case 4
            data.basic.nucleus = '7Li';
        case 5
            data.basic.nucleus = '13C';
        case 6
            data.basic.nucleus = '14N';
        case 7
            data.basic.nucleus = '15N';
        case 8
            data.basic.nucleus = '17O';
        case 9
            data.basic.nucleus = '19F';
        case 10
            data.basic.nucleus = '23Na';
        case 11
            data.basic.nucleus = '27Al';
        case 12
            data.basic.nucleus = '29Si';
        case 13
            data.basic.nucleus = '31P';
        case 14
            data.basic.nucleus = '57Fe';
        case 15
            data.basic.nucleus = '63Cu';
        case 16
            data.basic.nucleus = '67Zn';
        case 17
            data.basic.nucleus = '129Xe';
    end
    % update the gyromagnetic ration and Larmor frequency accordingly
    data.basic.gamma = getGyroRatio(data.basic.nucleus);
    data.basic.Omega0 = getOmega0(data.basic.gamma,data.basic.B0)/2/pi;
    
    % update the GUI data
    setappdata(fig,'data',data);
    
    % update the GUI
    set(gui.edit_handles.Omega0,'String',sprintf('%7.2f',data.basic.Omega0));
    set(gui.edit_handles.Gyro,'String',sprintf('%5.4e',data.basic.gamma));
    
    % if the pulse panel is activated
    switch data.basic.type
        case {'pulse','prepolpulse'}
            % adapt the B1 amplitude due to the changed gamma value
            onEditValue(gui.edit_handles.PulseTtau);
        otherwise
            % nothing to do
    end
    
else
    warndlg({'onPopupNucleus:','There is no figure with the BLOCHUS Tag open.'},...
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
