function onPopupPulseAxis(src,~)
%onPopupPulseAxis sets the orientation of the pulse axis
%
% Syntax:
%       onPopupPulseAxis(src)
%
% Inputs:
%       src - handle of the calling object
%
% Outputs:
%       none
%
% Example:
%        onPopupPulseAxis(src)
%
% Other m-files required:
%       getPulseParameters
%       plotBpulse
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
    
    % get the popup menu entry
    val = get(src,'Value');
    
    % set the corresponding pulse axis
    switch val
        case 1
            data.pulse.Axis = '+x';
        case 2
            data.pulse.Axis = '-x';
        case 3
            data.pulse.Axis = '+y';
        case 4
            data.pulse.Axis = '-y';
    end
    
    % update GUI data
    setappdata(fig,'data',data);
    % update pulse settings
    getPulseParameters(gui.figh);
    % get GUI data
    data = getappdata(fig,'data');
    % plot pulse
    plotBpulse(gui.figh);
    
else
    warndlg({'onPopupPulseAxis:','There is no figure with the BLOCHUS Tag open.'},...
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
