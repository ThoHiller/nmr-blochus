function onPushAxView(src,~)
%onPushAxView sets the view of the Bloch sphere plot to predefined sets
%
% Syntax:
%       onPushAxView(src)
%
% Inputs:
%       src - handle of the calling object
%
% Outputs:
%       none
%
% Example:
%       onPushAxView(src)
%
% Other m-files required:
%       none
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
    % get button tag
    tag = get(src,'String');
    
    switch tag
        case 'XZ'
            axes(gui.axes_handles.SphereL)
            view([0 0]);
            axes(gui.axes_handles.SphereR)
            view([0 0]);
        case 'YZ'
            axes(gui.axes_handles.SphereL)
            view([90 0]);
            axes(gui.axes_handles.SphereR)
            view([90 0]);
        case 'XY'
            axes(gui.axes_handles.SphereL)
            view([0 90]);
            axes(gui.axes_handles.SphereR)
            view([0 90]);
        case '3D'
            axes(gui.axes_handles.SphereL)
            view([-35 30]);
            axes(gui.axes_handles.SphereR)
            view([-35 30]);
    end
    
else
    warndlg({'onPushAxView:','There is no figure with the BLOCHUS Tag open.'},...
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
