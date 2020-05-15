function onMenuViewShow(src,~)
%onMenuViewShow handles the extra menu entries
%
% Syntax:
%       onMenuViewShow(src)
%
% Inputs:
%       src - handle of the calling object
%
% Outputs:
%       none
%
% Example:
%       onMenuViewShow(src,~)
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
    
    switch get(src,'Label')
        case 'Tooltips' % switch on/off Tooltips
            onoff = get(gui.menu.view_tooltips,'Checked');
            switchToolTips(gui,onoff);
            switch onoff
                case 'on'
                    set(gui.menu.view_tooltips,'Checked','off');
                    data.info.ToolTips = 0;
                case 'off'
                    set(gui.menu.view_tooltips,'Checked','on');
                    data.info.ToolTips = 1;
            end
            
        case 'Figure Toolbar' % switch on/off the default Figure Toolbar
            onoff = get(gui.menu.view_toolbar,'Checked');
            switch onoff
                case 'on'
                    set(gui.menu.view_toolbar,'Checked','off');
                    viewmenufcn('FigureToolbar');
                case 'off'
                    set(gui.menu.view_toolbar,'Checked','on');
                    viewmenufcn('FigureToolbar');
            end
    end
    
    % update GUI data
    setappdata(fig,'gui',gui);
    setappdata(fig,'data',data);
    % update status bar
    updateStatusInformation(fig);
    
else
    warndlg({'onMenuViewShow:','There is no figure with the BLOCHUS Tag open.'},...
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
