function gui = BLOCHUS_createMenus(gui)
%BLOCHUS_createMenus creates all GUI menus
%
% Syntax:
%       gui = NUCLEUSinv_createMenus(gui)
%
% Inputs:
%       gui - figure gui elements structure
%
% Outputs:
%       gui
%
% Example:
%       gui = BLOCHUS_createMenus(gui)
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

%% 1. File
gui.menu.file = uimenu(gui.figh,...
    'Label','File','Enable','off');

% 1.1 Import
gui.menu.file_import = uimenu(gui.menu.file,...
    'Label','Import','Enable','on');

% 1.1.1 BLOCHUS session file
gui.menu.file_import_session = uimenu(gui.menu.file_import,...
    'Label','Session file','Callback',@onMenuImport);

% 1.2 Export
gui.menu.file_export = uimenu(gui.menu.file,...
    'Label','Export');

% 1.2.1 BLOCHUS session file
gui.menu.file_export_session = uimenu(gui.menu.file_export,...
    'Label','Session file','Callback',@onMenuExport);

% 1.3 Restart
gui.menu.file_restart = uimenu(gui.menu.file,...
    'Label','Restart','Separator','on','Callback',@onMenuRestartQuit);

% 1.4 Quit
gui.menu.file_quit = uimenu(gui.menu.file,...
    'Label','Quit','Separator','on','Callback',@onMenuRestartQuit);

%% 2. Extras
gui.menu.view = uimenu(gui.figh,...
    'Label','View','Enable','off');
% 2.1 Tooltips (on/off)
gui.menu.view_tooltips = uimenu(gui.menu.view,...
    'Label','Tooltips','Checked','on','Callback',@onMenuViewShow);
% 2.2 Figure Toolbar
gui.menu.view_toolbar = uimenu(gui.menu.view,...
    'Label','Figure Toolbar','Callback',@onMenuViewShow);
% 2.3.1 Tigures
gui.menu.view_figures = uimenu(gui.menu.view,...
    'Label','Figures','Separator','on');
% 2.3.1.1 current
gui.menu.view_figures_current = uimenu(gui.menu.view_figures,...
    'Label','Current View','Callback',@onMenuViewFigure);
% 2.3.1.2 only magnetization
gui.menu.view_figures_mag = uimenu(gui.menu.view_figures,...
    'Label','Magnetization','Callback',@onMenuViewFigure);
% 2.3.1.3 only ramp
gui.menu.view_figures_ramp = uimenu(gui.menu.view_figures,...
    'Label','Switch-off Ramp','Enable','off','Callback',@onMenuViewFigure);
% 2.3.1.4 only pulse
gui.menu.view_figures_pulse = uimenu(gui.menu.view_figures,...
    'Label','Pulse','Enable','off','Callback',@onMenuViewFigure);

%% 3. Help
gui.menu.help = uimenu(gui.figh,...
    'Label','Help','Enable','off');

% 3.1 About
gui.menu.help_about = uimenu(gui.menu.help,...
    'Label','About','Callback',@onMenuHelp);

return

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
