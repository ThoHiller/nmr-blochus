function BLOCHUS
%BLOCHUS is a graphical user interface (GUI) to simulate NMR spin dynamics
%based on the Bloch equations
%
% Syntax:
%       BLOCHUS
%
% Inputs:
%       none
%
% Outputs:
%       none
%
% Example:
%       BLOCHUS
%
% Other m-files required:
%       BLOCHUS_createGUI.m
%       BLOCHUS_loadDefaults.m
%       calculateGuiOnMonitorPosition.m
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

%% GUI 'header' info and default GUI settings
myui.version = '0.1.2';
myui.date = '17.05.2020';
myui.author = 'Thomas Hiller';
myui.email = 'thomas.hiller[at]leibniz-liag.de';

myui.fontsize = 9;
myui.axfontsize = 11;
myui.linewidth = 2;
myui.color.basic = [143 188 143]./255;
myui.color.prepol = [222 184 135]./255;
myui.color.prepolB = [0.635 0.078 0.184];
myui.color.wait = [240 128 128]./255;
myui.color.pulse = [100 149 237]./255;

%% Default data settings
data = BLOCHUS_loadDefaults;

%% GUI initialization
gui.figh = figure('Name','BLOCHUS - BLOCHUniversalSimulator',...
    'NumberTitle','off','Tag','BLOCHUS','ToolBar','none','MenuBar','none',...
    'SizeChangedFcn',@onFigureSizeChange);

% position on screen
pos = BLOCHUS_setPositionOnScreen;
set(gui.figh,'Position',pos);

%% GUI data
gui.myui = myui;

% save the data struct within the GUI
setappdata(gui.figh,'data',data);
setappdata(gui.figh,'gui',gui);

%% Create GUI elements
BLOCHUS_createGUI(gui.figh,true);
% update status bar
updateStatusInformation(gui.figh);

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
