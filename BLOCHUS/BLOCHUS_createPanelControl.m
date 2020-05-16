function [gui,myui] = BLOCHUS_createPanelControl(gui,myui)
%BLOCHUS_createPanelControl creates "Control" panel
%
% Syntax:
%       [gui,myui] = BLOCHUS_createPanelControl(gui,myui)
%
% Inputs:
%       gui - figure gui elements structure
%       myui - individual GUI settings structure
%
% Outputs:
%       gui
%       myui
%
% Example:
%       [gui,myui] = BLOCHUS_createPanelControl(gui,myui)
%
% Other m-files required:
%       findjobj.m
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

gui.panels.Control.HBox = uix.HBox('Parent', gui.panels.Control.main,...
    'Spacing',3,'Padding',3);

tstr = 'Start calculation.';
gui.push_handles.Run = uicontrol('Style','pushbutton',...
    'Parent',gui.panels.Control.HBox,...
    'String','RUN',...
    'Tag','RUN',...
    'ToolTipString',tstr,...
    'FontSize',myui.fontsize,...
    'BackGroundColor','g',...
    'UserData',struct('Tooltipstr',tstr),...
    'Callback',@onPushRun);

tstr = 'Show animation of the latest result.';
gui.push_handles.Animate = uicontrol('Style','pushbutton',...
    'Parent',gui.panels.Control.HBox,...
    'String','ANIMATE',...
    'Tag','ANIMATE',...
    'ToolTipString',tstr,...
    'FontSize',myui.fontsize,...
    'UserData',struct('Tooltipstr',tstr),...
    'Enable','off',...
    'Callback',@onPushAnimate);

tstr = 'Reset axes grid to equally spaced 2x2.';
gui.push_handles.Grid = uicontrol('Style','pushbutton',...
    'Parent',gui.panels.Control.HBox,...
    'String','RESET GRID',...
    'Tag','GRID',...
    'ToolTipString',tstr,...
    'FontSize',myui.fontsize,...
    'UserData',struct('Tooltipstr',tstr),...
    'Callback',@onPushGrid);

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