function [gui,myui] = BLOCHUS_createStatusbar(gui,myui)
%BLOCHUS_createStatusbar creates the bottom status bar
%
% Syntax:
%       [gui,myui] = BLOCHUS_createStatusbar(gui,myui)
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
%       [gui,myui] = BLOCHUS_createStatusbar(gui,myui)
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

% create panels inside the bottom hbox to show persistent status
% information
gui.panels.status.main = uix.Panel('Parent',gui.bottom);
gui.panels.status.Timer = uix.Panel('Parent',gui.bottom);
gui.panels.status.PrePol = uix.Panel('Parent',gui.bottom);
gui.panels.status.Pulse = uix.Panel('Parent',gui.bottom);
gui.panels.status.MIDI = uix.Panel('Parent',gui.bottom);
gui.panels.status.Tooltips = uix.Panel('Parent',gui.bottom);
gui.panels.status.Version = uix.Panel('Parent',gui.bottom);

% adjust the panel widths
set(gui.bottom,'Widths',[330 -1 -1 -1 -1 -1 -1]);

gui.text_handles.Status = uicontrol('Style','Text',...
    'Parent',gui.panels.status.main,...
    'String','',...
    'HorizontalAlignment','left',...
    'FontSize',8);
gui.text_handles.TimerStat = uicontrol('Style','Text',...
    'Parent',gui.panels.status.Timer,...
    'String','Calc. Time: 0 s',...
    'HorizontalAlignment','left',...
    'FontSize',8);
gui.text_handles.PrePolStat = uicontrol('Style','Text',...
    'Parent',gui.panels.status.PrePol,...
    'String','Pre-Polarization: OFF',...
    'HorizontalAlignment','left',...
    'FontSize',8);
gui.text_handles.PulseStat = uicontrol('Style','Text',...
    'Parent',gui.panels.status.Pulse,...
    'String','Pulse: OFF',...
    'HorizontalAlignment','left',...
    'FontSize',8);
gui.text_handles.MIDIStat = uicontrol('Style','Text',...
    'Parent',gui.panels.status.MIDI,...
    'String','Discrete Pulse: OFF',...
    'HorizontalAlignment','left',...
    'FontSize',8);
gui.text_handles.TooltipsStat = uicontrol('Style','Text',...
    'Parent',gui.panels.status.Tooltips,...
    'String','Tooltips: ON',...
    'HorizontalAlignment','left',...
    'FontSize',8);
gui.text_handles.VersionStat = uicontrol('Style','Text',...
    'Parent',gui.panels.status.Version,...
    'String','Version: ',...
    'HorizontalAlignment','left',...
    'FontSize',8);

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
