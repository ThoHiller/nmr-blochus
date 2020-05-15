function [gui,myui] = BLOCHUS_createGridPlots(gui,myui)
%BLOCHUS_createGridPlots creates the "Plots" grid panel
%
% Syntax:
%       [gui,myui] = BLOCHUS_createGridPlots(gui,myui)
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
%       [gui,myui] = BLOCHUS_createGridPlots(gui,myui)
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

%% first create the individual tab panels
%(1,1) magnetization
gui.panels.Plot.Mag = uix.TabPanel('Parent',gui.right,...
    'BackgroundColor',myui.color.basic);
%(2,1) pre-polarization settings
gui.panels.Plot.PrePol = uix.TabPanel('Parent',gui.right,...
    'BackgroundColor',myui.color.prepol);
%(1,2) bloch sphere
gui.panels.Plot.Sphere = uix.TabPanel('Parent',gui.right,...
    'BackgroundColor',myui.color.basic);
%(2,2) pulse settings
gui.panels.Plot.Pulse = uix.TabPanel('Parent',gui.right,...
    'BackgroundColor',myui.color.pulse); 

%% magnetization
% because Matlab is buggy when resizing, put uicontainer inside the panels
plotXYZL = uicontainer('Parent',gui.panels.Plot.Mag);
plotXYZR = uicontainer('Parent',gui.panels.Plot.Mag);
plotXYZFFT = uicontainer('Parent',gui.panels.Plot.Mag);
gui.panels.Plot.Mag.TabTitles = {'Lab frame','Rot frame','FFT'};
gui.panels.Plot.Mag.TabWidth = 75;

gui.axes_handles.MagL = axes('Parent',plotXYZL,'Box','on');
gui.axes_handles.MagR = axes('Parent',plotXYZR,'Box','on');
gui.axes_handles.MagFFT = axes('Parent',plotXYZFFT,'Box','on');
clearSingleAxis(gui.axes_handles.MagL);
clearSingleAxis(gui.axes_handles.MagR);
clearSingleAxis(gui.axes_handles.MagFFT);

%% PrePol tab panel
plotBpre = uicontainer('Parent',gui.panels.Plot.PrePol);
plotalpha = uicontainer('Parent',gui.panels.Plot.PrePol);
plotdadt = uicontainer('Parent',gui.panels.Plot.PrePol);
plotwda = uicontainer('Parent',gui.panels.Plot.PrePol);
gui.panels.Plot.PrePol.TabTitles = {'Bp ramp',...
    [char(hex2dec('3B1')),'=',char(hex2dec('2221')),'B0 B',],...
    ['d',char(hex2dec('3B1')),' / dt'],...
    'adiab. cond.'};
gui.panels.Plot.PrePol.TabWidth = 75;

gui.axes_handles.Bpre = axes('Parent',plotBpre,'Box','on');
gui.axes_handles.alpha = axes('Parent',plotalpha,'Box','on');
gui.axes_handles.dadt = axes('Parent',plotdadt,'Box','on');
gui.axes_handles.wda = axes('Parent',plotwda,'Box','on');
clearSingleAxis(gui.axes_handles.Bpre);
clearSingleAxis(gui.axes_handles.alpha);
clearSingleAxis(gui.axes_handles.dadt);
clearSingleAxis(gui.axes_handles.wda);

%% 3D Bloch Sphere

% lab frame (L) and rot frame (R)
SphereBoxL = uix.VBox('Parent',gui.panels.Plot.Sphere,'Spacing',3,'Padding',3);
SphereBoxR = uix.VBox('Parent',gui.panels.Plot.Sphere,'Spacing',3,'Padding',3);
gui.panels.Plot.Sphere.TabTitles = {'Lab frame','Rot frame'};
gui.panels.Plot.Sphere.TabWidth = 75;

% add view buttons and axes to the panel
plotSphereL = uicontainer('Parent',SphereBoxL);
SphereButtonsL = uix.HButtonBox('Parent',SphereBoxL);
gui.push_handles.XZL = uicontrol('Parent',SphereButtonsL,...
    'String','XZ',...
    'Tag','XZ',...
    'FontSize',myui.fontsize,...
    'Callback',@onPushAxView);
gui.push_handles.YZL = uicontrol('Parent',SphereButtonsL,...
    'String','YZ',...
    'Tag','YZ',...
    'FontSize',myui.fontsize,...
    'Callback',@onPushAxView);
gui.push_handles.XYL = uicontrol('Parent',SphereButtonsL,...
    'String','XY',...
    'Tag','XY',...
    'FontSize',myui.fontsize,...
    'Callback',@onPushAxView);
gui.push_handles.DL = uicontrol('Parent',SphereButtonsL,...
    'String','3D',...
    'Tag','3D',...
    'FontSize',myui.fontsize,...
    'Callback',@onPushAxView);
set(SphereBoxL,'Heights',[-1 30]);
gui.axes_handles.SphereL = axes('Parent',plotSphereL,'Box','on');
view(3); axis equal;
clearSingleAxis(gui.axes_handles.SphereL);

% add view buttons and axes to the panel
plotSphereR = uicontainer('Parent',SphereBoxR);
SphereButtonsR = uix.HButtonBox('Parent',SphereBoxR);
gui.push_handles.XZR = uicontrol('Parent',SphereButtonsR,...
    'String','XZ',...
    'Tag','XZ',...
    'FontSize',myui.fontsize,...
    'Callback',@onPushAxView);
gui.push_handles.YZR = uicontrol('Parent',SphereButtonsR,...
    'String','YZ',...
    'Tag','YZ',...
    'FontSize',myui.fontsize,...
    'Callback',@onPushAxView);
gui.push_handles.XYR = uicontrol('Parent',SphereButtonsR,...
    'String','XY',...
    'Tag','XY',...
    'FontSize',myui.fontsize,...
    'Callback',@onPushAxView);
gui.push_handles.DR = uicontrol('Parent',SphereButtonsR,...
    'String','3D',...
    'Tag','3D',...
    'FontSize',myui.fontsize,...
    'Callback',@onPushAxView);
set(SphereBoxR,'Heights',[-1 30]);
gui.axes_handles.SphereR = axes('Parent',plotSphereR,'Box','on');
view(3); axis equal;
clearSingleAxis(gui.axes_handles.SphereR);

%% Pulse tab panel
plotPulseSetup = uix.HBox('Parent',gui.panels.Plot.Pulse);
plotPulseB = uicontainer('Parent',gui.panels.Plot.Pulse);
plotPulseFFT = uicontainer('Parent', gui.panels.Plot.Pulse);
gui.panels.Plot.Pulse.TabTitles = {'Pulse setup','Pulse','FFT'};
gui.panels.Plot.Pulse.TabWidth = 75;

plotPulseSetupF = uicontainer('Parent',plotPulseSetup);
plotPulseSetupI = uicontainer('Parent',plotPulseSetup);
gui.axes_handles.PulseSetupF = axes('Parent',plotPulseSetupF,'Box','on');
gui.axes_handles.PulseSetupI = axes('Parent',plotPulseSetupI,'Box','on');

gui.axes_handles.PulseB = axes('Parent',plotPulseB,'Box','on');
gui.axes_handles.PulseFFT = axes('Parent',plotPulseFFT,'Box','on');
clearSingleAxis(gui.axes_handles.PulseSetupF);
clearSingleAxis(gui.axes_handles.PulseSetupI);
clearSingleAxis(gui.axes_handles.PulseB);
clearSingleAxis(gui.axes_handles.PulseFFT);

% arrange the panels in a 2x2 grid
set(gui.right,'Widths',[-1 -1],'Heights',[-1 -1]);

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