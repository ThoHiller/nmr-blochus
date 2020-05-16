function [gui,myui] = BLOCHUS_createPanelBasic(data,gui,myui)
%BLOCHUS_createPanelBasic creates "Basic" settings panel
%
% Syntax:
%       [gui,myui] = BLOCHUS_createPanelBasic(gui,myui,data)
%
% Inputs:
%       data - figure data structure
%       gui - figure gui elements structure
%       myui - individual GUI settings structure
%
% Outputs:
%       gui
%       myui
%
% Example:
%       [gui,myui] = BLOCHUS_createPanelBasic(data,gui,myui)
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

%% create all boxes
gui.panels.Basic.VBox = uix.VBox('Parent', gui.panels.Basic.main,...
    'Spacing',3,'Padding',3);

setNuc = uix.HBox('Parent',gui.panels.Basic.VBox,'Spacing',3);
setB0 = uix.HBox('Parent',gui.panels.Basic.VBox,'Spacing',3);
setT1relax = uix.HBox('Parent',gui.panels.Basic.VBox,'Spacing',3);
setT2relax = uix.HBox('Parent',gui.panels.Basic.VBox,'Spacing',3);
setTsim = uix.HBox('Parent',gui.panels.Basic.VBox,'Spacing',3);
setM0 = uix.HBox('Parent',gui.panels.Basic.VBox,'Spacing',3);
setMinit = uix.HBox('Parent',gui.panels.Basic.VBox,'Spacing',3);

%% Nucleus & gyromagnetic ratio
gui.text_handles.Nuc = uicontrol('Style','Text',...
    'Parent',setNuc,...
    'String',['Nucleus | ',char(hex2dec('3B3')),' [rad/s/T]'],...
    'FontSize',myui.fontsize);
uix.Empty('Parent',setNuc);
tstr = 'Proton for spin simulation.';
gui.popup_handles.Nuc = uicontrol('Style', 'Popup',...
    'Parent',setNuc,...
    'String',{'1H','2H','3He','7Li','13C','14N','15N','17O','19F',...
    '23Na','27Al','29Si','31P','57Fe','63Cu','67Zn','129Xe'},...
    'Tag','Nuc',...
    'TooltipString',tstr,...
    'FontSize',myui.fontsize,...
    'UserData',struct('Tooltipstr',tstr),...
    'Callback',@onPopupNucleus);
gui.edit_handles.Gyro = uicontrol('Style','Edit',...
    'Parent',setNuc,...
    'String',sprintf('%5.4e',data.basic.gamma),...
    'Tag','Gyro',...
    'FontSize',myui.fontsize,...
    'Enable','off');
set(setNuc,'Widths',[120 30 -1 -1]);

%% B0 & omega0
gui.text_handles.B0 = uicontrol('Style','Text',...
    'Parent',setB0,...
    'String',['B0 [T] | ',char(hex2dec('3C9')),'0 [Hz]'],...
    'FontSize',myui.fontsize);
uix.Empty('Parent',setB0);
tstr = '<HTML><b>B</b><sub>0</sub> field strength in [T].<br>';
gui.edit_handles.B0 = uicontrol('Style','Edit',...
    'Parent', setB0,...
    'String',num2str(data.init.B0(1)),...
    'Tag','basic_B0',...
    'TooltipString',tstr,...
    'FontSize',myui.fontsize,...
    'UserData',struct('Tooltipstr',tstr,'defaults',data.init.B0),...    
    'Callback',@onEditValue);
tstr = ['<HTML>Larmor frequency <b>',char(hex2dec('3C9')),'</b><sub>0</sub> in [Hz].<br>'];
gui.edit_handles.Omega0 = uicontrol('Style','Edit',...
    'Parent',setB0,...
    'String',sprintf('%7.2f',data.init.Omega0(1)),...
    'Tag','basic_Omega0',...
    'TooltipString',tstr,...
    'FontSize',myui.fontsize,...
    'UserData',struct('Tooltipstr',tstr,'defaults',data.init.Omega0),... 
    'Callback',@onEditValue);
set(setB0,'Widths',[120 30 -1 -1]);

%% T1 relaxation time
gui.text_handles.T1relax = uicontrol('Style','Text',...
    'Parent',setT1relax,...
    'String','T1 [ms]',...
    'FontSize',myui.fontsize);
uix.Empty('Parent',setT1relax);
tstr = '<HTML>Longitudinal relaxation time T<sub>1</sub> in [ms].<br>';
gui.edit_handles.T1relax = uicontrol('Style','Edit',...
    'Parent',setT1relax,...
    'String',num2str(data.init.T1relax(1)),...
    'Tag','basic_T1relax',...
    'TooltipString',tstr,...
    'FontSize',myui.fontsize,...
    'UserData',struct('Tooltipstr',tstr,'defaults',data.init.T1relax),... 
    'Callback',@onEditValue);
set(setT1relax,'Widths',[120 30 -1]);

%% T2 relaxation time
gui.text_handles.T2relax = uicontrol('Style','Text',...
    'Parent',setT2relax,...
    'String','T2 [ms]',...
    'FontSize',myui.fontsize);
uix.Empty('Parent',setT2relax);
tstr = '<HTML>Transversal relaxation time T<sub>2</sub> in [ms].<br>';
gui.edit_handles.T2relax = uicontrol('Style','Edit',...
    'Parent',setT2relax,...
    'String',num2str(data.init.T2relax(1)),...
    'Tag','basic_T2relax',...
    'TooltipString',tstr,...
    'FontSize',myui.fontsize,...
    'UserData',struct('Tooltipstr',tstr,'defaults',data.init.T2relax),... 
    'Callback',@onEditValue);
set(setT2relax,'Widths',[120 30 -1]);

%% maximum simulation time
gui.text_handles.Tsim = uicontrol('Style','Text',...
    'Parent',setTsim,...
    'String','T sim. [ms]',...
    'FontSize',myui.fontsize);
uix.Empty('Parent',setTsim);
tstr = '<HTML>Total simulation time T<sub>sim</sub> in [ms].<br>';
gui.edit_handles.Tsim = uicontrol('Style','Edit',...
    'Parent',setTsim,...
    'String',num2str(data.init.Tsim(1)),...
    'Tag','basic_Tsim',...
    'TooltipString',tstr,...
    'FontSize',myui.fontsize,...
    'UserData',struct('Tooltipstr',tstr,'defaults',data.init.Tsim),... 
    'Callback',@onEditValue);
set(setTsim,'Widths',[120 30 -1]);

%% equilibrium magnetization M0
gui.text_handles.M0 = uicontrol('Style','Text',...
    'Parent',setM0,...
    'String','M0 [A/m]',...
    'FontSize',myui.fontsize);
uix.Empty('Parent',setM0);
gui.edit_handles.M0x = uicontrol('Style','Edit',...
    'Parent',setM0,...
    'String',num2str(data.basic.M0(1)),...
    'Tag','basic_M0x',...
    'FontSize',myui.fontsize,....
    'Enable','off');
gui.edit_handles.M0y = uicontrol('Style','Edit',...
    'Parent',setM0,...
    'String',num2str(data.basic.M0(2)),...
    'Tag','basic_M0y',...
    'FontSize',myui.fontsize,...
    'Enable','off');
gui.edit_handles.M0z = uicontrol('Style','Edit',...
    'Parent',setM0,...
    'String',num2str(data.basic.M0(3)),...
    'Tag','basic_M0z',...
    'FontSize',myui.fontsize,...
    'Enable','off');
set(setM0,'Widths',[120 30 -1 -1 -1]);

%% initial magnetization Minit
gui.text_handles.Minit = uicontrol('Style','Text',...
    'Parent',setMinit,...
    'String','Minit [A/m]',...
    'FontSize',myui.fontsize);
uix.Empty('Parent',setMinit);
tstr = 'X-comp. of initial magnetization M in [A/m]';
gui.edit_handles.Minitx = uicontrol('Style','Edit',...
    'Parent',setMinit,...
    'String',num2str(data.basic.Minit(1)),...
    'Tag','basic_Minitx',...
    'TooltipString',tstr,...
    'FontSize',myui.fontsize,...
    'UserData',struct('Tooltipstr',tstr,'defaults',[data.basic.Minit(1) -1e6 1e6]),...
    'Callback',@onEditValue);
tstr = 'Y-comp. of initial magnetization M in [A/m]';
gui.edit_handles.Minity = uicontrol('Style','Edit',...
    'Parent',setMinit,...
    'String',num2str(data.basic.Minit(2)),...
    'Tag','basic_Minity',...
    'TooltipString',tstr,...
    'FontSize',myui.fontsize,...
    'UserData',struct('Tooltipstr',tstr,'defaults',[data.basic.Minit(2) -1e6 1e6]),...
    'Callback',@onEditValue);
tstr = 'Z-comp. of initial magnetization M in [A/m]';
gui.edit_handles.Minitz = uicontrol('Style','Edit',...
    'Parent',setMinit,...
    'String',num2str(data.basic.Minit(3)),...
    'Tag','basic_Minitz',...
    'TooltipString',tstr,...
    'FontSize',myui.fontsize,...
    'UserData',struct('Tooltipstr',tstr,'defaults',[data.basic.Minit(3) -1e6 1e6]),...
    'Callback',@onEditValue);
set(setMinit,'Widths',[120 30 -1 -1 -1]);

% Java Hack to adjust the text fields vertical alignment
jh = findjobj(gui.text_handles.Nuc);
jh.setVerticalAlignment(javax.swing.JLabel.CENTER)
jh = findjobj(gui.text_handles.B0);
jh.setVerticalAlignment(javax.swing.JLabel.CENTER)
jh = findjobj(gui.text_handles.T1relax);
jh.setVerticalAlignment(javax.swing.JLabel.CENTER)
jh = findjobj(gui.text_handles.T2relax);
jh.setVerticalAlignment(javax.swing.JLabel.CENTER)
jh = findjobj(gui.text_handles.Tsim);
jh.setVerticalAlignment(javax.swing.JLabel.CENTER)
jh = findjobj(gui.text_handles.M0);
jh.setVerticalAlignment(javax.swing.JLabel.CENTER)
jh = findjobj(gui.text_handles.Minit);
jh.setVerticalAlignment(javax.swing.JLabel.CENTER)

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