function [gui,myui] = BLOCHUS_createPanelPrePol(data,gui,myui)
%BLOCHUS_createPanelPrePol creates "Pre-polarization" settings panel
%
% Syntax:
%       [gui,myui] = BLOCHUS_createPanelPrePol(data,gui,myui)
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
%       [gui,myui] = BLOCHUS_createPanelPrePol(data,gui,myui)
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
gui.panels.PrePol.VBox = uix.VBox('Parent', gui.panels.PrePol.main,...
    'Spacing',3,'Padding',3);

setCheck = uix.HBox('Parent',gui.panels.PrePol.VBox,'Spacing',3);
setPrePolRDS = uix.HBox('Parent',gui.panels.PrePol.VBox,'Spacing',3);
setPrePolRamp = uix.HBox('Parent',gui.panels.PrePol.VBox,'Spacing',3);
setPrePolFactor = uix.HBox('Parent',gui.panels.PrePol.VBox,'Spacing',3);
setPrePolTheta = uix.HBox('Parent',gui.panels.PrePol.VBox,'Spacing',3);
setPrePolPhi = uix.HBox('Parent',gui.panels.PrePol.VBox,'Spacing',3);
setPrePolSwitchFactor = uix.HBox('Parent',gui.panels.PrePol.VBox,'Spacing',3);
setPrePolTramp = uix.HBox('Parent',gui.panels.PrePol.VBox,'Spacing',3);
setPrePolTslope = uix.HBox( 'Parent', gui.panels.PrePol.VBox,'Spacing',3);

%%
uix.Empty('Parent',setCheck);
tstr = 'Activate pre-polarization switch-off ramp settings.';
gui.check_handles.PrePol = uicontrol('Style','Checkbox',...
    'Parent',setCheck,...
    'String','Use Pre-polarization switch-off',...
    'Value',0,...
    'Tag','checkPrePol',...
    'ToolTipString',tstr,...
    'FontSize',myui.fontsize,...
    'UserData',struct('Tooltipstr',tstr),...
    'Callback',@onCheckPrePolPulse);
uix.Empty('Parent',setCheck);
set(setCheck,'Widths',[-1 200 -1]);

%% relaxation during switch-off RDS switch
gui.text_handles.PrePolRDS = uicontrol('Style','Text',...
    'Parent',setPrePolRDS,...
    'String','Relax. during switch-off',...
    'FontSize',myui.fontsize);
uix.Empty('Parent',setPrePolRDS);
tstr = '<HTML>Switch on/off relaxation during switch-off - RDS.<br>';
gui.check_handles.PrePolRDS = uicontrol('Style','Checkbox',...
    'Parent',setPrePolRDS,...
    'String','RDS',...
    'Tag','PrePolRDS',...
    'ToolTipString',tstr,...
    'FontSize',myui.fontsize,...
    'UserData',struct('Tooltipstr',tstr),...
    'Enable','off',...
    'Callback',@onCheckRDP);
set(setPrePolRDS,'Widths',[140 10 -1]);

%% prepolarization switch-off ramp
gui.text_handles.PrePolRamp = uicontrol('Style','Text',...
    'Parent',setPrePolRamp,...
    'String','Ramp type',...
    'FontSize',myui.fontsize);
uix.Empty('Parent', setPrePolRamp);
tstr = ['<HTML>Choose between different switch-off ramps.<br><br>',...
    '<u>Available options:</u><br>',...
    '<b>exp.</b> Pure exponential.<br>',...
    '<b>linear & exp.</b> Lin + exp parts as introduced by Conradi et al. (2017).<br>',...
    '<b>half cos</b> 1+cos type ramp.<br>',...
    '<b>linear</b> Pure linear.<br>',...
    '<b>melton1995</b> Linear using switch-off <it>rates</it> after Melton et al. (1995).<br>',...
    '<b>MIDI</b> 1+cos type ramp with NMR Midi settings.<br><br>',...
    'Depending on the chosen method some edit fields have different meanings.<br><br>',...
    '<u>Default value:</u><br>',...
    '<b>exp</b><br>'];
gui.popup_handles.PrePolRamp = uicontrol('Style','Popup',...
    'Parent',setPrePolRamp,...
    'String',{'exp.','linear & exp.','half cos','linear','melton1995','MIDI'},...
    'Value',1,...
    'Tag','Ramp',...
    'TooltipString',tstr,...
    'FontSize',myui.fontsize,...
    'UserData',struct('Tooltipstr',tstr),...
    'Enable','off',...
    'Callback',@onPopupPrePolRamp);
set(setPrePolRamp,'Widths',[140 10 -1]);

%% pre-polarization factor
gui.text_handles.PrePolFactor = uicontrol('Style','Text',...
    'Parent',setPrePolFactor,...
    'String','PrePol. factor [B0]',...
    'FontSize',myui.fontsize);
uix.Empty('Parent', setPrePolFactor);
tstr = '<HTML>Pre-polarization factor in units of <b>B</b><sub>0</sub>.<br>';
gui.edit_handles.PrePolFactor = uicontrol('Style','Edit',...
    'Parent',setPrePolFactor,...
    'String',num2str(data.init.PrePolFactor(1)),...
    'Tag','prepol_Factor',...
    'TooltipString',tstr,...
    'FontSize',myui.fontsize,...
    'UserData',struct('Tooltipstr',tstr,'defaults',data.init.PrePolFactor),...
    'Enable','off',...
    'Callback',@onEditValue);
set(setPrePolFactor,'Widths',[140 10 -1]);

%% angle theta between B0 and Bpre
gui.text_handles.PrePolTheta = uicontrol('Style','Text',...
    'Parent',setPrePolTheta,...
    'String','Theta [deg]',...
    'FontSize',myui.fontsize);
uix.Empty('Parent', setPrePolTheta);
tstr = ['<HTML>Polar orientation of <b>B</b><sub>p</sub> in [deg].<br>',...
    'Rotation is around the y-axis [0 1 0].'];
gui.edit_handles.PrePolTheta = uicontrol('Style','Edit',...
    'Parent',setPrePolTheta,...
    'String',num2str(data.init.PrePolTheta(1)),...
    'Tag','prepol_Theta',...
    'TooltipString',tstr,...
    'FontSize',myui.fontsize,...
    'UserData',struct('Tooltipstr',tstr,'defaults',data.init.PrePolTheta),...
    'Enable','off',...
    'Callback',@onEditValue);
set(setPrePolTheta,'Widths',[140 10 -1]);

%% angle phi (azimuthal angle)
gui.text_handles.PrePolPhi = uicontrol('Style','Text',...
    'Parent',setPrePolPhi,...
    'String','Phi [deg]',...
    'FontSize',myui.fontsize);
uix.Empty('Parent', setPrePolPhi);
tstr = ['<HTML>Azimuthal orientation of <b>B</b><sub>p</sub> in [deg].<br>',...
    'Rotation is around the z-axis [0 0 1].'];
gui.edit_handles.PrePolPhi = uicontrol('Style','Edit',...
    'Parent',setPrePolPhi,...
    'String',num2str(data.init.PrePolPhi(1)),...
    'Tag','prepol_Phi',...
    'TooltipString',tstr,...
    'FontSize',myui.fontsize,...
    'UserData',struct('Tooltipstr',tstr,'defaults',data.init.PrePolPhi),...
    'Enable','off',...
    'Callback',@onEditValue);
set(setPrePolPhi,'Widths',[140 10 -1]);

%% pre-polarization switch factor for the 'linexp' and 'exp' ramps
gui.text_handles.PrePolSwitchFactor = uicontrol('Style','Text',...
    'Parent',setPrePolSwitchFactor,...
    'String','Switch factor [B0]',...
    'FontSize',myui.fontsize);
uix.Empty('Parent', setPrePolSwitchFactor);
tstr = ['<HTML><BODY>This field has different meanings depending on the chosen',...
    ' shut-off ramp:<br><br>',...
    '<u>exp:</u><br>',...
    'Switch factor <b>B</b>* in units of [<b>B</b><sub>0</sub>] as part of the',...
    ' exponential decay exp(-t <b>B</b><sub>max</sub> / <b>B</b>* / T<sub>slope</sub>)'...
    '.<br><br>',...
    '<u>linear & exp:</u><br>',...
    'Switch factor <b>B</b>* in units of [<b>B</b><sub>0</sub>] where the ramp',...
    ' changes from linear to exp.<br>',...
    '<u>half cos:</u><br>',...
    'Not needed.<br><br>',...
    '<u>linear:</u><br>',...
    'Not needed.<br><br>',...
    '<u>melton1995:</u><br>',...
    'Cutoff rate k / &#915, with k the pre-polarization factor and &#915=&omega<sub>',...
    '0</sub> / T<sub>ramp</sub><br><br>',...
    '<u>MIDI:</u><br>',...
    'Not needed.</BODY></HTML>'];
gui.edit_handles.PrePolSwitchFactor = uicontrol('Style','Edit',...
    'Parent',setPrePolSwitchFactor,...
    'String',num2str(data.init.PrePolSwitchFactor(1)),...
    'Tag','prepol_SwitchFactor',...
    'TooltipString',tstr,...
    'FontSize',myui.fontsize,...
    'UserData',struct('Tooltipstr',tstr,'defaults',data.init.PrePolSwitchFactor),...
    'Enable','off',...    
    'Callback',@onEditValue);
set(setPrePolSwitchFactor,'Widths',[140 10 -1]);

%% pre-polarization Tramp time
gui.text_handles.PrePolTramp = uicontrol('Style','Text',...
    'Parent',setPrePolTramp,...
    'String','T ramp [ms]',...
    'FontSize',myui.fontsize);
uix.Empty('Parent', setPrePolTramp);
tstr = '<HTML>Switch-off ramp time T<sub>ramp</sub> in [ms]';
gui.edit_handles.PrePolTramp = uicontrol('Style','Edit',...
    'Parent',setPrePolTramp,...
    'String',num2str(data.init.PrePolTramp(1)),...
    'Tag','prepol_Tramp',...
    'TooltipString',tstr,...
    'FontSize',myui.fontsize,...
    'UserData',struct('Tooltipstr',tstr,'defaults',data.init.PrePolTramp),...
    'Enable','off',...    
    'Callback',@onEditValue);
set(setPrePolTramp, 'Widths', [140 10 -1]);

%% pre-polarization Tslope parameter for e.g. the 'linexp' ramp
gui.text_handles.PrePolTslope = uicontrol('Style','Text',...
    'Parent',setPrePolTslope,...
    'String','T slope [ms]',...
    'FontSize',myui.fontsize);
uix.Empty('Parent', setPrePolTslope);
tstr = ['<HTML><BODY>This field has different meanings depending on the chosen',...
    ' switch-off ramp:<br><br>',...
    '<u>exp:</u><br>',...
    'T<sub>slope</sub> in [ms] as part of the exponential decay exp(-t <b>B</b>',...
    '<sub>max</sub> / <b>B</b>* / T<sub>slope</sub>).<br><br>',...
    '<u>linear & exp:</u><br>',...
    'Slope T<sub>slope</sub> in [ms] of the linear part of the switch-off ramp.<br>',...
    '<u>half cos:</u><br>',...
    'Not needed.<br><br>',...
    '<u>linear:</u><br>',...
    'Not needed.<br><br>',...
    '<u>melton1995:</u><br>',...
    'Not needed.<br><br>',...
    '<u>MIDI:</u><br>',...
    'Not needed.</BODY></HTML>'];
gui.edit_handles.PrePolTslope = uicontrol('Style','Edit',...
    'Parent',setPrePolTslope,...
    'String',num2str(data.init.PrePolTslope(1)),...
    'Tag','prepol_Tslope',...
    'TooltipString',tstr,...
    'FontSize',myui.fontsize,...
    'UserData',struct('Tooltipstr',tstr,'defaults',data.init.PrePolTslope),...
    'Enable','off',...    
    'Callback',@onEditValue);
set(setPrePolTslope,'Widths',[140 10 -1]);

%% Java Hack to adjust the text fields vertical alignment
jh = findjobj(gui.text_handles.PrePolRDS);
jh.setVerticalAlignment(javax.swing.JLabel.CENTER)
jh = findjobj(gui.text_handles.PrePolRamp);
jh.setVerticalAlignment(javax.swing.JLabel.CENTER)
jh = findjobj(gui.text_handles.PrePolFactor);
jh.setVerticalAlignment(javax.swing.JLabel.CENTER)
jh = findjobj(gui.text_handles.PrePolTheta);
jh.setVerticalAlignment(javax.swing.JLabel.CENTER)
jh = findjobj(gui.text_handles.PrePolPhi);
jh.setVerticalAlignment(javax.swing.JLabel.CENTER)
jh = findjobj(gui.text_handles.PrePolSwitchFactor);
jh.setVerticalAlignment(javax.swing.JLabel.CENTER)
jh = findjobj(gui.text_handles.PrePolTramp);
jh.setVerticalAlignment(javax.swing.JLabel.CENTER)
jh = findjobj(gui.text_handles.PrePolTslope);
jh.setVerticalAlignment(javax.swing.JLabel.CENTER)

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