function onMenuHelp(src,~)
%onMenuHelp shows the Help Information
%
% Syntax:
%       onMenuHelp(src)
%
% Inputs:
%       src - handle of the calling object
%
% Outputs:
%       none
%
% Example:
%       onMenuHelp(src)
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
    
    
    % BLOCHUS logo
    load('BLOCHUS_logo.mat','logo');
    
    % header info
    header{1,1}  = 'BLOCHUS:';
    header{end+1,1}  = ' ';
    header{end+1,1}  = ['author: ',gui.myui.author];
    header{end+1,1}  = ' ';
    header{end+1,1}  = ['version: ',gui.myui.version];
    header{end+1,1}  = ' ';
    header{end+1,1}  = ['date: ',gui.myui.date];
    header{end+1,1}  = ' ';
    
    % info text
    info{1,1} = ['BLOCHUS is a set of MATLAB tools, that allow some basic',...
        ' simulations of (S)NMR spin dynamics based on the Bloch equations.',...
        ' The Bloch equations are solved in the laboratory frame of reference',...
        ' with MATLABs built-in ode45 solver.',...
        ' Because it was developed within the scope of a near surface SNMR project,',...
        ' its main features are the simulation of',...
        ' (1) pre-polarization switch-off ramps and (2) excitation pulses.',...
        ' The main front-end to the underlying simulation tools is a',...
        ' graphical user interface (GUI)',...
        ' that allows playing around with the different features and helps',...
        ' to understand the basic concepts of (S)NMR spin dynamics.'];
    info{end+1,1} = ' ';
    info{end+1,1} = 'Have Fun!';
    
    % get BLOCHUS GUI position
    posf = get(fig,'Position');
    % default widht and height of About Figure
    ww = 560; hh = 420;
    xp = posf(1) + (posf(3)-ww)/2;
    yp = posf(2) + (posf(4)-hh)/2;
    % create Figure
    hf = figure('Name','About BLOCHUS',...
        'NumberTitle','off','Tag','Help','ToolBar','none','MenuBar','none',...
        'Resize','off','Position',[xp yp ww hh],'Visible','off');
    v1 = uix.VBox('Parent',hf,'Padding',10,'Spacing',10);
    
    % text area
    h1 = uix.VBox('Parent',v1);
    % button area
    h2 = uix.HBox('Parent',v1);
    set(v1,'Heights',[-1 30]);
    
    % text area
    h3 = uix.HBox('Parent',h1);
    % logo area
    h4 = uix.HBox('Parent',h1);
    set(h1,'Heights',[-1 -1]);
    
    % close button at the bottom
    uix.Empty('Parent',h2);
    p1 = uicontrol('Style','pushbutton','Parent',h2,'String','OK',...
        'FontSize',10,'Callback','closereq()');
    uix.Empty('Parent',h2);
    set(h2,'Widths',[-1 50 -1])
    
    % header
    uix.Empty('Parent',h3);
    t1 = uicontrol('Style','Text','Parent',h3,'String',header,...
        'FontSize',10,'HorizontalAlignment','left');
    % logo
    c1 = uicontainer('Parent',h3);
    ax1 = axes('Parent',c1);
    imshow(logo,'Parent',ax1);
    set(h3,'Widths',[50 -1 -1]);
    
    % info text
    uix.Empty('Parent',h4);
    t2 = uicontrol('Style','Text','Parent',h4,'String',info,...
        'FontSize',10,'HorizontalAlignment','left');
    uix.Empty('Parent',h4);
    set(h4,'Widths',[20 -1 20])
    
    % text hack
    jh = findjobj(t1);
    jh.setVerticalAlignment(javax.swing.JLabel.CENTER)
	
	% make the figure content visible
	set(hf,'Visible','on');
    
else
    warndlg({'onMenuHelp:','There is no figure with the BLOCHUS Tag open.'},...
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
