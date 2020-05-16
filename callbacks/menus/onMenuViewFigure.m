function onMenuViewFigure(src,~)
%onMenuViewFigure shows predefined figure layouts
%
% Syntax:
%       onMenuViewFigure(src)
%
% Inputs:
%       src - handle of the calling object
%
% Outputs:
%       none
%
% Example:
%       onMenuViewFigure(src)
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
    
    % get GUI position
    posf = get(fig,'Position');
    % opening the export figure
    expfig = figure;
    
    % due to the two axes on the pulse modulation panel, we need to know
    % before creating the new figure if this panel is active
    isdual = false;
    if get(gui.panels.Plot.Pulse,'Selection') == 1
        isdual = true;
    end
    
    % create the axes layout on the export figure and get the axes
    % positions
    switch get(src,'Label')
        case 'Current View'
            % we copy all visible axes in a 2x2 grid
            ax1 = subplot(2,2,1,'Parent',expfig);
            ax2 = subplot(2,2,2,'Parent',expfig);
            ax3 = subplot(2,2,3,'Parent',expfig);
            pos1 = get(ax1,'Position');
            pos2 = get(ax2,'Position');
            pos3 = get(ax3,'Position');
            delete(ax1);
            delete(ax2);
            delete(ax3);
            if isdual
                ax4a = subplot(2,4,7,'Parent',expfig);
                ax4b = subplot(2,4,8,'Parent',expfig);
                pos4a = get(ax4a,'Position');
                pos4b = get(ax4b,'Position');
                delete(ax4a);
                delete(ax4b);
            else
                ax4 = subplot(2,2,4,'Parent',expfig);
                pos4 = get(ax4,'Position');
                delete(ax4);
            end
            
        case {'Magnetization','Switch-off Ramp'}
            % we copy all visible axes in a 2x2 grid
            ax1 = subplot(2,2,1,'Parent',expfig);
            ax2 = subplot(2,2,2,'Parent',expfig);
            ax3 = subplot(2,2,3,'Parent',expfig);
            ax4 = subplot(2,2,4,'Parent',expfig);
            pos1 = get(ax1,'Position');
            pos2 = get(ax2,'Position');
            pos3 = get(ax3,'Position');
            pos4 = get(ax4,'Position');
            delete(ax1);
            delete(ax2);
            delete(ax3);
            delete(ax4);
            
        case 'Pulse'
            % we copy pulse parameter in a 3x2 grid
            ax1 = subplot(3,2,1,'Parent',expfig);
            ax2 = subplot(3,2,2,'Parent',expfig);
            ax3 = subplot(3,2,[3 4],'Parent',expfig);
            ax4 = subplot(3,2,[5 6],'Parent',expfig);
            pos1 = get(ax1,'Position');
            pos2 = get(ax2,'Position');
            pos3 = get(ax3,'Position');
            pos4 = get(ax4,'Position');
            delete(ax1);
            delete(ax2);
            delete(ax3);
            delete(ax4);
    end
    
    % copy the GUI axes to the export figure
    switch get(src,'Label')
        case 'Current View'
            mag = get(gui.panels.Plot.Mag,'Selection');
            sph = get(gui.panels.Plot.Sphere,'Selection');
            pre = get(gui.panels.Plot.PrePol,'Selection');
            pul = get(gui.panels.Plot.Pulse,'Selection');
            switch mag
                case 1
                    ax1 = copyobj(gui.axes_handles.MagL,expfig);
                case 2
                    ax1 = copyobj(gui.axes_handles.MagR,expfig);
                case 3
                    ax1 = copyobj(gui.axes_handles.MagFFT,expfig);
            end
            switch sph
                case 1
                    ax2 = copyobj(gui.axes_handles.SphereL,expfig);
                case 2
                    ax2 = copyobj(gui.axes_handles.SphereR,expfig);
            end
            switch pre
                case 1
                    ax3 = copyobj(gui.axes_handles.Bpre,expfig);
                case 2
                    ax3 = copyobj(gui.axes_handles.alpha,expfig);
                case 3
                    ax3 = copyobj(gui.axes_handles.dadt,expfig);
                case 4
                    ax3 = copyobj(gui.axes_handles.wda,expfig);
            end
            switch pul
                case 1
                    ax4a = copyobj(gui.axes_handles.PulseSetupF,expfig);
                    ax4b = copyobj(gui.axes_handles.PulseSetupI,expfig);
                case 2
                    ax4 = copyobj(gui.axes_handles.PulseB,expfig);
                case 3
                    ax4 = copyobj(gui.axes_handles.PulseFFT,expfig);
            end
            
            set(expfig,'Name','BLOCHUS: Current View');
            
        case 'Magnetization'
            ax1 = copyobj(gui.axes_handles.MagL,expfig);
            ax2 = copyobj(gui.axes_handles.SphereL,expfig);
            ax3 = copyobj(gui.axes_handles.MagR,expfig);
            ax4 = copyobj(gui.axes_handles.SphereR,expfig);
            
            set(expfig,'Name','BLOCHUS: Magnetization');
            
        case 'Switch-off Ramp'
            ax1 = copyobj(gui.axes_handles.Bpre,expfig);
            ax2 = copyobj(gui.axes_handles.alpha,expfig);
            ax3 = copyobj(gui.axes_handles.dadt,expfig);
            ax4 = copyobj(gui.axes_handles.wda,expfig);
            
            set(expfig,'Name','BLOCHUS: Switch-off Ramp');
            
        case 'Pulse'
            ax1 = copyobj(gui.axes_handles.PulseSetupF,expfig);
            ax2 = copyobj(gui.axes_handles.PulseSetupI,expfig);
            ax3 = copyobj(gui.axes_handles.PulseB,expfig);
            ax4 = copyobj(gui.axes_handles.PulseFFT,expfig);
            
            set(expfig,'Name','BLOCHUS: Pulse');
    end
    
    % adjust the axes positions
    switch get(src,'Label')
        case 'Current View'
            set(ax1,'Position',pos1);
            set(ax2,'Position',pos2);
            set(ax3,'Position',pos3);
            if isdual
                set(ax4a,'Position',pos4a);
                set(ax4b,'Position',pos4b);
            else
                set(ax4,'Position',pos4);
            end
            
        otherwise
            set(ax1,'Position',pos1);
            set(ax2,'Position',pos2);
            set(ax3,'Position',pos3);
            set(ax4,'Position',pos4);
    end
    
    % adjust the position of the export figure
    set(expfig,'Position',[posf(1)+300 posf(2) (posf(3)-300)*0.8 posf(4)*0.8]);
    
    % show legends
    switch get(src,'Label')
        case 'Current View'
            lgh1 = legend(ax1,'show');
            lgh3 = legend(ax3,'show');
            if isdual
                lgh4a = legend(ax4a,'show');
                lgh4b = legend(ax4b,'show');
            else
                lgh4 = legend(ax4,'show');
            end
            
        case 'Magnetization'
            lgh1 = legend(ax1,'show');
            lgh3 = legend(ax3,'show');
            
        case 'Switch-off Ramp'
            lgh1 = legend(ax1,'show');
            lgh2 = legend(ax2,'show');
            lgh4 = legend(ax4,'show');
            set(lgh4,'Interpreter','latex','Location','best');
            
        case 'Pulse'
            lgh1 = legend(ax1,'show');
            lgh2 = legend(ax2,'show');
            lgh3 = legend(ax3,'show');
            lgh4 = legend(ax4,'show');
    end
    
else
    warndlg({'onMenuViewFigure:','There is no figure with the BLOCHUS Tag open.'},...
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
