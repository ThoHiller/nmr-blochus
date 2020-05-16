function onPushAnimate(src,~)
%onPushAnimate animates the results of the latest simulation
%
% Syntax:
%       onPushAnimate(src)
%
% Inputs:
%       src - handle of the calling object
%
% Outputs:
%       none
%
% Example:
%       onPushAnimate(src)
%
% Other m-files required:
%       clearSingleAxis
%       plotBSphere
%       plotPulse
%       plotRamp
%       plotResults
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
    myui = gui.myui;
    
    % change the pushbutton color to indicate that an animation is running
    set(src,'BackGroundColor','r');
    % update the status info
    set(gui.text_handles.Status,'String','Animation running ...');
    pause(0.01);
    
    % ask the user for the number of frames to show    
    % default is 100
    prompt = {'\fontsize{9} Number of frames (default is 100):'};
    title = 'Number of Frames';
    dims = [1 50];
    definput = {'100'};
    opts.Interpreter = 'tex';
    answer = inputdlg(prompt,title,dims,definput);
    if ~isempty(answer)
        NoFrames = str2double(answer{1});
        if isnan(NoFrames) || NoFrames<1
            NoFrames = 100;
        end
    else
         NoFrames = 100;
    end    
    
    switch data.basic.type        
        case 'std'
            % the two important timings
            Tsim = data.basic.Tsim;
            
            % get the data to animate
            T = data.results.basic.T.*1e3;
            M = data.results.basic.M;
            Mrot = data.results.basic.Mrot;
            Mxy = sqrt(M(:,1).^2+M(:,2).^2);
            Mamp = sqrt(M(:,1).^2+M(:,2).^2+M(:,3).^2);
            
            % get the visible magnetization axis and data
            switch get(gui.panels.Plot.Mag,'Selection')
                case 1
                    ax1 = gui.axes_handles.MagL;
                    M1 = M;
                case 2
                    ax1 = gui.axes_handles.MagR;
                    M1 = Mrot;
                case 3
                    ax1 = gui.axes_handles.MagL;
                    % auto set FFT to lab frame
                    set(gui.panels.Plot.Mag,'Selection',1)
            end
            Mylim = get(ax1,'YLim');
            % get the visible Bloch sphere axis and data
            switch get(gui.panels.Plot.Sphere,'Selection')
                case 1
                    ax2 = gui.axes_handles.SphereL;
                    M2 = M;
                case 2
                    ax2 = gui.axes_handles.SphereR;
                    M2 = Mrot;
            end            
            Mview = get(ax2,'View');
            
			% time vector
            t = linspace(0,Tsim,NoFrames);
            
            % initialize all axes and animation lines
            axes(ax1);
            cla(ax1);
            hMx = animatedline('Color','r','LineWidth',2,'Parent',ax1);
            hMy = animatedline('Color','g','LineWidth',2,'Parent',ax1);
            hMz = animatedline('Color','b','LineWidth',2,'Parent',ax1);
            hMxy = animatedline('Color','m','LineWidth',2,'Parent',ax1);
            hMn = animatedline('Color','k','LineWidth',2,'LineStyle','--','Parent',ax1);
            set(get(ax1,'XLabel'),'String','t [ms]');
            set(get(ax1,'YLabel'),'String','magnetization M');
            grid(ax1,'on');
            set(ax1,'FontSize',myui.axfontsize);
            
            axes(ax2);
            cla(ax2);
            hold(ax2,'on');
            plotBSphere(18,18,ax2);
            h3dM = animatedline('Color',myui.color.basic,'LineWidth',2,'Parent',ax2);
            
            % loop over the data
            for nf = 1:NoFrames-1
                
                % get the index vector of the current data to plot
                indt = T>=t(nf) & T<t(nf+1);
                % current time vector
                x = T(indt);
                
                % get data to plot for M
                % left panel
                yMx = M1(indt,1);
                yMy = M1(indt,2);
                yMz = M1(indt,3);
                yMxy = Mxy(indt,1);
                yMn = Mamp(indt,1);

                % plot data into axes
                addpoints(hMx,x,yMx);
                addpoints(hMy,x,yMy);
                addpoints(hMz,x,yMz);
                addpoints(hMxy,x,yMxy);
                addpoints(hMn,x,yMn);
                
                legend(ax1,'x','M','z','Location','SouthWest');
                set(ax1,'XLim',[min(T) max(T)],'YLim',Mylim);
                lh = legend(ax1,'x','y','z','|xy|','|M|','Location','SouthWest');
                
                % Bloch sphere
                yMxx = M2(indt,1);
                yMyy = M2(indt,2);
                yMzz = M2(indt,3);
                
                addpoints(h3dM,yMxx,yMyy,yMzz);
                view(ax2,Mview);
                set(ax2,'Color','w','XColor','none','YColor','none','ZColor','none');
                axis(ax2,'equal');
                axis(ax2,'tight');
                
                % pause a bit to show the animation
                % drawnow limitrate works also but is a bit too fast
                pause(1/24);
                % drawnow limitrate;
            end
            hold(ax1,'off');
            hold(ax2,'off');
            
        case 'prepol'
            % get the results data
            basic = data.results.basic;
            prepol = data.results.prepol;
            
            % the two important timings
            Tsim = data.basic.Tsim;
            Tramp = data.prepol.Tramp;
            
            % get the data to animate
            T1 = basic.T.*1e3;
            M = basic.M;
            Mxy = sqrt(M(:,1).^2+M(:,2).^2);
            Mamp = sqrt(M(:,1).^2+M(:,2).^2+M(:,3).^2);
            M = M./data.basic.B0;
            Mxy = Mxy./data.basic.B0;
            Mamp = Mamp./data.basic.B0;
            % switch-off B-field trajectory
            Bamp = sqrt(prepol.Beff(:,1).^2+prepol.Beff(:,2).^2+prepol.Beff(:,3).^2);
            Beffn = prepol.Beff./Bamp;
            
            % set the M axes handles
            set(gui.panels.Plot.Mag,'Selection',1);
            set(gui.panels.Plot.Sphere,'Selection',1);
            ax1 = gui.axes_handles.MagL;
            ax2 = gui.axes_handles.SphereL;
            Mylim = get(ax1,'YLim');
            Mview = get(ax2,'View');
            
            % get the visible PrePol Ramp parameter axes and data
            % only this needs to be shown
            T2 = prepol.t.*1e3;
            switch get(gui.panels.Plot.PrePol,'Selection')
                case 1
                    ax3 = gui.axes_handles.Bpre;
                    Bval = prepol.Bp./data.basic.B0;
                    ylab = 'Bp [B0]';
                case 2
                    ax3 = gui.axes_handles.alpha;
                    Bval = rad2deg(prepol.alpha);
                    ylab = '\alpha [deg]';
                case 3
                    ax3 = gui.axes_handles.dadt;
                    Bval = prepol.dadt;
                    ylab = 'd\alpha / dt';
                case 4
                    ax3 = gui.axes_handles.wda;
                    Bval = prepol.dadt./prepol.omega(1:numel(prepol.dadt));
                    ylab = '(d\alpha/dt) / \gammaB';
            end
            Bylim = get(ax3,'YLim');

            % time vector
            dt = Tsim(end)/NoFrames;
            % including all relevant time markers
            t = unique([0:dt:Tsim Tramp]);
            % adapt the number frames
            NoFrames = numel(t);
            
            % initialize all axes and animation lines
            axes(ax1);
            cla(ax1);
            hMx = animatedline('Color','r','LineWidth',2,'Parent',ax1);
            hMy = animatedline('Color','g','LineWidth',2,'Parent',ax1);
            hMz = animatedline('Color','b','LineWidth',2,'Parent',ax1);
            hMxy = animatedline('Color','m','LineWidth',2,'Parent',ax1);
            hMn = animatedline('Color','k','LineWidth',2,'LineStyle','--','Parent',ax1);
            set(get(ax1,'XLabel'),'String','t [ms]');
            set(get(ax1,'YLabel'),'String','magnetization M');
            grid(ax1,'on');
            set(ax1,'FontSize',myui.axfontsize);
            
            axes(ax2);
            cla(ax2);
            hold(ax2,'on');
            plotBSphere(18,18,ax2);
            h3dMa = animatedline('Color',myui.color.prepol,'LineWidth',2,'Parent',ax2);
            h3dMb = animatedline('Color',myui.color.basic,'LineWidth',2,'Parent',ax2);
            h3dB = animatedline('Color',myui.color.prepolB,'LineWidth',2,'Parent',ax2);
            
            axes(ax3);
            clearSingleAxis(ax3);
            hB = animatedline('Color',myui.color.prepolB,'LineWidth',2,'Parent',ax3);
            set(get(ax3,'XLabel'),'String','t [ms]');
            set(get(ax3,'YLabel'),'String',ylab);
            grid(ax3,'on');
            set(ax3,'FontSize',myui.axfontsize);
            
            % loop over the data
            for nf = 1:NoFrames-1
                % get the index vector of the current data to plot
                indt = T1>=t(nf) & T1<t(nf+1);
                % get the current time vector for M
                x1 = T1(indt);
                % get data to plot for M
                yMx = M(indt,1);
                yMy = M(indt,2);
                yMz = M(indt,3);
                yMxy = Mxy(indt,1);
                yMn = Mamp(indt,1);
                
                % plot M into lab frame axes
                addpoints(hMx,x1,yMx);
                addpoints(hMy,x1,yMy);
                addpoints(hMz,x1,yMz);
                addpoints(hMxy,x1,yMxy);
                addpoints(hMn,x1,yMn);
                % M axis settings
                set(ax1,'XLim',[min(T1) max(T1)],'YLim',Mylim);
                lh = legend(ax1,'x','y','z','|xy|','|M|','Location','SouthWest');
                
                % as long as we are in the "ramp" phase gather the data of
                % the effective B-field to plot as a trace on the Bloch
                % sphere and the ramp parameter that is plotted into the
                % prepol-panel
                if t(nf) < Tramp
                    indt2 = T2>=t(nf) & T2<t(nf+1);
                    yBx = Beffn(indt2,1);
                    yBy = Beffn(indt2,2);
                    yBz = Beffn(indt2,3);
                    yBramp = Bval(indt2,1);
                    
                    % add the M data to the sphere (prepol color)
                    addpoints(h3dMa,yMx./Mamp(1),yMy./Mamp(1),yMz./Mamp(1));
                    % add the B-field trace to the sphere
                    addpoints(h3dB,yBx,yBy,yBz);
                else
                    % after the ramp is over
                    % add the M data to the sphere (basic color)
                    addpoints(h3dMb,yMx./Mamp(1),yMy./Mamp(1),yMz./Mamp(1));
                end
                % sphere axis settings
                view(ax2,Mview);
                set(ax2,'XLim',[-1.05 1.05],'YLim',[-1.05 1.05],'ZLim',[-1.05 1.05])
                set(ax2,'Color','w','XColor','none','YColor','none','ZColor','none');
                axis(ax2,'equal');
                axis(ax2,'tight');
                
                % plot the ramp parameter
                if t(nf) < Tramp
                    x2 = T2(indt2);
                    addpoints(hB,x2,yBramp);
                    set(ax3,'XLim',[min(T2) max(T2)],'YLim',Bylim);
                end
                
                % pause a bit to show the animation
                % drawnow limitrate works also but is a bit too fast
                pause(1/24);
                % drawnow limitrate;
                
            end
            hold(ax1,'off');
            hold(ax2,'off');
            hold(ax3,'off');
            
            % show the Ramp parameter again
            plotRamp(fig);
            
        case 'pulse'
            % get the results data
            basic = data.results.basic;
            pulse = data.results.pulse;
            
            % the two important timings
            Tsim = data.basic.Tsim;
            Ttau = data.pulse.Ttau;
            
            % get the data to animate
            T1 = basic.T.*1e3;
            M = basic.M;
            Mrot = basic.Mrot;
            Mxy = sqrt(M(:,1).^2+M(:,2).^2);
            Mamp = sqrt(M(:,1).^2+M(:,2).^2+M(:,3).^2);
            
            % get the visible magnetization axis and data
            switch get(gui.panels.Plot.Mag,'Selection')
                case 1
                    ax1 = gui.axes_handles.MagL;
                    M1 = M;
                case 2
                    ax1 = gui.axes_handles.MagR;
                    M1 = Mrot;
                case 3
                    ax1 = gui.axes_handles.MagL;
                    % auto set FFT to lab frame
                    set(gui.panels.Plot.Mag,'Selection',1)
            end
            Mylim = get(ax1,'YLim');
            % get the visible Bloch sphere axis and data
            switch get(gui.panels.Plot.Sphere,'Selection')
                case 1
                    ax2 = gui.axes_handles.SphereL;
                    M2 = M;
                case 2
                    ax2 = gui.axes_handles.SphereR;
                    M2 = Mrot;
            end
            Mview = get(ax2,'View');
            
            % get the visible Pulse parameter axes and data
            % only this needs to be shown
            T2 = pulse.t;
            % switch for df/I modulation axes
            showdual = false;
            switch get(gui.panels.Plot.Pulse,'Selection')
                case 1
                    ax3a = gui.axes_handles.PulseSetupF;
                    ax3b = gui.axes_handles.PulseSetupI;
                    P1 = [pulse.df pulse.I];
                    ylaba = 'df [Hz]';
                    ylabb = 'I [A]';
                    Pylima = get(ax3a,'YLim');
                    Pylimb = get(ax3b,'YLim');
                    showdual = true;
                case {2,3}
                    ax3 = gui.axes_handles.PulseB;
                    P1 = pulse.Bxy./data.basic.B0;
                    ylab = 'B_1 [B_0]';
                    % if FFT is visible, set axes to pulse amplitude
                    set(gui.panels.Plot.Pulse,'Selection',2)
                    Pylim = get(ax3,'YLim');
            end
            
            % time vector
            dt = Tsim(end)/NoFrames;
            % including all relevant time markers
            t = unique([0:dt:Tsim Ttau]);
            % adapt the number frames
            NoFrames = numel(t);
            
            % initialize all axes and animation lines
            axes(ax1);
            cla(ax1);
            hMx = animatedline('Color','r','LineWidth',2,'Parent',ax1);
            hMy = animatedline('Color','g','LineWidth',2,'Parent',ax1);
            hMz = animatedline('Color','b','LineWidth',2,'Parent',ax1);
            hMxy = animatedline('Color','m','LineWidth',2,'Parent',ax1);
            hMn = animatedline('Color','k','LineWidth',2,'LineStyle','--','Parent',ax1);
            set(get(ax1,'XLabel'),'String','t [ms]');
            set(get(ax1,'YLabel'),'String','magnetization M');
            grid(ax1,'on');
            set(ax1,'FontSize',myui.axfontsize);
            
            axes(ax2);
            cla(ax2);
            hold(ax2,'on');
            plotBSphere(18,18,ax2);
            h3dMa = animatedline('Color',myui.color.pulse,'LineWidth',2,'Parent',ax2);
            h3dMb = animatedline('Color',myui.color.basic,'LineWidth',2,'Parent',ax2);
            
            if showdual
                clearSingleAxis(ax3a);
                hPx = animatedline('Color',myui.color.pulse,'LineWidth',2,'Parent',ax3a);
                set(get(ax3a,'XLabel'),'String','t [ms]');
                set(get(ax3a,'YLabel'),'String',ylaba);
                grid(ax3a,'on');
                set(ax3a,'FontSize',myui.axfontsize);
                
                clearSingleAxis(ax3b);
                hPy = animatedline('Color',myui.color.pulse,'LineWidth',2,'Parent',ax3b);
                set(get(ax3b,'XLabel'),'String','t [ms]');
                set(get(ax3b,'YLabel'),'String',ylabb);
                grid(ax3b,'on');
                set(ax3b,'FontSize',myui.axfontsize);
            else
                clearSingleAxis(ax3);
                hPx = animatedline('Color','r','LineWidth',2,'Parent',ax3);
                hPy = animatedline('Color','g','LineWidth',2,'Parent',ax3);
                set(get(ax3,'XLabel'),'String','t [ms]');
                set(get(ax3,'YLabel'),'String',ylab);
                grid(ax3,'on');
                set(ax3,'FontSize',myui.axfontsize);
            end
            
            % loop over the data
            for nf = 1:NoFrames-1
                % get the index vector of the current data to plot
                indt = T1>=t(nf) & T1<t(nf+1);
                % get the current time vector for M
                x1 = T1(indt);
                % get data to plot for M
                yMx = M1(indt,1);
                yMy = M1(indt,2);
                yMz = M1(indt,3);
                yMxx = M2(indt,1);
                yMyy = M2(indt,2);
                yMzz = M2(indt,3);
                yMxy = Mxy(indt,1);
                yMn = Mamp(indt,1);
                
                % plot M into lab frame axes
                addpoints(hMx,x1,yMx);
                addpoints(hMy,x1,yMy);
                addpoints(hMz,x1,yMz);
                addpoints(hMxy,x1,yMxy);
                addpoints(hMn,x1,yMn);
                % M axis settings
                set(ax1,'XLim',[min(T1) max(T1)],'YLim',Mylim);
                lh = legend(ax1,'x','y','z','|xy|','|M|','Location','SouthWest');
                
                % as long as we are in the "pulse" phase
                if t(nf) < Ttau
                    % add the M data to the sphere (pulse color)
                    addpoints(h3dMa,yMxx,yMyy,yMzz);
                else
                    % after the ramp is over
                    % add the M data to the sphere (basic color)
                    addpoints(h3dMb,yMxx,yMyy,yMzz);
                end
                % sphere axis settings
                view(ax2,Mview);
                set(ax2,'XLim',[-1.05 1.05],'YLim',[-1.05 1.05],'ZLim',[-1.05 1.05])
                set(ax2,'Color','w','XColor','none','YColor','none','ZColor','none');
                axis(ax2,'equal');
                axis(ax2,'tight');
                
                % plot the Pulse parameter
                if t(nf) < Ttau
                    indt2 = T2>=t(nf) & T2<t(nf+1); % t2?
                    x2 = T2(indt2);
                    yBx = P1(indt2,1);
                    yBy = P1(indt2,2);
                    
                    addpoints(hPx,x2,yBx);
                    addpoints(hPy,x2,yBy);
                    if showdual
                        set(ax3a,'XLim',[min(T2) max(T2)],'YLim',Pylima);
                        lh = legend(ax3a,'df mod.','Location','SouthEast');
                        
                        set(ax3b,'XLim',[min(T2) max(T2)],'YLim',Pylimb);
                        lh = legend(ax3b,'I mod.','Location','SouthEast');
                    else
                        set(ax3,'XLim',[min(T2) max(T2)],'YLim',Pylim);
                        lh = legend(ax3,'x','y','Location','SouthWest');
                    end
                end
                
                % pause a bit to show the animation
                % drawnow limitrate works also but is a bit too fast
                pause(1/24);
                % drawnow limitrate;
                
            end
            hold(ax1,'off');
            hold(ax2,'off');
            if showdual
                hold(ax3a,'off');
                hold(ax3b,'off');
            else
                hold(ax3,'off');
            end
            % show the Pulse parameter again
            plotPulse(fig);
            
        case 'prepolpulse'
            % get the results data
            basic = data.results.basic;
            prepol = data.results.prepol;
            pulse = data.results.pulse;
            
            % the important timings
            Tsim = data.basic.Tsim;
            Tramp = data.prepol.Tramp;
            Ttau = data.pulse.Ttau;
            Twait = data.pulse.Twait;
            Trelax = Tsim-Tramp-Twait-Ttau;
            
            % get the data to animate
            T1 = basic.T.*1e3;
            M = basic.M;
            Mrot = basic.Mrot;
            Mxy = sqrt(M(:,1).^2+M(:,2).^2);
            Mamp = sqrt(M(:,1).^2+M(:,2).^2+M(:,3).^2);
            % normalize the magnetization
            M = M./data.basic.B0;
            Mrot = Mrot./data.basic.B0;
            Mxy = Mxy./data.basic.B0;
            Mamp = Mamp./data.basic.B0;
            % switch-off B-field trajectory
            Bamp = sqrt(prepol.Beff(:,1).^2+prepol.Beff(:,2).^2+prepol.Beff(:,3).^2);
            Beffn = prepol.Beff./Bamp;
            
            % get the visible magnetization axis and data
            switch get(gui.panels.Plot.Mag,'Selection')
                case 1
                    ax1 = gui.axes_handles.MagL;
                    M1 = M;
                case 2
                    ax1 = gui.axes_handles.MagR;
                    M1 = Mrot;
                case 3
                    ax1 = gui.axes_handles.MagL;
                    % auto set FFT to lab frame
                    set(gui.panels.Plot.Mag,'Selection',1)
            end
            Mylim = get(ax1,'YLim');            
            % get the visible Bloch sphere axis and data
            switch get(gui.panels.Plot.Sphere,'Selection')
                case 1
                    ax2 = gui.axes_handles.SphereL;
                    M2 = M;
                case 2
                    ax2 = gui.axes_handles.SphereR;
                    M2 = Mrot;
            end
            Mview = get(ax2,'View');
            
            % get the visible PrePol Ramp parameter axis and data
            T2 = prepol.t.*1e3;
            switch get(gui.panels.Plot.PrePol,'Selection')
                case 1
                    ax3 = gui.axes_handles.Bpre;
                    Bval = prepol.Bp./data.basic.B0;
                    ylab = 'Bp [B0]';
                case 2
                    ax3 = gui.axes_handles.alpha;
                    Bval = rad2deg(prepol.alpha);
                    ylab = '\alpha [deg]';
                case 3
                    ax3 = gui.axes_handles.dadt;
                    Bval = prepol.dadt;
                    ylab = 'd\alpha / dt';
                case 4
                    ax3 = gui.axes_handles.wda;
                    Bval = prepol.dadt./prepol.omega(1:numel(prepol.dadt));
                    Bval = [Bval;Bval(end)];
                    ylab = '(d\alpha/dt) / \gammaB';
            end
            Bylim = get(ax3,'YLim');
            
            % get the visible Pulse parameter axes and data
            T3 = pulse.t;
            % switch for df/I modulation axes
            showdual = false;
            switch get(gui.panels.Plot.Pulse,'Selection')
                case 1
                    ax4a = gui.axes_handles.PulseSetupF;
                    ax4b = gui.axes_handles.PulseSetupI;
                    P1 = [pulse.df pulse.I];
                    ylaba = 'df [Hz]';
                    ylabb = 'I [A]';
                    Pylima = get(ax4a,'YLim');
                    Pylimb = get(ax4b,'YLim');
                    showdual = true;
                case {2,3}
                    ax4 = gui.axes_handles.PulseB;
                    P1 = pulse.Bxy./data.basic.B0;
                    ylab = 'B_1 [B_0]';
                    % if FFT is visible, set axes to pulse amplitude
                    set(gui.panels.Plot.Pulse,'Selection',2)
                    Pylim = get(ax4,'YLim');
            end

            % time vector
            dt = Tsim(end)/NoFrames;
            % including all relevant time markers
            t = unique([0:dt:Tsim Tramp Tramp+Twait Tramp+Twait+Ttau]);
            % adapt the number frames
            NoFrames = numel(t);

            % initialize all axes and animation lines
            axes(ax1);
            cla(ax1);
            hMx = animatedline('Color','r','LineWidth',2,'Parent',ax1);
            hMy = animatedline('Color','g','LineWidth',2,'Parent',ax1);
            hMz = animatedline('Color','b','LineWidth',2,'Parent',ax1);
            hMxy = animatedline('Color','m','LineWidth',2,'Parent',ax1);
            hMn = animatedline('Color','k','LineWidth',2,'LineStyle','--','Parent',ax1);
            set(get(ax1,'XLabel'),'String','t [ms]');
            set(get(ax1,'YLabel'),'String','magnetization M');
            grid(ax1,'on');
            set(ax1,'FontSize',myui.axfontsize);
            
            axes(ax2);
            cla(ax2);
            hold(ax2,'on');
            plotBSphere(18,18,ax2);
            h3dMa = animatedline('Color',myui.color.prepol,'LineWidth',2,'Parent',ax2);
            h3dMb = animatedline('Color',myui.color.wait,'LineWidth',2,'Parent',ax2);
            h3dMc = animatedline('Color',myui.color.pulse,'LineWidth',2,'Parent',ax2);
            h3dMd = animatedline('Color',myui.color.basic,'LineWidth',2,'Parent',ax2);
            h3dB = animatedline('Color',myui.color.prepolB,'LineWidth',2,'Parent',ax2);
            
            axes(ax3);
            clearSingleAxis(ax3);
            hB = animatedline('Color',myui.color.prepolB,'LineWidth',2,'Parent',ax3);
            set(get(ax3,'XLabel'),'String','t [ms]');
            set(get(ax3,'YLabel'),'String',ylab);
            grid(ax3,'on');
            set(ax3,'FontSize',myui.axfontsize);
            
            if showdual
                clearSingleAxis(ax4a);
                hPx = animatedline('Color',myui.color.pulse,'LineWidth',2,'Parent',ax4a);
                set(get(ax4a,'XLabel'),'String','t [ms]');
                set(get(ax4a,'YLabel'),'String',ylaba);
                grid(ax4a,'on');
                set(ax4a,'FontSize',myui.axfontsize);
                
                clearSingleAxis(ax4b);
                hPy = animatedline('Color',myui.color.pulse,'LineWidth',2,'Parent',ax4b);
                set(get(ax4b,'XLabel'),'String','t [ms]');
                set(get(ax4b,'YLabel'),'String',ylabb);
                grid(ax4b,'on');
                set(ax4b,'FontSize',myui.axfontsize);
            else
                clearSingleAxis(ax4);
                hPx = animatedline('Color','r','LineWidth',2,'Parent',ax4);
                hPy = animatedline('Color','g','LineWidth',2,'Parent',ax4);
                set(get(ax4,'XLabel'),'String','t [ms]');
                set(get(ax4,'YLabel'),'String',ylab);
                grid(ax4,'on');
                set(ax4,'FontSize',myui.axfontsize);
            end
            
            % loop over the data
            for nf = 1:NoFrames-1
                % get the index vector of the current M data to plot
                indt = T1>=t(nf) & T1<t(nf+1);
                % get the current time vector for M
                x1 = T1(indt);
                % get data to plot for M
                yMx = M1(indt,1);
                yMy = M1(indt,2);
                yMz = M1(indt,3);
                yMxx = M2(indt,1)./Mamp(1);
                yMyy = M2(indt,2)./Mamp(1);
                yMzz = M2(indt,3)./Mamp(1);
                yMxy = Mxy(indt,1);
                yMn = Mamp(indt,1);
                
                % plot M into the left panel
                addpoints(hMx,x1,yMx);
                addpoints(hMy,x1,yMy);
                addpoints(hMz,x1,yMz);
                addpoints(hMxy,x1,yMxy);
                addpoints(hMn,x1,yMn);
                % M axis settings
                set(ax1,'XLim',[min(T1) max(T1)],'YLim',Mylim);
                lh = legend(ax1,'x','y','z','|xy|','|M|','Location','SouthWest');
                
                % plot into the Bloch sphere
                % as long as we are in the "pulse" phase
                if t(nf) < Tramp
                    indt2 = T2>t(nf) & T2<=t(nf+1); % t2?
                    yBx = Beffn(indt2,1);
                    yBy = Beffn(indt2,2);
                    yBz = Beffn(indt2,3);
                    yBramp = Bval(indt2,1);
                    
                    % add the M data to the sphere (prepol color)
                    addpoints(h3dMa,yMxx,yMyy,yMzz);
                    % add the B-field trace to the sphere
                    addpoints(h3dB,yBx,yBy,yBz);
                elseif t(nf)>Tramp && t(nf) <= Tramp+Twait
                    % add the M data to the sphere (wait color)
                    addpoints(h3dMb,yMxx,yMyy,yMzz);
                elseif t(nf)>Tramp+Twait && t(nf) <= Tsim-Trelax
                    % add the M data to the sphere (pulse color)
                    addpoints(h3dMc,yMxx,yMyy,yMzz);
                elseif t(nf) > Tsim-Trelax
                    % after the pulse is over
                    % add the M data to the sphere (basic color)
                    addpoints(h3dMd,yMxx,yMyy,yMzz);
                end
                % sphere axis settings
                view(ax2,Mview);
                set(ax2,'XLim',[-1.05 1.05],'YLim',[-1.05 1.05],'ZLim',[-1.05 1.05])
                set(ax2,'Color','w','XColor','none','YColor','none','ZColor','none');
                axis(ax2,'equal');
                axis(ax2,'tight');
                
                % plot the ramp parameter
                if t(nf) < Tramp
                    x2 = T2(indt2);
                    addpoints(hB,x2,yBramp);
                    set(ax3,'XLim',[min(T2) max(T2)],'YLim',Bylim);
                end
                
                % plot the Pulse parameter
                if t(nf)>=Tramp+Twait && t(nf) < Tsim-Trelax
                    indt3 = T3+Tramp+Twait>=t(nf) & T3+Tramp+Twait<t(nf+1);
                    x3 = T3(indt3);
                    yBx = P1(indt3,1);
                    yBy = P1(indt3,2);
                    
                    addpoints(hPx,x3,yBx);
                    addpoints(hPy,x3,yBy);
                    if showdual
                        set(ax4a,'XLim',[min(T3) max(T3)],'YLim',Pylima);
                        lh = legend(ax4a,'df mod.','Location','SouthEast');
                        
                        set(ax4b,'XLim',[min(T3) max(T3)],'YLim',Pylimb);
                        lh = legend(ax4b,'I mod.','Location','SouthEast');
                    else
                        set(ax4,'XLim',[min(T3) max(T3)],'YLim',Pylim);
                        lh = legend(ax4,'x','y','Location','SouthWest');
                    end
                end
                
                % pause a bit to show the animation
                % drawnow limitrate works also but is a bit too fast
                pause(1/24);
                % drawnow limitrate;
                
            end
            hold(ax1,'off');
            hold(ax2,'off');
            hold(ax3,'off');
            if showdual
                hold(ax4a,'off');
                hold(ax4b,'off');
            else
                hold(ax4,'off');
            end
            % show the Ramp parameter again
            plotRamp(fig);
            % show the Pulse parameter again
            plotPulse(fig);    
    end
    
    % reset the pushbutton color
    set(src,'BackGroundColor',[.94 .94 .94]);
    % update the status info
    set(gui.text_handles.Status,'String','Animation running ... finished.');
    % done
    pause(0.01);
    % show the results parameter again
    plotResults(fig);
    
else
    warndlg({'onPushAnimate:','There is no figure with the BLOCHUS Tag open.'},...
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
