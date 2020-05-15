function onPushAnimate(src,~)
fig  = findobj('Tag','BLOCHUS');
gui  = getappdata(fig,'gui');
data = getappdata(fig,'data');

% change the pushbutton color to indicate that an animation is running
set(src,'BackGroundColor','r');
pause(0.01);

NoFrames = 100;
FPS = 24;

newblue = [0 0.4470 0.7410];
neworange = [0.8500 0.3250 0.0980];

switch data.basic.type
    
    case 'std'
        set(gui.text_handles.Status,'String','Animation running ...');
        pause(0.01);
        Mamp = zeros(numel(data.TT),1);
        for i=1:numel(data.TT)
            Mamp(i,1) = norm(data.MM(i,:)); % amplitude of M
        end
        
        dt = data.TT(end)/NoFrames;
        time = 0:dt:data.TT(end);
        
        axes(gui.axes_handles.MagL);
        cla(gui.axes_handles.MagL);
        h1 = animatedline('Color','r','LineWidth',2,'Parent',gui.axXYZL);
        h2 = animatedline('Color','g','LineWidth',2,'Parent',gui.axXYZL);
        h3 = animatedline('Color','b','LineWidth',2,'Parent',gui.axXYZL);
        h4 = animatedline('Color','k','LineWidth',2,'LineStyle','--','Parent',gui.axXYZL);
        
        axes(gui.axes_handles.MagR);
        cla(gui.axes_handles.MagR);
        h7 = animatedline('Color','r','LineWidth',2,'Parent',gui.axXYZR);
        h8 = animatedline('Color','g','LineWidth',2,'Parent',gui.axXYZR);
        h9 = animatedline('Color','b','LineWidth',2,'Parent',gui.axXYZR);
        h10 = animatedline('Color','k','LineWidth',2,'LineStyle','--','Parent',gui.axXYZR);
        
        axes(gui.axSphereL);
        cla(gui.axSphereL);
        hold(gui.axSphereL,'on');
        plotBSphere(18,18,gui.axSphereL);
        h5 = animatedline('Color',myui.color.basic,'LineWidth',2,'Parent',gui.axSphereL);
        
        axes(gui.axSphereR);
        cla(gui.axSphereR);
        hold(gui.axSphereR,'on');
        plotBSphere(18,18,gui.axSphereR); 
        h6 = animatedline('Color',myui.color.basic,'LineWidth',2,'Parent',gui.axSphereR);
        
        for nf = 1:NoFrames-1
            
            x1 = data.TT(data.TT>=time(nf) & data.TT<time(nf+1)).*1e3;
            y1 = data.MM(data.TT>=time(nf) & data.TT<time(nf+1),1);
            y2 = data.MM(data.TT>=time(nf) & data.TT<time(nf+1),2);
            y3 = data.MM(data.TT>=time(nf) & data.TT<time(nf+1),3);
            y1R = data.MMrot(data.TT>=time(nf) & data.TT<time(nf+1),1);
            y2R = data.MMrot(data.TT>=time(nf) & data.TT<time(nf+1),2);
            y3R = data.MMrot(data.TT>=time(nf) & data.TT<time(nf+1),3);
            y4 = Mamp(data.TT>=time(nf) & data.TT<time(nf+1))./Mamp(1);
            addpoints(h1,x1,y1);
            addpoints(h2,x1,y2);
            addpoints(h3,x1,y3);
            addpoints(h4,x1,y4);
            
            legend(gui.axes_handles.MagL,'M_x','M_y','M_z','Location','SouthEast');
            set(gui.axes_handles.MagL,'XLim',[min(data.TT.*1e3) max(data.TT.*1e3)],...
                'XTick',linspace(min(data.TT.*1e3),max(data.TT.*1e3),5))
            set(get(gui.axes_handles.MagL,'XLabel'),'String','t [ms]');
            set(get(gui.axes_handles.MagL,'YLabel'),'String','M');
            set(get(gui.axes_handles.MagL,'Title'),'String','magnetization components - lab frame');
            
            
            addpoints(h7,x1,y1R);
            addpoints(h8,x1,y2R);
            addpoints(h9,x1,y3R);
            addpoints(h10,x1,y4);
            
            legend(gui.axes_handles.MagR,'M_x','M_y','M_z','Location','SouthEast');
            set(gui.axXYZR,'XLim',[min(data.TT.*1e3) max(data.TT.*1e3)],...
                'XTick',linspace(min(data.TT.*1e3),max(data.TT.*1e3),5))
            set(get(gui.axes_handles.MagR,'XLabel'),'String','t [ms]');
            set(get(gui.axes_handles.MagR,'YLabel'),'String','M');
            set(get(gui.axes_handles.MagR,'Title'),'String','magnetization components - rot frame');
            
            addpoints(h5,y1,y2,y3);
            view([134 30]);
            set(gui.axSphereL,'Color','w','XColor','none','YColor','none','ZColor','none');
            axis equal tight
            
            addpoints(h6,y1R,y2R,y3R);
            set(gui.axSphereR,'Color','w','XColor','none','YColor','none','ZColor','none');
            axis equal tight
            
            pause(1/FPS);
        end
        hold(gui.axes_handles.MagL,'off');
        hold(gui.axes_handles.MagR,'off');
        hold(gui.axSphereL,'off');
        hold(gui.axSphereR,'off');
        set(gui.text_handles.Status,'String','Animation running ... done.');
        pause(0.01);
        
    case 'prepol'
        set(gui.text_handles.Status,'String','Animation running ...');
        pause(0.01);
        Mamp = zeros(numel(data.TT),1);
        MM = zeros(numel(data.TT),3);
        for i=1:numel(data.TT)
            Mamp(i,1) = norm(data.MM(i,:)); % amplitude of M
            MM(i,:) = data.MM(i,:)./Mamp(1); % normalize to initial M
        end
        Beffn = zeros(numel(data.results.prepol.t),3);
        for i=1:numel(data.results.prepol.t)
            Beffn(i,:) = data.results.prepol.Beff(i,:)./norm(data.results.prepol.Beff(i,:));
        end
        
        dt = data.TT(end)/NoFrames;
        time = 0:dt:data.TT(end);
        
        cla(gui.axXYZL);
        axes(gui.axXYZL);
        h1 = animatedline('Color','r','LineWidth',2,'Parent',gui.axXYZL);
        h2 = animatedline('Color','g','LineWidth',2,'Parent',gui.axXYZL);
        h3 = animatedline('Color','b','LineWidth',2,'Parent',gui.axXYZL);
        h4 = animatedline('Color','k','LineWidth',2,'LineStyle','--','Parent',gui.axXYZL);
        if data.prepol.Tramp < data.basic.Tsim
            h11 = animatedline('Color',[0.8 0.8 0.8],'LineWidth',1,'Parent',gui.axXYZL);
        end
        
        cla(gui.axSphereL);
        axes(gui.axSphereL);
        plotBlochSphere(30,30); hold on;
        h5a = animatedline('Color',myui.color.prepol,'LineWidth',2,'Parent',gui.axSphereL);
        h5b = animatedline('Color',myui.color.basic,'LineWidth',2,'Parent',gui.axSphereL);
        h6 = animatedline('Color',myui.color.prepolB,'LineWidth',2,'Parent',gui.axSphereL);
        
        cla(gui.axBpre);
        axes(gui.axBpre);
        h10 = animatedline('Color',myui.color.prepolB,'LineWidth',2,'Parent',gui.axBpre);
        y10a = data.results.prepol.Bp./data.basic.B0;
        
        cla(gui.axalpha);
        axes(gui.axalpha);
        h7 = animatedline('Color',myui.color.prepolB,'LineWidth',2,'Parent',gui.axalpha);
        
        cla(gui.axes_handles.dadt);
        axes(gui.axes_handles.dadt);
        h8 = animatedline('Color',myui.color.prepolB,'LineWidth',2,'Parent',gui.axes_handles.dadt);
        y8a = data.results.prepol.dadt;
        tt = data.results.prepol.t(1:end-1);
        
        cla(gui.axes_handles.wda);
        axes(gui.axes_handles.wda);
        h9 = animatedline('Color',myui.color.prepolB,'LineWidth',2,'Parent',gui.axes_handles.wda);
        y9a = data.results.prepol.omega(1:numel(data.results.prepol.dadt))./data.results.prepol.dadt;
        
        for nf = 1:NoFrames-1
            
            x1 = data.TT(data.TT>=time(nf) & data.TT<time(nf+1)).*1e3;
            y1 = MM(data.TT>=time(nf) & data.TT<time(nf+1),1);
            y2 = MM(data.TT>=time(nf) & data.TT<time(nf+1),2);
            y3 = MM(data.TT>=time(nf) & data.TT<time(nf+1),3);
            
            if time(nf)<data.prepol.Tramp/1e3
                y1a = MM(data.TT>=time(nf) & data.TT<time(nf+1),1);
                y2a = MM(data.TT>=time(nf) & data.TT<time(nf+1),2);
                y3a = MM(data.TT>=time(nf) & data.TT<time(nf+1),3);
            else
                y1b = MM(data.TT>=time(nf) & data.TT<time(nf+1),1);
                y2b = MM(data.TT>=time(nf) & data.TT<time(nf+1),2);
                y3b = MM(data.TT>=time(nf) & data.TT<time(nf+1),3);
            end
            
            y4 = Mamp(data.TT>=time(nf) & data.TT<time(nf+1))./Mamp(1);
            addpoints(h1,x1,y1);
            addpoints(h2,x1,y2);
            addpoints(h3,x1,y3);
            addpoints(h4,x1,y4);
            if data.prepol.Tramp < data.basic.Tsim
                addpoints(h11,[data.prepol.Tramp data.prepol.Tramp],[min(min(MM)) max(max(MM))]);
            end
            
            legend(gui.axXYZL,'M_x','M_y','M_z','Location','SouthEast');
            set(gui.axXYZL,'XLim',[min(data.TT.*1e3) max(data.TT.*1e3)],...
                'XTick',linspace(min(data.TT.*1e3),max(data.TT.*1e3),5),...
                'YLim',[min(min(MM)) max(max(MM))])
            set(get(gui.axXYZL,'XLabel'),'String','t [ms]');
            set(get(gui.axXYZL,'YLabel'),'String','M');
            set(get(gui.axXYZL,'Title'),'String',['projection M \rightarrow B_0 p = ',sprintf('%4.3f',data.results.prepol.p)]);
            
            y6(:,1) = Beffn(data.results.prepol.t./1e3>=time(nf) & data.results.prepol.t./1e3<time(nf+1),1);
            y6(:,2) = Beffn(data.results.prepol.t./1e3>=time(nf) & data.results.prepol.t./1e3<time(nf+1),2);
            y6(:,3) = Beffn(data.results.prepol.t./1e3>=time(nf) & data.results.prepol.t./1e3<time(nf+1),3);
            if time(nf)<data.prepol.Tramp/1e3
                addpoints(h5a,y1a,y2a,y3a);
            else
                addpoints(h5b,y1b,y2b,y3b);
            end
            addpoints(h6,y6(:,1),y6(:,2),y6(:,3));
            set(gui.axSphereL,'XGrid','on','YGrid','on','ZGrid','on','Box','on');
            axis(gui.axSphereL,'equal','tight'); view(3);
            set(get(gui.axSphereL,'XLabel'),'String','X');
            set(get(gui.axSphereL,'YLabel'),'String','Y');
            set(get(gui.axSphereL,'ZLabel'),'String','Z');
            set(get(gui.axSphereL,'Title'),'String','"Bloch sphere"');
            
            if data.results.prepol.t(end)/1e3 >= time(nf)
                x2 = data.results.prepol.t(data.results.prepol.t./1e3>=time(nf) & data.results.prepol.t./1e3<time(nf+1));
                x3 = tt(tt./1e3>=time(nf) & tt./1e3<time(nf+1));
                
                y10 = y10a(data.results.prepol.t./1e3>=time(nf) & data.results.prepol.t./1e3<time(nf+1));
                addpoints(h10,x2,y10);
                set(gui.axBpre,'XLim',[min(data.results.prepol.t) max(data.results.prepol.t)],'YLim',[0 data.prepol.Factor]);
                set(get(gui.axBpre,'XLabel'),'String','t [ms]');
                set(get(gui.axBpre,'YLabel'),'String','B_{pre} [B_0]');
                set(get(gui.axBpre,'Title'),'String','B_{pre} switch-off ramp');
                
                y7 = data.results.prepol.alpha(data.results.prepol.t./1e3>=time(nf) & data.results.prepol.t./1e3<time(nf+1));
                addpoints(h7,x2,y7);
                set(gui.axalpha,'XLim',[min(data.results.prepol.t) max(data.results.prepol.t)],'YLim',[min(data.results.prepol.alpha) max(data.results.prepol.alpha)]);
                set(get(gui.axalpha,'XLabel'),'String','t [ms]');
                set(get(gui.axalpha,'YLabel'),'String','\alpha [deg]');
                set(get(gui.axalpha,'Title'),'String','angle \alpha between B_0 and B_{eff}');
                
                y8 = y8a(tt./1e3>=time(nf) & tt./1e3<time(nf+1));
                addpoints(h8,x3,y8);
                set(gui.axes_handles.dadt,'XLim',[min(data.results.prepol.t) max(data.results.prepol.t)],'YLim',[min(data.results.prepol.dadt) max(data.results.prepol.dadt)]);
                set(get(gui.axes_handles.dadt,'XLabel'),'String','t [ms]');
                set(get(gui.axes_handles.dadt,'YLabel'),'String','d\alpha / dt');
                set(get(gui.axes_handles.dadt,'Title'),'String','time derivative of \alpha');
                
                y9 = y9a(tt./1e3>=time(nf) & tt./1e3<time(nf+1));
                addpoints(h9,x3,y9);
                set(gui.axes_handles.wda,'XLim',[min(data.results.prepol.t) max(data.results.prepol.t)],'YScale','log','YLim',[min(y9a) max(y9a)]);
                set(get(gui.axes_handles.wda,'XLabel'),'String','t [ms]');
                set(get(gui.axes_handles.wda,'YLabel'),'String','\omega / d\alpha');
                set(get(gui.axes_handles.wda,'Title'),'String','adiabatic criterion(?)');
            end
            clear y1 y2 y3 y4 y6 y7 y8 y9
            pause(1/FPS);
        end
        hold(gui.axXYZL,'off');
        hold(gui.axSphereL,'off');
        set(gui.text_handles.Status,'String','Animation running ... done.');
        pause(0.01);
        
    case 'pulse'
        mgh = msgbox('Not yet implemented.','BLOCHUS info');
        
    case 'prepolpulse'
        mgh = msgbox('Not yet implemented.','BLOCHUS info');
end


% change the pushbutton color back to green
% set(src,'BackGroundColor','g');
set(src,'BackGroundColor',[.94 .94 .94]);
plotResults(fig);
end % onPushAnimate
