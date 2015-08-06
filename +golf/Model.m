classdef Model < models.BaseModel
% Model: 
%
% @docupdate
%
% @author Daniel Wirtz @date 2015-08-06
%
% @new{0,7,dw,2015-08-06} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties
        % Golfballradius in Meter gut: 0.015 
        rad_ball = 0.021; % m (42mm Durchmesser von Wikipedia)
        
        rad_hole = 0.054;
        
        % Geländeauflösung
        detail = 0.02;
        
        vmin = 0.03;
        
        Green;
        Hole;
        Ball;
    end
    
    properties(Access=private)
        HX_;
        HY_;
        HZ_;
    end
    
    properties(Access=private, Transient)
        tx_ball;
        tx_rasen;
    end
    
    methods
        function this = Model(greenfile)
            if nargin < 1
                greenfile = 'bsp4';
            end
            this.tx_ball = imread(fullfile(fileparts(mfilename('fullpath')),'tx_golfball.png'));
%             this.tx_ball = imread(fullfile(fileparts(mfilename('fullpath')),'tx_golfball_falschrum.png'));
            this.tx_rasen = imread(fullfile(fileparts(mfilename('fullpath')),'tx_rasen.png'));
            
            this.T = 20;
            this.dt = 0.04;
            
            this.System = models.golf.System(this);
            
            this.ODESolver = solvers.MLode15i;
            
            this.loadGreen(greenfile);
        end
        
        function [t, x, ctime] = computeTrajectory(this, mu, inputidx)
            if isa(this.ODESolver,'solvers.MLWrapper')
                opt = this.ODESolver.odeopts;
                opt.OutputFcn = @this.ODE_Callback;
                this.ODESolver.odeopts = opt;
            else
                error('Not working with non-MatLab solver wrappers.');
            end
            [t, x, ctime] = computeTrajectory@models.BaseModel(this, mu, inputidx);
        end
        
        function status = ODE_Callback(this, t, y, flag)
            status = 0;
            if isempty(flag)
                % For cases with more than one computed state
                y = y(:,end);
                
                x = y(1:2);
                v = y(3:4);
                
                % Area check
                area = this.Green.area;
                status = status || ...
                    ~(x(1)<=area(2) && x(1)>=area(1)) && (x(2)<=area(4) && x(2)>=area(3));
                
                % Prüft, ob der Ball "über dem Loch schwebt" und langsam genug ist, um
                % herein gefallen zu sein
                sys = this.System;
                mu = sys.mu;
                status = status || ...
                    ((norm(v) < sys.maxvforhole) && (norm(x-mu(5:6))<this.rad_hole));
                
                % Prüft, ob der Ball noch schnell genug ist.
                % Dabei muss beachtet werden, dass der Ball sich nicht an einem Hang
                % befindet, er konnte dort während einer Umkehr sehr langsam werden.
                % Dies indiziert der Gradient, der normiert auf ebener Fläche klein ist.
                gradient = this.gradgreen(x);
                status = status || ...
                    ~((norm(v) > this.vmin) || (norm(gradient) > 0.040));
            end
            if status
                a = 5;
            end
        end
        
        function loadGreen(this,filename)
            %LOADGREEN Loads a stored golf setting
            if nargin < 2
                filename = fullfile(fileparts(mfilename('fullpath')),'default.mat');
            elseif exist(filename,'file') ~= 2
                filename_ = fullfile(fileparts(mfilename('fullpath')),filename);
                if exist(filename_,'file') == 2
                    filename = filename_;
                else
                    filename_ = fullfile(fileparts(mfilename('fullpath')),[filename '.mat']);
                    if exist(filename_,'file') == 2
                        filename = filename_;
                    else
                        error('Invalid file name: %s',filename);
                    end
                end
            end
            
            s = load(filename);
            g = struct;
            a = s.('area');
            g.area = a;
            g.hills = s.('hills');
            
            % Set params to default parameter values
            h = s.('hole');
            this.DefaultMu(5:6) = h(1:2);
            p = s.('params');
            this.DefaultMu(1:4) = p(:);
            this.Green = g;
            
            % Das Loch grafisch erzeugen
            phi = linspace(0,2*pi,50);
            theta = linspace(0,pi,50);
            [PHI,THETA]=meshgrid(phi,theta);
            this.HX_ = sin(THETA).*cos(PHI);
            this.HY_ = sin(THETA).*sin(PHI);
        end
        
        function [ z ] = green(this,x,y)
            % Green - Erzeugt eine Hügellandschaft mit gegebenen Parametern
            hills = this.Green.hills;
            z = 0; 
            for numhill = 1:size(hills,1)
                rad = -log(0.05)/hills(numhill,4)^2;
                z = z + hills(numhill,3) * exp(-((x-hills(numhill,1)).^2 +(y-hills(numhill,2)).^2)*rad);
            end
        end
        
        function [ dx ] = gradgreen( this, x ) % x \in R^2
            % Gradient von g
            h = sqrt(eps);
            hlp = this.green(x(1),x(2));
            dx = 1/h*[this.green(x(1)+h,x(2))-hlp ; this.green(x(1),x(2)+h)-hlp];
        end
        
        function map = golfcolormap(~)
            steps = 200;
            % Erste Farbe: Schwarz für das Loch!
            map([1 2 3]) = [0 0 0];
            % Zweite Farbe: Verlauf für das Rasengrün!
            map = [map; [zeros(steps+1,1) (0:1/steps:1)' zeros(steps+1,1)] ];
            % Zweite Farbe: Verlauf für den Golfball!
            %map = [map; [(0:1/steps:1)' (0:1/steps:1)' (0:1/steps:1)']];
        end
        
        function ax = initPlot(this, pm)
            % Plotgreen: Zeichnet das Putting-Green und das zugehörige Loch
            
            mu = this.System.mu;
            
            %% Init hole plotting
            HX = mu(5) + this.rad_hole*this.HX_;
            HY = mu(6) + this.rad_hole*this.HY_;
            HZ = this.green(HX,HY);
            % Loch als horizontale Kreisfläche
%             HZ = ones(size(HX))*hz);

            %% 
            ax = pm.nextPlot('golf','Trace of golf ball','[x] = m','[y] = m');
            % ,'Position',[100 100 640 480]
            set(gcf,'Name','Modellierung','Renderer','opengl');
            
            % Grünfläche zeichnen
            a = this.Green.area;
            [x,y] = meshgrid(a(1):this.detail:a(2),...
                a(3):this.detail:a(4));
            z = this.green(x,y);
            hGreen = surfl(ax,x,y,z);
            colormap(this.golfcolormap);
            hold(ax,'on');

            % Lichteffekte
            set(hGreen,'AmbientStrength',.75);
            light('Position',[0 0 5],'Style','local');
            set(hGreen,'FaceLighting','phong','AmbientStrength',1)

            % Rasentextur - Sehr langsam! 
%             hGreenTex = surfl(ax,x,y,z+0.005);
%             set(hGreenTex,'CDATA',this.tx_rasen,'FaceColor','texturemap','FaceAlpha',.3);

            % Loch ins Grün einzeichnen (ein bisschen anheben, damit die grafiken sich
             % nicht überschneiden ..
             hHole = surfl(ax,HX,HY,HZ+0.01);
             % Lochfarbe schwarz setzen
             set(hHole,'CData',0,'FaceColor','texturemap');

             shading flat;
             daspect(ax,[1 1 1]);
        end
        
        function plot(this, t, y, pm)
            %             global rad_ball BX BY BZ t hFigure sol;% hBall;
            if nargin < 4
                pm = PlotManager;
                pm.LeaveOpen = true;
            end
            
            ax = this.initPlot(pm);
            withline = true;
            
            %             makevideo = false;
            %             video_folder = 'Video\';
            %
            %             if makevideo
            %                 set(hFigure,'DoubleBuffer','on');                                    % fps anpassen!
            %                 mov = avifile(strcat(video_folder,'output.avi'),'compression','None','fps',round(1/t));
            %             end
            
            oldgr = this.green(y(1,1),y(2,1));
            
            % Ballumfang berechnen
            ball_circ = 2*pi*this.rad_ball;
            % Startoffset bestimmen
            oldoffset = [-this.gradgreen([y(1,1) y(2,1)]); 1];
            oldoffset = this.rad_ball * oldoffset / norm(oldoffset);
            
            % Ball am anfang einmal zeichnen
            [BX,BY,BZ] = sphere(50);
            hBall = surfl(ax,this.rad_ball*BX + y(1,1)+oldoffset(1),...
                this.rad_ball * BY + y(2,1)+oldoffset(2),...
                this.rad_ball* BZ + oldgr +oldoffset(3));
            % s = size(get(hBall,'CDATA'));
            % [xg,yg] = meshgrid(-1:2/(s(1)-1):1);
            % col = .5*exp(-(xg.^2+yg.^2));
            % col = min(min(get(hBall,'CDATA'))) + col*(max(max(get(hBall,'CDATA')))/max(max(col)));
            % Textur für den Ball
            set(hBall,'CDATA',this.tx_ball,'FaceColor','texturemap');
            shading flat;
            zoom reset;
            
            dt = this.dt;
            for index = 2:size(y,2)
                
                % Einmal die Höhe berechnen
                gr = this.green(y(1,index),y(2,index));
                
                % Verschiebung bestimmen
                rel_move = [y(1,index)-y(1,index-1); y(2,index)-y(2,index-1); gr-oldgr];
                
                % Richtung des Balles
                % Die Rotationsschse ist senkrecht zur Bewegungsrichtung
                rot_axis = [-rel_move(2) rel_move(1) 0];
                
                if withline
                    plot3([y(1,index-1)+0.01 y(1,index)+0.01], [y(2,index-1)+0.01 y(2,index)+0.01], [oldgr+0.01 gr+0.01], 'r');
                end
                %             g = 9.81;
                %             gradi = grad(sol(:,index));
                %             nor = [gradi; -1];
                %             nor = nor/norm(nor);
                %             fn = nor * -g*nor(3);
                %             fb = [0; 0; -g] - fn;
                % Gradienten
                %             plot3([sol(1,index) sol(1,index)+gradi(1)], [sol(2,index) sol(2,index)+gradi(2)], [gr gr], 'r');
                % Normale
                %             plot3([sol(1,index) sol(1,index)+nor(1)], [sol(2,index) sol(2,index)+nor(2)], [gr gr+nor(3)], 'r');
                % Normalenkraft
                %             plot3([sol(1,index) sol(1,index)+fn(1)], [sol(2,index) sol(2,index)+fn(2)], [gr gr+fn(3)], 'b');
                % Hangabtriebskraft
                %             plot3([sol(1,index) sol(1,index)+fb(1)], [sol(2,index) sol(2,index)+fb(2)], [gr gr+fb(3)], 'y');
                % Zeichnet die Drehachsen
                %             plot3([sol(1,index) rot_axis(1)+sol(1,index)],[sol(2,index) rot_axis(2)+sol(2,index)],[gr gr]+.5,'w');
                
                % Korrekte Positionierung zum Untergrund, der Ball
                offset = [-this.gradgreen([y(1,index) y(2,index)]); 1];
                offset = this.rad_ball * offset / norm(offset);
                
                % Korrigierten Verschiebungsvektor bestimmen
                % (Altes Offset weg, verschieben, neues Offset drauf
                absmove = rel_move-oldoffset+offset;
                
                % Verschiebung der Daten vornehmen
                set(hBall,'XData',get(hBall,'XData')+absmove(1));
                set(hBall,'YData',get(hBall,'YData')+absmove(2));
                set(hBall,'ZData',get(hBall,'ZData')+absmove(3));
                
                % Entspricht Drehung des Balles um Winkel alpha
                alpha = 360 * norm(rot_axis)/ball_circ;
                % Zentrum der Drehung festlegen
                center = [y(1,index)+offset(1) y(2,index)+offset(2) gr+offset(3)];
                
                % Ball drehen!
                rotate(hBall,rot_axis,alpha,center);
                
                oldgr = gr;
                oldoffset = offset;
                
                % Beschriftung
                title(strcat('Zeit: ',num2str(t(index)),'s, v=',num2str(norm((y(:,index)-y(:,index-1))*1/dt)),'m/s'));
                
%                 if makevideo
%                     mov = addframe(mov,getframe(hFigure));
%                     % Erzeugt einzelne Dateien
%                     % print(hFigure,'-djpeg95',strcat(video_folder,'frame',num2str(v),'.jpg'));
%                 else
                    pause(dt);
%                 end

                if ~ishandle(ax)
                    return;
                end
            end
            hold(ax,'off');
%             if makevideo
%                 close(mov);
%             end
        end
    end
    
    methods(Static)
        function res = test_Golf
             m = models.golf.Model('4_2');
             [t,y] = m.simulate;
             m = models.golf.Model('4_3');
             [t,y] = m.simulate;
             m = models.golf.Model('4_4');
             [t,y] = m.simulate;
             m = models.golf.Model('4_5');
             [t,y] = m.simulate;
             
             m = models.golf.Model('bsp1');
             [t,y] = m.simulate;
             m = models.golf.Model('bsp2');
             [t,y] = m.simulate;
             m = models.golf.Model('bsp3');
             [t,y] = m.simulate;
             m = models.golf.Model('bsp4');
             [t,y] = m.simulate;
             
             m = models.golf.Model;
             [t,y] = m.simulate;
             m.plot(t,y);
             res = true;
        end
    end
end