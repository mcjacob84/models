classdef HillGen < handle
% HillGen: 
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
        % Einstellungen für zufällige Hügel
        % Angabe eines Wertes > 0 erstellt zufällige hügel der angegebenen
        % Anzahl mit zufallswerten, die durch die angaben max_.. gegeben sind
        hillcount = 0; % hills = 0 stellt auf manuelle Definition der Hügel um!
        max_hillheight = .5; % Höhe in Meter
        max_hillradius = 3; % Maximaler Radius in Meter
        min_hillradius = .5; % Minimaler Radius in Meter
    end
    
    methods
        function this = HillGen
        end
        
        function generate(this)
            % Orte der Hügel
            [hxpos,hypos] = equidisthills(hillcount);
            hills = [hxpos hypos];
            randomizehills; % Fügt automatisch die beiden Spalten mit Höhe+Radius ein, und plottet
        end
        
        function [h,r] = generatehills(hills,maxheight,maxrad,minrad)
            % Höhen, zwischen -1 und 1
               h = maxheight*2*(rand(hills,1)-0.5);
            % Radien der Hügel
               r = (maxrad-minrad)*rand(hills,1)+minrad;
            % Erzeugt zufällige Hügel
        end
        
        function randomize(this)
            [h,r] = generatehills(size(hills,1),max_hillheight,max_hillradius,min_hillradius);
            hills(:,3) = h; hills(:,4) = r;
            % Lochdaten anpassen
            HZ = green(HX,HY);
            plotgreen();
            % Erzeugt neue höhen für das Gras
        end
        
        function [ xpos,ypos ] = equidisthills( hills )
            global area;
            xhills = round(sqrt((area(2)-area(1))/(area(4)-area(3))*hills));
            yhills = round(hills/xhills);
            xdist = (area(2)-area(1))/(xhills+1);
            ydist = (area(4)-area(3))/(yhills+1);
            % Initialisieren
            xpos = zeros(hills,1); ypos = xpos;
            % Positionieren
            for hill = 0:hills-1
                xpos(hill+1) = area(1) + (mod(hill,xhills)+1)*xdist;
                ypos(hill+1) = area(3) + (floor(hill/xhills)+1)*ydist;
            end
            % Erzeugt eine äquidistante Verteilung von Hügeln im Gelände
        end
    end
    
end