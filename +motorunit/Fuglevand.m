classdef Fuglevand < handle
% Fuglevand: 
%
% @docupdate
%
% @author Daniel Wirtz @date 2014-05-23
%
% @new{0,7,dw,2014-05-23} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties
        fsamp = 1000; % sampling rate [samples / s]
        len = 10;     % length of stimulation in [s]
        
        RR = 30;      % range of recruitment threshold values
        fmin = 8;     % min discharge rate
        fmax = 35;    % max discharge rate
        
        sd_ISI = 0.1; % Standard dev. of ISI [Standard dev. del Jitter medio]
        
        RP = 100;  % range of twitch forces RP
        RT = 3;    % range of contraction times RT
        
        lengthWin = 1000; % duration of a single twitch in [ms]
        
        Resolution = 100;
    end
    
    properties(SetAccess=private)
        pende;
    end
    
    properties(Access=private)
        r;
        fSeed;
    end
    
    properties(Dependent)
        Seed;
    end
    
    methods
        function this = Fuglevand
            %Defines the frequency change will know how to increase the level of the contraction frequency becomes maximum when the force reaches to 60% more?
            %the minimum force recruitment
            this.pende = (this.fmax-this.fmin)/60; % gain g_e in the paper
            
            this.Seed = 1;
        end
        
        function forces = getForces(this, mu, T, dt)
            % Obtain the modeled force response for given parameter and
            % times
            %
            % Parameters:
            % mu: The parameter containing fibre type and exitation level
            % @type rowvec<double> @default [0.1 5]
            % T: The final time `T` @type double @default 2000
            % dt: The time-step `\Delta t` @type double @default 1
            
            if nargin < 4
                dt = 1;
                if nargin < 3
                    T = 2000;
                    if nargin < 2
                        mu = [.1 5];
                    end
                end
            end
            
            % Create desired time grid before messing with T
            t = 0:dt:T;
                       
            % Add window of single peak to simulation time
            T = T + this.lengthWin;
            
            % buffer as the start time is inconsistent
            buffer = 300; %[ms]
            T = T + buffer;
            
            firing = zeros(1,T);
            forces = firing;
            
            k = 2;
            gain = 1;
            
            % recruitment threshold excitation
            rte = exp(log(this.Resolution)*mu(1))*this.RR/this.Resolution;
            exDrive = mu(2)*this.RR/10;
            
            if rte <= exDrive*1.001
            
                % peak twitch force
                P = exp(log(this.Resolution)*mu(1))*this.RP/this.Resolution;
                Tfac = 90*(1./P)^(log(this.RT)/log(this.RP));

                % instantaneous discharge frequency
                finst = this.pende*(exDrive-rte)+this.fmin;
                finst = min(finst, this.fmax);

                % previous firing time
                % rand is a pseudorandom value drawn from the standard uniform distribution on the open interval (0,1).
                istprec = this.r.rand/finst;

                istatt = 0;  % current firing time
                while(istatt+0.02 < T/1000)
                    % variability due to the jitter
                    % standard normal distri with mean 0 and standard deviation sd_ISI
                    % rand is a pseudorandom value drawn from the standard normal distribution
                    var_jitter = this.sd_ISI*this.r.randn*0;

                    % moment of activation (time)
                    istatt = istprec+(1+var_jitter)/finst;
                    % check that the firing is not too close to the previous one
                    if(istatt < istprec+0.02)
                        istatt = istprec+0.02;
                    end
                    firing(round(istatt*1000)) = 1;
                    Ti_ISI = Tfac/(istatt-istprec)/1000;
                    if(Ti_ISI < 0.4)
                        gain(k) = 1.0;%#ok
                    else
                        gain(k) = (1-exp(-2*Ti_ISI^3))/Ti_ISI/((1-exp(-2*0.4^3))/0.4);%#ok
                    end
                    k=k+1;
                    istprec = istatt;
                end % while
                hlp = [0:round(this.lengthWin/1000*this.fsamp)]/this.fsamp/(Tfac/1000);
                winCoeff = P.*(hlp.*exp(1-hlp));
                winCoeff = [0 winCoeff(winCoeff>1.e-4)];

                spikes = find(firing==1);
                for i=1:size(spikes,2)
                    if(spikes(i)+size(winCoeff,2) <= length(forces))
                        idx = spikes(i) + (1:size(winCoeff,2));
                        forces(idx) = forces(idx)+gain(i)*winCoeff;
                    end
                end
            end
            % Subtract window of single peak from force
            forces = forces(1:end-this.lengthWin);
            
            startpos = find(forces,1)-10;
            if startpos > 0
                forces(:,1:startpos) = [];
                forces(:,end-buffer+startpos+1:end) = [];
            end
            
            % Interpolate to desired time grid
            forces = interp1(1:length(forces),forces,t,'cubic');
        end
        
        function set.Seed(this, value)
            this.r = RandStream('mt19937ar','Seed',value);
            this.fSeed = value;
        end
        
        function value = get.Seed(this)
            value = this.fSeed;
        end
    end
    
    methods(Static)
        function test_ForceResponse
            f = models.motorunit.Fuglevand;
            [x,y] = meshgrid(0:.05:1,0:.4:10);
            p = [x(:) y(:)]';
            force = zeros(1,numel(x));
            for k=1:size(p,2)
                fmu = f.getForces(p(:,k));
                force(k) = max(fmu);
            end
            force = reshape(force,size(x,1),[]);
            surf(x,y,force,'FaceColor','interp','EdgeColor','interp');
        end
    end
    
end