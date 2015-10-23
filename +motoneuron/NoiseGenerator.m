classdef NoiseGenerator < handle
% NoiseGenerator: Models the noise input for the motoneuron model as done
% in the original script by Francesco
%
% @author Daniel Wirtz @date 2014-05-21
%
% @new{0,7,dw,2014-05-21} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties
        RandSeed = 100000;
        
        % Set this to true to disable noisy signal output.
        %
        % This will cause the getInput function to simply return the mean
        % values of the respective noises for each time t.
        %
        % @type logical @default false
        DisableNoise = false;
    end
    
    properties(SetAccess=private)
        totalNoise;
        baseNoise;
        a;
        b;
        AP;
        factor;
        baseMean;
        totalMean;
    end
    
    methods
        function this = NoiseGenerator
            % Loads the experimental data for motoneuron noise input
            %
            % The noise data is available on a millisecond sampling interval.
            %
            % Relevant fields:
            % - AP: Mean current
            % - thetaP: Variance of the independent noise
            % - LOWPASSP: Lowpass bandwith of the independent noise
            % - scaleP: total standard deviation, "AP/scale"
            % - noiseP: The noise data. "Common noise (alpha+beta)"
            mc = metaclass(this);
            p = load(fullfile(fileparts(which(mc.Name)),'neuro_input'));
            this.factor = sqrt(p.thetaP)/p.scaleP;
            this.AP = p.AP;
            [this.b, this.a] = butter(2,p.LOWPASSP*2/1000,'low');
            % Base noise
            this.baseNoise = p.noiseP;
            this.baseMean = mean(p.noiseP);
        end
        
        function setFibreType(this, mu_fibretype)
            rs = RandStream('mt19937ar','Seed',round(this.RandSeed*mu_fibretype*100));
            noiseSI = filter(this.b,this.a,rs.randn(1,length(this.baseNoise)));
            % Independent noise
            indepNoise = this.AP*this.factor/std(noiseSI)*noiseSI;
            this.totalNoise = this.baseNoise + indepNoise;
            this.totalMean = this.baseMean + mean(indepNoise);
        end
        
        function u = getInput(this, t)
            % The raw noise data is available in a sampling interval of one millisecond.
            % As the global model unit is also milliseconds, `t=1` equals an elapsed time of 
            % one millisecond.
            %
            % As theoretically the finite noise samples are ending at some stage, we put them
            % in an infinite loop via mod(t,numsamples).
            %
            % Return values:
            % u: A 2xlength(t) vector containing the mean current as first
            % component and the total (fibre-type dependent) noise as
            % second component.
            
            if this.DisableNoise
                u = zeros(2,length(t));
                u(1,:) = this.AP; 
                u(2,:) = this.totalMean;
            else
                pos = ceil(mod(t+eps,length(this.totalNoise)));
                % total noise. the mean current factor is build into the affine input
                % mapping B.
                u = [ones(size(pos))*this.AP
                     this.totalNoise(pos)];
            end
        end
    end
    
end