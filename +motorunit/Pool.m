classdef Pool < handle
% Pool: A motorunit pool with either Fuglevand or Shorten model for force
% generation.
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
        FibreTypes;
        Version;
        %ShortenInput = @(t)3;
    end
    
    properties(SetAccess=private)
        fugle;
        shorten;
        alpha;
    end
    
    methods
        
        function this = Pool(version)
            if nargin < 1
                version = 1;
            end
            this.Version = version;
            this.fugle = models.motorunit.Fuglevand;
            % We use the Neumann shorten version in the motorunit pools
            this.shorten = models.motorunit.Shorten('SarcoVersion',2);
        end
        
        function prepareSimulation(this, model, meancurrent)
            ft = this.FibreTypes;
            nf = length(ft);
            a = zeros(nf,model.T+300);
            this.fugle.Seed = 1;
            if this.Version == 1
                for fidx = 1:nf
                    a(fidx,:) = this.fugle.getForces([ft(fidx); meancurrent],T+300);
                end
                startpos = find(sum(a,1),1)-10;
                if startpos > 0
                    a(:,1:startpos) = [];
                    a(:,end-300+startpos+1:end) = [];
                end
            elseif this.Version == 2
                m = this.shorten;
                m.T = model.T;
                m.dt = model.dt;
                %m.System.Inputs{1} = this.ShortenInput;
                for fidx = 1:nf
                    [~,y] = this.shorten.simulate([ft(fidx); meancurrent],1);
                    a(fidx,:) = y(2,:);
                end
            end
            this.alpha = a;
            mval = max(a(:));
            if mval > 0
                this.alpha = a/mval;
            end
        end
        
        function alpha = getActivation(this, t)
            a = this.alpha;
            alpha = a(:,min(size(a,2),round(t+1)));
        end
    end
    
end