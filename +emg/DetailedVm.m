classdef DetailedVm < handle
    %VMGENERATOR A separate class used to compute detailed Vm signals for
    % defined parameters, to be used within the models.emg suite
    %
    
    properties(SetAccess=private)
        submesh_idx;
        % musclefibre model for full simulations
        musclefibremodel;
        mus_precomp;
    end
    
    methods
        
        function this = DetailedVm(len, dims, sv)
            dx = models.musclefibre.Model.dxDefault;
            N = round(len/dx);
            % Coarse index within fine resolution
            this.submesh_idx = round(linspace(1,N,dims));
            mf = models.musclefibre.Model(...
                'SarcoVersion',sv,...
                'N',N,'dx',dx,...
                'DynamicIC',true,'SPM',false,'OutputScaling',false,...
                'Spindle',false,'Noise',true,...
                'JunctionN',1);
            mf.EnableTrajectoryCaching = true;
            this.musclefibremodel = mf;
            
            % Hack for 
            mudir = fileparts(which('models.emg.Model'));
            s = load(fullfile(mudir,'data','mus.mat'));
            this.mus_precomp = s.mus;
        end
        
        function sig = computeSignal(this, t, mu)
            pos = Utils.findVecInMatrix(mu,this.mus_precomp);
            if pos > 0
                s = load(sprintf('/data/local/musclefibre/Vm_%d.mat',pos));
                if isequal(s.t,0:dt:T)
                    fine_signal = s.Vm;
                end
            else
                mf = this.musclefibremodel;
                mf.T = t(end); % infer from passed parameters
                mf.dt = t(2)-t(1);
                [~,~,~,x] = mf.simulate(mu,1);
                % Find positions of Vm on each sarcomere
                moto_off = mf.System.dm;
                % Each first sarcomere model dimension contains the current
                % Vm value
                pos = moto_off + (1:mf.System.dsa:size(x,1)-moto_off);
                fine_signal = x(pos,:);
            end
            % Only return the effectively needed signals from the
            % finer submodel simulation.
            sig = fine_signal(this.submesh_idx,:);
        end
    end
    
end

