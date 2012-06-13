classdef Tests
% Tests: Some tests and simulation settings for Burger's equation models.
%
%
%
% @author Daniel Wirtz @date 2012-06-13
%
% @new{0,6,dw,2012-06-13} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    methods(Static)
        function m = test_Burgers
            m = models.burgers.Burgers;
            
            %% Sampling - manual
            s = sampling.ManualSampler;
            s.Samples = logspace(log10(0.04),log10(0.08),150);
            m.Sampler = s;
            
            %% Approx
            a = approx.KernelApprox;
            a.Kernel = kernels.GaussKernel;
            a.ParamKernel = kernels.GaussKernel;
            
            al = approx.algorithms.VectorialKernelOMP;
            al.MaxGFactor = [1.5 0 1.5];
            al.MinGFactor = [.2 0 .6];
            al.gameps = 1e-4;
            al.MaxExpansionSize = 300;
            al.MaxAbsErrFactor = 1e-5;
            al.MaxRelErr = 1e-3;
            al.NumGammas = 40;
            a.Algorithm = al;
            m.Approx = a;
            
            m.offlineGenerations;
            save test_Burgers;
        end
        
        function m = test_Burgers_DEIM(dim, version)
            if nargin < 2
                version = 1;
                if nargin < 1
                    dim = 2000;
                end
            end
            m = models.burgers.Burgers(dim, version);
            
            if dim < 1000
                m.Data = data.MemoryModelData;
            else
                m.Data = data.FileModelData;
            end
            
            %% Sampling - manual
            s = sampling.ManualSampler;
            s.Samples = logspace(log10(0.04),log10(0.08),100);
            m.Sampler = s;
            
            %% Approx
            a = approx.DEIM;
            a.MaxOrder = 80;
            m.Approx = a;
            
            offline_times = m.offlineGenerations;
            gitbranch = KerMor.getGitBranch;
            
            clear a s;
            d = fullfile(KerMor.App.DataStoreDirectory,'test_Burgers_DEIM');
            mkdir(d);
            oldd = pwd;
            cd(d);
            eval(sprintf('save test_Burgers_DEIM_d%d_v%d',dim,version));
            cd(oldd);
        end
        
        function m = test_Burgers_DEIM_B(dim)
            if nargin < 1
                dim = 200;
            end
            m = models.burgers.Burgers(dim, 2);
            m.PlotAzEl = [-130 26];
            
            if dim <= 2000
                m.Data = data.MemoryModelData;
            else
                m.Data = data.FileModelData(m);
            end
            
            %% Sampling - manual
            s = sampling.ManualSampler;
            s.Samples = logspace(log10(0.005),log10(0.02),100);
            m.Sampler = s;
            
            %% Approx
            a = approx.DEIM;
            a.MaxOrder = 80;
            m.Approx = a;
            
            s = m.System;
            s.x0 = dscomponents.ConstInitialValue(zeros(dim,1));
            
            x = linspace(m.Omega(1),m.Omega(2),dim+2);
            x = x(2:end-1);
            pos1 = logical((x >= .1) .* (x <= .3));
            pos2 = logical((x >= .6) .* (x <= .7));
            %% Not that exciting set of inputs
%             s.Inputs{1} = @(t)2*sin(2*t*pi);
%             B = zeros(dim,1);
%             B(pos1) = 4*exp(-((x(pos1)-.2)/.03).^2);
%             B(pos2) = 1;
%             
            s.Inputs{1} = @(t)[sin(2*t*pi); (t>.2)*(t<.4)];
            B = zeros(dim,2);
            B(pos1,1) = 4*exp(-((x(pos1)-.2)/.03).^2);
            B(pos2,2) = 4;
            
            s.B = dscomponents.LinearInputConv(B);

            m.TrainingInputs = 1;
            
            offline_times = m.offlineGenerations;
            gitbranch = KerMor.getGitBranch;
            
            clear a s;
            d = fullfile(KerMor.App.DataStoreDirectory,'test_Burgers_DEIM_B');
            mkdir(d);
            oldd = pwd;
            cd(d);
            eval(sprintf('save test_Burgers_DEIM_B_d%d',dim));
            cd(oldd);
        end
    end
    
end