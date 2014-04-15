classdef Tests
% Tests: Some test settings regarding the PCD model simulations.
%
%
%
% @author Daniel Wirtz @date 2012-06-11
%
% @new{0,6,dw,2012-06-11} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing
    
    methods(Static)
        
        function m = demo_SpatialDependency
            m = models.pcdi.PCDIModel;
            m.System.h = 2e-7;
            m.T = 100;
            m.dt = 5;
            mu = [.5 .6 .5 2.5 .5]';
            
            % Plot the currently used coefficient
            m.System.plotDiffusionCoeff(mu(5));
            
            % Plot range of possible coefficients
            m.System.plotDiffusionCoeff(0:.1:1);
            
            [t,y,~,x] = m.simulate(mu);
            
            % Normal output
            m.plot(t,y);
            
            % Plot c3 concentration at center line
            m.plotState(t,x,1);
            
            % Plot top view over time
            m.plotState(t,x,2);
            
            % Plot augmented concentration differences
            m.plotState(t,x,3);
            
        end
    end   
end