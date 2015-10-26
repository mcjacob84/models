classdef EMG < matlab.unittest.TestCase
    %EMG Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(TestParameter)
        shapes = {'precomp' 'actual' 'rosen'};
        FT = {'precomp' 'actual'};
        sv = {1 2};
    end
    
    methods(Test)
        function ModelOptions(~, shapes, FT, sv)
            dim = [20 10 7 3];
            types = [0 rand(1,3) 1];
            m = models.emg.Model('SarcoVersion',sv,...
                'Shapes',shapes,...
                'MUTypes',types,...
                'FiringTimes',FT,...
                'Dim',[dim(1:2) sum(dim(3:4))],'Geo',.1*dim);
            m.T = 20;
            m.simulate;
        end
        
        function Plots(~)
            m = models.emg.Model('Shapes','precomp','FiringTimes','precomp');
            m.T = 10;
            m.dt = .1;
            [t,y] = m.simulate;
            m.plot(t,y);
            m.plot1D(t,y);
            m.plotMU3d(1);
            m.plotMU3d;
            m.plotMuscleConfiguration;
            m.plotPrecompAPShapes;
        end
        
        function FullSignal(~)
            m = models.emg.Model('Shapes','full','MUTypes',[0 1]);
            m.T = 10;
            m.dt = .1;
            m.simulate;
        end
    end
    
end

