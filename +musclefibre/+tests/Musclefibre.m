classdef Musclefibre < matlab.unittest.TestCase
    % Tests for musclefibre model
    
    properties(TestParameter)
        %sv = {1 2};
        sv = {2};
        ic = {true false};
        spm = {true false};
        os = {true false};
        ft = {0 1};
    end
    
    methods(Test)
        function MusclefibreVersions(~,sv,ic,spm,os,ft)
            m = models.musclefibre.Model('N',4,...
                'SarcoVersion',sv,...
                'DynamicIC',ic,...
                'SPM',spm,...
                'OutputScaling',os);
            m.T = 100;
            m.dt = .01;
            [t,y] = m.simulate([ft; 9]);
            m.plot(t,y);
            drawnow;
        end
        
        function LargeFibre(~)
            m = models.musclefibre.Model('N',100);
            m.T = 100;
            [t,y] = m.simulate;
            m.plot(t,y);
        end
    end
end

