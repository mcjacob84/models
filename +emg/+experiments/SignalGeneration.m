%% Init
%setpref('MWHR15','base','/home/dwirtz/papers/MWHR15 - EMG Signal generation');
base = getpref('MWHR15','base');

clear classes;%#ok
args = {'Shapes','precomp','FiringTimes','precomp','MUTypes',.6};
%args = [args {'Dim',[150 2 3],'Geo',[15 3 2 1]}];
mw = models.emg.Model(args{:});
mo = models.emg.Model(args{:},'DynAmpPS',false);
mw.T = 500;
mo.T = 500;
mw.dt = .1;
mo.dt = .1;

mw.System.prepareSimulation(mw.DefaultMu,mw.DefaultInput);
mo.System.prepareSimulation(mo.DefaultMu,mo.DefaultInput);

%%
mw.plotAPShapes;

%%
[t,yw] = mw.simulate;
[t,yo] = mo.simulate;
steps = 1:20:length(t);
%%
mw.plot(t(steps),yw(:,steps));
%%
mo.plot(t(steps),yo(:,steps));

%%
t = mw.Times;
Rw = mw.System.f.getVm(t);
Ro = mo.System.f.getVm(t);

%%
pm = PlotManager(false,1,2);
pm.LeaveOpen = true;
ax = pm.nextPlot('nodyn','Original Vm','t [ms]','Vm [mV]');
%[~,idx] = max(max(Ro,[],2));
idx = 1:size(Ro,1);
plot(ax,t,Ro(idx,:));
ax = pm.nextPlot('dyn','Dynamic Vm','t [ms]','Vm [mV]');
plot(ax,t,Rw(idx,:));
%%
pm.savePlots(base,'Format',{'pdf','jpg'});

%%
pm = PlotManager(false,1,2);
ax = pm.nextPlot('diff_rhs','Difference in Vm RHS signals','t [ms]','Difference (L2-space) [mV]');
plot(ax,t,Norm.L2(Rw-Ro));
ax = pm.nextPlot('diff_vm','Difference in effective Vm signals','t [ms]','Difference (L2-space) [mV]');
plot(ax,t,Norm.L2(yw-yo));
pm.done;