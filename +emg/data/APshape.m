m = models.musclefibre.Model(500,5);
m.T = 50;
m.dt = .1;
dx = 0.01;
MU = 0:0.1:1;
U = [2,2,2,2,2,2,2.5,3,3,4,4];
APshape=cell(11,1);
Vmidx = 18:58:29017;
for i = 1:11
    [t,y] = m.simulate([MU(i);U(i)],1);
    APstart = 30;
    tdx = find((y(18+(APstart-1)*58,:)> -80).*(1:length(t)>5),1,'first');
    y = y(Vmidx,tdx); % V_m at some time
    APend = find((y<-80).*(1:500 > APstart)',1,'first')-1;
    APshape{i} = y(APstart:APend);
end

save APshape APshape MU dx