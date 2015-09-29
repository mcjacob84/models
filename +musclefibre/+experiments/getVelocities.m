function [T,V,tgrid,Vinterp,Vpoly] = getVelocities(t, Vms, distances, tgrid, minV)

% Defaults
if nargin < 5
    minV = -20; %[mV]
    if nargin < 4
        tgrid = linspace(0,t(end),1000);
    end
end

nruns = length(Vms);
if isscalar(distances)
    distances = ones(1,nruns)*distances;
end
V = cell(1,nruns);
T = cell(1,nruns);
tsteps = length(tgrid);
Vinterp = zeros(nruns,tsteps);
Vpoly = zeros(nruns,tsteps);
pi = ProcessIndicator('Extracting velocities for %d runs',nruns,false,nruns);
for idx = 1:nruns
    Vm = Vms{idx};
    
    % Find all locations with negative derivative
    negpos = find(diff(Vm(1,:))<0);
    % Find the positions where the negative derivative starts
    firstpos = [1 find(diff(negpos) > 1)+1];
    % Check that the found locations are peaks and not in the lower noise
    ispeak = Vm(1,negpos(firstpos)) > minV;
    % Remove unwanted
    firstpos(~ispeak) = [];
    peaks_junction = negpos(firstpos);
    % Same for end peaks
    negpos = find(diff(Vm(2,:))<0);
    firstpos = [1 find(diff(negpos) > 1)+1];
    ispeak = Vm(2,negpos(firstpos)) > minV;
    firstpos(~ispeak) = [];
    peaks_end = negpos(firstpos);
    
    %     plot(t,Vm(2,:),'r',t(negpos),Vm(2,negpos),'bx')
    %     plot(t,Vm(2,:),'r',t(negpos(firstpos)),Vm(2,negpos(firstpos)),'bx')
    % Trivial case
    ndiff = length(peaks_junction) - length(peaks_end);
    if ndiff > 0
        plot(t,Vm(1,:),'r',t(peaks_junction),Vm(1,peaks_junction),'bx',...
            t,Vm(2,:),'b',t(peaks_end),Vm(2,peaks_end),'rx')
        peaks_junction(end-ndiff+1:end) = [];
        fprintf(2,'Check plot for correct removal of last junction peak!\n');
    end
    
    % Get the time instances the peak was recorded
    peaktimes = [];
    peaktimes(1,:) = t(peaks_junction);
    peaktimes(2,:) = t(peaks_end);
    % Compute velocities
    tdiff = diff(peaktimes,1);
    V{idx} = 10*distances(idx)./tdiff; % [cm/ms = 0.01m/0.001s = .1m/s]
    T{idx} = peaktimes(1,:);
    if size(peaktimes,2) > 1
        Vi = interp1(T{idx},V{idx},tgrid,'cubic');
        Vinterp(idx,:) = Vi;
        c = polyfit(T{idx},V{idx},3);
        Vpoly(idx,:) = polyval(c,tgrid);
    end
    pi.step;
end
pi.stop;