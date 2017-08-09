% Attempt at creating a map of clusters

%% Try 1
% Base cluster # on the smallest largest eigenvalue gap
N = 35;M=35;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% As i forgot to add sigVal/TauVal in previous parameterAnalysis code,
% this'll have to do
sigP = [0 10 35];
TauP = [0.01 1 35];

sigVal = sigP(1) + (sigP(2) - sigP(1))./sigP(3).*[1:sigP(3)];
TauVal = TauP(1) + (TauP(2) - TauP(1))./TauP(3).*[1:TauP(3)];

[Tau,sig] = meshgrid(TauVal,sigVal);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

map(1).sig  = sig;
map(1).Tau  = Tau;
map(1).Z    = zeros(N,M);

map(2).sig  = sig;
map(2).Tau  = Tau;
map(2).Z    = zeros(N,M);

map(3).sig  = sig;
map(3).Tau  = Tau;
map(3).Z    = zeros(N,M);

% Results sometimes has 40+ clusters, add a filter?

for ii=1:N
    for jj=1:M
        idx = out1(ii,jj).c.nc>20;
        out1(ii,jj).c.nc(idx)=[];
        out1(ii,jj).c.gc(idx)=[];
        [v,i] = min(out1(ii,jj).c.gc);
        map(1).Z(ii,jj) = out1(ii,jj).c.nc(i);
        
        idx = out2(ii,jj).c.nc>20;
        out2(ii,jj).c.nc(idx)=[];
        out2(ii,jj).c.gc(idx)=[];
        [v,i] = min(out2(ii,jj).c.gc);
        map(2).Z(ii,jj) = out2(ii,jj).c.nc(i);
        
        idx = out3(ii,jj).c.nc>20;
        out3(ii,jj).c.nc(idx)=[];
        out3(ii,jj).c.gc(idx)=[];
        [v,i] = min(out3(ii,jj).c.gc);
        map(3).Z(ii,jj) = out3(ii,jj).c.nc(i);
    end
end
