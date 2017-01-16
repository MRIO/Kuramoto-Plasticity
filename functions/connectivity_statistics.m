function out = connectivity_statistics(in)

W = in.W;




stats.stdW = std(W(W~=0));
out.stats.connections = sum(W>0);
out.stats.meanweight  = mean(W(find(W>0)));

if exist('clustering_coef_wd')
    out.stats.clustercoeff.wd = clustering_coef_wd(W);
    out.stats.clustercoeff.bu = clustering_coef_bu(W);
else
    out.stats.clustercoeff = 'WARNING: did not find connectivity toolbox';
end

