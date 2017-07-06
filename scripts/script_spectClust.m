x = sin(out2.state(:,9001:end));
cc = partialcorr(x', mean(x)'); % partial corherelation w.r.t the mean

%%
subplot(331)
plot(x'), axis tight
title('Data')
subplot(332)
pcolor(cc), shading flat
title('Partial correlation')
subplot(333)
hist(cc(:),100), axis tight
title('Partial correlation')

%%

rho = (cc + 1)/2; % normalize rho to get an affinity matrix 0<=rho<=1
rho = (rho+rho')/2; % rho had better be symmetric
c = spect_clust(rho, 20);

%%

subplot(334)
plot(c.nc,c.gc,'-o'), axis tight
xlabel('Number of clusters')
ylabel('Eigenvalue gap')

subplot(335)
[~, ic] = sort(c.label(:,2));
pcolor(cc(ic,ic)), shading flat
title('Partial correlation')

subplot(336)
%fshow(c.label(:,2)) ?
plot(c.label(:,2), 'x'); %replacement of fshow function?

%% Check with information criteria

subplot(337)
[bic, aic] = baic(x, c);
plot([bic aic]), axis tight % These criterion give different results.. 
legend('BIC','AIC')
subplot(338)
[~,ic] = sort(c.label(:,10)); % Let's go for n = 10 as BIC suggests..
pcolor(cc(ic,ic)), shading flat % For some reasons, it doesn't seem to work..
subplot(339)
%fshow(c.label(:,10)) ?
plot(c.d(1:10),'*')
title('Eigenvalues')


