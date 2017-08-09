%% NOTES
% 
%plot(out2.phase')
%plot(diff(out2.phase)')
%figure
%subplot(2,1,1)
%plot(mod(out2.phase',2*pi))
%subplot(2,1,2)
%plot(mod(diff(out2.phase)',2*pi))

%%

clear x cc rho c
clf

%x = sin(out2.state(:,9001:end));
%x = diff(out2.phase(:,200:end)');
x = out1.state(:,9001:end);
cc = calcSimilarity(x);

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

c = spect_clust(cc,20);

%%

subplot(334)
plot(c.nc,c.gc,'-o'), axis tight
xlabel('Number of clusters')
ylabel('Eigenvalue gap')

subplot(335)
[~, ic] = sort(c.label(:,4));
pcolor(cc(ic,ic)), shading flat
title('Partial correlation')

subplot(336)
%fshow(c.label(:,2)) ?
plot(c.label(:,4), 'x'); %replacement of fshow function?

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


