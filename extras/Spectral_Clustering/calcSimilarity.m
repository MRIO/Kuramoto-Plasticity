function similarity = calcSimilarity(phases)
% function similarity = calcSimilarity(phases)
% Simple similarity between two time series
% phases = N*M (N signals, M timesteps)
similarity = zeros(size(phases,1));
for ii = 1:size(phases,1)
    for jj = 1:size(phases,1)
        similarity(ii,jj)=sum(abs(bsxfun(@minus,phases(ii,:),phases(jj,:))));
    end
end

m = max(max(similarity));
if m~=0
    similarity = 1-similarity./max(max(similarity));
end
end