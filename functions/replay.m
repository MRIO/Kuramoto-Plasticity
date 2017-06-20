function output = replay(PP, range)
% Short function to replay the kuramotoSheet results
if range~='all'
    PP = PP(:,(end-range):end);
end

N=10;M=10;

f = figure(100);

[XX, YY] = meshgrid([1:M],[1:N]);
XX = XX(:); YY = YY(:);

for t = 2:size(PP,2)
    if ~ishghandle(f)
        break;
    end
    SS = reshape(PP(:,t),N,M);
    cla
    % imagesc(reshape(theta_t(:,t),N,M))
    mesh(SS); hold on
    scatter3(XX, YY ,PP(:,t),60,'filled')
    title('phase')
    caxis([0 2*pi])
    zlim([-3 3])
    
    drawnow
end
end