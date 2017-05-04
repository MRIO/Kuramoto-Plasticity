%% 0. Connectivity test
oscillators = ones(1,100)*10*2*pi;
oscillators = [];

plasticity = {'test' 1 10*pi};
%plasticity = {'seliger' 10 100};
plasticity2 = {'test' 1 0.5};

time = 1;
steps = 1000;
out1 = kuramotoSheet([10 10],1,'plotme', 0,'connectivity', 'all to all', 'oscillators',oscillators', 'time', time, 'dt', time/steps, 'plasticity', plasticity)
%out2 = kuramotoSheet([10 10],1,'plotme', 0,'connectivity', 'all to all', 'oscillators',oscillators', 'time', time, 'dt', time/steps, 'plasticity', plasticity2)

close all

% oscmean = mean(out1.oscillators);
% oscstd = std(out1.oscillators);
% find( (out1.oscillators > (oscmean+2.5*oscstd)) | (out1.oscillators < (oscmean-2.5*oscstd)))

figure('Position', [350 300 1250 500])
subplot(1,2,1)
imagesc(out1.adjacency{1})
colorbar
subplot(1,2,2)
imagesc(out1.connectivity)
colorbar

%% 0.0 Timelapse Adjacency
close all
f = figure
for i=1:steps
    if ~ishghandle(f)
        break;
    end
    
    cla;
    imagesc(out1.adjacency{i});
    colorbar
    
    val = num2str(round(i/steps*100)) ;
    title([val '%'])
    drawnow;
end

%% 0.1 HeatMap Cluster
%Continue;

%Improve, see log

cgobj = clustergram(flipud(out1.connectivity), 'cluster', 'row');
clusterlabels = cellfun(@str2num, cgobj.ColumnLabels);
clustered_state = out1.state(clusterlabels,:);

%% 0.1.1 Replay data
clear idx XX YY SS MOV PP
    PP = clustered_state;

    f = figure(100);
    makemovie = 0;
    
    idx = ones(1,N*M);
    	
	[XX, YY] = meshgrid([1:M],[1:N]);
	XX = XX(:); YY = YY(:);
	

% 	if exist('linspecer')
% 		LSpec = linspecer(length(unique(idx)))
% 		set(a(1), 'colororder', LSpec);
%     else
% 		set(a(1), 'colororder', jet(length(unique(idx))));
%     end
		
		for t = 2:size(PP,2)
            if ~ishghandle(f)
                break;
            end
			SS = reshape(PP(:,t),N,M);
			cla
			% imagesc(reshape(theta_t(:,t),N,M))
			mesh(SS); hold on
			scatter3(XX, YY ,PP(:,t),60,LSpec(idx,:),'filled')
			title('phase')
			caxis([0 2*pi])
			zlim([-3 3])

			drawnow
			if makemovie
                % TODO first frame of movie missing t=1???
				MOV(t-1) = getframe(f);
			end

        end

%% 2 Temp adjacency matrix
figure
N = 4; M = 4; radius = 1;
[X Y] = meshgrid([1:N],[1:M]);
X = X(:); Y = Y(:); 

% # compute adjacency matrix
connectivity = squareform( pdist([X Y], 'euclidean') <= radius );

imagesc(connectivity)
disp('a')
        
%% 3 Temp

connectivity = fliplr(checkerboard(50,1,1) > 0.5);

%%
a = -1; b = -0.0001;

sign(b)>sign(a)