%% 0. Connectivity test
checker = fliplr(checkerboard(50,1,1) > 0.5);
%oscillators2 = (double(1:100<51)+10)*2*pi;
oscillators2 = double(1:100<51)*pi;
oscillators = ones(1,100)*10*2*pi;

plasticity = {'STDP' 1 pi};
plasticity = {'seliger' 10 100};
plasticity2 = {'STDP' 1 0.5};

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
imagesc(out1.connectivity)
colorbar
subplot(1,2,2)
imagesc(out2.connectivity)
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
initial_condition = out1.phase(clusterlabels,1)

out_cluster = kuramotoSheet([10 10],1,'plotme', 1,'init_cond', initial_condition, 'connectivity', 'all to all', 'oscillators',oscillators', 'time', time, 'dt', time/steps, 'plasticity', plasticity)
%% 0.2 Rate of change
d_adjacency = zeros(1,steps);
for i=2:steps
    d_adjacency(i) = sum(sum(abs(out2.adjacency{i}-out2.adjacency{i-1})));
end

%% 1 Temp adjacency matrix
figure
N = 4; M = 4; radius = 1;
[X Y] = meshgrid([1:N],[1:M]);
X = X(:); Y = Y(:); 

% # compute adjacency matrix
connectivity = squareform( pdist([X Y], 'euclidean') <= radius );

imagesc(connectivity)
disp('a')
        
%% 2 Temp

connectivity = fliplr(checkerboard(50,1,1) > 0.5);