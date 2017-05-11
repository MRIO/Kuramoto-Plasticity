%% 0. Connectivity test
oscillators = ones(1,100)*10*2*pi;
oscillators = [];

plasticity = {'STDP' 1 [1 1] [60 60]}; %try different parameters STDP
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

%% - Initial Conditions test
N=10;M=10;
oscillators = rand(N*M,1)*2*pi*10; %use default, oscillators = [];
% oscillators{2} = ones(N*M,1)*2*pi*10;
% oscillators{3} = (reshape(fliplr(checkerboard(5,1,1) > 0.5),100,1)*10+5)*2*pi;
% oscillators{4} = zeros(100,1);
% oscillators{4}(51:end,1)=5*2*pi;
% oscillators{4}=oscillators{4}+10*2*pi;


init_cond{1} = rand(N*M,1)*2*pi;
init_cond{2} = zeros(N*M,1);
init_cond{3} = (reshape(fliplr(checkerboard(5,1,1) > 0.5),100,1))*pi+pi/2;
init_cond{4} = zeros(N*M,1);
init_cond{4}(51:end,1)=pi;
init_cond{4} = init_cond{4}+pi/2;

connectivity = zeros(N*M);

plasticity = {'seliger' 10 100};

for i=1:4
    out(i) = kuramotoSheet([N M],1,'plotme',0,'connectivity', connectivity,'oscillators',oscillators,'init_cond',init_cond{i},'plasticity',plasticity);
end

figure
for i=1:4
    subplot(2,2,i)
    imagesc(out(i).connectivity);
    colorbar
end

%% - Training Test
N=10;M=10;
time = 1;

%oscillators = 7*2*pi;

oscillators = [];
omega_mean = 7; omega_std = 2;
init_cond = rand(N*M,1)*2*pi;
connectivity = zeros(N*M);

plasticity = {'seliger' 10 100};

training = [100 7*2*pi];
training_time = 0.5;

% Signal %
training_signal = zeros(10);
switch 1
    case 1
training_signal(:,1:3)=pi;
training_signal(:,7:10)=4*pi/3;
    case 2
training_signal(6:10,1:5)=pi/2;
training_signal(1:5,6:10)=pi;
training_signal(6:10,6:10)=3*pi/2;
    case 3
        training_signal = fliplr(checkerboard(5,1,1) > 0.5)*pi;
end
training_signal = reshape(training_signal,100,1);
%%%%%%%%%%

out = kuramotoSheet([N M],1,'plotme',0,'time',time,'connectivity', connectivity,'oscillators',oscillators,'omega_mean',omega_mean,'omega_std',omega_std,'init_cond',init_cond,'plasticity',plasticity,'training', training, 'training_signal', training_signal,'training_time', training_time);
%out1 = kuramotoSheet([N M],1,'plotme',0,'connectivity', connectivity,'oscillators',oscillators,'init_cond',init_cond,'plasticity',plasticity);

figure
imagesc(out.connectivity);
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
%Improve?

cgobj = clustergram(flipud(out.connectivity), 'cluster', 'row');
clusterlabels = cellfun(@str2num, cgobj.ColumnLabels);
clustered_state = out.state(clusterlabels,:);

%% 0.1.1 Replay data
clear idx XX YY SS MOV PP
N=10;M=10;
    %PP = clustered_state;
    PP = out.state;

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
			scatter3(XX, YY ,PP(:,t),60,'filled')
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

connectivity = fliplr(checkerboard(5,1,1) > 0.5);
imagesc(connectivity);

%%
a = -1; b = -0.0001;

sign(b)>sign(a)