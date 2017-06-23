%% 0. Connectivity test
%oscillators = ones(100,1)*10*2*pi;
oscillators = [];
oscillators = (randn(4,1)*1+7)*2*pi;

plasticity = {'STDP' 0.1 [1 -1] [20 20]}; %try different parameters STDP
plasticity2 = {'STDP' 0.1 [1 -1] [20 20]}; %try different parameters STDP
%plasticity = {'seliger' 10 100};

% scaling 10
%sigmoid = [1.25 5];
%sigmoid2 = [1.25 5];

% scaling 1
sigmoid = [10 0.5];
sigmoid2 = [10 0.5];
init_scaling = 0;
decay = 0.1;
decay2 = 0;

time = 1;
steps = 1000;
out1 = kuramotoSheet([2 2],1,'plotme', 0,'connectivity', 'all to all', 'oscillators',oscillators, 'time', time, 'dt', time/steps, 'plasticity', plasticity,'sigmoid',sigmoid, 'init_scaling',init_scaling, 'decay', decay);
out2 = kuramotoSheet([2 2],1,'plotme', 0,'connectivity', 'all to all', 'oscillators',oscillators, 'time', time, 'dt', time/steps, 'plasticity', plasticity2,'sigmoid',sigmoid2, 'init_scaling',init_scaling, 'decay', decay2);

close all

% oscmean = mean(out1.oscillators);
% oscstd = std(out1.oscillators);
% find( (out1.oscillators > (oscmean+2.5*oscstd)) | (out1.oscillators < (oscmean-2.5*oscstd)))

figure('Position', [350 300 1250 500])
subplot(1,2,1)
imagesc(out1.adjacency{end})
colorbar
subplot(1,2,2)
imagesc(out2.adjacency{end})
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

%plasticity = {'seliger' 10 100};
plasticity = {'STDP' 1 [1 1] [10 10]}; %try different parameters STDP

training = [100 7*2*pi];
training_time = 0.5;

% Signal %
training_signal = zeros(10);
switch 4
    case 1
        training_signal(:,1:3)=pi;
        training_signal(:,7:10)=4*pi/3;
    case 2
        training_signal(6:10,1:5)=pi/2;
        training_signal(1:5,6:10)=pi;
        training_signal(6:10,6:10)=3*pi/2;
    case 3
        training_signal = fliplr(checkerboard(5,1,1) > 0.5)*pi;
    case 4
        training_signal(:,1:5) = pi;
end
training_signal = reshape(training_signal,100,1);
%%%%%%%%%%

out2 = kuramotoSheet([N M],1,'plotme',0,'time',time,'connectivity', connectivity,'oscillators',oscillators,'omega_mean',omega_mean,'omega_std',omega_std,'init_cond',init_cond,'plasticity',plasticity,'training', training, 'training_signal', training_signal,'training_time', training_time);
%out1 = kuramotoSheet([N M],1,'plotme',0,'connectivity', connectivity,'oscillators',oscillators,'init_cond',init_cond,'plasticity',plasticity);

figure
imagesc(out2.connectivity);
colorbar
%% - Training STDP test
disp('old crap')
return
%oscillators = ones(100,1)*10*2*pi;
oscillators = [];
oscillators = (randn(N*M,1)*1+7)*2*pi;

plasticity = {'STDP' 1 [1 0] [1 1]}; %try different parameters STDP
%plasticity = {'seliger' 10 100};

sigmoid = [1.25 5];
sigmoid2 = [1.25 5];
init_scaling = 0;

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
    case 4
        training_signal(:,1:5) = pi;
end
training_signal = reshape(training_signal,100,1);
%%%%%%%%%%

time = 2;
steps = 2000;
out1 = kuramotoSheet([10 10],10,'plotme', 0,'connectivity', 'all to all', 'oscillators',oscillators, 'time', time, 'dt', time/steps, 'plasticity', plasticity,'sigmoid',sigmoid, 'init_scaling',init_scaling,'training', training, 'training_signal', training_signal,'training_time', training_time);

close all

% oscmean = mean(out1.oscillators);
% oscstd = std(out1.oscillators);
% find( (out1.oscillators > (oscmean+2.5*oscstd)) | (out1.oscillators < (oscmean-2.5*oscstd)))

figure(100)
imagesc(out1.adjacency{end})
colorbar

%% - Training STDP decay test
%oscillators = ones(100,1)*10*2*pi;
%oscillators = [];
oscillators = (randn(100,1)*0.05+10)*2*pi;

scaling = 4*pi;
scaling2 = 4*pi;

% scaling 4*pi, tau 002 latest optimal values for 3cluster 0625

plasticity = {'STDP' 0.1 [1 1] [0.002 0.002]}; %try different parameters STDP
plasticity2 = {'STDP' 0.1 [1 0] [0.002 0.002]}; %try different parameters STDP

sigmoid = [10 0.5];
sigmoid2 = [10 0.5];
init_scaling = 0;
decay = 0;
decay2 = 0;

training = [0 10*2*pi];
training_time = 0.5;

% Signal %
training_signal = zeros(10);
switch 2
    case 1 % 3cluster
        training_signal(:,1:3)=2/3*pi;
        training_signal(:,7:10)=4/3*pi;
    case 2 % 4cluster
        training_signal(6:10,1:5)=pi/2;
        training_signal(1:5,6:10)=pi;
        training_signal(6:10,6:10)=3*pi/2;
    case 3
        training_signal = fliplr(checkerboard(5,1,1) > 0.5)*pi;
    case 4 % 2cluster
        training_signal(:,1:5) = pi;
        
end
training_signal = reshape(training_signal,100,1);
%%%%%%%%%%


time = 10;
steps = 10000;
out1 = kuramotoSheet([10 10],scaling,'plotme', 0,'connectivity', 'all to all', 'oscillators',oscillators, 'time', time, 'dt', time/steps, 'plasticity', plasticity,'sigmoid',sigmoid, 'init_scaling',init_scaling, 'decay', decay,'training', training, 'training_signal', training_signal,'training_time', training_time);
out2 = kuramotoSheet([10 10],scaling2,'plotme', 0,'connectivity', 'all to all', 'oscillators',oscillators, 'time', time, 'dt', time/steps, 'plasticity', plasticity2,'sigmoid',sigmoid2, 'init_scaling',init_scaling, 'decay', decay2,'training', training, 'training_signal', training_signal,'training_time', training_time);

close all

% oscmean = mean(out1.oscillators);
% oscstd = std(out1.oscillators);
% find( (out1.oscillators > (oscmean+2.5*oscstd)) | (out1.oscillators < (oscmean-2.5*oscstd)))

figure('Position', [350 300 1250 500])
subplot(1,2,1)
imagesc(out1.adjacency{end})
colorbar
subplot(1,2,2)
imagesc(out2.adjacency{end})
colorbar

%% STDP June 20+
%parameters to vary: scaling, oscillators (distribution parametrize),
%connectivity, radius, plasticity

sigmoid = [10 0.5];
sigmoid2 = [10 0.5];
% No decay? Maybe for STDP [1 0]
decay = 0;
decay2 = 0;

% No training
training = [0 10*2*pi];
training_time = 0.5;
training_signal = [];

init_scaling = 1; % Initial connectivity used in these tests

% Change distribution (uniform? bi/multi-modal gaussian? See notes), see
% effect on cluster formation

switch 'uniform'
    case 'normal'
        oscillators = (randn(100,1)*0.05+10);
    case 'uniform'
        interval = [5 15];
        oscillators = interval(1) + (interval(2)-interval(1)).*rand(100,1);
    case 'multimodal normal'
        %%% create multimodal
        %mu = [7.5 10 12.5];
        %sigma = 0.05;
        %oscillators = normpdf(
end
oscillators = oscillators*2*pi;

connectivity = 'euclidean';
radius = 3;

% scaling 4*pi, tau 002 latest optimal values for 3cluster 0625
scaling = 4*pi;
scaling2 = 4*pi;

stepval = 0.501;

% Try varying STDP functions/parameters and effect on clusters
plasticity = {'STDP' 0.2 [1 1] [0.002 0.002]};
plasticity2 = {'STDP' 0.2 [1 1] [0.002 0.002]};

time = 10;
steps = 10000;
out1 = kuramotoSheet([10 10],scaling,'plotme', 0,'connectivity', connectivity, 'radius', radius, 'oscillators',oscillators, 'time', time, 'dt', time/steps, 'plasticity', plasticity,'sigmoid',sigmoid, 'init_scaling',init_scaling, 'decay', decay,'training', training, 'training_signal', training_signal,'training_time', training_time, 'stepval',stepval);
out2 = kuramotoSheet([10 10],scaling2,'plotme', 0,'connectivity', 'all to all', 'radius', radius, 'oscillators',oscillators, 'time', time, 'dt', time/steps, 'plasticity', plasticity2,'sigmoid',sigmoid2, 'init_scaling',init_scaling, 'decay', decay2,'training', training, 'training_signal', training_signal,'training_time', training_time, 'stepval',stepval);

close all

if 0
    figure('Position', [50 100 1250 500])
else
    figure('Position', [350 300 1250 500])
end
subplot(1,2,1)
imagesc(out1.adjacency{end})
colorbar
subplot(1,2,2)
imagesc(out2.adjacency{end})
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
%Improve/replace

cgobj = clustergram(flipud(out.connectivity), 'cluster', 'row');
clusterlabels = cellfun(@str2num, cgobj.ColumnLabels);
clustered_state = out.state(clusterlabels,:);

%% 0.1.1 Replay data
clear idx XX YY SS MOV PP
N=10;M=10;
%PP = out1.state(:,9001:end);
%PP = test_STDP_Random05_training{2}.state;
%PP = test_STDP_Random05_4cluster{3}.state;

PP = out1.state(:,9001:end);
%PP = result.state;
%PP = clustered_state;

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

%% STDP Function plots
% Sigmoid
clear
close all
figure
a = 10;
c = 0.5;
t = linspace(-1,2,1000);

y = sigmf(t, [a c]);

plot(t,y)
xlabel('Original connectivity')
ylabel('Adjusted connectivity')
title('Sigmoid Connectivity Adjustment')

% STDP function
figure
dTimePos = 1;
dTimeNeg = 1;
A = [0.5 0.5];
tau = [0.002 0.002];
t = linspace(-0.02,0.02,1000);
ts = t>0;

y = ts.*A(1).*exp(- t./tau(1)) -~ts.*A(2).*exp( t./tau(2));

plot(t,y)
xlabel('Time between spikes (s)')
ylabel('Weight change')
title('STDP Weight Change Function')

%% Generate parameter test
%return
clear

oscillators = [];
sigmoid = [1.25 5];
init_scaling = 0;

time = 10;
steps = 10000;

tau_range = [0.1 0.2 0.5 1 5 10 20 50 100 250];
tic
for j=1:10
    parfor i=1:10
        
        plasticity = {'STDP' 1 [1 0] [tau_range(i) tau_range(j)]};
        output = kuramotoSheet([10 10],10,'plotme', 0,'connectivity', 'all to all', 'oscillators',oscillators, 'time', time, 'dt', time/steps, 'plasticity', plasticity,'sigmoid',sigmoid, 'init_scaling',init_scaling);
        results_parameters{i,j}= output.parameters;
        results_adjacency{i,j} = output.adjacency{end};
        
        output = [];
        
    end
end
toc

%% Collect

results{1} = results_parameters;
results{2} = results_adjacency;

%% Remove Unnecessary data from result
result = rmfield(out1,{'adjacency';'orderparameter';'phase';'meanphase'})

%% Plot parameter test
close all

figure
for i=1:10
    subplot(3,4,i)
    imagesc(results_adjacency{5,i})
    colorbar
end