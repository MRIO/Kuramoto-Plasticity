clear

% kuramotoSheet Profiling
N = 50;
M = 50;

sigmoid = [10 0.5];
sigmoid2 = [10 0.5];
% No decay? Maybe for STDP [1 0]
decay = 0;
decay2 = 0;

% No training
training = [0 10*2*pi];
training_time = 0.5;
training_signal = [];

init_scaling = 0.5; % Initial connectivity used in these tests

switch 'uniform'
    case 'normal'
        oscillators = (randn(N*M,1)*0.05+10);
    case 'uniform'
        interval = [9.5 10.5];
        oscillators = interval(1) + (interval(2)-interval(1)).*rand(N*M,1);
    case 'multimodal normal'
        mu = [9 11 13];
        sigma = 0.05;
        oscillators = trimodalGaussian(N,M,mu,sigma);
end
oscillators = oscillators*2*pi;

connectivity = 'euclidean';
radius = 2;

scaling = 4*pi;
scaling2 = 4*pi;

stepval = 0.501;

% Try varying STDP functions/parameters and effect on clusters
plasticity = {'STDP' 0.3 [1 1] [0.002 0.002]};

GPU = 1;

time = 10;
steps = 1000;
tic
out1 = kuramotoSheet([N M],scaling,'record_adjacency',0,'GPU',GPU, 'plotme', 0,'connectivity', connectivity, 'radius', radius, 'oscillators',oscillators, 'time', time, 'dt', time/steps, 'plasticity', plasticity,'sigmoid',sigmoid, 'init_scaling',init_scaling, 'decay', decay,'training', training, 'training_signal', training_signal,'training_time', training_time, 'stepval',stepval);
toc