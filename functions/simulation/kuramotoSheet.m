function out = kuramotoSheet(varargin)
% function out = kuramotoSheet([N,M], K, 'Param', 'Value')]
%
% Calculates the sheet of kuramoto oscillators with different connectivity schemes.
% User can select the oscillator frequencies, initial conditions, and
% time base (dt for the forward euler solver)
% The script also outputs the synchronization (kuramoto parameter), i.e., the
% centroid of the phase of the group of oscillators.
%
% inputs:
% 	networksize, [N M]
% 	scaling for coupling - K
% outputs:
% 	out.state = sin ( theta_t) ;
% 	out.phase = phase;
% 	out.parameters = parameters
% 	out.oscillators = intrinsic frequency of oscillators;
% 	out.orderparameter = synchrony measure (kuramoto parameter)
% 	out.meanphase = mean phase of oscillators at t
% 	out.seed = random seed
%
% parameter value pairs:
%
% ('radius', 4)
% ('dt', 1e-3)
% ('simtime',1) % in seconds
% ('omega_mean', 10)
% ('plotme', 1)
% ('noise', 0)
% ('connectivity', [])
% connectivity  = 'inverse dist';
% connectivity  = 'euclidean';
% connectivity = 'inverse dist';
% connectivity = 'chebychev';
%
% kuramotoSheet([50 50],105)
%
% author: m@negrello.org
% all rights to kuramoto and mathworks, all wrongs are mine ;D

anim = 1; makemovie = 0;

% [=================================================================]
%  parse inputs
% [=================================================================]

p = inputParser;
p.addRequired('networksize')
p.addRequired('scaling')

p.addParameter('radius', 3)
p.addParameter('dt', 1e-3)
p.addParameter('time',1) % in seconds
p.addParameter('omega_mean', 7) % in Hz
p.addParameter('init_cond', []) % in Hz
p.addParameter('omega_std', 2) % in Hz
p.addParameter('oscillators', [])
p.addParameter('plotme', 1)
p.addParameter('connectivity', 'euclidean')  % adjacency matrix
p.addParameter('seed', 0)
p.addParameter('plasticity', {0 1 1}, @iscell) % See implementation

%%% Needs testing
p.addParameter('training',[0 10*2*pi]) %(1) - Training amplitude, (2) - training frequency
p.addParameter('training_signal',[]) % Training Phases
p.addParameter('training_time',0.25) % in seconds
p.addParameter('sigmoid',[0 5]); % [1.25 5] gives range between 0 and 10
p.addParameter('init_scaling', 1); %Scales the initial connectivity, different from 'scaling' which only affects the connectivity after sigmoid has been applied.
p.addParameter('decay',0) % decay/second
p.addParameter('stepval',0.20001);
p.addParameter('record_adjacency',1);
p.addParameter('GPU',0);

p.parse(varargin{:});

netsize = p.Results.networksize;
scaling = p.Results.scaling;
radius  = p.Results.radius;
connectivity = p.Results.connectivity;
dt = p.Results.dt;
simtime = p.Results.time;
omega_mean = p.Results.omega_mean;
omega_std = p.Results.omega_std;
plotme = p.Results.plotme;
seed = p.Results.seed;
init_cond = p.Results.init_cond;
oscillators = p.Results.oscillators;
plasticity = p.Results.plasticity;
plasticity_type = plasticity{1};
plasticity{1} = [];
plasticity_data = cell2mat(plasticity);
training = p.Results.training;
training_signal = p.Results.training_signal;
training_time = p.Results.training_time;
sigmoid = p.Results.sigmoid;
init_scaling = p.Results.init_scaling;
decay = p.Results.decay;
stepval = p.Results.stepval;
record_adjacency = p.Results.record_adjacency;
gpu = p.Results.GPU;

N = netsize(1);
M = netsize(2);
NO = prod(netsize);
idx = ones(1,NO);

% [=================================================================]
%  randomize oscillator intrinsic frequencies
% [=================================================================]

rng(seed,'twister')
if isvector(oscillators)
    omega_i = oscillators;
else
    omega_i = (randn(N*M,1)*omega_std+omega_mean)*2*pi;
end

scale_to_intrinsic_freq = 0;

% [=================================================================]
%  connectivity
% [=================================================================]

if ischar(connectivity)
    switch connectivity
        case 'all to all'
            connectivity = ones(N*M) - eye(N*M);
            
        case 'null'
            connectivity = zeros(N*M);
            
        case 'chebychev'
            
            [X Y] = meshgrid([1:N],[1:M]);
            X = X(:); Y = Y(:);
            
            % # compute adjacency matrix
            connectivity = squareform( pdist([X Y], 'chebychev') <= radius );
            
        case 'euclidean'
            
            [X Y] = meshgrid([1:N],[1:M]);
            X = X(:); Y = Y(:);
            
            % # compute adjacency matrix
            connectivity = squareform( pdist([X Y], 'euclidean') <= radius );
            
            
        case 'inverse dist'
            depth   = [1:N];
            breadth = [1:M];
            
            [X Y] = meshgrid(depth,breadth);
            X = X(:); Y = Y(:);
            
            % # compute adjacency matrix
            connectivity = 1./squareform( pdist([X Y], 'euclidean') );
            
            
            
        case 'random'
            % W = (ones(N*M)-eye(N*M) ) .* rand(N*M);
    end
end

connectivity = connectivity.*init_scaling;

% ensure that there are no self connections
connectivity(logical(eye(M*N))) = 0;

% [=================================================================]
%  Training/Input
% [=================================================================]

if ~isvector(training_signal)
    training_signal = zeros(N*M,1);
end

% [=================================================================]
%  randomize initial condition
% [=================================================================]
%initial condition (initial phase)
if isempty(init_cond)
    theta_t(:,1) = rand(N*M,1)*2*pi;
else
    % disp('using initial condition (assuming between 0 and 2pi)')
    theta_t(:,1) = init_cond;
end

if gpu
    theta_t = gpuArray(theta_t);
    connectivity = gpuArray(connectivity);
    
    %%% New
    scaling = gpuArray(scaling);
    omega_i = gpuArray(omega_i);
    dt = gpuArray(dt);
    training = gpuArray(training);
    training_signal = gpuArray(training_signal);
    training_time = gpuArray(training_time);
    plasticity_data = gpuArray(plasticity_data);
    simtime = gpuArray(simtime);
    
    spikeTime = zeros(N*M,1,'gpuArray');
    requiresUpdate = zeros(N*M,'gpuArray');
    PP = zeros(size(theta_t),'gpuArray');
    dw = zeros(1,simtime/dt,'gpuArray');
    dwabs = zeros(1,simtime/dt,'gpuArray');
else
    spikeTime = zeros(N*M,1);
    requiresUpdate = zeros(N*M);
    PP = zeros(size(theta_t));
    dw = zeros(1,simtime/dt);
    dwabs = zeros(1,simtime/dt);
end

% [=================================================================]
%  simulate
% [=================================================================]

if gpu %Find alternative to IF statement
    sConnectivity = arrayfun(@sigmoidFun,connectivity,sigmoid(1),sigmoid(2),stepval);
else
    sConnectivity = sigmoidFun(connectivity,sigmoid(1),sigmoid(2),stepval);
end

if record_adjacency
    dt2=gather(dt);
    simtime2=gather(simtime);
    adjacency = cell(1,simtime2/dt2);
    adjacency{1} = sConnectivity;
end

if plotme; f = figure(100); a(1) = subplot(121);a(2) = subplot(122); end
for t = 2:simtime/dt
    disp(t)
    phasedifferences = bsxfun(@minus, theta_t(:,t-1)',theta_t(:,t-1));
    
    phasedifferences_W = sConnectivity.*scaling.*sin(phasedifferences);
    
    summed_sin_diffs = mean(phasedifferences_W,2);
    
    
    if gpu % i know...
        theta_t(:,t) = arrayfun(@computeKuramoto,theta_t(:,t-1), omega_i, summed_sin_diffs, training_time, training(1),training(2), training_signal,t,dt);
    else
        theta_t(:,t) = theta_t(:,t-1) + dt.*( omega_i + summed_sin_diffs + ((t.*dt)<training_time).*training(1).*sin(theta_t(:,t-1)-training(2).*dt.*t-training_signal(:,1)) );
    end
    
    PP(:,t) = sin(mod(theta_t(:,t),2.*pi));
    
    switch plasticity_type
        case 'seliger' %{2} - epsilon, {3} - alpha
            % Currently ignores spatial distance between oscillators, should
            % it?
            connectivity = connectivity + dt.*plasticity_data(1).* ...
                ( plasticity_data(2) .* cos(phasedifferences) - connectivity );
            
        case 'STDP' % {2} - delta-adjustment: factor affecting weight change, {3} - [A1;A2], {4} - [tau1;tau2]
            %%% Does not yet account for sign(0)=0, also does not account for bsxfun result deltaTime=0 (as 0 is used to ignore elements)
            
            % Determines upward zero crossover (= spike) of each oscillator
            % and calculates time difference between all spikes
            upwardZeroCross = sign(mod(theta_t(:,t),2.*pi)-pi) > sign(mod(theta_t(:,t-1),2.*pi)-pi);
            spikeTime(upwardZeroCross) = t.*dt;
            
            deltaTime = bsxfun(@minus, spikeTime,spikeTime'); %deltaTime = t_i - t_j
            
            % Spiked oscillators are memorized and all updated oscillator
            % couples are determined
            requiresUpdate(:,upwardZeroCross) = 1;
            updateMatrix = requiresUpdate.*requiresUpdate';
            updateMatrix(logical(eye(N.*M))) = 0;
            requiresUpdate = requiresUpdate - updateMatrix;
            
            % Actual implementation of Spike timing-dependent plasticity
            % function, see http://www.scholarpedia.org/article/Spike-timing_dependent_plasticity#Basic_STDP_Model
            deltaTime = deltaTime.*updateMatrix;
            dTimePos = deltaTime > 0;
            dTimeNeg = deltaTime < 0;
            
            if gpu
                W = arrayfun(@computeWeights,dTimePos,dTimeNeg,plasticity_data(2),plasticity_data(3), ...
                    plasticity_data(4),plasticity_data(5),deltaTime);
            else
            W = dTimePos.*plasticity_data(2).*exp( -deltaTime./plasticity_data(4)) ...
                -dTimeNeg.*plasticity_data(3).*exp( deltaTime./plasticity_data(5));
            end
            
            
            
            dw(t) = sum(sum(W));
            dwabs(t) = sum(sum(abs(W)));
            
            % Update connectivity
            connectivity = connectivity + plasticity_data(1).*W;
            
        case 'test'
            connectivity = connectivity + dt*plasticity_data(1)*((abs(phasedifferences)<=plasticity_data(2)).*phasedifferences - connectivity);
            
        case 'null'
        otherwise
            disp('Running without a plasticity rule.')
            plasticity_type = 'null';
    end
    
    connectivity = connectivity - decay.*dt;
    connectivity(connectivity<0)=0;
    
    if gpu %Find alternative to IF statement
        sConnectivity = arrayfun(@sigmoidFun,connectivity,sigmoid(1),sigmoid(2),stepval);
    else
        sConnectivity = sigmoidFun(connectivity,sigmoid(1),sigmoid(2),stepval);
    end
    
    if record_adjacency
        adjacency{t} = sConnectivity;
    end
    
end


% [=================================================================]
%  plots
% [=================================================================]



if plotme
    
    ffff = figure
    
    
    subplot(2,2,1)
    plot(linspace(0,simtime, simtime*dt^-1), mod(theta_t,2*pi)')
    ylabel('phase (theta)')
    subplot(2,2,2)
    imagesc(connectivity), colorbar
    title('connectivity')
    subplot(2,2,3)
    plot(linspace(0,simtime, simtime*dt^-1), sin(theta_t'))
    ylabel('sin(theta)')
    xlabel('seconds')
    subplot(2,2,4)
    hist(omega_i/(2*pi),N)
    ylabel('#')
    xlabel('frequency (Hz)')
    
    figure(f)
    axes(a(2))
    line(repmat(linspace(0,simtime, simtime*dt^-1), length(unique(idx)), 1)', [abs(k)]')
    
    
    title('kuramoto parameter')
    xlabel('time (ms)')
    ylim([0 1])
    
    [XX YY] = meshgrid([1:M],[1:N]);
    XX = XX(:); YY = YY(:);
    
    
    if exist('linspecer','file')
        LSpec = linspecer(length(unique(idx)))
        set(ffff, 'colormap',   LSpec);
        set(a(1), 'colororder', LSpec);
        set(a(2), 'colororder', LSpec);
    else
        set(ffff, 'colormap',   jet(length(unique(idx))));
        set(a(1), 'colororder', jet(length(unique(idx))));
        set(a(2), 'colororder', jet(length(unique(idx))));
    end
    
    
    while anim
        
        for t = 2:simtime/dt
            if ~ishghandle(a(1))
                break;
            end
            SS = reshape(PP(:,t),N,M);
            axes(a(1)); 			cla
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
        
        if ~ishghandle(a(1))
            break;
        end
        anim = input(['repeat? ; 0 - no, 1 - yes \n'  ])
    end
    
end


% [================================================]
%  Write out
% [================================================]


if gpu
    theta_t = gather(theta_t);
    connectivity = gather(connectivity);
    
    PP = gather(PP);
    omega_i = gather(omega_i);
    dw = gather(dw);
    dwabs = gather(dwabs);
end

out.state = PP;
out.phase = theta_t;
out.parameters = p.Results;
out.oscillators = omega_i/(2*pi);
out.connectivity = connectivity;
if record_adjacency
    out.adjacency = adjacency;
end
out.seed = seed;
if makemovie && plotme ; out.movie = MOV; end
out.dw = dw;
out.dwabs = dwabs;
% out.all =  sin(mod(theta_t,2*pi));
end

function result = sigmoidFun(x,a,c,stepval)
result = (1./(1+exp(-a.*(x-c)))) .* (x>stepval);
% sigmoid .* heaviside
end

function result = computeKuramoto(prevTheta, omega_i, summed_sin_diffs, training_time, training1,training2, training_signal,t,dt)
result = prevTheta + dt.*( omega_i + summed_sin_diffs + ((t.*dt)<training_time).*training1.*sin(prevTheta-training2.*dt.*t-training_signal) );
end

function W = computeWeights(dTimePos,dTimeNeg,pData2,pData3,pData4,pData5,deltaTime)
    W = dTimePos.*pData2.*exp( -deltaTime./pData4) -dTimeNeg.*pData3.*exp( deltaTime./pData5); 
end