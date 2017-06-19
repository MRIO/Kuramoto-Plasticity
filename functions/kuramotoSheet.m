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

%gpu = 0;

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
	p.addParameter('noise', 0) 
	p.addParameter('connectivity', 'euclidean')  % adjacency matrix
	p.addParameter('clusterize', [0 0 0 0],@isvector)
	p.addParameter('seed', 0)
    p.addParameter('plasticity', {0 1 1}, @iscell) % See implementation
    
    %%% Needs testing
    p.addParameter('training',[0 10*2*pi]) %(1) - Training amplitude, (2) - training frequency
    p.addParameter('training_signal',[]) % Training Phases
    p.addParameter('training_time',0.25) % in seconds
    p.addParameter('sigmoid',[0 5]); % [1.25 5] gives range between 0 and 10
    p.addParameter('init_scaling', 10); %Scales the initial connectivity, different from 'scaling' which only affects the connectivity after sigmoid has been applied.
    p.addParameter('decay',0) % decay/second

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
	noise = p.Results.noise;
	clusterize = p.Results.clusterize;
	seed = p.Results.seed;
	init_cond = p.Results.init_cond;
	oscillators = p.Results.oscillators;
    plasticity = p.Results.plasticity;
    training = p.Results.training;
    training_signal = p.Results.training_signal;
    training_time = p.Results.training_time;
    sigmoid = p.Results.sigmoid;
    init_scaling = p.Results.init_scaling;
    decay = p.Results.decay;

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


% [=================================================================]
%  this
% [=================================================================]


if clusterize(1)
	W = connectivity;

    groupsize = clusterize(2);

    k = round(NO/groupsize);
    [idx] = kmeans([X Y], k);

    ProbCluster = clusterize(3);
    ProbOriginal = clusterize(4);

    cW = zeros(NO);
    for ii = 1:NO
        for jj = 1:NO
            if idx(ii)==idx(jj)
                cW(ii,jj) = 1;
            end
        end
    end

    W = cW.* ( (rand(NO)+(eye(NO))) <= ProbCluster) + W.* ( rand(NO) <= ProbOriginal) ;
    

	    out.stats.clusters = idx;


connectivity = W;
end


% ensure that there are no self connections
connectivity(logical(eye(M*N))) = 0;

% [=================================================================]
%  Training/Input
% [=================================================================]

if ~isvector(training_signal)
    training_signal = zeros(N*M,1);
end


% [=================================================================]
%  Scale coupling parameter?
% [=================================================================]

% Hu, X., Boccaletti, S., Huang, W., Zhang, X., Liu, Z., Guan, S., & Lai, C.-H. (2014). Exact solution for first-order synchronization transition in a generalized Kuramoto model. Scientific Reports, 4, 7262–6. http://doi.org/10.1038/srep07262

% Do not use scaling combined with plasticity? Or incorporate it
% in the plasticity function?
% if scale_to_intrinsic_freq
% 	connectivity = bsxfun(@times, omega_i, connectivity) * scaling / NO;
% else
% 	connectivity = scaling*connectivity;
% end


% [=================================================================]
%  create noise
% [=================================================================]

ou_noise = zeros(NO, simtime/dt);
if noise

	mixalpha = .3;
	sig = .05;
	th = 10;
	mu = 0;

	for t = 2:simtime/dt
		ou_noise(:,t) =  ou_noise(:,t-1) +th*(mu-ou_noise(:,t-1))*dt + ...
	        (1-mixalpha)*sig*sqrt(dt)*randn(NO,1) + ...
	         mixalpha*sig*sqrt(dt)*randn*ones(NO,1);
	end
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

k = zeros(1,simtime*(1/dt)); 
MP = zeros(simtime*(1/dt));

% if gpu
% 	theta_t = gpuArray(theta_t);
% 	connectivity = gpuArray(connectivity);
% 	ou_noise = gpuArray(ou_noise);
% 	k = gpuArray(k);
% end

% [=================================================================]
%  simulate
% [=================================================================]

sConnectivity = sigmoidConnectivity(connectivity, sigmoid(1), sigmoid(2));
adjacency{1} = sConnectivity;
spikeTime = zeros(N*M,1); requiresUpdate = zeros(N*M);

if plotme; f = figure(100); a(1) = subplot(121);a(2) = subplot(122); end
for t = 2:simtime/dt

	phasedifferences = bsxfun(@minus, theta_t(:,t-1)',theta_t(:,t-1));

	phasedifferences_W = sConnectivity.*scaling.*sin(phasedifferences);
	
	summed_sin_diffs = mean(phasedifferences_W,2); %ignore self?

    %%% Training requires testing
	theta_t(:,t) = theta_t(:,t-1) + dt*( omega_i + summed_sin_diffs + (t*dt<training_time)*training(1)*sin(theta_t(:,t-1)-training(2)*dt*t-training_signal(:,1)) ) + ou_noise(:,t);
    

	PP(:,t) = sin(mod(theta_t(:,t),2*pi));
    
    switch plasticity{1}
        case 'seliger' %{2} - epsilon, {3} - alpha
            % Currently ignores spatial distance between oscillators, should
            % it?
            connectivity = connectivity + dt*plasticity{2}* ...
                ( plasticity{3} * cos(phasedifferences) - connectivity );
            
        case 'STDP' % {2} - delta-adjustment: factor affecting weight change, {3} - [A1;A2], {4} - [tau1;tau2]
            %%% Does not yet account for sign(0)=0, also does not account for bsxfun result deltaTime=0 (as 0 is used to ignore elements)
            
            % Determines upward zero crossover (= spike) of each oscillator
            % and calculates time difference between all spikes
            upwardZeroCross = sign(mod(theta_t(:,t),2*pi)-pi) > sign(mod(theta_t(:,t-1),2*pi)-pi);
            spikeTime(upwardZeroCross) = t.*dt;
            
            deltaTime = bsxfun(@minus, spikeTime,spikeTime'); %deltaTime = t_i - t_j
            
            % Spiked oscillators are memorized and all updated oscillator
            % couples are determined
            requiresUpdate(:,upwardZeroCross) = 1;
            updateMatrix = requiresUpdate.*requiresUpdate';
            updateMatrix(logical(eye(N*M))) = 0;
            requiresUpdate = requiresUpdate - updateMatrix;
            
            % Actual implementation of Spike timing-dependent plasticity
            % function, see http://www.scholarpedia.org/article/Spike-timing_dependent_plasticity#Basic_STDP_Model
            deltaTime = deltaTime.*updateMatrix;
            dTimePos = deltaTime > 0;
            dTimeNeg = deltaTime < 0;
            
            W = dTimePos.*plasticity{3}(1).*exp(- deltaTime./plasticity{4}(1)) ...
                -dTimeNeg.*plasticity{3}(2).*exp( deltaTime./plasticity{4}(2));

            % Update connectivity
            connectivity = connectivity + plasticity{2}.*W;
            
            case 'test'
                connectivity = connectivity + dt*plasticity{2}*((abs(phasedifferences)<=plasticity{3}).*phasedifferences - connectivity);
            
        case 'null'
        otherwise
            disp('Running without a plasticity rule.')
            plasticity{1} = 'null';
    end
    
    connectivity = connectivity - decay.*dt;
    connectivity(connectivity<0)=0;
    sConnectivity = sigmoidConnectivity(connectivity,sigmoid(1),sigmoid(2));
    adjacency{t} = sConnectivity;
	% [=================================================================]
	%  order parameter
	% [=================================================================]
	if ~clusterize(1)
		GP = theta_t(:,t);
		MP = circ_mean(GP)+pi;
		
		k(t) = mean( exp(1i*(bsxfun(@minus, GP, MP))));
	else
		for ui = unique(idx)'
			GP = theta_t(find(idx==ui),t);
			MP = circ_mean(GP)+pi;
			
			k(ui,t) = mean( exp(1i*(bsxfun(@minus, GP, MP))));
		end
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
	

	if exist('linspecer')
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


% if gpu
% 	theta_t = gather(theta_t);
% 	connectivity = gather(connectivity);
% 	ou_noise = gather(ou_noise);
% 	k = gather(k);
% end

out.state = PP;
out.phase = theta_t;
out.parameters = p.Results;
out.oscillators = omega_i/(2*pi);
out.orderparameter = abs(k);
out.connectivity = connectivity;
out.adjacency = adjacency;
out.meanphase = MP;
out.seed = seed;
 if makemovie && plotme ; out.movie = MOV; end
% out.all =  sin(mod(theta_t,2*pi));
end

function W = sigmoidConnectivity(connectivity, a, c)
    if a==0
        W=connectivity;
    else
        W = sigmf(connectivity, [a c]);
    end
    
%    W(connectivity<0)=0; % lower bound (actual) connectivity before this function
%    instead.
end