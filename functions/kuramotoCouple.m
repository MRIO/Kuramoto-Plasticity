function out = kuramotoCouple(varargin)
    % Reduced variant of kuramotoSheet.m; a simple model of two coupled
    % oscillators
    %
    % Inputs:
    %   scaling for coupling - K
    % Outputs:

    anim = 1;
    idx = [1 1];
% [=================================================================]
%  parse inputs
% [=================================================================]

    p = inputParser;
    p.addRequired('scaling');
    p.addRequired('plotResults');

    p.addParameter('dt', 1e-3)
    p.addParameter('time',1) % in seconds
    p.addParameter('omega_mean', 7) % in Hz
    p.addParameter('omega_std', 2) % in Hz
    p.addParameter('init_cond', []) % in Hz
    p.addParameter('oscillators', [])
    p.addParameter('connectivity', []);
    p.addParameter('plasticity', [0 1 1]); % TODO, test and change default values

    p.parse(varargin{:});

    scaling = p.Results.scaling;
    plotResults = p.Results.plotResults;
    dt = p.Results.dt;
    simtime = p.Results.time;
    omega_mean = p.Results.omega_mean;
    omega_std = p.Results.omega_std;
    init_cond = p.Results.init_cond;
    oscillators = p.Results.oscillators;
    connectivity{1} = p.Results.connectivity;
    plasticity = p.Results.plasticity;

% [=================================================================]
%  randomize oscillator intrinsic frequencies
% [=================================================================]

    rng(0,'twister')
    if oscillators
        omega_i = oscillators;
    else
        omega_i = (randn(2,1)*omega_std+omega_mean)*2*pi;
    end

% [=================================================================]
%  connectivity
% [=================================================================]

    if isempty(connectivity{1})
        connectivity{1} = (ones(2) - eye(2))*scaling;
    else
        connectivity{1} = scaling*connectivity{1};
    end

% [=================================================================]
%  randomize initial condition
% [=================================================================]
    %initial condition (initial phase)
    if isempty(init_cond)
        theta_t(:,1) = rand(2,1)*2*pi;
    else
        % disp('using initial condition (assuming between 0 and 2pi)')
        theta_t(:,1) = init_cond;
    end

    k = zeros(1,simtime*(1/dt));
    PP = zeros(2,simtime/dt);

% [=================================================================]
%  simulate
% [=================================================================]


    for t = 2:simtime/dt

        phasedifferences = bsxfun(@minus, theta_t(:,t-1)',theta_t(:,t-1));

        phasedifferences_W = connectivity{t-1}.*sin(phasedifferences);

        summed_sin_diffs = mean(phasedifferences_W,2);

        theta_t(:,t) = theta_t(:,t-1) + dt*( omega_i + summed_sin_diffs  );

        PP(:,t) = sin(mod(theta_t(:,t),2*pi));
        
        if plasticity(1)
             connectivity{t} = connectivity{t-1} + dt*plasticity(2)* ...
            ( plasticity(3) * cos(phasedifferences) - connectivity{t-1} );
        else
            connectivity{t} = connectivity{t-1};
        end

        % [=================================================================]
        %  order parameter
        % [=================================================================]
        GP = theta_t(:,t);
        MP = circ_mean(GP)+pi;

        k(t) = mean( exp(1i*(bsxfun(@minus, GP, MP))));
    end

    writeOutput()

% [=================================================================]
%  plots
% [=================================================================]
        if ~plotResults
            return;
        end
        
        f = figure(100);
        a(1) = subplot(121);
        a(2) = subplot(122);

        ffff = figure;

        subplot(2,2,1)
        plot(linspace(0,simtime, simtime*dt^-1), mod(theta_t,2*pi)')
        ylabel('phase (theta)')
        subplot(2,2,2)
        imagesc(connectivity{1}), colorbar %% TODO, plot connectivity over time? movie?
        title('Initial connectivity')
        subplot(2,2,3)
        plot(linspace(0,simtime, simtime*dt^-1), sin(theta_t'))
        ylabel('sin(theta)')
        xlabel('seconds')
        subplot(2,2,4)
        hist(omega_i/(2*pi),2)
        ylabel('#')
        xlabel('frequency (Hz)')

    figure(f)
        axes(a(2))
            line(repmat(linspace(0,simtime, simtime*dt^-1), length(unique(idx)), 1)', abs(k)')
            title('kuramoto parameter')
            xlabel('time (ms)')
            ylim([0 1])

        axes(a(1));
            while anim

                for t = 2:simtime/dt
                    if ~ishghandle(a(1))
                        return
                    end
                    cla

                    scatter([-1 1], PP(:,t), 60, 'filled');
                    title('phase')
                    axis([-2 2 -2 2])

                    drawnow
                end
                anim = input('repeat? ; 0 - no, 1 - yes \n');
            end

    function writeOutput()
        out.state = PP;
        out.phase = theta_t;
        out.parameters = p.Results;
        out.oscillators = omega_i/(2*pi);
        out.orderparameter = abs(k);
        out.meanphase = MP;
        out.connectivity = connectivity;
        out.timesteps = simtime/dt;
    end
end
