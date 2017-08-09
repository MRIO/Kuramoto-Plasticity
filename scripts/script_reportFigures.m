%% 0.1.1 Replay data
clear idx XX YY SS MOV PP
N=10;M=10;
PP = out1.state(:,4551:4552);
%PP = out1.state;

f = figure(101);
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

%colormap(f,[[1 0 0];[0 0 1]])
%scatter3(XX, YY ,PP(:,t),60,labs+2,'filled')
% use colormap/labs to mark clusters (retrieve labs from spect_cluster
% output.c.labels(:,#of_clusters)

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

%% figstab
figure
subplot(2,2,1)
plot(out1.state(:,1:600)');
xlabel('Time (ms)')
title('STDP')
subplot(2,2,2)
plot(out1.state(:,4401:5000)');
set(gca,'XTick',0:200:600)
set(gca,'XTickLabel',4400:200:5000)
xlabel('Time (ms)')
title('STDP')
subplot(2,2,3)
plot(out2.state(:,1:600)');
xlabel('Time (ms)')
title('No Plasticity')
subplot(2,2,4)
plot(out2.state(:,4401:5000)');
set(gca,'XTick',0:200:600)
set(gca,'XTickLabel',4400:200:5000)
xlabel('Time (ms)')
title('No Plasticity')

%% figstab 2
close all
figure
subplot(1,2,1)
N=10;M=10;
t=4550;
PP = out1.state; %#1

SS = reshape(PP(:,t),N,M);
    mesh(SS); hold on
    scatter3(XX, YY ,PP(:,t),60,'filled')
    
    title('STDP (t = 4550ms)')
    zlabel('Phase')
    caxis([0 2*pi])
    zlim([-3 3])

subplot(1,2,2)
PP = out2.state; %#2

SS = reshape(PP(:,t),N,M);
    mesh(SS); hold on
    scatter3(XX, YY ,PP(:,t),60,'filled')
    
    title('No Plasticity (t = 4550ms)')
    zlabel('Phase')
    caxis([0 2*pi])
    zlim([-3 3])
    
figure
N=10;M=10;
t=75;
PP = out1.state; %#1

SS = reshape(PP(:,t),N,M);
    mesh(SS); hold on
    scatter3(XX, YY ,PP(:,t),60,'filled')
    
    title('Training (t = 75ms)')
    zlabel('Phase')
    caxis([0 2*pi])
    zlim([-3 3])

%% STDP Function plots
% Sigmoid
clear
close all
figure
a = 10;
c = 0.5;
t = linspace(-1,2,1000);

y = (1./(1+exp(-a.*(t-c)))) .* (t>0.501);



plot(t,y)
xlabel('Original connectivity')
ylabel('Adjusted connectivity')
title('Sigmoid Connectivity Adjustment')

% STDP function
figure
subplot(1,3,1)
dTimePos = 1;
dTimeNeg = 1;
A = [1 1];
tau = [0.002 0.002];
t = linspace(-0.02,0.02,1000);
ts = t>0;

y = ts.*A(1).*exp(- t./tau(1)) -~ts.*A(2).*exp( t./tau(2));

plot(t,y)
xlabel('Time between spikes (s)')
ylabel('Weight change')
title('STDP Function (Tau = 0.002)')

subplot(1,3,2)
dTimePos = 1;
dTimeNeg = 1;
A = [1 1];
tau = [0.0005 0.0005];
t = linspace(-0.02,0.02,1000);
ts = t>0;

y = ts.*A(1).*exp(- t./tau(1)) -~ts.*A(2).*exp( t./tau(2));

plot(t,y)
xlabel('Time between spikes (s)')
ylabel('Weight change')
title('STDP Function (Tau = 0.0005)')

subplot(1,3,3)
dTimePos = 1;
dTimeNeg = 1;
A = [1 1];
tau = [0.008 0.008];
t = linspace(-0.02,0.02,1000);
ts = t>0;

y = ts.*A(1).*exp(- t./tau(1)) -~ts.*A(2).*exp( t./tau(2));

plot(t,y)
xlabel('Time between spikes (s)')
ylabel('Weight change')
title('STDP Function (Tau = 0.008)')

figure
t=linspace(0,20,10000);
y=mod(t,2*pi)-pi;
hold on
plot(t,y,'b')
plot(t,zeros(1,length(t)),'black')
plot([pi 3*pi 5*pi],[0 0 0],'red.','MarkerSize',25)
hold off
xlabel('Time')
ylabel('Phase')
title('Spike Definition')

%%
close all
clear A
figure

for ii=1:3
subplot(1,3,ii)
A{1}(:) = [1 1];
A{2}(:) = [1 -1];
A{3}(:) = [1 0];
tau = [0.002 0.002];
t = linspace(-0.02,0.02,1000);
ts = t>0;

y = ts.*A{ii}(1).*exp(- t./tau(1)) -~ts.*A{ii}(2).*exp( t./tau(2));
hold on

plot(t,zeros(1,length(t)),'black')
plot(t,y, 'b')
plot([0 0],[-1 1], 'black')

hold off
axis([-0.02 0.02 -1 1])
%xlabel('Time between spikes (s)')
%ylabel('Weight change')

end