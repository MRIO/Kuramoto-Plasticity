%% Input
clear

plotResult = 0;
alpha = 2.5;
epsilon = 1;
connectivity = 2;
oscillators = [11;10]; % in hz

%% Model
result = kuramotoCouple(1,plotResult,'time', 2, 'connectivity', (ones(2)-eye(2))*connectivity,'plasticity',[1 epsilon alpha], 'oscillators', oscillators*2*pi)

% %fuck it
% size=size(result.connectivity,2);
% K = zeros(1,size);
% for i=1:size
%     K(i)=result.connectivity{i}(2);
% end

K = result.Kval;
size = size(K,2);

%phi = result.phase(1,:) - result.phase(2,:);
phi = result.Dphase;
DOmega = result.oscillators(1)-result.oscillators(2);
%DOmega = result.Domega;

K_tilde = K/DOmega;
alpha_tilde = alpha/DOmega;
epsilon_tilde = epsilon/DOmega;

%% Plot
f = figure(1)
clf


t=linspace(-2*pi,4*pi,1000);

y=1./sin(t); %phi nullcline
y2=alpha_tilde*cos(t); %K_tilde nullcline

hold on
legend('Sim result','phi nullcline','K nullcline') %fix
xlabel('phi')
ylabel('K tilde')
axis([0 2*pi -3 3])


for i=2:size
    if ~ishghandle(f)
        hold off
        close
        return
    end
    
    cla
    
    plot(phi(1:i),K_tilde(1:i))
    %scatter(mod(phi(i),2*pi),K_tilde(i))
    plot(t,y) %phi null
    plot(t,y2) %K_tilde null
    
    drawnow
end

hold off

