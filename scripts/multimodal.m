clear
clf
clc
close

sigma = 0.5;
mu = [4 7];
range = [0 12 100000]; % min; max; steps
data = rand(1,100000);

%% RNG
% add multiple gaussians and normalize
bi = @(x) (normpdf(x,mu(1),sigma) + normpdf(x,mu(2),sigma))./2;


integral(bi,range(1),range(2)) % check normalization
t = linspace(range(1),range(2),range(3));
y = bi(t); % pdf
y2 = cumsum(y).*((range(2)-range(1))./range(3)); % cdf, riemann sum
y2(end) % check cdf, should be 1

idx =(y2<0.0001|y2>0.9999);
y2(idx) = [];
t2 = t;
t2(idx) = [];

Vq = interp1(y2,t2,data,'linear'); %To pull values from (inverted) cdf

%% Plot
figure('Position', [75 500 1750 400])
subplot(1,4,1)
plot(t,y)
title('PDF')
subplot(1,4,2)
hold on
plot(t2,y2, 'blue')
%plot(Vq,data,'x')
hold off
title('CDF')
subplot(1,4,3)
hold on
plot(y2,t2,'blue')
%plot(data,Vq,'xblue')
hold off
title('CDF inverted')
subplot(1,4,4)
histogram(Vq,100,'Normalization','probability')
title('Random sample')