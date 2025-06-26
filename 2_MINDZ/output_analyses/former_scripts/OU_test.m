% OU returns a value following an Ornstein-Uhlenbeck stochastic process.
%
% http://planetmath.org/encyclopedia/OrnsteinUhlenbeckProcess.html
%
% MapsF 2011


% Test the independence of the O-U process to the value of the timestep.


%--- Initialization

damp = 5e-4;          % Rate of mean reversion

sigma = 0.1;          % Volatility or average magnitude of the random fluctuations

mu = 5;              % Long-term mean of the process

dt = [1800 14400];    % Time-steps

tmax = 5*86400/dt(1); % # of timesteps

var = nan(tmax,2);    % Variable; here ~ depth

var(1,:) = 100;

bottom = 300;


%--- Computation

for t = 2:tmax

    M  = mu+(var(t-1,:)-mu).*exp(-damp.*dt);

    SD = sqrt( (0.5*sigma/damp)*(1-exp(-2*damp.*dt)) );

    dummy(t,:) = M + SD.*randn(1,2);

    % correct for negative values = above surface!
    var(t,:) = max(0.001, bottom-abs(bottom-abs(dummy(t,:))) );

end


%--- Plot

figure


subplot(2,1,1)

[n,xn] = hist(var(:,1),300);

bar(xn,n,'edgecolor','k','facecolor','k')

axis([0 100 0 30])

hold on

[n,xn] = hist(var(1:end/dt(2)*dt(1),2),10);

bar(xn,n,'edgecolor','r','facecolor','none','linewidth',2)

set(gca,'fontsize',16)

title('Depth frequency distribution','fontsize',18)

xlabel('Depth')

ylabel('Count')


subplot(2,1,2)

plot([1:tmax]*dt(1)/86400,var(:,1),'k.')

axis([0 0.1*tmax*dt(1)/86400 0 150])

hold on

plot([1:dt(2)/dt(1):tmax]*dt(1)/86400,var(1:end/dt(2)*dt(1),2),'rs')

set(gca,'fontsize',16)

title('Depth time evolution','fontsize',18)

%xlabel('Time-step')

datetick('x','HH:MM','keeplimits')

ylabel('Depth')

legend('dt=10min','dt=4h')

legend('boxoff')

%print('-depsc2','-tiff','-painters','Fig_OU_test01.eps')


figure

[h,stats] = cdfplot(var(:,1))

set(h,'color','k')

hold on

[h,stats] = cdfplot(var(:,2))

set(h,'color','r')

xlabel('Depth (m)')

legend('dt=10min','dt=4h')

legend('boxoff')

%print('-depsc2','-tiff','-painters','Fig_OU_test02.eps')
