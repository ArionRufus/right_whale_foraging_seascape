function [z] = OU_test_module(z,sigma,mu,dt)
% OU returns a value following an Ornstein-Uhlenbeck stochastic process.
%
% http://planetmath.org/encyclopedia/OrnsteinUhlenbeckProcess.html
%
% MapsF 2011

%--- Initialization

% z = the variable, in our case depth
% sigma = volatility or average magnitude of the random fluctuations
% mu = long-term mean of the process or in our case the targetted depth
% dt = time-step in s

damp = 5e-4;          % Rate of mean reversion = parameterize the "swimming speed"

%--- Computation

M  = mu+(z-mu).*exp(-damp.*dt);

SD = sigma; %sqrt( (0.5.*sigma.^2./damp)*(1-exp(-2.*damp.*dt)) );

z = M + SD.*randn(size(SD));

end