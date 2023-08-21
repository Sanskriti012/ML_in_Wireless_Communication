%% Question 3
%% In this question we generate t for given values of parameter and for different values of noice varience.
clc;clear;close all;
rng(1);
% specify parameters
N = 20; % number of observations
M = 40; % number of features
D0 = 7; % number of non-zero entities in weight vector
w = sprand(40,1,0,0.175);
vr_dB = [-20, -15, -10, -5, 0]; % noise variance in dB (vr is varience)

% our t is in form as, t = Phi*w + n.

% now we initialize our variables to store results
Phi_all = cell(length(vr_dB), 1);
w_all = cell(length(vr_dB), 1);
n_all = cell(length(vr_dB), 1);
t_all = cell(length(vr_dB), 1);

% generate data for each noise variance that are given to us.
for i = 1:length(vr_dB)
    sigma = 10^(vr_dB(i)/10); % convert dB to linear scale
    % generate design/dictionary matrix
    Phi = randn(N, M);
    Phi_all{i} = Phi;
    
    % generate sparse weight vector
    w = zeros(M, 1);
    idx = randperm(M, D0); % randomly select D0 indices
    w(idx) = randn(D0, 1); % set the selected indices to be random Gaussian values
    w_all{i} = w;
    
    % generate noise and observations
    n = sigma*randn(N, 1); % generate noise
    n_all{i} = n;
    t = Phi*w + n; % generate observations
    t_all{i} = t;

    figure(i);
    hold off
    plot(Phi,t,'k.','markersize',10);
    xlabel('phi');
    ylabel('t');
end
