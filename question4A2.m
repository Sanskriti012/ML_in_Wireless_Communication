%% Question4
%% In this question we have to apply SBL for regression from [R1] to get the maximum apoaterior estimate of the weight vector w.

% clc; clear all; close all;
rng(1)
% Set hyperparameters
alpha = 1e-6; % Prior strength
beta = 1e-6; % Noise precision

% Initialize variables
N = 20;
M = 40;
D0 = 7;
vr_dB = [-20, -15, -10, -5, 0];
alpha = ones(M,1)*100;

% initialize some matrix
NSME=zeros(5,1);
w = zeros(M,1); 
w_hats = ones(M,5);
gamma=zeros(M,1);

%generate design/dictionary matrix
X = randn(N, M);

alpha_new=alpha;

for j = 1:length(vr_dB)
    sigma = 10^(vr_dB(j)/10); % convert dB to linear scale
    % generate design/dictionary matrix
    Phi = randn(N, M);
    Phi_all{j} = Phi;
    
    % generate sparse weight vector
    w = zeros(M, 1);
    idx = randperm(M, D0); % randomly select D0 indices
    w(idx) = randn(D0, 1); % set the selected indices to be random Gaussian values
    w_all{j} = w;
    
    % generate noise and observations
    n = sigma*randn(N, 1); % generate noise
    n_all{j} = n;
    t = Phi*w + n; % generate observations
    t_all{j} = t;

% Run SBL algorithm
max_iterations = 100;
w_hat_prev= zeros(M,1);

for i = 1:max_iterations
    % Compute posterior precision matrix
    S = diag(alpha) + beta*(X'*X);
    
    % Compute posterior mean
    w_hat= beta*(S\X')*t;
    
    % Compute posterior covariance
    Sigma = inv(S);
  % computing alpha_new parameter. 
  for p=1:M
      gamma(p) = 1- alpha(p)*Sigma(p,p)
      alpha_new(p) = gamma(p)/(w_hat(p))^2;
  end 
    
    % Check convergence
    if norm(w_hat-w_hat_prev)^2 < (1e-3)* norm(w_hat_prev)^2
        break;
    end
    w_hat_prev = w_hat;
    alpha = alpha_new;
end
% creating w_hats matrix
for k = 1:M
    w_hats(k,j)= w_hat(k);
end
end