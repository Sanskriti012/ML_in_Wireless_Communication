%% Question 2
%% In this question we  have to generate systhetic data 
clc;clear;close all;
rng(1);
%% first part where we have to generate NXM matrix.
% specify parameters.
N = 100; % number of rows
M = 50; % number of columns
%generate design/dictionary matrix
Phi = randn(N, M);

%% second part where we have to generte MX1 sparse weight vector
% specify parameter
D0 = 10; % desired number of non-zero entries
% generate sparse weight vector
w = zeros(M, 1);
idx = randperm(M, D0); % randomly select D0 indices
w(idx) = randn(D0, 1); % set the selected indices to be random Gaussian values

%% third part where we have to generate noice entries and observations t.
% specify parameter
sigma = 0.5; % variance of the noice
% generate noise and observations of the data.
n = sigma*randn(N, 1); % generate noise
t = Phi*w + n; % generate observations

display(t)