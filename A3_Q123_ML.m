clc;clear all;close all;
rng(1);
%% Q1- Generation of synthetic data set:
% Set the parameters
mu1 = [3; 3]; % mean
sigma1 = [1 0; 0 2]; % covarience
mu2 = [1; -3]; % mean
sigma2 = [2 0; 0 1]; % covarience

% Mixing probabilities for the two Gaussians
pi1 = 0.8;
pi2 = 0.2;

N = 500;% Number of data points

% Generate data points from the gaussian
X = [mvnrnd(mu1', sigma1, round(pi1*N)); mvnrnd(mu2', sigma2, round(pi2*N))]';

% Plot the generated dataset
figure(1);
plot(X(1,:),X(2,:),'k.','markersize',10) %plot of synthetic data in cluster 1
xlabel('x1');
ylabel('x2');
title('Synthetic dataset');
hold on;

%% Q2 and 3 
% Now we apply Expectation-Maximisation (EM) algorithm to the above dataset
% Iteratively perform the Step 1 and Step 2, till the lower bound converges.

% Number of clusters
K = 2;

% Initialize the parameters
pi = [0.5, 0.5];
mu = [1, -1; 1, 1];
sigma(:,:,1) = [1 0; 0 1];
sigma(:,:,2) = [1 0; 0 1];

% Define the tolerance for convergence
tol = 1e-6;

% Initialize the log-likelihood
log_likelihood = -inf;

while true
    % Expectation step
    
    % Compute the responsibilities
    for k = 1:K
        numerator(:,k) = pi(k) * mvnpdf(X', mu(:,k)', sigma(:,:,k));
    end
    denominator = sum(numerator, 2);
    q = numerator ./ denominator;
    
    % Compute the log-likelihood 
    prev_log_likelihood = log_likelihood;
    log_likelihood = sum(log(denominator));
    
    % Check for convergence
    if abs(log_likelihood - prev_log_likelihood) < tol
        break;
    end
    
    % Maximization step
    
    % Update the mixing probabilities
    pi = sum(q, 1) / size(X, 2);
    
    % Update the means
    for k = 1:K
        mu(:,k) = (X * q(:,k)) / sum(q(:,k));
    end
    
    % Update the covariance matrices
    for k = 1:K
        X_centered = X - mu(:,k);
        sigma(:,:,k) = (X_centered * diag(q(:,k)) * X_centered') / sum(q(:,k));
    end
end

% Assign data points to clusters
[~, C] = max(q, [], 2);

% Display results
fprintf('Mixing Coefficients:\n');
disp(pi);
fprintf('mean 1:\n');
disp(mu(:,1));
fprintf('Covariance: 1 \n');
disp(sigma(:, :, 1));
fprintf(' Mean: 2 \n');
disp(mu(:,2));
fprintf('Covariance 2 :\n');
disp(sigma(:, :, 2));


% Create grid for contour plot
x = linspace(-3, 8, 100);  
y = linspace(-6, 9, 100);  
[X0, Y0] = meshgrid(x, y); % create meshgrid for x and y

% Calculate Gaussian distribution
F = mvnpdf([X0(:), Y0(:)], mu(:,1)', sigma(:, :, 1)); 
G = mvnpdf([X0(:), Y0(:)], mu(:,2)', sigma(:, :, 2)); 

% Reshape data for contour plot
 U= reshape(F, length(y), length(x)); 
V = reshape(G, length(y), length(x)); 

% Plot the Contour map and Data points
figure(2);
scatter(X(1,:), X(2,:),10,C,'filled' );
hold on;
contour(X0, Y0, U); 
hold on
contour(X0, Y0, V);
hold on
xlabel('x1');
ylabel('x2');
legend({'cluster 1','cluster 2'})
title('EM Algorithm for Gaussian Mixture Model');