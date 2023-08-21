%% Predictive variance example
clc; clear; close all;
rng(1);
%% Generate a dataset from the given model
N=100; % Number of data points
% Generate uniformly distributed values of x between -5 and 5
x=10*(sort(rand(N,1))-0.5);
% Define  true parameters
w0 = 0; w1 = 1; w2 = -1; w3 = 5;
% Define t
t = w0 + w1*x + w2*(x.^2) + w3*(x.^3);
% Add some noise to our t which is normally distributed
noise_var = 300;
t = t + sqrt(noise_var)*randn(N,1);

%% plot the data
figure(1);
hold off
plot(x,t,'k.', 'markersize',10);
xlabel('x');
ylabel('t');
hold on
title('Data');

%% Predictive error bar for linear(first order) model
testx =[-5:0.1:5]';
order1 = 1;
i = 1:length(order1);
X = [];
testX = [];
for k = 0:order1(i);
        X = [X x.^k];
        testX = [testX testx.^k];
end 
% write the value of w_hat as w and sigma square
w = inv(X'*X)*X'*t;
sigmasquare = (1/N)*(t'*t - t'*X*w);
%calculate the mean and variance 
lin_mean = testX*w;
lin_var = sigmasquare * diag(testX*inv(X'*X)*testX');
% Plot the data and predictions for linear(first order) model
figure(2);
hold off
plot(x,t,'k.','markersize',10);
xlabel('x');
ylabel('t');
hold on
% plot the error bar for linear model
errorbar(testx,lin_mean,lin_var,'r')
title('linear model');

%% Predictive error bar for cubic (third order) model
testx =[-5:0.1:5]';
order2 = 3;
i = 1:length(order2);
X = [];
testX = [];
for k = 0:order2(i);
        X = [X x.^k];
        testX = [testX testx.^k];
end 
% write the value of w_hat as w and sigma square
w = inv(X'*X)*X'*t;
sigmasquare = (1/N)*(t'*t - t'*X*w);
%calculate the mean and variance 
cub_mean = testX*w;
cub_var = sigmasquare * diag(testX*inv(X'*X)*testX');
% Plot the data and predictions for cubic(third oder) model
figure(3);
hold off
plot(x,t,'k.','markersize',10);
xlabel('x');
ylabel('t');
hold on
% plot the error bar for cubic model
errorbar(testx,cub_mean,cub_var,'g')
title('cubic model');

%% Predictive error bar for sixth order model
testx =[-5:0.1:5]';
order3 = 6;
i = 1:length(order3);
X = [];
testX = [];
for k = 0:order3(i);
        X = [X x.^k];
        testX = [testX testx.^k];
end 
% write the value of w_hat as w and sigma square
w = inv(X'*X)*X'*t;
sigmasquare = (1/N)*(t'*t - t'*X*w);
%calculate the mean and variance 
six_mean = testX*w;
six_var = sigmasquare * diag(testX*inv(X'*X)*testX');
% Plot the data and predictions for sixth oder model
figure(4);
hold off
plot(x,t,'k.','markersize',10);
xlabel('x');
ylabel('t');
hold on
% plot the error bar for sixth order model
errorbar(testx,six_mean,six_var,'b')
title('sixth oder model');

%% combined plot of all models
figure(5);
hold off
plot(x,t,'k.','markersize',10);
xlabel('x');
ylabel('t');
hold on
errorbar(testx,lin_mean,lin_var,'r')
errorbar(testx,cub_mean,cub_var,'g')
errorbar(testx,six_mean,six_var,'b')
legend('data','linear','cubic','sixth-oder');
title('combined plot');