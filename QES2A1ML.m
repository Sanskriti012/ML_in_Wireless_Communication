%% Regularized least squares:
clc;clear all;close all;
rng(1);
%% Generate a synthetic dataset from a linear function
N = 6; %Number of data points

% Generate random x values between 0 and 1
%x = 1*(sort(rand(N,1))+0);
x = [0:0.2:1]';
%% Define the true parameters
w_0 = -3;
w_1 = 2;

%% Define t
y = w_0 + w_1*x;

%% add noice to our t
noicevar = 3;
t = y +sqrt(noicevar)*randn(size(x));

%% plot the data
plot(x,t,'k.','markersize',25);

%% Build up the data so that it has up to fifth order terms
testx = [0:0.01:1]';
X = [];
testX = [];
for k = 0:5
    X = [X x.^k];
    testX = [testX testx.^k];
end

%% Fit the model with different values of the regularization parameter
% $$lambda$$
lam = [0 1e-6 1e-2 1e-1];

%% for first value of lambda = 0
    lambda = lam(1);
    w=size(x,1);
    w=(X'*X + w*lambda*eye(size(X,2)))\(X'*t);
    figure(1);hold on
    plot(x,t,'k.','markersize',20);
    hold on
    plot(testx,testX*w,'r','linewidth',2)
    xlim([-0.1 1.1])
    xlabel('$x$','interpreter','latex','fontsize',15);
    ylabel('$f(x)$','interpreter','latex','fontsize',15);
    title('combine plot','Interpreter','latex','fontsize',20)
    
    %% for second value of lambda = 1e-6
    i=1:length(lam)
    lambda = lam(2);
    w=size(x,1);
    w=(X'*X + w*lambda*eye(size(X,2)))\(X'*t);
    figure(1);hold on
    plot(x,t,'k.','markersize',20);
    hold on
    plot(testx,testX*w,'b','linewidth',2)
    xlim([-0.1 1.1])
    xlabel('$x$','interpreter','latex','fontsize',15);
    ylabel('$f(x)$','interpreter','latex','fontsize',15);
   
    %% for third value of lambda = 1e-2
    i=1:length(lam)
    lambda = lam(3);
    w=size(x,1);
    w=(X'*X + w*lambda*eye(size(X,2)))\(X'*t);
    figure(1);hold on
    plot(x,t,'k.','markersize',20);
    hold on
    plot(testx,testX*w,'g','linewidth',2)
    xlim([-0.1 1.1])
    xlabel('$x$','interpreter','latex','fontsize',15);
    ylabel('$f(x)$','interpreter','latex','fontsize',15);

    %% for forth value of lambda = 1e-1
    i=1:length(lam)
    lambda = lam(4);
    w=size(x,1);
    w=(X'*X + w*lambda*eye(size(X,2)))\(X'*t);
    figure(1);hold on
    plot(x,t,'k.','markersize',20);
    hold on
    plot(testx,testX*w,'y','linewidth',2)
    xlim([-0.1 1.1])
    xlabel('$x$','interpreter','latex','fontsize',15);
    ylabel('$f(x)$','interpreter','latex','fontsize',15);

 legend('','data','lambda=0','','lambda=1e-6','','lambda=0.01','','lambda=0.1','')