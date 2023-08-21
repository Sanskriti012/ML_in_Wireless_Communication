%% pdf of two dimensional gaussian random vector
%% gaussian_surface
%% Define the parameters for (a) part
mu = [2; 1]';
sigma = [1, 0; 0, 1];

% Create a grid of points for evaluating the PDF
x1 = -2:0.1:5;
x2 = -2:0.1:5;
[x1, x2] = meshgrid(x1,x2);
x = [x1(:)-mu(1) x2(:)-mu(2)];

% Evaluate the PDF
t=(1/sqrt(2*pi));
t=t./sqrt(det(sigma));
pdf = t*exp(-0.5*diag(x*inv(sigma)*x'));

% Reshape the PDF into a grid
pdf = reshape(pdf, length(x2),length(x1));

% Plot the PDF(Probability density)
figure;
surf(x1, x2, pdf);
xlabel('x_1');
ylabel('x_2');
zlabel('PDF');
% Giving title to the plot
title('PDF of two-dimensional Gaussian (a)');

% Plot the contour plot
figure;
contour(x1, x2, pdf);
xlabel('x_1');
ylabel('x_2');
title('Contour plot of two-dimensional Gaussian (a)');

%% Define the parameters for (b) part
mu = [2; 1]';
sigma = [1, 0.8; 0.8, 1];

% Create a grid of points for evaluating the PDF
[x1, x2] = meshgrid(-2:0.1:5, -2:0.1:5);
x = [x1(:)-mu(1) x2(:)-mu(2)];

% Evaluate the PDF
t=(1/sqrt(2*pi));
t=t./sqrt(det(sigma));
pdf = t*exp(-0.5*diag(x*inv(sigma)*x'));

% Reshape the PDF into a grid
pdf = reshape(pdf, size(x1));

% Plot the PDF
figure;
surf(x1, x2, pdf);
caxis([min(pdf(:))-.5*range(pdf(:)), max(pdf(:))]);
xlabel('x_1');
ylabel('x_2');
zlabel('PDF');
% Giving title to the plot
title('PDF of two-dimensional Gaussian (b)');

% Plot the contour plot
figure;
contour(x1, x2, pdf);
xlabel('x_1');
ylabel('x_2');
title('Contour plot of two-dimensional Gaussian (b)');
