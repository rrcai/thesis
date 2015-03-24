% Calculates mean square error and mean absolute error as a function of the
% number of data points available in prior belief of cross-covariance
% matrix.

clear all;
close all;

n_x = 2;
n_y = 1;
sigma_x = 1;
sigma_y = 1;
sigma_X = sigma_x*eye(n_x);
sigma_Y = sigma_y*eye(n_y);
sigma_XY = [0.1; 0.2];

% Choose beta.
beta = 100;
% Choose m (number of times to sample the Wishart distribution).
m = 100;
% Choose variety of n values (number of data points available).
allN = 5*(100:1000);

% Set up array to store values.
MSE_data = zeros(size(allN,2),3);
MSE_data(:,1) = allN;
MSE_errdata = zeros(size(allN,2),3);
MSE_errdata(:,1) = allN;

[TrueA, Truebeta_crit] = gib_optimize(sigma_X,sigma_Y,sigma_XY,beta);
trueVal = Truebeta_crit(1,1);

for n = 1:size(allN,2);
    Data = sample_wishart(sigma_X, sigma_Y, sigma_XY, beta, allN(n), m);
    
    % Plot different values of W, A, and beta_crit
    allW = cell2mat(Data(:,1));
    allA = cell2mat(Data(:,2));
    allBeta_crit = cell2mat(Data(:,3));
    
    indices = (((1:(size(allBeta_crit,1)/2)).*2)-1)';
    thisBetacritData = allBeta_crit(indices,1);
    
    this_abs = abs(thisBetacritData-trueVal)./abs(trueVal);
    assert(max(this_abs)<5);
    MSE_data(n,2) = sum(this_abs)/m;
    MSE_errdata(n,2) = std(this_abs)/sqrt(m);
    
    this_abs = (thisBetacritData-trueVal).^2;
    MSE_data(n,3) = sum(this_abs)/m;
    MSE_errdata(n,3) = std(this_abs)/sqrt(m);
end

% Plot this data.

figure;
errorbar(MSE_data(:,1),MSE_data(:,2),MSE_errdata(:,2));
axis([0 1.05*allN(end) 0 1.05*max(MSE_data(:,2))]);
hold on;
% plot([allN(1) allN(end)],[TrueA(1,k) TrueA(1,k)]);

str = sprintf('Error of estimates of {\\beta_1}^c with %d samples',m);
title(str);

str = sprintf('percent error of estimated value of {\\beta_1}^c');
ylabel(str);
xlabel('number of initial data points given');
str = sprintf('beta_perc_%dsamples.png',m);
% print('-dpng', str);
hold off;