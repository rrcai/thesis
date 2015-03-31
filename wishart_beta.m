% Calculates mean square error and mean absolute error as a function of the
% number of data points available in prior belief of cross-covariance
% matrix.

clear all;
close all;

tic;

n_x = 2;
n_y = 1;
sigma_x = 1;
sigma_y = 1;
sigma_X = sigma_x*eye(n_x);
sigma_Y = sigma_y*eye(n_y);
sigma_XY = [0.1; 0.2];

% sigma_X = [1 0.5; 0.5 2];

% Choose beta.
beta = 100;
% Choose m (number of times to sample the Wishart distribution).
m = 500;
% Choose variety of n values (number of data points available).
allN = 10*(50:500);

% Set up array to store values.
MSE_data = zeros(size(allN,2),4);
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
    MSE_data(n,2) = sum(this_abs)/m;
    MSE_errdata(n,2) = std(this_abs)/sqrt(m);
    
    this_mse = (thisBetacritData-trueVal).^2;
    MSE_data(n,3) = sum(this_mse)/m;
    MSE_errdata(n,3) = std(this_mse)/sqrt(m);
    
    overage = logical(thisBetacritData > trueVal);
    MSE_data(n,4) = sum(overage)/size(thisBetacritData,1);
    n
end

toc;

% Plot this data.
figure;
h = errorbar(MSE_data(:,1),MSE_data(:,2),MSE_errdata(:,2));
% set(get(h,'Parent'), 'XScale', 'log');
% axis tight;
str = sprintf('Error of estimates of {\\beta_1}^c with %d samples',m);
title(str);
str = sprintf('percent error of estimated value of {\\beta_1}^c');
ylabel(str);
xlabel('number of initial data points given');
str = sprintf('beta_perc_%dsamples_%dminN_%dmaxN_sym.png',m,min(allN),max(allN));
print('-dpng', str);

% Plot this data.
figure;
h = errorbar(MSE_data(:,1),MSE_data(:,3),MSE_errdata(:,3));
hold on;
% set(get(h,'Parent'), 'XScale', 'log');
% axis tight;
plot([allN(1) allN(end)],[Truebeta_crit(1,1) Truebeta_crit(1,1)]);
legend('mean square error', 'true {\beta_1}^c value');
str = sprintf('Error of estimates of {\\beta_1}^c with %d samples',m);
title(str);
str = sprintf('mean square error of estimated value of {\\beta_1}^c');
ylabel(str);
xlabel('number of initial data points given');
str = sprintf('beta_mse_%dsamples_%dminN_%dmaxN_sym.png',m,min(allN),max(allN));
print('-dpng', str);
hold off;

% How frequently is beta over or under the true value?
figure;
scatter(MSE_data(:,1),MSE_data(:,4));
hold on;
h = lsline;
set(h,'LineWidth',3);
plot([allN(1) allN(end)],[0.5 0.5],'LineWidth',3);
legend('frequency of larger values', 'regression line','50%');
str = sprintf('Percent of {\\beta_1}^c larger than true {\\beta_1}^c with %d samples',m);
title(str);
str = sprintf('percent of estimated {\\beta_1}^c larger than true value');
ylabel(str);
xlabel('number of initial data points given');
str = sprintf('beta_overage_%dsamples_%dminN_%dmaxN_sym.png',m,min(allN),max(allN));
print('-dpng', str);
hold off;

testlogic = logical(MSE_data(:,4)<0.5);
testlogic = sum(testlogic)/size(testlogic,1);