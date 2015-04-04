% Script that generates plots.

close all;
clear all;

%% Generate plots for information curve. 
tic;

% Non-symmetric Sigma_X
n_x = 2;
n_y = 1;
sigma_X = [1 0.5; 0.5 2];
sigma_y = 1;
sigma_Y = sigma_y*eye(n_y);
sigma_XY = [0.1; 0.2];

% Choose beta.
allBeta = 10*(1:500);
% Choose m (number of times to sample the Wishart distribution).
m = 500;
% Choose one n value (number of data points available).
n = 1000;

[IC_data_nosym, IC_errdata_nosym] = plotCostFunc(sigma_X,sigma_Y, ...
    sigma_XY, allBeta, m, n );

% Symmetric Sigma_X
sigma_x = 1;
sigma_X = sigma_x*eye(n_x);
[IC_data_sym, IC_errdata_sym] = plotCostFunc(sigma_X,sigma_Y, ...
    sigma_XY, allBeta, m, n );

% Plot "scatterplot" information curve for the non-sym case and sym case. 
figure;
hold on;
h1 = plot(IC_data_nosym(:,2), IC_data_nosym(:,3),'k');
h2 = errorbar(IC_data_nosym(:,4),IC_data_nosym(:,5),IC_errdata_nosym(:,3),'.g');
h3 = herrorbar(IC_data_nosym(:,4),IC_data_nosym(:,5),IC_errdata_nosym(:,2),'.g');

h4 = plot(IC_data_sym(:,2), IC_data_sym(:,3),'b');
h5 = errorbar(IC_data_sym(:,4),IC_data_sym(:,5),IC_errdata_sym(:,3),'.m');
h6 = herrorbar(IC_data_sym(:,4),IC_data_sym(:,5),IC_errdata_sym(:,2),'.m');
legend([h1 h2 h4 h5], {'true information curve for non-symmetric {\Sigma}_x', ...
    'estimated information curve for non-symmetric {\Sigma}_x',...
    'true information curve for symmetric {\Sigma}_x', ...
    'estimated information curve for symmetric {\Sigma}_x'});
str = 'Information curve error bars';
title(str);
ylabel('I(T;Y)');
xlabel('I(X;T)');
str = sprintf('%s_%dn_%dm.png','infoCurve_both',n,m);
print('-dpng', str);
hold off;

% Should be roughly 210 sec (3.5 min)
timeElapsed = toc;

% Plot for just symmetric
figure;
hold on;
h4 = plot(IC_data_sym(:,2), IC_data_sym(:,3),'b');
h5 = errorbar(IC_data_sym(:,4),IC_data_sym(:,5),IC_errdata_sym(:,3),'.m');
h6 = herrorbar(IC_data_sym(:,4),IC_data_sym(:,5),IC_errdata_sym(:,2),'.m');
legend([h4 h5], {'true information curve', 'estimated information curve'});
str = 'Information curve error bars for symmetric {\Sigma}_x';
title(str);
ylabel('I(T;Y)');
xlabel('I(X;T)');
str = sprintf('%s_%dn_%dm.png','infoCurve_sym',n,m);
print('-dpng', str);
hold off;

% Plot for just non-symmetric
figure;
hold on;
h1 = plot(IC_data_nosym(:,2), IC_data_nosym(:,3),'k');
h2 = errorbar(IC_data_nosym(:,4),IC_data_nosym(:,5),IC_errdata_nosym(:,3),'.g');
h3 = herrorbar(IC_data_nosym(:,4),IC_data_nosym(:,5),IC_errdata_nosym(:,2),'.g');
legend([h1 h2], {'true information curve', 'estimated information curve'});
str = 'Information curve error bars for non-symmetric {\Sigma}_x';
title(str);
ylabel('I(T;Y)');
xlabel('I(X;T)');
str = sprintf('%s_%dn_%dm.png','infoCurve_nosym',n,m);
print('-dpng', str);
hold off;

%% Generate plots for beta critical error analysis.
%% Error in beta, log plot

% Log plot 
beta = 100;
% Choose m (number of times to sample the Wishart distribution).
m = 500;
% Choose variety of n values (number of data points available).
allN = 100*(10:100);

sigma_X = [1 0.5; 0.5 2];
sigma_y = 1;
sigma_Y = sigma_y*eye(n_y);
sigma_XY = [0.1; 0.2];

[Beta_data_nosym, Beta_errdata_nosym] = plotBetaError(sigma_X,sigma_Y, ...
    sigma_XY, beta, m, allN );

[TrueA, Truebeta_crit_nosym] = gib_optimize(sigma_X,sigma_Y,sigma_XY,beta);

sigma_X = sigma_x*eye(n_x);
[Beta_data_sym, Beta_errdata_sym] = plotBetaError(sigma_X,sigma_Y, ...
    sigma_XY, beta, m, allN );

[TrueA, Truebeta_crit_sym] = gib_optimize(sigma_X,sigma_Y,sigma_XY,beta);

% Plot percent error on log scale. 
figure;
h1 = errorbar(Beta_data_sym(1:end,1),Beta_data_sym(1:end,2),...
    Beta_errdata_sym(1:end,2));
hold on;
h2 = errorbar(Beta_data_nosym(1:end,1),Beta_data_nosym(1:end,2),...
    Beta_errdata_nosym(1:end,2));
% ylim([0 max(MSE_data(91:end,2))*1.2]);
% set(get(h1,'Parent'), 'XScale', 'log');
% axis tight;
legend('error for symmetric {\Sigma}_x','error for non-symmetric {\Sigma}_x');
str = sprintf('Error of estimates of {\\beta_1}^c with %d samples',m);
title(str);
str = sprintf('percent error of estimated value of {\\beta_1}^c');
ylabel(str);
xlabel('number of initial data points given');
str = sprintf('beta_perc_%dsamples_%dminN_%dmaxN.png',m,min(allN),max(allN));
print('-dpng', str);

% Plot mean square error on log scale. 
figure;
h1 = errorbar(Beta_data_sym(:,1),Beta_data_sym(:,3),Beta_errdata_sym(:,3));
hold on;
plot([allN(1) allN(end)],[Truebeta_crit_sym(1,1) Truebeta_crit_sym(1,1)]);
h2 = errorbar(Beta_data_nosym(:,1),Beta_data_nosym(:,3),Beta_errdata_nosym(:,3));
set(get(h1,'Parent'), 'XScale', 'log');
axis tight;
plot([allN(1) allN(end)],[Truebeta_crit_nosym(1,1) Truebeta_crit_nosym(1,1)]);
legend('mean square error for symmetric {\Sigma}_x', ...
    'true {\beta_1}^c value for symmetric {\Sigma}_x',...
    'mean square error for non-symmetric {\Sigma}_x',...
    'true {\beta_1}^c value for non-symmetric {\Sigma}_x');
str = sprintf('Error of estimates of {\\beta_1}^c with %d samples',m);
title(str);
str = sprintf('mean square error of estimated value of {\\beta_1}^c');
ylabel(str);
xlabel('number of initial data points given');
str = sprintf('beta_mse_%dsamples_%dminN_%dmaxN.png',m,min(allN),max(allN));
print('-dpng', str);
hold off;

% How frequently is beta over or under the true value?
figure;
h1 = scatter(Beta_data_sym(:,1),Beta_data_sym(:,4));
hold on;
h2 = scatter(Beta_data_nosym(:,1),Beta_data_nosym(:,4));
plot([allN(1) allN(end)],[0.5 0.5],'LineWidth',3);
set(get(h1,'Parent'), 'XScale', 'log');
axis tight;
legend('frequency of overage for symmetric {\Sigma}_x', ...
    'frequency of overage for non-symmetric {\Sigma}_x',...
    '50%');
str = sprintf('Percent of {\\beta_1}^c larger than true {\\beta_1}^c with %d samples',m);
title(str);
str = sprintf('percent of estimated {\\beta_1}^c larger than true value');
ylabel(str);
xlabel('number of initial data points given');
str = sprintf('beta_overage_%dsamples_%dminN_%dmaxN.png',m,min(allN),max(allN));
print('-dpng', str);
hold off;

%% Beta analysis plots, linear plot

% Log plot 
beta = 100;
% Choose m (number of times to sample the Wishart distribution).
m = 500;
% Choose variety of n values (number of data points available).
allN = 10*(50:500);

sigma_X = [1 0.5; 0.5 2];
sigma_y = 1;
sigma_Y = sigma_y*eye(n_y);
sigma_XY = [0.1; 0.2];

[Beta_data_nosym, Beta_errdata_nosym] = plotBetaError(sigma_X,sigma_Y, ...
    sigma_XY, beta, m, allN );

[TrueA, Truebeta_crit_nosym] = gib_optimize(sigma_X,sigma_Y,sigma_XY,beta);

sigma_X = sigma_x*eye(n_x);
[Beta_data_sym, Beta_errdata_sym] = plotBetaError(sigma_X,sigma_Y, ...
    sigma_XY, beta, m, allN );

[TrueA, Truebeta_crit_sym] = gib_optimize(sigma_X,sigma_Y,sigma_XY,beta);

% Plot percent error on log scale. 
figure;
h1 = errorbar(Beta_data_sym(1:end,1),Beta_data_sym(1:end,2),...
    Beta_errdata_sym(1:end,2));
hold on;
h2 = errorbar(Beta_data_nosym(1:end,1),Beta_data_nosym(1:end,2),...
    Beta_errdata_nosym(1:end,2));
% ylim([0 max(MSE_data(91:end,2))*1.2]);
legend('error for symmetric {\Sigma}_x','error for non-symmetric {\Sigma}_x');
str = sprintf('Error of estimates of {\\beta_1}^c with %d samples',m);
title(str);
str = sprintf('percent error of estimated value of {\\beta_1}^c');
ylabel(str);
xlabel('number of initial data points given');
str = sprintf('beta_perc_%dsamples_%dminN_%dmaxN.png',m,min(allN),max(allN));
print('-dpng', str);

% Plot mean square error on log scale. 
figure;
h1 = errorbar(Beta_data_sym(:,1),Beta_data_sym(:,3),Beta_errdata_sym(:,3));
hold on;
plot([allN(1) allN(end)],[Truebeta_crit_sym(1,1) Truebeta_crit_sym(1,1)]);
h2 = errorbar(Beta_data_nosym(:,1),Beta_data_nosym(:,3),Beta_errdata_nosym(:,3));
plot([allN(1) allN(end)],[Truebeta_crit_nosym(1,1) Truebeta_crit_nosym(1,1)]);
legend('mean square error for symmetric {\Sigma}_x', ...
    'true {\beta_1}^c value for symmetric {\Sigma}_x',...
    'mean square error for non-symmetric {\Sigma}_x',...
    'true {\beta_1}^c value for non-symmetric {\Sigma}_x');
str = sprintf('Error of estimates of {\\beta_1}^c with %d samples',m);
title(str);
str = sprintf('mean square error of estimated value of {\\beta_1}^c');
ylabel(str);
xlabel('number of initial data points given');
str = sprintf('beta_mse_%dsamples_%dminN_%dmaxN.png',m,min(allN),max(allN));
print('-dpng', str);
hold off;

% How frequently is beta over or under the true value?
figure;
h1 = scatter(Beta_data_sym(:,1),Beta_data_sym(:,4),'b','*');
hold on;
h2 = polyfit(Beta_data_sym(:,1),Beta_data_sym(:,4),1);
h2 = polyval(h2,Beta_data_sym(:,1));
plot(Beta_data_sym(:,1),h2);
h3 = scatter(Beta_data_nosym(:,1),Beta_data_nosym(:,4),'k','x');
h4 = polyfit(Beta_data_nosym(:,1),Beta_data_nosym(:,4),1);
h4 = polyval(h4,Beta_data_nosym(:,1));
plot(Beta_data_nosym(:,1),h4);
h5 = plot([allN(1) allN(end)],[0.5 0.5],'LineWidth',2);
legend('frequency of overage for symmetric {\Sigma}_x', ...
    'regression line for symmetric {\Sigma}_x',...
    'frequency of overage for non-symmetric {\Sigma}_x',...
    'regression line for non-symmetric {\Sigma}_x',...
    '50%','Location','southeast');
str = sprintf('Percent of {\\beta_1}^c larger than true {\\beta_1}^c with %d samples',m);
title(str);
str = sprintf('percent of estimated {\\beta_1}^c larger than true value');
ylabel(str);
xlabel('number of initial data points given');
str = sprintf('beta_overage_%dsamples_%dminN_%dmaxN.png',m,min(allN),max(allN));
print('-dpng', str);
hold off;

%% Generate plots for analysis of error in A and degenerate cases.  
%% Error in A, log plot

% Log plot 
beta = 100;
% Choose m (number of times to sample the Wishart distribution).
m = 500;
% Choose variety of n values (number of data points available).
allN = 100*(10:10000);

sigma_X = [1 0.5; 0.5 2];
[ MSE_data_nosym, MSE_errdata_nosym, ABS_data_nosym, ABS_errdata_nosym ] = plotAError( sigma_X, ...
    sigma_Y, sigma_XY, beta, m, allN);

[TrueA_nosym, Truebeta_crit_nosym] = gib_optimize(sigma_X,sigma_Y,sigma_XY,beta);

sigma_X = sigma_x*eye(n_x);
[ MSE_data_sym, MSE_errdata_sym, ABS_data_sym, ABS_errdata_sym ] = plotAError( sigma_X, ...
    sigma_Y, sigma_XY, beta, m, allN);

[TrueA_sym, Truebeta_crit_sym] = gib_optimize(sigma_X,sigma_Y,sigma_XY,beta);

% Plot how often betacrit is greater than beta.
figure;
semilogx(MSE_data_nosym(:,1), MSE_data_nosym(:,end));
hold on;
semilogx(MSE_data_sym(:,1), MSE_data_sym(:,end));
set(get(h,'Parent'), 'XScale', 'log');
str = 'Frequency of {\beta}_{1}^{c} greater than \beta';
title(str);
ylabel('percent of degenerate solutions');
xlabel('number of initial data points given');
str = sprintf('beta_failure_%dsamples_largeN.png',m);
print('-dpng', str);
hold off;

% Plot MSE error in elements of A, log.
for k = 1:2
    figure;
    h1 = errorbar(MSE_data_sym(:,1),MSE_data_sym(:,k+1),MSE_errdata_sym(:,k+1));
    hold on;
    h2 = plot([allN(1) allN(end)],[abs(TrueA_sym(1,k)) abs(TrueA_sym(1,k))]);
    h3 = errorbar(MSE_data_nosym(:,1), MSE_data_nosym(:,k+1),MSE_errdata_nosym(:,k+1));
    h4 = plot([allN(1) allN(end)],[abs(TrueA_nosym(1,k)) abs(TrueA_nosym(1,k))]);
    legend('error for symmetric {\Sigma}_x', 'absolute value of true value for symmetric',...
        'error for non-symmetric {\Sigma}_x', 'absolute value of true value for non-symmetric');
    xlim([0 1.05*allN(end)]);
    set(get(h1,'Parent'), 'XScale', 'log');
                
    if(k==1)
        str = sprintf('Error of estimates of %dst element of projection matrix with %d samples',k,m);
    elseif(k==2)
        str = sprintf('Error of estimates of %dnd element of projection matrix with %d samples',k,m);
    elseif(k==3)
        str = sprintf('Error of estimates of %drd element of projection matrix with %d samples',k,m);
    else
        str = sprintf('Error of estimates of %dth element of projection matrix with %d samples',k,m);
    end
    title(str);
    
    if(k==1)
        str = sprintf('mean square error of estimated value of %dst element of A', k);
    elseif(k==2)
        str = sprintf('mean square error of estimated value of %dnd element of A', k);
    elseif(k==3)
        str = sprintf('mean square error of estimated value of %drd element of A', k);
    else
        str = sprintf('mean square error of estimated value of %dth element of A', k);
    end
    ylabel(str);
        
    xlabel('number of initial data points given');
    str = sprintf('mse_filt_%dth_%dsamples_%dminN_%dmaxN.png',k,m,min(allN),allN(end));
    print('-dpng', str);
    
    hold off;
end

% Plot percent error in elements of A, log. 
for k = 1:2
    figure;
    h = errorbar(ABS_data_sym(:,1),ABS_data_sym(:,k+1),ABS_errdata_sym(:,k+1));
    hold on;
    h2 = errorbar(ABS_data_nosym(:,1),ABS_data_nosym(:,k+1),ABS_errdata_nosym(:,k+1));
    xlim([0 1.05*allN(end)]);
    legend('error for symmetric {\Sigma}_x', 'error for non-symmetric {\Sigma}_x');
    set(get(h,'Parent'), 'XScale', 'log');

    if(k==1)
        str = sprintf('Error of estimates of %dst element of projection matrix with %d samples',k,m);
    elseif(k==2)
        str = sprintf('Error of estimates of %dnd element of projection matrix with %d samples',k,m);
    elseif(k==3)
        str = sprintf('Error of estimates of %drd element of projection matrix with %d samples',k,m);
    else
        str = sprintf('Error of estimates of %dth element of projection matrix with %d samples',k,m);
    end
    title(str);
    
    if (k==1)
        str = sprintf('percent error of estimated value of %dst element of A', k);
    elseif(k==2)
        str = sprintf('percent error of estimated value of %dnd element of A', k);
    elseif(k==3)
        str = sprintf('percent error of estimated value of %drd element of A', k);
    else
        str = sprintf('percent error of estimated value of %dth element of A', k);
    end
    ylabel(str);
    
    xlabel('number of initial data points given');
    str = sprintf('perc_filt_%dth_%dsamples_%dminN_%dmaxN.png',k,m,min(allN),allN(end));
    print('-dpng', str);
    
    hold off;
end

%% Error in A, linear plot

% Pick beta. 
beta = 100;
% Choose m (number of times to sample the Wishart distribution).
m = 500;
% Choose variety of n values (number of data points available).
allN = 10*(5:500);

sigma_X = [1 0.5; 0.5 2];
[ MSE_data_nosym, MSE_errdata_nosym, ABS_data_nosym, ABS_errdata_nosym ] = plotAError( sigma_X, ...
    sigma_Y, sigma_XY, beta, m, allN);

[TrueA_nosym, Truebeta_crit_nosym] = gib_optimize(sigma_X,sigma_Y,sigma_XY,beta);

sigma_X = sigma_x*eye(n_x);
[ MSE_data_sym, MSE_errdata_sym, ABS_data_sym, ABS_errdata_sym ] = plotAError( sigma_X, ...
    sigma_Y, sigma_XY, beta, m, allN);

[TrueA_sym, Truebeta_crit_sym] = gib_optimize(sigma_X,sigma_Y,sigma_XY,beta);

% Plot how often betacrit is greater than beta.
figure;
plot(MSE_data_nosym(:,1), MSE_data_nosym(:,end));
hold on;
plot(MSE_data_sym(:,1), MSE_data_sym(:,end));
legend('non-symmetric {\Sigma}_x', 'symmetric {\Sigma}_x');
str = 'Frequency of {\beta}_{1}^{c} greater than \beta';
title(str);
ylabel('percent of degenerate solutions');
xlabel('number of initial data points given');
str = sprintf('beta_failure_%dsamples_smallN.png',m);
print('-dpng', str);
hold off;

% Plot MSE error in elements of A, linear.
for k = 1:2
    figure;
    h1 = errorbar(MSE_data_sym(:,1),MSE_data_sym(:,k+1),MSE_errdata_sym(:,k+1));
    hold on;
    h2 = plot([allN(1) allN(end)],[abs(TrueA_sym(1,k)) abs(TrueA_sym(1,k))]);
    h3 = errorbar(MSE_data_nosym(:,1), MSE_data_nosym(:,k+1),MSE_errdata_nosym(:,k+1));
    h4 = plot([allN(1) allN(end)],[abs(TrueA_nosym(1,k)) abs(TrueA_nosym(1,k))]);
    legend('error for symmetric {\Sigma}_x', 'absolute value of true value for symmetric',...
        'error for non-symmetric {\Sigma}_x', 'absolute value of true value for non-symmetric');
    xlim([0 1.05*allN(end)]);
%    set(get(h1,'Parent'), 'XScale', 'log');
                
    if(k==1)
        str = sprintf('Error of estimates of %dst element of projection matrix with %d samples',k,m);
    elseif(k==2)
        str = sprintf('Error of estimates of %dnd element of projection matrix with %d samples',k,m);
    elseif(k==3)
        str = sprintf('Error of estimates of %drd element of projection matrix with %d samples',k,m);
    else
        str = sprintf('Error of estimates of %dth element of projection matrix with %d samples',k,m);
    end
    title(str);
    
    if(k==1)
        str = sprintf('mean square error of estimated value of %dst element of A', k);
    elseif(k==2)
        str = sprintf('mean square error of estimated value of %dnd element of A', k);
    elseif(k==3)
        str = sprintf('mean square error of estimated value of %drd element of A', k);
    else
        str = sprintf('mean square error of estimated value of %dth element of A', k);
    end
    ylabel(str);
        
    xlabel('number of initial data points given');
    str = sprintf('mse_filt_%dth_%dsamples_%dminN_%dmaxN.png',k,m,min(allN),allN(end));
    print('-dpng', str);
    
    hold off;
end

% Plot percent error in elements of A, linear. 
for k = 1:2
    figure;
    h = errorbar(ABS_data_sym(:,1),ABS_data_sym(:,k+1),ABS_errdata_sym(:,k+1));
    hold on;
    h2 = errorbar(ABS_data_nosym(:,1),ABS_data_nosym(:,k+1),ABS_errdata_nosym(:,k+1));
    xlim([0 1.05*allN(end)]);
    legend('error for symmetric {\Sigma}_x', 'error for non-symmetric {\Sigma}_x');
%     set(get(h,'Parent'), 'XScale', 'log');

    if(k==1)
        str = sprintf('Error of estimates of %dst element of projection matrix with %d samples',k,m);
    elseif(k==2)
        str = sprintf('Error of estimates of %dnd element of projection matrix with %d samples',k,m);
    elseif(k==3)
        str = sprintf('Error of estimates of %drd element of projection matrix with %d samples',k,m);
    else
        str = sprintf('Error of estimates of %dth element of projection matrix with %d samples',k,m);
    end
    title(str);
    
    if (k==1)
        str = sprintf('percent error of estimated value of %dst element of A', k);
    elseif(k==2)
        str = sprintf('percent error of estimated value of %dnd element of A', k);
    elseif(k==3)
        str = sprintf('percent error of estimated value of %drd element of A', k);
    else
        str = sprintf('percent error of estimated value of %dth element of A', k);
    end
    ylabel(str);
    
    xlabel('number of initial data points given');
    str = sprintf('perc_filt_%dth_%dsamples_%dminN_%dmaxN.png',k,m,min(allN),allN(end));
    print('-dpng', str);
    hold off;
end