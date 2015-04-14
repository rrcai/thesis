% Script that generates plots.

close all;
clear all;

%% Generate plots for information curve. 
tic;

% Non-diagonal Sigma_X
n_x = 2;
n_y = 1;
sigma_X = [1 0.5; 0.5 2];
sigma_y = 1;
sigma_Y = sigma_y*eye(n_y);
sigma_XY = [0.1; 0.2];

% Choose beta.
allBeta = 1*(1:500);

% Choose m (number of times to sample the Wishart distribution).
m = 500;
% Choose one n value (number of data points available).
n = 5000;

disp('running IC nodiag\n');
[IC_data_nosym, IC_errdata_nosym] = plotCostFunc(sigma_X,sigma_Y, ...
    sigma_XY, allBeta, m, n );

% diagonal Sigma_X
sigma_x = 1;
sigma_X = sigma_x*eye(n_x);
disp('running IC diag\n');
[IC_data_sym, IC_errdata_sym] = plotCostFunc(sigma_X,sigma_Y, ...
    sigma_XY, allBeta, m, n );

IC_perc_error_sym = zeros(size(IC_data_sym,1),3);
IC_perc_error_nosym = zeros(size(IC_data_nosym,1),3);

IC_perc_error_sym(:,1) = allBeta;
IC_perc_error_nosym(:,1) = allBeta;

IC_perc_error_sym(:,2) = abs((IC_data_sym(:,4)-IC_data_sym(:,2))./IC_data_sym(:,2));
IC_perc_error_sym(:,3) = abs((IC_data_sym(:,5)-IC_data_sym(:,3))./IC_data_sym(:,3));

IC_perc_error_nosym(:,2) = abs((IC_data_nosym(:,4)-IC_data_nosym(:,2))./IC_data_sym(:,2));
IC_perc_error_nosym(:,3) = abs((IC_data_nosym(:,5)-IC_data_nosym(:,3))./IC_data_nosym(:,3));

% Plot percent errors.
figure;
hold on;
h1 = scatter(IC_perc_error_sym(21:end,1),IC_perc_error_sym(21:end,2)); %past
h3 = plot([allBeta(21) allBeta(end)], [mean(IC_perc_error_sym(21:end,2)) mean(IC_perc_error_sym(21:end,2))]);
h2 = scatter(IC_perc_error_sym(21:end,1),IC_perc_error_sym(21:end,3)); %future
h4 = plot([allBeta(21) allBeta(end)], [mean(IC_perc_error_sym(21:end,3)) mean(IC_perc_error_sym(21:end,3))]);
legend([h1 h2 h3 h4], {'percent error in past information','percent error in future information',...
    'mean percent error in past information','mean percent error in future information'});
xlabel('{\beta}');
ylabel('percent error');
title('Percent error of past and future information for diagonal {\Sigma}_x');
hold off;
filestr = 'percErrInfoCurveDiag';
print('-dpng',filestr);

figure;
hold on;
h1 = scatter(IC_perc_error_nosym(44:end,1),IC_perc_error_nosym(44:end,2)); %past
h2 = scatter(IC_perc_error_nosym(44:end,1),IC_perc_error_nosym(44:end,3)); %future
legend([h1 h2], {'percent error in past information','percent error in future information'});
xlabel('{\beta}');
ylabel('percent error');
title('Percent error of past and future information for non-diagonal {\Sigma}_x');
hold off;
filestr = 'percErrInfoCurveNoDiag';
print('-dpng',filestr);

% Histogram of percent errors? 
figure; 
h1 = histogram(IC_perc_error_sym(21:end,2),40);
h = findobj(gca,'Type','patch');
set(h,'FaceColor','r','EdgeColor','w','facealpha',0.75);
hold on;
h3 = histogram(IC_perc_error_nosym(50:end,2),40);
h = findobj(gca,'Type','patch');
set(h,'facealpha',0.75);
ylabel('frequency of percent error');
xlabelstr = 'percent error in past information';
xlabel(xlabelstr)
legend([h1 h3], {'diagonal {\Sigma}_x','non-diagonal {\Sigma}_x'});
titlestr = sprintf('Frequency of percent error in past information');
title(titlestr);
str = sprintf('pastInfoErr_hist_%dsamples.png',m);
print('-dpng', str);
hold off;

figure; 
h1 = histogram(IC_perc_error_sym(21:end,3),30);
h = findobj(gca,'Type','patch');
set(h,'FaceColor','r','EdgeColor','w','facealpha',0.75);
hold on;
h3 = histogram(IC_perc_error_nosym(50:end,3),40);
h = findobj(gca,'Type','patch');
set(h,'facealpha',0.75);
ylabel('frequency of percent error');
xlabelstr = 'percent error in future information';
xlabel(xlabelstr)
legend([h1 h3], {'diagonal {\Sigma}_x','non-diagonal {\Sigma}_x'});
titlestr = sprintf('Frequency of percent error in future information');
title(titlestr);
str = sprintf('futureInfoErr_hist_%dsamples.png',m);
print('-dpng', str);
hold off;

% Plot "scatterplot" information curve for the non-sym case and sym case. 
figure;
hold on;
h1 = plot(IC_data_nosym(:,2), IC_data_nosym(:,3),'k');
h2 = scatter(IC_data_nosym(:,4),IC_data_nosym(:,5));
h4 = plot(IC_data_sym(:,2), IC_data_sym(:,3),'b');
h5 = scatter(IC_data_sym(:,4),IC_data_sym(:,5));
legend([h1 h2 h4 h5], {'true information curve for non-diagonal {\Sigma}_x', ...
    'estimated information curve for non-diagonal {\Sigma}_x',...
    'true information curve for diagonal {\Sigma}_x', ...
    'estimated information curve for diagonal {\Sigma}_x'});
str = 'Information curve error bars';
yl = get(gca, 'ylim');
yl(1) = 0;
ylim(yl);
xl = get(gca, 'xlim');
xl(1) = 0;
xlim(xl);
title(str);
ylabel('I(T;Y)');
xlabel('I(X;T)');
str = sprintf('%s_%dminBeta_%dmaxBeta_%dn_%dm.png','infoCurve_both',min(allBeta),...
    max(allBeta),n,m);
print('-dpng', str);
hold off;

% Plot for just symmetric
figure;
hold on;
h4 = plot(IC_data_sym(:,2), IC_data_sym(:,3),'b');
h5 = scatter(IC_data_sym(:,4),IC_data_sym(:,5));
legend([h4 h5], {'true information curve', 'estimated information curve'});
str = 'Information curve error bars for diagonal {\Sigma}_x';
title(str);
ylabel('I(T;Y)');
xlabel('I(X;T)');
yl = get(gca, 'ylim');
yl(1) = 0;
ylim(yl);
xl = get(gca, 'xlim');
xl(1) = 0;
xlim(xl);
str = sprintf('%s_%dminBeta_%dmaxBeta_%dn_%dm.png','infoCurve_diag',min(allBeta),...
    max(allBeta),n,m);
print('-dpng', str);
hold off;

% Plot for just symmetric with bars 
figure;
hold on;
h4 = plot(IC_data_sym(:,2), IC_data_sym(:,3),'b');
h5 = scatter(IC_data_sym(:,4),IC_data_sym(:,5),'.m');
testx = [IC_data_sym(:,2) IC_data_sym(:,4)];
testy = [IC_data_sym(:,3) IC_data_sym(:,5)];
plot(testx', testy');
yl = get(gca, 'ylim');
yl(1) = 0;
ylim(yl);
xl = get(gca, 'xlim');
xl(1) = 0;
xlim(xl);
legend([h4 h5], {'true information curve', 'estimated information curve'});
str = 'Information curve error bars for diagonal {\Sigma}_x';
title(str);
ylabel('I(T;Y)');
xlabel('I(X;T)');
str = sprintf('%s_withBars_%dminBeta_%dmaxBeta_%dn_%dm.png','infoCurve_diag',...
    min(allBeta),max(allBeta),n,m);
print('-dpng', str);
hold off;

% Plot for just non-symmetric
figure;
hold on;
h4 = plot(IC_data_nosym(:,2), IC_data_nosym(:,3),'b');
h5 = scatter(IC_data_nosym(:,4),IC_data_nosym(:,5));
legend([h4 h5], {'true information curve', 'estimated information curve'});
str = 'Information curve error bars for non-diagonal {\Sigma}_x';
title(str);
ylabel('I(T;Y)');
xlabel('I(X;T)');
yl = get(gca, 'ylim');
yl(1) = 0;
ylim(yl);
xl = get(gca, 'xlim');
xl(1) = 0;
xlim(xl);
str = sprintf('%s_%dminBeta_%dmaxBeta_%dn_%dm.png','infoCurve_nodiag',min(allBeta),...
    max(allBeta),n,m);
print('-dpng', str);
hold off;

% Plot for just non-diagonal with bars
figure;
hold on;
h1 = plot(IC_data_nosym(:,2), IC_data_nosym(:,3),'k');
h2 = scatter(IC_data_nosym(:,4),IC_data_nosym(:,5),'.m');
testx = [IC_data_nosym(:,2) IC_data_nosym(:,4)];
testy = [IC_data_nosym(:,3) IC_data_nosym(:,5)];
h3 = plot(testx', testy');
legend([h1 h2], {'true information curve', 'estimated information curve'});
str = 'Information curve error bars for non-diagonal {\Sigma}_x';
title(str);
ylabel('I(T;Y)');
xlabel('I(X;T)');
yl = get(gca, 'ylim');
yl(1) = 0;
ylim(yl);
xl = get(gca, 'xlim');
xl(1) = 0;
xlim(xl);
str = sprintf('%s_withBars_%dminBeta_%dmaxBeta_%dn_%dm.png','infoCurve_nodiag',...
    min(allBeta),max(allBeta),n,m);
print('-dpng', str);
hold off;

%% Generate plots for beta critical error analysis.
%% Error in beta, log plot
tic;
% Log plot 
beta = 100;
% Choose m (number of times to sample the Wishart distribution).
m = 500;
% Choose variety of n values (number of data points available).
allN = logspace(3,6,900);

n_x = 2;
n_y = 1;

sigma_X = [1 0.5; 0.5 2];
sigma_y = 1;
sigma_x = 1;
sigma_Y = sigma_y*eye(n_y);
sigma_XY = [0.1; 0.2];

disp('running beta error no sym');
[Beta_data_nosym, Beta_errdata_nosym] = plotBetaError(sigma_X,sigma_Y, ...
    sigma_XY, beta, m, allN );

[TrueA, Truebeta_crit_nosym] = gib_optimize(sigma_X,sigma_Y,sigma_XY,beta);

disp('running beta error sym');
sigma_X = sigma_x*eye(n_x);
[Beta_data_sym, Beta_errdata_sym] = plotBetaError(sigma_X,sigma_Y, ...
    sigma_XY, beta, m, allN );

[TrueA, Truebeta_crit_sym] = gib_optimize(sigma_X,sigma_Y,sigma_XY,beta);

BetaLog_TE = toc;

% Plot percent error on log scale. 
figure;
h1 = errorbar(Beta_data_sym(1:end,1),Beta_data_sym(1:end,2),...
    Beta_errdata_sym(1:end,2));
hold on;
h2 = errorbar(Beta_data_nosym(1:end,1),Beta_data_nosym(1:end,2),...
    Beta_errdata_nosym(1:end,2));
set(get(h1,'Parent'), 'XScale', 'log');
axis tight;
legend('error for diagonal {\Sigma}_x','error for non-diagonal {\Sigma}_x');
str = sprintf('Error of estimates of {\\beta_1}^c with %d samples',m);
title(str);
str = sprintf('percent error of estimated value of {\\beta_1}^c');
ylabel(str);
xlabel('number of initial data points given');
str = sprintf('beta_perc_%dsamples_%dminN_%dmaxN.png',m,min(allN),max(allN));
print('-dpng', str);
hold off;

% Plot mean square error on log scale. Also do separate sym and nosym
% plots.
figure;
h1 = errorbar(Beta_data_sym(:,1),Beta_data_sym(:,3),Beta_errdata_sym(:,3));
hold on;
plot([allN(1) allN(end)],[Truebeta_crit_sym(1,1) Truebeta_crit_sym(1,1)]);
h2 = errorbar(Beta_data_nosym(:,1),Beta_data_nosym(:,3),Beta_errdata_nosym(:,3));
set(get(h1,'Parent'), 'XScale', 'log');
axis tight;
plot([allN(1) allN(end)],[Truebeta_crit_nosym(1,1) Truebeta_crit_nosym(1,1)]);
legend('mean square error for diagonal {\Sigma}_x', ...
    'true {\beta_1}^c value for diagonal {\Sigma}_x',...
    'mean square error for non-diagonal {\Sigma}_x',...
    'true {\beta_1}^c value for non-diagonal {\Sigma}_x');
str = sprintf('Error of estimates of {\\beta_1}^c with %d samples',m);
title(str);
str = sprintf('mean square error of estimated value of {\\beta_1}^c');
ylabel(str);
xlabel('number of initial data points given');
str = sprintf('beta_mse_%dsamples_%dminN_%dmaxN.png',m,min(allN),max(allN));
print('-dpng', str);
hold off;

% Symmetric plot for error in beta
figure;
h1 = errorbar(Beta_data_sym(:,1),Beta_data_sym(:,3),Beta_errdata_sym(:,3));
hold on;
plot([allN(1) allN(end)],[Truebeta_crit_sym(1,1) Truebeta_crit_sym(1,1)]);
set(get(h1,'Parent'), 'XScale', 'log');
axis tight;
legend('mean square error for diagonal {\Sigma}_x', ...
    'true {\beta_1}^c value for diagonal {\Sigma}_x');
str = sprintf('Error of estimates of {\\beta_1}^c with %d samples',m);
title(str);
str = sprintf('mean square error of estimated value of {\\beta_1}^c');
ylabel(str);
xlabel('number of initial data points given');
str = sprintf('beta_mse_%dsamples_%dminN_%dmaxN_diag.png',m,min(allN),max(allN));
print('-dpng', str);
hold off;

% Non symmetric plot of mean square error for beta
figure;
hold on;
h2 = errorbar(Beta_data_nosym(:,1),Beta_data_nosym(:,3),Beta_errdata_nosym(:,3));
set(get(h2,'Parent'), 'XScale', 'log');
axis tight;
plot([allN(1) allN(end)],[Truebeta_crit_nosym(1,1) Truebeta_crit_nosym(1,1)]);
legend('mean square error for non-diagonal {\Sigma}_x',...
    'true {\beta_1}^c value for non-diagonal {\Sigma}_x');
str = sprintf('Error of estimates of {\\beta_1}^c with %d samples',m);
title(str);
str = sprintf('mean square error of estimated value of {\\beta_1}^c');
ylabel(str);
xlabel('number of initial data points given');
str = sprintf('beta_mse_%dsamples_%dminN_%dmaxN_nosym.png',m,min(allN),max(allN));
print('-dpng', str);
hold off;

% How frequently is beta over or under the true value? make histogram
figure;
h1 = scatter(Beta_data_sym(:,1),Beta_data_sym(:,4));
hold on;
h2 = scatter(Beta_data_nosym(:,1),Beta_data_nosym(:,4));
plot([allN(1) allN(end)],[0.5 0.5],'LineWidth',3);
set(get(h1,'Parent'), 'XScale', 'log');
axis tight;
legend('frequency of overage for diagonal {\Sigma}_x', ...
    'frequency of overage for non-diagonal {\Sigma}_x',...
    '50%');
str = sprintf('Percent of {\\beta_1}^c larger than true {\\beta_1}^c with %d samples',m);
title(str);
str = sprintf('percent of estimated {\\beta_1}^c larger than true value');
ylabel(str);
xlabel('number of initial data points given');
str = sprintf('beta_overage_%dsamples_%dminN_%dmaxN.png',m,min(allN),max(allN));
print('-dpng', str);
hold off;

% Histogram of overages
figure; 
h1 = histogram(Beta_data_sym(:,4),20);
h = findobj(gca,'Type','patch');
set(h,'FaceColor','r','EdgeColor','w','facealpha',0.75);
hold on;
h2 = line([0.5 0.5], [0 0.2*size(allN,2)]);
h3 = histogram(Beta_data_nosym(:,4),20);
h = findobj(gca,'Type','patch');
set(h,'facealpha',0.75);
ylabel('frequency of overage in estimates');
xlabelstr = sprintf('frequency across %d samples',m);
xlabel(xlabelstr)
legend([h1 h3], {'diagonal {\Sigma}_x','non-diagonal {\Sigma}_x'});
titlestr = sprintf('Frequency of overage in errors of {\\beta_1}^c across n from %d to %d', min(allN),max(allN));
title(titlestr);
str = sprintf('beta_overage_hist_%dsamples_%dminN_%dmaxN.png',m,min(allN),max(allN));
print('-dpng', str);
hold off;

percOver50symBig = logical(Beta_data_sym(:,4)>0.5);
percOver50symBig = sum(percOver50symBig)/size(percOver50symBig,1);

percOver50nosymBig = logical(Beta_data_nosym(:,4)>0.5);
percOver50nosymBig = sum(percOver50nosymBig)/size(percOver50nosymBig,1);

%% Beta analysis plots, linear plot

tic;

% Log plot 
beta = 100;
% Choose m (number of times to sample the Wishart distribution).
m = 500;
% Choose variety of n values (number of data points available).
allN = 10*(100:1000);

sigma_X = [1 0.5; 0.5 2];
sigma_y = 1;
sigma_Y = sigma_y*eye(n_y);
sigma_XY = [0.1; 0.2];

disp('beta linear error no diag');
[Beta_data_nosym, Beta_errdata_nosym] = plotBetaError(sigma_X,sigma_Y, ...
    sigma_XY, beta, m, allN );

[TrueA, Truebeta_crit_nosym] = gib_optimize(sigma_X,sigma_Y,sigma_XY,beta);

disp('beta linear error diag');
sigma_X = sigma_x*eye(n_x);
[Beta_data_sym, Beta_errdata_sym] = plotBetaError(sigma_X,sigma_Y, ...
    sigma_XY, beta, m, allN );

[TrueA, Truebeta_crit_sym] = gib_optimize(sigma_X,sigma_Y,sigma_XY,beta);

BetaLinear_TE = toc;

% Plot percent error on linear scale. 
figure;
h1 = errorbar(Beta_data_sym(1:end,1),Beta_data_sym(1:end,2),...
    Beta_errdata_sym(1:end,2));
hold on;
h2 = errorbar(Beta_data_nosym(1:end,1),Beta_data_nosym(1:end,2),...
    Beta_errdata_nosym(1:end,2));
ylim([0 max(Beta_data_nosym(1:end,2))*1.2]);
legend('error for diagonal {\Sigma}_x','error for non-diagonal {\Sigma}_x');
str = sprintf('Error of estimates of {\\beta_1}^c with %d samples',m);
title(str);
str = sprintf('percent error of estimated value of {\\beta_1}^c');
ylabel(str);
xlabel('number of initial data points given');
str = sprintf('beta_perc_%dsamples_%dminN_%dmaxN.png',m,min(allN),max(allN));
print('-dpng', str);
hold off;

% Plot mean square error on linear scale. 
figure;
h1 = errorbar(Beta_data_sym(:,1),Beta_data_sym(:,3),Beta_errdata_sym(:,3));
hold on;
plot([allN(1) allN(end)],[Truebeta_crit_sym(1,1) Truebeta_crit_sym(1,1)]);
h2 = errorbar(Beta_data_nosym(:,1),Beta_data_nosym(:,3),Beta_errdata_nosym(:,3));
plot([allN(1) allN(end)],[Truebeta_crit_nosym(1,1) Truebeta_crit_nosym(1,1)]);
legend('mean square error for diagonal {\Sigma}_x', ...
    'true {\beta_1}^c value for diagonal {\Sigma}_x',...
    'mean square error for non-diagonal {\Sigma}_x',...
    'true {\beta_1}^c value for non-diagonal {\Sigma}_x');
str = sprintf('Error of estimates of {\\beta_1}^c with %d samples',m);
title(str);
str = sprintf('mean square error of estimated value of {\\beta_1}^c');
ylabel(str);
xlabel('number of initial data points given');
str = sprintf('beta_mse_%dsamples_%dminN_%dmaxN.png',m,min(allN),max(allN));
print('-dpng', str);
hold off;

% mean square error, linear, symmetric
figure;
h1 = errorbar(Beta_data_sym(:,1),Beta_data_sym(:,3),Beta_errdata_sym(:,3));
hold on;
plot([allN(1) allN(end)],[Truebeta_crit_sym(1,1) Truebeta_crit_sym(1,1)]);
legend('mean square error for symmetric {\Sigma}_x', ...
    'true {\beta_1}^c value for symmetric {\Sigma}_x');
str = sprintf('Error of estimates of {\\beta_1}^c with %d samples',m);
title(str);
str = sprintf('mean square error of estimated value of {\\beta_1}^c');
ylabel(str);
xlabel('number of initial data points given');
str = sprintf('beta_mse_%dsamples_%dminN_%dmaxN_sym.png',m,min(allN),max(allN));
print('-dpng', str);
hold off;

% mean square error, non-symmetric, linear
figure;
h2 = errorbar(Beta_data_nosym(:,1),Beta_data_nosym(:,3),Beta_errdata_nosym(:,3));
hold on;
plot([allN(1) allN(end)],[Truebeta_crit_nosym(1,1) Truebeta_crit_nosym(1,1)]);
legend('mean square error for non-symmetric {\Sigma}_x',...
    'true {\beta_1}^c value for non-symmetric {\Sigma}_x');
str = sprintf('Error of estimates of {\\beta_1}^c with %d samples',m);
title(str);
str = sprintf('mean square error of estimated value of {\\beta_1}^c');
ylabel(str);
xlabel('number of initial data points given');
str = sprintf('beta_mse_%dsamples_%dminN_%dmaxN_nosym.png',m,min(allN),max(allN));
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
legend('frequency of overage for diagonal {\Sigma}_x', ...
    'regression line for diagonal {\Sigma}_x',...
    'frequency of overage for non-diagonal {\Sigma}_x',...
    'regression line for non-diagonal {\Sigma}_x',...
    '50%','Location','southeast');
str = sprintf('Percent of {\\beta_1}^c larger than true {\\beta_1}^c with %d samples',m);
title(str);
str = sprintf('percent of estimated {\\beta_1}^c larger than true value');
ylabel(str);
xlabel('number of initial data points given');
str = sprintf('beta_overage_%dsamples_%dminN_%dmaxN.png',m,min(allN),max(allN));
print('-dpng', str);
hold off;

% Histogram of overages
figure; 
h1 = histogram(Beta_data_sym(:,4),20);
h = findobj(gca,'Type','patch');
set(h,'FaceColor','r','EdgeColor','w','facealpha',0.75);
hold on;
line([0.5 0.5], [0 0.2*size(allN,2)]);
h2 = histogram(Beta_data_nosym(:,4),20);
h = findobj(gca,'Type','patch');
set(h,'facealpha',0.75);
ylabel('frequency of overage in estimates');
xlabelstr = sprintf('frequency across %d samples',m);
xlabel(xlabelstr);
legend([h1 h2],{'diagonal {\Sigma}_x','non-diagonal {\Sigma}_x'});
titlestr = sprintf('Frequency of overage in errors of {\\beta_1}^c across n from %d to %d', min(allN), max(allN));
title(titlestr);
str = sprintf('beta_overage_%dsamples_%dminN_%dmaxN.png',m,min(allN),max(allN));
print('-dpng', str);
hold off;

percOver50symSmall = logical(Beta_data_sym(:,4)>0.5);
percOver50symSmall = sum(percOver50symSmall)/size(percOver50symSmall,1);

percOver50nosymSmall = logical(Beta_data_nosym(:,4)>0.5);
percOver50nosymSmall = sum(percOver50nosymSmall)/size(percOver50nosymSmall,1);

%% Generate plots for analysis of error in A and degenerate cases.  
%% Error in A, log plot
tic;

% Log plot 
beta = 100;
% Choose m (number of times to sample the Wishart distribution).
m = 500;
% Choose variety of n values (number of data points available).
allN = logspace(3,6,900);

sigma_X = [1 0.5; 0.5 2];
disp('A error analysis log nodiag\n');
[ MSE_data_nosym, MSE_errdata_nosym, ABS_data_nosym, ABS_errdata_nosym ] = plotAError( sigma_X, ...
    sigma_Y, sigma_XY, beta, m, allN);

[TrueA_nosym, Truebeta_crit_nosym] = gib_optimize(sigma_X,sigma_Y,sigma_XY,beta);

sigma_X = sigma_x*eye(n_x);
disp('A error analysis log diag\n');
[ MSE_data_sym, MSE_errdata_sym, ABS_data_sym, ABS_errdata_sym ] = plotAError( sigma_X, ...
    sigma_Y, sigma_XY, beta, m, allN);

[TrueA_sym, Truebeta_crit_sym] = gib_optimize(sigma_X,sigma_Y,sigma_XY,beta);

ALog_TE = toc;

% Plot how often betacrit is greater than beta.
figure;
h = semilogx(MSE_data_nosym(:,1), MSE_data_nosym(:,end));
hold on;
semilogx(MSE_data_sym(:,1), MSE_data_sym(:,end));
set(get(h,'Parent'), 'XScale', 'log');
str = 'Frequency of {\beta}_{1}^{c} greater than \beta';
title(str);
ylabel('percent of degenerate solutions');
xlabel('number of initial data points given');
legend('non-diagonal {\Sigma}_x','diagonal {\Sigma}_x');
str = sprintf('beta_failure_%dsamples_largeN.png',m);
print('-dpng', str);
hold off;

% Also do separate plots for sym and nosym. 
figure;
h = semilogx(MSE_data_nosym(:,1), MSE_data_nosym(:,end));
set(get(h,'Parent'), 'XScale', 'log');
str = 'Frequency of {\beta}_{1}^{c} greater than \beta for non-diagonal {\Sigma}_x';
title(str);
ylabel('percent of degenerate solutions');
xlabel('number of initial data points given');
str = sprintf('beta_failure_%dsamples_largeN_nodiag.png',m);
print('-dpng', str);
hold off;

figure;
h = semilogx(MSE_data_sym(:,1), MSE_data_sym(:,end));
set(get(h,'Parent'), 'XScale', 'log');
str = 'Frequency of {\beta}_{1}^{c} greater than \beta for diagonal {\Sigma}_x';
title(str);
ylabel('percent of degenerate solutions');
xlabel('number of initial data points given');
str = sprintf('beta_failure_%dsamples_largeN_diag.png',m);
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
    legend('error for diagonal {\Sigma}_x', 'absolute value of true value for diagonal',...
        'error for non-diagonal {\Sigma}_x', 'absolute value of true value for non-diagonal');
    xlim([0 1.05*allN(end)]);
    set(get(h1,'Parent'), 'XScale', 'log');
                
    if(k==1)
        titlestr = sprintf('Error of estimates of %dst element of projection matrix with %d samples',k,m);
    elseif(k==2)
        titlestr = sprintf('Error of estimates of %dnd element of projection matrix with %d samples',k,m);
    elseif(k==3)
        titlestr = sprintf('Error of estimates of %drd element of projection matrix with %d samples',k,m);
    else
        titlestr = sprintf('Error of estimates of %dth element of projection matrix with %d samples',k,m);
    end
    title(titlestr);
    
    if(k==1)
        ystr = sprintf('mean square error of estimated value of %dst element of A', k);
    elseif(k==2)
        ystr = sprintf('mean square error of estimated value of %dnd element of A', k);
    elseif(k==3)
        ystr = sprintf('mean square error of estimated value of %drd element of A', k);
    else
        ystr = sprintf('mean square error of estimated value of %dth element of A', k);
    end
    ylabel(ystr);    
    xlabel('number of initial data points given');
    filestr = sprintf('mse_filt_%dth_%dsamples_%dminN_%dmaxN.png',k,m,min(allN),allN(end));
    print('-dpng', filestr);
    hold off;
    
    % symmetric plot
    figure;
    h1 = errorbar(MSE_data_sym(:,1),MSE_data_sym(:,k+1),MSE_errdata_sym(:,k+1));
    hold on;
    h2 = plot([allN(1) allN(end)],[abs(TrueA_sym(1,k)) abs(TrueA_sym(1,k))]);
    legend('error for diagonal {\Sigma}_x', 'absolute value of true value for diagonal');
    xlim([0 1.05*allN(end)]);
    set(get(h1,'Parent'), 'XScale', 'log');
    title(titlestr);
    ylabel(ystr);
    xlabel('number of initial data points given');
    filestr = sprintf('mse_filt_%dth_%dsamples_%dminN_%dmaxN_diag.png',k,m,min(allN),allN(end));
    print('-dpng', filestr);
    hold off;
    
    % non-symmetric plot
    figure;
    h3 = errorbar(MSE_data_nosym(:,1), MSE_data_nosym(:,k+1),MSE_errdata_nosym(:,k+1));
    hold on;
    h4 = plot([allN(1) allN(end)],[abs(TrueA_nosym(1,k)) abs(TrueA_nosym(1,k))]);
    legend('error for non-symmetric {\Sigma}_x', 'absolute value of true value for non-symmetric');
    xlim([0 1.05*allN(end)]);
    set(get(h1,'Parent'), 'XScale', 'log');
    title(titlestr);
    ylabel(ystr);
    xlabel('number of initial data points given');
    filestr = sprintf('mse_filt_%dth_%dsamples_%dminN_%dmaxN_nosym.png',k,m,min(allN),allN(end));
    print('-dpng', filestr);
    hold off;
end

% Plot percent error in elements of A, log. 
for k = 1:2
    figure;
    h = errorbar(ABS_data_sym(:,1),ABS_data_sym(:,k+1),ABS_errdata_sym(:,k+1));
    hold on;
    h2 = errorbar(ABS_data_nosym(:,1),ABS_data_nosym(:,k+1),ABS_errdata_nosym(:,k+1));
    axis tight;
    legend('error for diagonal {\Sigma}_x', 'error for non-diagonal {\Sigma}_x');
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

tic;

% Pick beta. 
beta = 100;
% Choose m (number of times to sample the Wishart distribution).
m = 500;
% Choose variety of n values (number of data points available).
allN = 10*(100:500);

sigma_X = [1 0.5; 0.5 2];
disp('A linear error no diag\n');
[ MSE_data_nosym, MSE_errdata_nosym, ABS_data_nosym, ABS_errdata_nosym ] = plotAError( sigma_X, ...
    sigma_Y, sigma_XY, beta, m, allN);

[TrueA_nosym, Truebeta_crit_nosym] = gib_optimize(sigma_X,sigma_Y,sigma_XY,beta);

disp('A linear error diag\n');
sigma_X = sigma_x*eye(n_x);
[ MSE_data_sym, MSE_errdata_sym, ABS_data_sym, ABS_errdata_sym ] = plotAError( sigma_X, ...
    sigma_Y, sigma_XY, beta, m, allN);

[TrueA_sym, Truebeta_crit_sym] = gib_optimize(sigma_X,sigma_Y,sigma_XY,beta);

ALinear_TE = toc;

% Plot how often betacrit is greater than beta.
figure;
plot(MSE_data_nosym(:,1), MSE_data_nosym(:,end));
hold on;
plot(MSE_data_sym(:,1), MSE_data_sym(:,end));
legend('non-diagonal {\Sigma}_x', 'diagonal {\Sigma}_x');
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
    legend('error for diagonal {\Sigma}_x', 'absolute value of true value for symmetric',...
        'error for non-diagonal {\Sigma}_x', 'absolute value of true value for non-symmetric');
    xlim([0 1.05*allN(end)]);
                
    if(k==1)
        titlestr = sprintf('Error of estimates of %dst element of projection matrix with %d samples',k,m);
    elseif(k==2)
        titlestr = sprintf('Error of estimates of %dnd element of projection matrix with %d samples',k,m);
    elseif(k==3)
        titlestr = sprintf('Error of estimates of %drd element of projection matrix with %d samples',k,m);
    else
        titlestr = sprintf('Error of estimates of %dth element of projection matrix with %d samples',k,m);
    end
    title(titlestr);
    
    if(k==1)
        ystr = sprintf('mean square error of estimated value of %dst element of A', k);
    elseif(k==2)
        ystr = sprintf('mean square error of estimated value of %dnd element of A', k);
    elseif(k==3)
        ystr = sprintf('mean square error of estimated value of %drd element of A', k);
    else
        ystr = sprintf('mean square error of estimated value of %dth element of A', k);
    end
    ylabel(ystr);
        
    xlabel('number of initial data points given');
    filestr = sprintf('mse_filt_%dth_%dsamples_%dminN_%dmaxN.png',k,m,min(allN),allN(end));
    print('-dpng', filestr);
    hold off;
    
    % make symmetric linear A plot
    figure;
    h1 = errorbar(MSE_data_sym(:,1),MSE_data_sym(:,k+1),MSE_errdata_sym(:,k+1));
    hold on;
    h2 = plot([allN(1) allN(end)],[abs(TrueA_sym(1,k)) abs(TrueA_sym(1,k))]);
    legend('error for symmetric {\Sigma}_x', 'absolute value of true value for symmetric');
    xlim([0 1.05*allN(end)]);
    title(titlestr);
    ylabel(ystr);
    xlabel('number of initial data points given');
    filestr = sprintf('mse_filt_%dth_%dsamples_%dminN_%dmaxN_sym.png',k,m,min(allN),allN(end));
    print('-dpng', filestr);
    
    % non-symmetric plot
    figure;
    h3 = errorbar(MSE_data_nosym(:,1), MSE_data_nosym(:,k+1),MSE_errdata_nosym(:,k+1));
    hold on;
    h4 = plot([allN(1) allN(end)],[abs(TrueA_nosym(1,k)) abs(TrueA_nosym(1,k))]);
    legend('error for non-symmetric {\Sigma}_x', 'absolute value of true value for non-symmetric');
    xlim([0 1.05*allN(end)]);   
    title(titlestr);
    ylabel(ystr);
    xlabel('number of initial data points given');
    filestr = sprintf('mse_filt_%dth_%dsamples_%dminN_%dmaxN_nosym.png',k,m,min(allN),allN(end));
    print('-dpng', filestr);
        
end

% Plot percent error in elements of A, linear. 
for k = 1:2
    figure;
    h = errorbar(ABS_data_sym(:,1),ABS_data_sym(:,k+1),ABS_errdata_sym(:,k+1));
    hold on;
    h2 = errorbar(ABS_data_nosym(:,1),ABS_data_nosym(:,k+1),ABS_errdata_nosym(:,k+1));
    legend('error for diagonal {\Sigma}_x', 'error for non-diagonal {\Sigma}_x');
    ylim([0 max(ABS_data_nosym(:,1))*1.1]);
    
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

%% Generate error convergence plots.

%% symmetric sigma_X
tic;

sigma_X = sigma_x*eye(n_x);
sigma_Y = sigma_y*eye(n_y);
sigma_XY = [0.1; 0.2];

% Choose beta.
beta = 100;
% Choose m (number of times to sample the Wishart distribution).
m = 500;
% Choose variety of n values (number of data points available).
allN = 10*(100:1000);

disp('convergence data');
[ ConvLin_data, ConvQuad_data, FiltConvLin_data1, ...
    FiltConvQuad_data1, FiltConvLin_data2, FiltConvQuad_data2 ] = ...
    plotErrorConv( sigma_X, sigma_Y, sigma_XY, ...
    beta, m, allN);

ConvSym_TE = toc;

% Plot filtered values for linear convergence rates. 
figure;
scatter(FiltConvLin_data1(2:end,1), FiltConvLin_data1(2:end,2));
hold on;
scatter(FiltConvLin_data2(2:end,1), FiltConvLin_data2(2:end,2));
legend('A_{11}','A_{12}');
title('Filtered sequence of linear convergence rates for A_{11} and A_{12} for diagonal {\Sigma}_x');
xlabel('number of data points observed');
ylabel('(linear) convergence ratio');
str = sprintf('conv_filt_lin_%dsamples_sym.png',m);
print('-dpng', str);
hold off;

% Plot filtered values for quadratic convergence rates. 
figure;
scatter(FiltConvQuad_data1(2:end,1), FiltConvQuad_data1(2:end,2));
hold on;
scatter(FiltConvQuad_data2(2:end,1), FiltConvQuad_data2(2:end,2));
legend('A_{11}','A_{12}');
title('Filtered sequence of quadratic convergence rates for A_{11} and A_{12} for diagonal {\Sigma}_x');
xlabel('number of data points observed');
ylabel('(quadratic) convergence ratio');
str = sprintf('conv_filt_quad_%dsamples_sym.png',m);
print('-dpng', str);
hold off;

figure;
scatter(ConvLin_data(2:end,1), ConvLin_data(2:end,2));
hold on;
scatter(ConvLin_data(2:end,1), ConvLin_data(2:end,3));
legend('A_{11}','A_{12}');
title('Sequence of linear convergence rates for A_{11} and A_{12} for diagonal {\Sigma}_x');
xlabel('number of data points observed');
ylabel('(linear) convergence ratio');
str = sprintf('conv_lin_11_%dsamples_sym.png',m);
print('-dpng', str);
hold off;

figure;
scatter(ConvQuad_data(2:end,1), ConvQuad_data(2:end,2));
hold on;
scatter(ConvQuad_data(2:end,1), ConvQuad_data(2:end,3));
title('Sequence of linear convergence rates for A_{11} and A_{12} for diagonal {\Sigma}_x');
xlabel('number of data points observed');
ylabel('(quadratic) convergence ratio');
legend('A_{11}','A_{12}');
str = sprintf('conv_quad_%dsamples_sym.png',m);
print('-dpng', str);
hold off;

% non-symmetric sigma_X, convergence
sigma_X = [1 0.5; 0.5 2];

[ ConvLin_data, ConvQuad_data, FiltConvLin_data1, ...
    FiltConvQuad_data1, FiltConvLin_data2, FiltConvQuad_data2 ] = ...
    plotErrorConv( sigma_X, sigma_Y, sigma_XY, ...
    beta, m, allN);

ConvNoSym_TE = toc;

% Plot filtered values for linear convergence rates. 
figure;
scatter(FiltConvLin_data1(2:end,1), FiltConvLin_data1(2:end,2));
hold on;
scatter(FiltConvLin_data2(2:end,1), FiltConvLin_data2(2:end,2));
legend('A_{11}','A_{12}');
title('Filtered sequence of linear convergence rates for A_{11} and A_{12} for non-diagonal {\Sigma}_x');
xlabel('number of data points observed');
ylabel('(linear) convergence ratio');
str = sprintf('conv_filt_lin_%dsamples_nosym.png',m);
print('-dpng', str);
hold off;

% Plot filtered values for quadratic convergence rates. 
figure;
scatter(FiltConvQuad_data1(2:end,1), FiltConvQuad_data1(2:end,2));
hold on;
scatter(FiltConvQuad_data2(2:end,1), FiltConvQuad_data2(2:end,2));
legend('A_{11}','A_{12}');
title('Filtered sequence of quadratic convergence rates for A_{11} and A_{12} for non-diagonal {\Sigma}_x');
xlabel('number of data points observed');
ylabel('(quadratic) convergence ratio');
str = sprintf('conv_filt_quad_%dsamples_nosym.png',m);
print('-dpng', str);
hold off;

figure;
scatter(ConvLin_data(2:end,1), ConvLin_data(2:end,2));
hold on;
scatter(ConvLin_data(2:end,1), ConvLin_data(2:end,3));
legend('A_{11}','A_{12}');
title('Sequence of linear convergence rates for A_{11} and A_{12} for non-diagonal {\Sigma}_x');
xlabel('number of data points observed');
ylabel('(linear) convergence ratio');
str = sprintf('conv_lin_11_%dsamples_nosym.png',m);
print('-dpng', str);
hold off;

figure;
scatter(ConvQuad_data(2:end,1), ConvQuad_data(2:end,2));
hold on;
scatter(ConvQuad_data(2:end,1), ConvQuad_data(2:end,3));
title('Sequence of linear convergence rates for A_{11} and A_{12} for non-diagonal {\Sigma}_x');
xlabel('number of data points observed');
ylabel('(quadratic) convergence ratio');
legend('A_{11}','A_{12}');
str = sprintf('conv_quad_%dsamples_nosym.png',m);
print('-dpng', str);
hold off;

%% Time elapsed
fileID = fopen('timing.txt','w');
fprintf(fileID,'information curve: %f seconds\n', IC_TE);
fprintf(fileID,'beta error on log scale: %f seconds\n',BetaLog_TE);
fprintf(fileID,'beta error on linear scale: %f seconds\n',BetaLinear_TE);
fprintf(fileID,'A error on log scale: %f seconds\n', ALog_TE);
fprintf(fileID,'A error on linear scale: %f seconds\n', ALinear_TE);
fprintf(fileID,'Symmetric convergence: %f seconds\n', ConvSym_TE);
fprintf(fileID,'Asymmetric convergence: %f seconds\n', ConvNoSym_TE);
fprintf(fileID,'percent over 50 in symmetric, big: %f\n', percOver50symBig);
fprintf(fileID,'percent over 50 in non-symmetric, big: %f\n', percOver50nosymBig);
fprintf(fileID,'percent over 50 in symmetric, small: %f\n', percOver50symSmall);
fprintf(fileID,'percent over 50 in non-symmetric, small: %f\n', percOver50nosymSmall);
fclose(fileID);