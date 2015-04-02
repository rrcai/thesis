clear all;
close all;

n_x = 2;
n_y = 1;
sigma_X = [1 0.5; 0.5 2];
sigma_y = 1;
sigma_Y = sigma_y*eye(n_y);
sigma_XY = [0.1; 0.2];

% Choose beta.
allBeta = 10*(0.1:1000);
% Choose m (number of times to sample the Wishart distribution).
m = 500;
% Choose one n value (number of data points available).
n = 10000;

% Set up array to store values. 
% Column 1: beta
% Column 2: true past info
% Column 3: true future info
% Column 4: mean past info
% Column 5: mean future info
IC_data = zeros(size(allBeta,2), 5);
IC_data(:,1) = allBeta;

% Column 1: beta
% Column 2: stderr of past info
% Column 3: stderr of future info
IC_errdata = zeros(size(allBeta,2), 3);
IC_errdata(:,1) = allBeta;

MSE_data = zeros(size(allBeta,2), (n_x)^2+2);
MSE_data(:,1) = allBeta;
MSE_errdata = zeros(size(allBeta,2), (n_x)^2+1);
MSE_errdata(:,1) = allBeta;

ABS_data = zeros(size(allBeta,2), (n_x)^2+2);
ABS_data(:,1) = allBeta;
ABS_errdata = zeros(size(allBeta,2), (n_x)^2+1);
ABS_errdata(:,1) = allBeta;

for beta = 1:size(allBeta,2);
    Data = sample_wishart(sigma_X, sigma_Y, sigma_XY, allBeta(beta), n, m);

    % Find values of W, A, and beta crit.
    allW = cell2mat(Data(:,1));
    allA = cell2mat(Data(:,2));
    allBeta_crit = cell2mat(Data(:,3));
    
    indices = (((1:(size(allBeta_crit,1)/2)).*2)-1)';
    indices2 = (2*(1:size(allBeta_crit,1)/2))';
    thisBetacritData = [allBeta_crit(indices,1) allBeta_crit(indices2,2)];
    trueBetacritVals = thisBetacritData(1,:);
    trueEvals = 1-1./trueBetacritVals;
    expBetacritVals = thisBetacritData(2:end,:);
    expEvals = 1-1./expBetacritVals;
    
    % eigenvalues
    findIndex = @(x) min(abs(x-beta));
    [min_diff, nBeta] = cellfun(findIndex, num2cell(expBetacritVals,2), 'UniformOutput', false);
    nBeta = cell2mat(nBeta);
    testlogical = beta < expBetacritVals(nBeta,2);
    nBeta(testlogical) = nBeta(testlogical) - 1;

    % true past info eigenvalues
    [min_diff, trueNBeta] = findIndex(trueBetacritVals);
    if beta < trueBetacritVals(trueNBeta,2)
        trueNBeta = trueNBeta - 1;
    end

    % true past and future info values
    if (trueNBeta == 0)
        IC_data(beta,2) = 0;
        IC_data(beta,3) = IC_data(beta,2)-0;
    elseif (trueNBeta == 1)
        IC_data(beta,2) = (1/2)*log((beta-1)*(1-trueEvals(trueNBeta,1))/trueEvals(trueNBeta,1));
        IC_data(beta,3) = IC_data(beta,2)-(1/2)*log(beta*(1-trueEvals(trueNBeta,1)));
    elseif (trueNBeta == 2)
        IC_data(beta,2) = (1/2)*(log((beta-1)*(1-trueEvals(trueNBeta,1))/trueEvals(trueNBeta,1))...
            +log((beta-1)*(1-trueEvals(1,1))/trueEvals(1,1)));
        IC_data(beta,3) = IC_data(beta,2)-(1/2)*(log(beta*(1-trueEvals(trueNBeta,1)))...
            +log(beta*(1-trueEvals(1,1))));
    end
    
    % mean past and future info
    findInfo = @(x) infoCalculator(x, beta, expEvals);
    [expPastInfo, expFutureInfo] = cellfun(findInfo, num2cell(nBeta,2), 'UniformOutput', false);
    
    IC_data(beta, 4) = mean(expPastInfo);
    IC_data(beta, 5) = mean(expFutureInfo);
    
    IC_errdata(beta,2) = std(expPastInfo)/sqrt(size(expPastInfo,1));
    IC_errdata(beta,3) = std(expFutureInfo)/sqrt(size(expFutureInfo,1));
    
end

% Plot how often beta is less than beta crit.
figure;
plot(IC_data(:,1), MSE_data(:,end));

str = 'Frequency of {\beta}_{1}^{c} less than \beta';
title(str);

ylabel('percent of degenerate solutions');
xlabel('number of initial data points given');
str = sprintf('beta_failure_%dsamples_nosym_%dmaxN.png',m,allN(end));
print('-dpng', str);
        
% Plot this data. 
for k = 1:2
    figure;
%     [ax, h1, h2] = plotyy(MSE_data(:,1), MSE_data(:,k+1), MSE_data(:,1), MSE_data(:,end), ...
%         @(x,y) errorbar(x,y,MSE_errdata(:,k+1)), @plot);
    errorbar(MSE_data(:,1),MSE_data(:,k+1),MSE_errdata(:,k+1));
    hold on;
    plot([allN(1) allN(end)],[abs(TrueA(1,k)) abs(TrueA(1,k))]);
%     set(ax(2), 'YLim', [0 0.3]);
    legend('mean square error', 'absolute value of true value');
    xlim([0 1.05*allN(end)]);
                
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
%     ylabel(ax(2),'percent of degenerate solutions');
        
    xlabel('number of initial data points given');
    str = sprintf('mse_filtered_%dth_%dsamples_%dmaxN_nosym.png',k,m,allN(end));
    print('-dpng', str);
    
    hold off;
end

% Plot this data. 
for k = 1:2
    figure;
%     [ax, h1, h2] = plotyy(ABS_data(:,1), ABS_data(:,k+1), ABS_data(:,1), ABS_data(:,end), ...
%         @(x,y) errorbar(x,y,ABS_errdata(:,k+1)), @plot);
    errorbar(ABS_data(:,1),ABS_data(:,k+1),ABS_errdata(:,k+1));
%     set(ax(2), 'YLim', [0 0.3]);
    xlim([0 1.05*allN(end)]);

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
%     ylabel(ax(2),'percent of degenerate solutions');
    
    xlabel('number of initial data points given');
    str = sprintf('perc_filtered_%dth_%dsamples_%dmaxN_nosym.png',k,m,allN(end));
    print('-dpng', str);
end