clear all;
close all;

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

for betaIndex = 1:size(allBeta,2);
    
    beta = allBeta(betaIndex);
    Data = sample_wishart(sigma_X, sigma_Y, sigma_XY, beta, n, m);

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
    findIndex = @(x) min(abs(x-allBeta(betaIndex)));
    [min_diff, nBeta] = cellfun(findIndex, num2cell(expBetacritVals,2), 'UniformOutput', false);
    nBeta = cell2mat(nBeta);
    testidx = sub2ind(size(expBetacritVals), 1:size(expBetacritVals,1), nBeta');
    testlogical = beta < expBetacritVals(testidx);
    nBeta(testlogical) = nBeta(testlogical) - 1;

    % true past info eigenvalues
    [min_diff, trueNBeta] = findIndex(trueBetacritVals);
    if beta < trueBetacritVals(1,trueNBeta);
        trueNBeta = trueNBeta - 1;
    end

    % true past and future info values
    if (trueNBeta == 0)
        IC_data(betaIndex,2) = 0;
        IC_data(betaIndex,3) = IC_data(betaIndex,2)-0;
    elseif (trueNBeta == 1)
        IC_data(betaIndex,2) = (1/2)*log((beta-1)*(1-trueEvals(trueNBeta,1))/trueEvals(trueNBeta,1));
        IC_data(betaIndex,3) = IC_data(betaIndex,2)-(1/2)*log(beta*(1-trueEvals(trueNBeta,1)));
    elseif (trueNBeta == 2)
        IC_data(betaIndex,2) = (1/2)*(log((beta-1)*(1-trueEvals(trueNBeta,1))/trueEvals(trueNBeta,1))...
            +log((beta-1)*(1-trueEvals(1,1))/trueEvals(1,1)));
        IC_data(betaIndex,3) = IC_data(betaIndex,2)-(1/2)*(log(beta*(1-trueEvals(trueNBeta,1)))...
            +log(beta*(1-trueEvals(1,1))));
    end
    
    % mean past and future info
    findInfo = @(x) infoCalculator(x, beta, expEvals);
    [expPastInfo, expFutureInfo] = cellfun(findInfo, num2cell(nBeta,2), 'UniformOutput', false);
    expPastInfo = cell2mat(expPastInfo);
    expFutureInfo = cell2mat(expFutureInfo);
    
    IC_data(betaIndex, 4) = mean(expPastInfo);
    IC_data(betaIndex, 5) = mean(expFutureInfo);
    
    IC_errdata(betaIndex,2) = std(expPastInfo)/sqrt(size(expPastInfo,1));
    IC_errdata(betaIndex,3) = std(expFutureInfo)/sqrt(size(expFutureInfo,1));
    
end

% Plot true information curve. 
figure;
plot(IC_data(:,2), IC_data(:,3));
hold on;
errorbar(IC_data(:,4),IC_data(:,5),IC_errdata(:,3),'.k');
herrorbar(IC_data(:,4),IC_data(:,5),IC_errdata(:,2),'.k');
legend('true information curve', 'estimated information curve');
str = 'Information curve error bars';
title(str);
ylabel('I(T;Y)');
xlabel('I(X;T)');
str = sprintf('infoCurve_%dn_%dm_nosym.png',n,m);
print('-dpng', str);

hold off;