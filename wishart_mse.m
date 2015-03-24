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
m = 1000;
% Choose variety of n values (number of data points available).
allN = 5*(5:1000);

% Set up array to store values. 
MSE_data = zeros(size(allN,2), (n_x)^2+2);
MSE_data(:,1) = allN;
MSE_errdata = zeros(size(allN,2), (n_x)^2+1);
MSE_errdata(:,1) = allN;

ABS_data = zeros(size(allN,2), (n_x)^2+2);
ABS_data(:,1) = allN;
ABS_errdata = zeros(size(allN,2), (n_x)^2+1);
ABS_errdata(:,1) = allN;

[TrueA, Truebeta_crit] = gib_optimize(sigma_X,sigma_Y,sigma_XY,beta);
trueVal = Truebeta_crit(1,1);

for n = 1:size(allN,2);
    Data = sample_wishart(sigma_X, sigma_Y, sigma_XY, beta, allN(n), m);

    % Find values of W, A, and beta. 
    allW = cell2mat(Data(:,1));
    allA = cell2mat(Data(:,2));
    allBeta_crit = cell2mat(Data(:,3));
    
    indices = (((1:(size(allBeta_crit,1)/2)).*2)-1)';
    thisBetacritData = allBeta_crit(indices,1);
    trueVal = Truebeta_crit(1,1);
    
    % A will always be n_x by n_x. 
    for k = 1:4
        % For this element of A and this value of n, find the mean 
        % square error. 
        row = -rem(ceil(k/2),2)+2*(1:(m+1));
        col = rem(k,2);
        if (col == 0)
            col = 2;
        end

        thisAData = allA(row,col);
        thisWData = allW(row,col);
        % Remove the first row, since it is the true value. 
        trueVal = thisAData(1);
        expVal = thisAData(2:end);
        
        indices = logical(thisBetacritData < beta);
        
        % How often does this happen? Store this data in another column. 
        MSE_data(n,end) = 1-(sum(indices(2:end))/(size(indices,1)-1));
        ABS_data(n,end) = 1-(sum(indices(2:end))/(size(indices,1)-1));
        
        this_mse = (thisAData(indices)-trueVal).^2;
        MSE_data(n,k+1) = sum(this_mse(2:end))/size(this_mse,1);
        MSE_errdata(n,k+1) = std(this_mse(2:end))/sqrt(size(this_mse,1));
        
        this_mse = abs(expVal-trueVal)./abs(trueVal);
        ABS_data(n,k+1) = sum(this_mse(2:end))/size(this_mse,1);
        ABS_errdata(n,k+1) = std(this_mse(2:end))/sqrt(size(this_mse,1));
    end
end

% Plot how often beta is less than beta crit.
% figure;
% plot(MSE_data(:,1), MSE_data(:,end));
% 
% str = 'Frequency of {\beta}_{1}^{c} less than \beta';
% title(str);
% 
% ylabel('percent of degenerate solutions');
% xlabel('number of initial data points given');
% str = sprintf('beta_failure_%dsamples.png',m);
% print('-dpng', str);
        
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
    str = sprintf('mse_filtered_%dth_%dsamples_%dmaxN.png',k,m,allN(end));
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
    str = sprintf('perc_filtered_%dth_%dsamples_%dmaxN.png',k,m,allN(end));
    print('-dpng', str);
end

