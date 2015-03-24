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
allN = 5*(5:100);

% Set up array to store values. 
MSE_data = zeros(size(allN,2), (n_x)^2+1);
MSE_data(:,1) = allN;
MSE_errdata = zeros(size(allN,2), (n_x)^2+1);
MSE_errdata(:,1) = allN;

[TrueA, Truebeta_crit] = gib_optimize(sigma_X,sigma_Y,sigma_XY,beta);

for n = 1:size(allN,2);
    Data = sample_wishart(sigma_X, sigma_Y, sigma_XY, beta, allN(n), m);

    % Plot different values of W, A, and beta_crit
    allW = cell2mat(Data(:,1));
    allA = cell2mat(Data(:,2));
    
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
        
        % Plot this information in a histogram
        % OR
        % Calculate mean square error and error bars
        % this_mse = (expVal-trueVal).^2;
        this_abs = abs(expVal-trueVal)./abs(trueVal);
        MSE_data(n,k+1) = sum(this_abs)/m;
        MSE_errdata(n,k+1) = std(this_abs)/sqrt(m);
    end
end

% Plot this data. 
for k = 1:2
    figure;
    errorbar(MSE_data(:,1),MSE_data(:,k+1),MSE_errdata(:,k+1));
    axis([0 1.05*allN(end) 0 1.05*max(MSE_data(:,k+1))]);
    hold on;
    % plot([allN(1) allN(end)],[TrueA(1,k) TrueA(1,k)]);
    
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
    str = sprintf('perc_%dth_%dsamples.png',k,m);
    print('-dpng', str);
    
    hold off;
end


