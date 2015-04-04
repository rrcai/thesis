function [ MSE_data, MSE_errdata, ABS_data, ABS_errdata ] = plotAError( sigma_X, ...
    sigma_Y, sigma_XY, beta, m, allN)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% Calculates mean square error and mean absolute error in elements of A 
% as a function of the
% number of data points available in prior belief of covariance matrices.
% matrix. 

n_x = size(sigma_X, 2);

% Set up arrays to store values. 
MSE_data = zeros(size(allN,2), (n_x)^2+2);
MSE_data(:,1) = allN;
MSE_errdata = zeros(size(allN,2), (n_x)^2+1);
MSE_errdata(:,1) = allN;

ABS_data = zeros(size(allN,2), (n_x)^2+2);
ABS_data(:,1) = allN;
ABS_errdata = zeros(size(allN,2), (n_x)^2+1);
ABS_errdata(:,1) = allN;

for n = 1:size(allN,2);
    Data = sample_wishart(sigma_X, sigma_Y, sigma_XY, beta, allN(n), m);

    % Find values of W, A, and beta. 
    allA = cell2mat(Data(:,2));
    allBeta_crit = cell2mat(Data(:,3));
    
    indices = (((1:(size(allBeta_crit,1)/2)).*2)-1)';
    thisBetacritData = allBeta_crit(indices,1);
    
    indices = logical(thisBetacritData < beta);
    % How often does this happen? Store this data in another column. 
    MSE_data(n,end) = 1-(sum(indices(2:end))/(size(indices,1)-1));
    
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
        % Remove the first row, since it is the true value. 
        trueVal = thisAData(1);
        expVal = thisAData(2:end);
        
        this_mse = (thisAData(indices)-trueVal).^2;
        MSE_data(n,k+1) = sum(this_mse(2:end))/size(this_mse,1);
        MSE_errdata(n,k+1) = std(this_mse(2:end))/sqrt(size(this_mse,1));

        this_mse = abs(expVal-trueVal)./abs(trueVal);
        ABS_data(n,k+1) = sum(this_mse(2:end))/size(this_mse,1);
        ABS_errdata(n,k+1) = std(this_mse(2:end))/sqrt(size(this_mse,1));
    end

end