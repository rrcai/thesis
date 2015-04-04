function [ MSE_data, MSE_errdata ] = plotBetaError( sigma_X, sigma_Y, sigma_XY, beta, m, allN )
%plotBetaError finds the error in beta critical values for different 
% degrees of freedom in the Wishart distribution. 

% Set up array to store values.
% Column 1: value of n
% Column 2: percent error
% Column 3: mean square error
% Column 4: overage/underage
MSE_data = zeros(size(allN,2),4);
MSE_data(:,1) = allN;
MSE_errdata = zeros(size(allN,2),3);
MSE_errdata(:,1) = allN;

[TrueA, Truebeta_crit] = gib_optimize(sigma_X,sigma_Y,sigma_XY,beta);
trueVal = Truebeta_crit(1,1);

for n = 1:size(allN,2);
    Data = sample_wishart(sigma_X, sigma_Y, sigma_XY, beta, allN(n), m);
    
    % Plot different values of W, A, and beta_crit
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
end

end

