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
m = 500;
% Choose variety of n values (number of data points available).
allN = 5*(5:1000);

% Set up array to store values. 
Conv_data = zeros(size(allN,2), (n_x)^2+1);
Conv_data(:,1) = allN;
Conv_errdata = zeros(size(allN,2), (n_x)^2+1);
Conv_errdata(:,1) = allN;

[TrueA, Truebeta_crit] = gib_optimize(sigma_X,sigma_Y,sigma_XY,beta);

for n = 1:size(allN,2);
    Data = sample_wishart(sigma_X, sigma_Y, sigma_XY, beta, allN(n), m);

    % Plot different values of W, A, and beta_crit
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
        
        this_mse = (thisAData(indices)-trueVal).^2;
        Conv_data(n,k+1) = sum(this_mse(2:end))/size(this_mse,1);
        Conv_errdata(n,k+1) = std(this_mse(2:end))/sqrt(size(this_mse,1));
    end
end

MSE_conv = Conv_data(2:end,2:3)./(Conv_data(1:end-1,2:3));
scatter(Conv_data(2:5:end,1), MSE_conv(1:5:end,1));
lsline;
title('Sequence of Q-rates of convergence for A_{11}');
xlabel('number of data points observed');
ylabel('(linear) convergence ratio');
str = sprintf('conv_lin_11_%dsamples_sym.png',m);
print('-dpng', str);

figure;
scatter(Conv_data(2:5:end,1), MSE_conv(1:5:end,2));
lsline;
% repn = repmat(MSE_data(2:10:end,1),2,1);
% repconv = [MSE_conv(1:10:end,1); MSE_conv(1:10:end,2)];
% k1 = convhull(MSE_data(2:end,1),MSE_conv(1:end,1));
% k2 = convhull(MSE_data(2:end,1),MSE_conv(1:end,2));
% plot(MSE_data(k1(2:end/2),1),MSE_conv(k1(2:end/2),1),'r-');
% plot(MSE_data(k2(2:end/2),1),MSE_conv(k2(2:end/2),2),'g-');
title('Sequence of Q-rates of convergence for A_{12}');
xlabel('number of data points observed');
ylabel('(linear) convergence ratio');
str = sprintf('conv_lin_12_%dsamples_sym.png',m);
print('-dpng', str);

figure;
MSE_conv2 = Conv_data(2:end,2:3)./((Conv_data(1:end-1,2:3)).^2);
scatter(Conv_data(2:5:end,1), MSE_conv2(1:5:end,1));
hold on;
scatter(Conv_data(2:5:end,1), MSE_conv2(1:5:end,2));
title('Sequence of Q-rates of convergence');
xlabel('number of data points observed');
ylabel('(quadratic) convergence ratio');
legend('A_{11}','A_{12}');
str = sprintf('conv_quad_%dsamples_sym.png',m);
print('-dpng', str);