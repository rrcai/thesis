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
allN = 10*(50:1000);

% Set up array to store values. 
ConvLin_data = zeros(size(allN,2)-1, (n_x)+1);
ConvLin_data(:,1) = allN(1:end-1);
ConvQuad_data = zeros(size(allN,2)-1, (n_x)+1);
ConvQuad_data(:,1) = allN(1:end-1);
EstA_data = zeros(size(allN,2), n_x+1);

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
    for k = 1:2
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
        
        EstA_data(n,k+1) = sum(thisAData(indices))/sum(indices);
   end
end

ConvLin_data(:,2) = (EstA_data(2:end,2)-TrueA(1,1))./(EstA_data(1:end-1,2)-TrueA(1,1));
ConvLin_data(:,3) = (EstA_data(2:end,3)-TrueA(1,2))./(EstA_data(1:end-1,3)-TrueA(1,2));
ConvQuad_data(:,2) = (EstA_data(2:end,2)-TrueA(1,1))./((EstA_data(1:end-1,2)-TrueA(1,1)).^2);
ConvQuad_data(:,3) = (EstA_data(2:end,3)-TrueA(1,2))./((EstA_data(1:end-1,3)-TrueA(1,2)).^2);
scatter(ConvLin_data(2:5:end,1), ConvLin_data(2:5:end,2));
hold on;
h = lsline;
scatter(ConvLin_data(2:5:end,1), ConvLin_data(2:5:end,3));
% plot([min(ConvLin_data(:,1)) max(ConvLin_data(:,1))],[TrueA(1,1) TrueA(1,1)]);
legend('scatterplot of q-ratios for A_{11}','scatterplot of q-ratios for A_{12}');
title('Sequence of Q-rates of convergence for A_{11} and A_{12}');
xlabel('number of data points observed');
ylabel('(linear) convergence ratio');
str = sprintf('conv_lin_11_%dsamples_sym.png',m);
print('-dpng', str);
hold off;

% figure;
% scatter(ConvLin_data(2:5:end,1), ConvLin_data(2:5:end,3));
% hold on;
% h = lsline;
% set(h,'LineWidth',3);
% % plot([min(ConvLin_data(:,1)) max(ConvLin_data(:,1))],[TrueA(1,2) TrueA(1,2)]);
% legend('scatterplot of q-ratios','linear regression','true value of A_{12}');
% % repn = repmat(MSE_data(2:10:end,1),2,1);
% % repconv = [MSE_conv(1:10:end,1); MSE_conv(1:10:end,2)];
% % k1 = convhull(MSE_data(2:end,1),MSE_conv(1:end,1));
% % k2 = convhull(MSE_data(2:end,1),MSE_conv(1:end,2));
% % plot(MSE_data(k1(2:end/2),1),MSE_conv(k1(2:end/2),1),'r-');
% % plot(MSE_data(k2(2:end/2),1),MSE_conv(k2(2:end/2),2),'g-');
% title('Sequence of linear Q-rates of convergence for A_{12}');
% xlabel('number of data points observed');
% ylabel('(linear) convergence ratio');
% str = sprintf('conv_lin_12_%dsamples_sym.png',m);
% print('-dpng', str);

figure;
scatter(ConvQuad_data(2:5:end,1), ConvQuad_data(2:5:end,2));
hold on;
% plot([min(ConvLin_data(:,1)) max(ConvLin_data(:,1))],[TrueA(1,1) TrueA(1,1)]);
scatter(ConvQuad_data(2:5:end,1), ConvQuad_data(2:5:end,3));
title('Sequence of quadratic Q-rate of convergence for A_{11} and A_{12}');
xlabel('number of data points observed');
ylabel('(quadratic) convergence ratio');
legend('A_{11}','A_{12}');
str = sprintf('conv_quad_%dsamples_sym.png',m);
print('-dpng', str);
hold off;