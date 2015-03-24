function [ Data ] = sample_wishart( sigma_X, sigma_Y, sigma_XY, beta, ...
    numDataPointsObserved, numTrials )
%sample_wishart: Calculates resulting projection matrix for each 
% sample of covariance matrix. 

n_x = size(sigma_X,1);
n_y = size(sigma_Y,1);

% Choose n (number of data points observed).
n = numDataPointsObserved;
% Choose m (number of times to sample the Wishart distribution).
m = numTrials;

Sigma = [sigma_X, sigma_XY; sigma_XY.', sigma_Y];
[TrueA, Truebeta_crit] = gib_optimize(sigma_X,sigma_Y,sigma_XY,beta);
Data = cell(m+1, 3);
Data{1,1} = Sigma;
Data{1,2} = TrueA;
Data{1,3} = Truebeta_crit;

% Calculate the output matrix A for each sampled covariance matrix. 
for i = 1:m
    W = wishrnd(Sigma, n)/n;
    [A, beta_crit] = gib_optimize(W(1:n_x,1:n_x),...
        W(n_x+1:n_x+n_y,n_x+1:n_x+n_y), W(1:n_x,n_x+1:n_x+n_y),...
        beta);
    if(sign(A(1,1))~=sign(TrueA(1,1)))
        A = -1*A;
    end
    Data{i+1,1} = W;
    Data{i+1,2} = A;
    Data{i+1,3} = beta_crit;
end

% Plot different values of W and A 
allW = cell2mat(Data(:,1));
allA = cell2mat(Data(:,2));

end

