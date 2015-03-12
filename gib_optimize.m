function [ A, beta_crit ] = gib_optimize( sigma_X, sigma_Y, sigma_XY, beta )
% gib_optimize: find projection matrix A for X
%   Implements Chechik's algorithm to find projection matrix A that 
%   compresses X into T
% Author: Regina Cai

% index of size doesn't matter because sigma_X is square.
n_x = size(sigma_X,1);

% Find left eigenvectors of the relevant matrix. 
sigma_XcondY = sigma_X-sigma_XY*inv(sigma_Y)*sigma_XY';
B = sigma_XcondY/sigma_X;
% Requires MATLAB R2014b
[V,eval,evec] = eig(B);

% Find critical beta values
[sortEvec, sortEval] = sortem(evec,eval);
beta_crit = 1./(1.-sortEval);

% Find between which beta critical values beta falls
[min_diff, array_position] = min(abs(diag(beta_crit)-beta));
% If beta is larger than largest critical beta, use all eigenvalues and
% eigenvectors
if (beta > beta_crit(array_position,array_position))
    max_index = array_position;
else
    max_index = array_position-1;
end

% Construct A
A = zeros(n_x,n_x);
for i=1:max_index
    % calculate r
    r = sortEvec(:,i)'*sigma_X*sortEvec(:,i);
    % calculate alpha
    alpha = sqrt((beta*(1-sortEval(i,i))-1)/sortEval(i,i)*r);
    A(:,i) = alpha.*sortEvec(:,i);
end
A = A';

end

