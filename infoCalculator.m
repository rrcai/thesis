function [ expPastInfo, expFutureInfo ] = infoCalculator( nBeta, beta, ...
    expEvals )
% Calculates the past and future information values based on GIB. 
    
    if (nBeta == 0)
        expPastInfo = 0;
        expFutureInfo = expPastInfo-0;
    elseif (nBeta == 1)
        expPastInfo = (1/2)*log((beta-1)*(1-expEvals(nBeta,1))/expEvals(nBeta,1));
        expFutureInfo = expPastInfo-(1/2)*log(beta*(1-expEvals(nBeta,1)));
    elseif (nBeta == 2)
        expPastInfo = (1/2)*(log((beta-1)*(1-expEvals(nBeta,1))/expEvals(nBeta,1))...
            +log((beta-1)*(1-expEvals(1,1))/expEvals(1,1)));
        expFutureInfo = expPastInfo-(1/2)*(log(beta*(1-expEvals(nBeta,1)))...
            +log(beta*(1-expEvals(1,1))));
    end

end

