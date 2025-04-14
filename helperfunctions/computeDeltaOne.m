function output = computeDeltaOne(etaChar,dcChar, dcProb, randomSwapFlag)
%computes the value of deltaone.

dcMin = dcProb*(1-dcChar);
dcMax = dcProb*(1+dcChar);

etaRatio = (1-etaChar)/(1+etaChar);


if randomSwapFlag == 0

    zeroPhotonPart = (1 - (1-(1-dcMin)^2) / (1-(1-dcMax)^2) ) * ( (dcMax*(2-dcMin)) / (1-(1-dcMin)^2));
    OtherPart = 4*(1 - sqrt(1 - (1-dcMin)^2 *(1-etaRatio)));
    
    
    output = max(zeroPhotonPart, OtherPart);
    return;

elseif randomSwapFlag == 1
    output =  4*(1 - sqrt(1 - (1-dcMin)^2 *(1-etaRatio)^2 / 2));
    return;
else
    warning("random swap flag not set correctly!!!");
    return;

end
