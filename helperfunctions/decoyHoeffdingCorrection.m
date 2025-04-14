function output = decoyHoeffdingCorrection(n , logEpsilon)
%Computes the finite-size correction term due to Hoeffding
arguments        
        n (1,1) double {mustBeNonnegative}
        logEpsilon (1,1) double {mustBeNegative}
    end

    output = sqrt( (n/2)*(log(2) - 2*log(2)*logEpsilon)  );
    %We need natural log of Epsilon. We are supplied with log base 2. 
end



