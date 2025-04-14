function output = gammaSerf(nX,nK, logEpsilon)
%Computes the finite-size correction term due to Serfling. 

arguments        
        nX (1,1) double {mustBePositive}
        nK (1,1) double {mustBePositive}
        logEpsilon (1,1) double {mustBeNegative}
    end

    % Define the f_serf function as a local function
    function f = fSerf(m, n)
        f = (n * m^2) / ((n + m) * (m + 1));
    end

    % Compute f_serf for nX and nK
    fSerfValue = fSerf(nX, nK);
    
    % Compute gamma_serf using the f_serf value and epsilon
    output = sqrt( -log(2)*logEpsilon /  fSerfValue)  ;
    %We need natural log of Epsilon. We are supplied with log base 2. 
end



