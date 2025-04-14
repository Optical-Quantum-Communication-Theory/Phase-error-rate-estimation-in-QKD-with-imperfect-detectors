function c = gammaBin(n,delta,logEpsilon)
% Computes c such that
% \sum_{k = nc+ndelta)^n (n choose k) delta^k (1-delta)^(n-k) \leq
% epsilon^2

% Basically c such that prob of getting greater than delta+c in the
% binomial distribution(n,delta) is smaller than epsilon^2

    
 arguments
        n (1,1) double {mustBePositive}
        delta (1,1) double {mustBeNonnegative}
        logEpsilon (1,1) double {mustBeNegative}
    end
    n = floor(n); % we need n to be an integer. 
    epsilon = 2^(logEpsilon);

    if epsilon^2 > 1e-16 % no-finite precision issues.
    
        % Define the function of ce to find its root
        func = @(c) binocdf(ceil(n*c + n*delta) - 1, n, delta, 'upper') - epsilon^2;
        
        % Initial guess for c
        cInitialGuess = 0;
        
        % Use fzero to find the root of the function
        options = optimset('TolX',1e-9); % Set tolerance for more accurate result
        c = fzero(func, cInitialGuess, options);
        
    
        % We want to make sure that the c is valid.
        StepSize = 1e-9;
    
        while binocdf(ceil(n*c+n*delta)-1,n,delta,'upper') > epsilon^2
            c = c+StepSize;
        end
    
    
        % Ensure c is within the valid range [0, 1]
        c = min(max(c, 0), 1);
    
    else % we use hoeffdings which can directly work with logepsilon

        c = sqrt( - log(2) *logEpsilon / n  );
    end


    

end



