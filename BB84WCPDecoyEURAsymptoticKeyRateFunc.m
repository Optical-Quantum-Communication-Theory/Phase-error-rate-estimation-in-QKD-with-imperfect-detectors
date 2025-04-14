function [keyRate, modParser] = BB84WCPDecoyEURAsymptoticKeyRateFunc(params,options,mathSolverFunc,debugInfo)
% BB84WCPDecoyEURAsymptoticKeyRateFunc A key rate function in the asymptotic regime for a WCP BB84 protocol
% using decoy states using the EUR proof technique for imperfect detectors.
%
%
% Input parameters:
% * fEC: error correction effiency. If set to 1, we are correcting at the
%   Shannon limit. 
% * expectationsConditional: The conditional expectations (as an array)
%   from Alice and Bob's measurements that line up with it's corresponding
%   observable in observablesJoint. These values should be between 0 and 1,
% * DeltaLoss : Parameter determining error in characterizing of the loss in
%   detectors
% * DeltaDarkCounts : Parameter determining error in characterization of
%   dark counts in detectors.
% Outputs:
% * keyrate: Key rate of the QKD protocol.
% Options:
% * verboseLevel: (global option) See makeGlobalOptionsParser for details.
% * decoyTolerance (1e-14): Tolerance on decoy analysis linear program.
%   must be greater than or equal to 0.
% * decoySolver ("SDPT3"): Solver to use for decoy analysis.
% * decoyPrecision ("high"): CVX precision to use for decoy analysis.
% 
% * decoyPhotonCutOff (10): Photon number cut off for decoy analysis. Must
%   be a positive integer.
% DebugInfo:
% 
arguments
    params (1,1) struct
    options (1,1) struct
    mathSolverFunc (1,1) function_handle
    debugInfo (1,1) DebugInfo
end





%% modParser
modParser = moduleParser(mfilename);

modParser.addRequiredParam("pT",@(x) mustBeInRange(x,0,1));

modParser.addRequiredParam("expectationsConditional",@(x) eachRowMustBeAProbDist(x));
modParser.addRequiredParam("decoys",@(x) allCells(x,@(y) y>=0));
modParser.addRequiredParam("decoyProbs",@mustBeProbDistCell); %,must sum to 1. 
% modParser.addAdditionalConstraint(@isEqualSize,["observablesJoint","expectationsConditional"]);
modParser.addRequiredParam("probSignalsA",@mustBeProbDist);

modParser.addRequiredParam("fEC", @(x) mustBeGreaterThanOrEqual(x,1));

modParser.addRequiredParam("DeltaLoss", @(x) x>=0 & x<=1 );
modParser.addRequiredParam("DeltaDarkCounts", @(x) x>=0 & x<=1 );

modParser.addRequiredParam("etaDet", @(x) x>=0 & x<=1 );
modParser.addRequiredParam("dcRate", @(x) x>=0 & x<=1 );

%finite Size stuff is not necessary;

modParser.parse(params);

params = modParser.Results;



%% simple setup (Not needed since we do not do numerical optimization)
%debugMathSolver = debugInfo.addLeaves("mathSolver");
%mathSolverInput = struct(); 


%% Postprocessing step
numDecoys = numel(params.decoys);
decoyProbs = cell2mat(params.decoyProbs);
decoys = cell2mat(params.decoys); %decoy intensity.

for i = 1:numDecoys
   FrequencyCondIntensity(:,:,i) = ( diag(params.probSignalsA)*params.expectationsConditional(:,:,i) ); %Frequency of events conditioned on intensity..
end

JointFreqofEvents = FrequencyCondIntensity.*reshape(decoyProbs,1,1,3);
%Alice is ordered HVAD Bob isordered  noclick- H -V -doubleclicks- noclicck
%-A -D-doubleclick


%These are actually frequencies, we still call them nX for easy comparison
%with finite size code.

nXListFreqFreq = zeros(numDecoys,1);
nXneqListFreqFreq = zeros(numDecoys,1);
nZListFreq = zeros(numDecoys,1);
nZneqListFreq = zeros(numDecoys,1);


for i = 1:numDecoys
    nZListFreq(i) = sum(JointFreqofEvents(1:2,2:4,i),'all');
    nXListFreq(i) = sum(JointFreqofEvents(3:4,6:8,i),'all');
   
    nXneqListFreq(i) = JointFreqofEvents(3,7)+ JointFreqofEvents(4,6) ...  %single licks
                  +0.5*(JointFreqofEvents(3,8)+JointFreqofEvents(4,8)); %double clicks
    nZneqListFreq(i) =  JointFreqofEvents(1,3)+ JointFreqofEvents(2,2) ...  %single licks
                  +0.5*(JointFreqofEvents(1,4)+JointFreqofEvents(2,4)) ; %doube clicks

end
nKListFreq = (1-params.pT)*nZListFreq; %number of rounds used for key generation.



%% Error correction
% We use a simply error-correction formula. 
eZobs = sum(nZneqListFreq,'all')/ sum(nZListFreq,'all'); %error rate observed in Z basis..
deltaLeak = params.fEC*sum(nKListFreq,'all')*binaryEntropy(eZobs);
debugInfo.storeInfo("deltaLeak",deltaLeak);






%% Decoy annalysis

%setup variables
nXFreq = sum(nXListFreq);
nKFreq = sum(nKListFreq);
nXneqFreq = sum(nXneqListFreq);

% tau n = probability n photon being emitted;
tauZero = computeTau(0, decoys, decoyProbs);
tauOne = computeTau(1, decoys, decoyProbs);


factor3 = exp(-decoys(3)) / decoyProbs(3);
factor2 = exp(-decoys(2)) / decoyProbs(2);
factor1 = exp(-decoys(1)) / decoyProbs(1);


% lower bound on zero-photon nZ
nKzeroLower = tauZero * ( decoys(2)*factor3*nKListFreq(3) - decoys(3)*factor2*nKListFreq(2) ) / (decoys(2) - decoys(3)) ;

% lower bound on zero-photon nX
nXzeroLower = tauZero * ( decoys(2)* nXListFreq(3) - decoys(3)*nXListFreq(2) ) / (decoys(2) - decoys(3)) ;


%Lower bound on single-photon nZ component
numerK = tauOne*decoys(1)* (factor2*nKListFreq(2) - factor3*nKListFreq(3) - (decoys(2)^2 - decoys(3)^2) / decoys(1)^2 * (factor1*nKListFreq(1) - nKzeroLower / tauZero)); 
denomK = decoys(1)*(decoys(2) - decoys(3)) - decoys(2)^2 + decoys(3)^2;
nKoneLower = numerK / denomK;

%Lower bound on single-photon nX component
numerX = tauOne*decoys(1)* (factor2*nXListFreq(2) - factor3*nXListFreq(3) - (decoys(2)^2 - decoys(3)^2) / decoys(1)^2 * (factor1*nXListFreq(1) - nXzeroLower / tauZero)); 
denomX = decoys(1)*(decoys(2) - decoys(3)) - decoys(2)^2 + decoys(3)^2;
nXoneLower = numerX / denomX;


% upper bound on zero-photon nXneq
nXneqoneUpper = tauOne * ( factor2*nXneqListFreq(2) - factor3*nXneqListFreq(3)) / (decoys(2) - decoys(3)) ;


% upper bound on zero-photon eXobs
eXobsoneUpper = nXneqoneUpper / nXoneLower;





% Just write keyrate expression here...   


deltaone = computeDeltaOne(params.DeltaLoss, params.DeltaDarkCounts , params.dcRate);
deltatwo = computeDeltaTwo(params.DeltaLoss, params.DeltaDarkCounts , params.dcRate);



gammaserf = 0;
gammaone = 0;
gammatwo = 0;

      



eXphaseBound = (eXobsoneUpper + gammaserf + deltaone + gammaone)/(1-gammatwo-deltatwo);

    
keyRate = nKoneLower*(1- binaryEntropy(eXphaseBound)) - deltaLeak;



end




%%VALIDATION FUNCTIONS

function eachRowMustBeAProbDist(expectationsConditional)

% get the dimensions of the conditional expectations. Then based on that
% pick a strategy to handle it
dimExpCon = size(expectationsConditional);

errorID ="BasicBB84WCPDecoyKeyRateFunc:InvalidRowsAreNotProbDists";
errorTXT = "A row in the conditional distribution is not a valid probability distribution.";

if numel(dimExpCon) == 2 % Matlab's minimum number of dimensions is 2.
    % The array is 2d and the slicing is easy
    for index = 1:dimExpCon(1)
        if~isProbDist(expectationsConditional(index,:))
           throwAsCaller(MException(errorID,errorTXT));
        end
    end
else
    % We have some tricky slicing to do for 3 plus dimensions.
    % We need to index the first dimension and the combination of
    % dimensions 3 and up. The second dimension will just use :.
    maskedDims = [dimExpCon(1),prod(dimExpCon(3:end))];

    for index = 1:prod(maskedDims)
        vecIndex = ind2subPlus(maskedDims,index);
        if ~isProbDist(expectationsConditional(vecIndex(1),:,vecIndex(2)))
            throwAsCaller(MException(errorID,errorTXT));
        end
    end
end
end


function mustBeProbDistCell(input)
mustBeProbDist([input{:}])
end

function output = computeTau(n, intList, probList)
output = 0;

for i = 1 : numel(intList)
    output = output + probList(i)*exp(-intList(i)) * (intList(i))^n * probList(i) / factorial(n);
end

end


