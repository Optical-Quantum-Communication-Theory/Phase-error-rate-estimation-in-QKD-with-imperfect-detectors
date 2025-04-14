function [keyRate, modParser] = BB84WCPDecoyEURKeyRateFunc(params,options,mathSolverFunc,debugInfo)
% BB84WCPDecoyEURKeyRateFunc A key rate function for a WCP BB84 protocol
% using decoy states using the EUR proof technique for imperfect detectors.
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

%% options parser
optionsParser = makeGlobalOptionsParser(mfilename);
optionsParser.addOptionalParam("decoyTolerance",1e-14,@(x) x>=0);
optionsParser.addOptionalParam("decoySolver","SDPT3");
optionsParser.addOptionalParam("decoyPrecision","high");
optionsParser.addOptionalParam("decoyForceSep",false, @islogical);
optionsParser.addOptionalParam("decoyPhotonCutOff",10,@(x)mustBeInteger(x));
optionsParser.addAdditionalConstraint(@(x) x>0,"decoyPhotonCutOff");
optionsParser.parse(options);
options = optionsParser.Results;


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
modParser.addOptionalParam("DeltaDarkCounts", 0 ,   @(x) x>=0 & x<=1 ); % Default value is overwritten.
modParser.addOptionalParam("OverwriteDeltaDarkCounts", 0, @(x) x==0 || x==1); 
modParser.addRequiredParam("randomSwapFlag", @(x) x==0 || x == 1);

modParser.addRequiredParam("etaDet", @(x) x>=0 & x<=1 );
modParser.addRequiredParam("dcRate", @(x) x>=0 & x<=1 );

%finite Size stuff;


modParser.addRequiredParam("N", @(x) mustBeGreaterThan(x, 0));
modParser.addRequiredParam("logEpsilon");

modParser.parse(params);

params = modParser.Results;


logEpsilon = params.logEpsilon;


if params.OverwriteDeltaDarkCounts == 1
    params.DeltaDarkCounts = params.DeltaLoss;
end


%% simple setup (Not needed since we do not do numerical optimization)
%debugMathSolver = debugInfo.addLeaves("mathSolver");
%mathSolverInput = struct(); 


%% Postprocessing step
numDecoys = numel(params.decoys);
decoyProbs = cell2mat(params.decoyProbs);
decoys = cell2mat(params.decoys); %decoy intensity.

for i = 1:numDecoys
    NumberCondIntensity(:,:,i) = params.N*( diag(params.probSignalsA)*params.expectationsConditional(:,:,i) ); %number of events conditioned on intensity..
end

JointNumberofEvents = NumberCondIntensity.*reshape(decoyProbs,1,1,3);
%Alice is ordered HVAD Bob isordered  noclick- H -V -doubleclicks- noclicck
%-A -D-doubleclick


nXList = zeros(numDecoys,1);
nXneqList = zeros(numDecoys,1);
nZList = zeros(numDecoys,1);
nZneqList = zeros(numDecoys,1);


for i = 1:numDecoys
    nZList(i) = sum(JointNumberofEvents(1:2,2:4,i),'all');
    nXList(i) = sum(JointNumberofEvents(3:4,6:8,i),'all');
   
    nXneqList(i) = JointNumberofEvents(3,7,i)+ JointNumberofEvents(4,6,i) ...  %single licks
                  +0.5*(JointNumberofEvents(3,8,i)+JointNumberofEvents(4,8,i)); %double clicks
    nZneqList(i) =  JointNumberofEvents(1,3,i)+ JointNumberofEvents(2,2,i) ...  %single licks
                  +0.5*(JointNumberofEvents(1,4,i)+JointNumberofEvents(2,4,i)) ; %doube clicks

end
nKList = (1-params.pT)*nZList; %number of rounds used for key generation.



%% Error correction
% We use a simply error-correction formula. 
eZobs = sum(nZneqList,'all')/ sum(nZList,'all'); %error rate observed in Z basis..
deltaLeak = params.fEC*sum(nKList,'all')*binaryEntropy(eZobs);
debugInfo.storeInfo("deltaLeak",deltaLeak);






%% Decoy annalysis

%setup variables
nX = sum(nXList);
nK = sum(nKList);
nXneq = sum(nXneqList);

% tau n = probability n photon being emitted;
tauZero = computeTau(0, decoys, decoyProbs);
tauOne = computeTau(1, decoys, decoyProbs);

nKmu3Minus = exp( +decoys(3)) / decoyProbs(3) * ( nKList(3) - decoyHoeffdingCorrection(nK, logEpsilon.decoy) );
nKmu2Minus = exp( +decoys(2)) / decoyProbs(2) * ( nKList(2) - decoyHoeffdingCorrection(nK, logEpsilon.decoy) );
nKmu1Minus = exp( +decoys(1)) / decoyProbs(1) * ( nKList(1) - decoyHoeffdingCorrection(nK, logEpsilon.decoy) );

nKmu3Plus = exp( +decoys(3)) / decoyProbs(3) * ( nKList(3) + decoyHoeffdingCorrection(nK, logEpsilon.decoy) );
nKmu2Plus = exp( +decoys(2)) / decoyProbs(2) * ( nKList(2) + decoyHoeffdingCorrection(nK, logEpsilon.decoy) );
nKmu1Plus = exp( +decoys(1)) / decoyProbs(1) * ( nKList(1) + decoyHoeffdingCorrection(nK, logEpsilon.decoy) );


nXmu3Minus = exp( +decoys(3)) / decoyProbs(3) * ( nXList(3) - decoyHoeffdingCorrection(nX, logEpsilon.decoy) );
nXmu2Minus = exp( +decoys(2)) / decoyProbs(2) * ( nXList(2) - decoyHoeffdingCorrection(nX, logEpsilon.decoy) );
nXmu1Minus = exp( +decoys(1)) / decoyProbs(1) * ( nXList(1) - decoyHoeffdingCorrection(nX, logEpsilon.decoy) );

nXmu3Plus = exp( +decoys(3)) / decoyProbs(3) * ( nXList(3) + decoyHoeffdingCorrection(nX, logEpsilon.decoy) );
nXmu2Plus = exp( +decoys(2)) / decoyProbs(2) * ( nXList(2) + decoyHoeffdingCorrection(nX, logEpsilon.decoy) );
nXmu1Plus = exp( +decoys(1)) / decoyProbs(1) * ( nXList(1) + decoyHoeffdingCorrection(nX, logEpsilon.decoy) );


nXneqmu3Minus = exp( +decoys(3)) / decoyProbs(3) * ( nXneqList(3) - decoyHoeffdingCorrection(nXneq, logEpsilon.decoy) );
nXneqmu2Minus = exp( +decoys(2)) / decoyProbs(2) * ( nXneqList(2) - decoyHoeffdingCorrection(nXneq, logEpsilon.decoy) );
nXneqmu1Minus = exp( +decoys(1)) / decoyProbs(1) * ( nXneqList(1) - decoyHoeffdingCorrection(nXneq, logEpsilon.decoy) );

nXneqmu3Plus = exp( +decoys(3)) / decoyProbs(3) * ( nXneqList(3) + decoyHoeffdingCorrection(nXneq, logEpsilon.decoy) );
nXneqmu2Plus = exp( +decoys(2)) / decoyProbs(2) * ( nXneqList(2) + decoyHoeffdingCorrection(nXneq, logEpsilon.decoy) );
nXneqmu1Plus = exp( +decoys(1)) / decoyProbs(1) * ( nXneqList(1) + decoyHoeffdingCorrection(nXneq, logEpsilon.decoy) );









% lower bound on zero-photon nZ
nKzeroLower = tauZero * ( decoys(2)*nKmu3Minus - decoys(3)*nKmu2Plus ) / (decoys(2) - decoys(3)) ;
nKzeroLower = max(nKzeroLower,0);

% lower bound on zero-photon nX
nXzeroLower = tauZero * ( decoys(2)*nXmu3Minus - decoys(3)*nXmu2Plus ) / (decoys(2) - decoys(3)) ;
nXzeroLower = max(nXzeroLower,0);

%Lower bound on single-photon nZ component
numerK = tauOne*decoys(1)* (nKmu2Minus - nKmu3Plus - (decoys(2)^2 - decoys(3)^2) / decoys(1)^2 * (nKmu1Plus - nKzeroLower / tauZero)); 
denomK = decoys(1)*(decoys(2) - decoys(3)) - decoys(2)^2 + decoys(3)^2;
nKoneLower = numerK / denomK;
nKoneLower = max(nKoneLower,0);

%Lower bound on single-photon nX component
numerX = tauOne*decoys(1)* (nXmu2Minus - nXmu3Plus - (decoys(2)^2 - decoys(3)^2) / decoys(1)^2 * (nXmu1Plus - nXzeroLower / tauZero)); 
denomX = decoys(1)*(decoys(2) - decoys(3)) - decoys(2)^2 + decoys(3)^2;
nXoneLower = numerX / denomX;
nXoneLower = max(nXoneLower,0);


% upper bound on zero-photon nXneq
nXneqoneUpper = tauOne * ( nXneqmu2Plus - nXneqmu3Minus) / (decoys(2) - decoys(3)) ;
% upper bound on zero-photon eXobs
eXobsoneUpper = nXneqoneUpper / nXoneLower;


eXobsoneUpper = max( min(eXobsoneUpper,1)    ,0);


% Just write keyrate expression here...   


deltaone = computeDeltaOne(params.DeltaLoss, params.DeltaDarkCounts , params.dcRate, params.randomSwapFlag);
deltatwo = computeDeltaTwo(params.DeltaLoss, params.DeltaDarkCounts , params.dcRate, params.randomSwapFlag);


if nXoneLower == 0 || nKoneLower == 0 || eXobsoneUpper == 1
    keyRate = 0;
    return;
end



gammaserf = gammaSerf(nXoneLower,nKoneLower , logEpsilon.ATa);
gammaone = gammaBin(nKoneLower, deltaone, logEpsilon.ATb);
gammatwo = gammaBin(nKoneLower, deltatwo, logEpsilon.ATc);

      



eXphaseBound = (eXobsoneUpper + gammaserf + deltaone + gammaone)/(1-gammatwo-deltatwo);

    
length = nKoneLower*(1- binaryEntropyModified(eXphaseBound)) - deltaLeak +(2+2*logEpsilon.PA) + (logEpsilon.EV-1);

%disp(length);
%disp(nK);

keyRate = length / params.N;

if ~isreal(keyRate)
    keyRate = 0;
end


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
    output = output + probList(i)*exp(-intList(i)) * (intList(i))^n  / factorial(n);
end

end


