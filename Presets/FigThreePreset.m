function qkdInput = FigThreePreset(InputValues)
% BasicBB84WCPDecoyPreset A preset for BB84 using decoy states using EUR as the proof
% technique. 

qkdInput = QKDSolverInput();

%% Parameters


qkdInput.addScanParameter("DeltaLoss", num2cell( linspace(0.00, 0.50, 25)) )
qkdInput.addFixedParameter("OverwriteDeltaDarkCounts",1);

%DeltaDarkCounts is set to DeltaLoss in KeyRateFunc.


qkdInput.addFixedParameter("eta", 10.^(-25/10) ); %Loss is 25db
% qkdInput.addFixedParameter("transmittance",1);

% decoys should be in one group, which can be created with these lines:

mu1 = 0.9;
mu2 = 0.1;
mu3= 0.00;

if (mu1 <= mu2+mu3 || mu2 <= mu3)
    warning('Intensities do not satisfy mu1 >= mu2 + mu3 or mu2 > mu3. Key rates are NOT VALID');
end

qkdInput.addFixedParameter("GROUP_decoys_1", mu1); %signal intensity
qkdInput.addFixedParameter("GROUP_decoys_2", mu2); % decoy intensity 1
qkdInput.addFixedParameter("GROUP_decoys_3", mu3); % decoy intensity 2 

qkdInput.addFixedParameter("GROUP_decoyProbs_1",1/3);
qkdInput.addFixedParameter("GROUP_decoyProbs_2",1/3);
qkdInput.addFixedParameter("GROUP_decoyProbs_3",1/3); %these are decoyProbs for X basis rounds

qkdInput.addFixedParameter("randomSwapFlag", InputValues(2)); 


qkdInput.addFixedParameter("pz",0.5);
qkdInput.addFixedParameter("pT",0.05); % probability that Z basis rounds used for testing. 

qkdInput.addFixedParameter("fEC",1.16);
qkdInput.addFixedParameter("misalignmentAngle", 2*pi/180); 

dcRate = 1e-6;
qkdInput.addFixedParameter("dcRate", dcRate);
etaDet = 0.7;
qkdInput.addFixedParameter("etaDet", dcRate);


% assign value if loss and dark count rate to each detector.
qkdInput.addFixedParameter("darkCountRate", [dcRate,dcRate,dcRate,dcRate]);
qkdInput.addFixedParameter("detectorEfficiency",[etaDet,etaDet,etaDet,etaDet]); 

% add epsilons   %parametes below are unused if we compute asymptotic key
% rates

logEpsilon.decoy = log2(1e-12); %for decoy analysis
logEpsilon.ATa = log2(1e-12); %for phase-error estimate on the single photon space. 
logEpsilon.ATb = log2(1e-12);
logEpsilon.ATc = log2(1e-12);
%logEpsilon.decoyvac = log2(1e-8); %for the 1/2 error rate in vacuum constraints.
logEpsilon.PA = log2(1e-12);
logEpsilon.EV = log2(1e-12); 

qkdInput.addFixedParameter("logEpsilon",logEpsilon);

qkdInput.addFixedParameter("N",InputValues(1));


% We actually do not need a description function since we do not do
% numerical key rate optimizations, but we include a trivial decription
% function here for compatibility with the software framework.
descriptionModule = QKDDescriptionModule(@BB84WCPDecoyEURDescriptionFunc);
qkdInput.setDescriptionModule(descriptionModule);

% channel model.
channelModule = QKDChannelModule(@ActiveBB84WCPDecoyChannelFunc);
qkdInput.setChannelModule(channelModule);

% Key rate module performs squashing and decoy analysis.
keyRateOptions = struct();


keyMod = QKDKeyRateModule(@BB84WCPDecoyEURKeyRateFunc, keyRateOptions);
qkdInput.setKeyRateModule(keyMod);

    


optimizerMod = QKDOptimizerModule(@coordinateDescentFunc,struct("verboseLevel",0),struct("verboseLevel",0));
qkdInput.setOptimizerModule(optimizerMod);
%qkdInput.addOptimizeParameter("GROUP_decoys_1",struct("lowerBound",0.2,"initVal",0.9,"upperBound",1));

% math solver options % THESE DO NOT MATTER SINCE we do not do numerics.

mathSolverOptions = struct();
mathSolverMod = QKDMathSolverModule(@FW2StepSolver,mathSolverOptions,mathSolverOptions);
qkdInput.setMathSolverModule(mathSolverMod);

% global options
qkdInput.setGlobalOptions(struct("errorHandling",3,"verboseLevel",1,"cvxSolver","mosek"));
