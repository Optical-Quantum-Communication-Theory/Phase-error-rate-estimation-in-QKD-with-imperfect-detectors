
function [newParams, modParser]= ActiveBB84WCPDecoyChannelFunc(params,options,debugInfo)
% ActiveBB84WCPDecoyChannelFunc A channel function for active BB84 using WCP
% states, supporting decoy intensities. Given a collection of decoy
% intensities, this channel produces a group of 4x8 tables of
% expectations, one for each decoy intensity, which are the conditional
% probability for each of Bob's 8 detector patterns given Alice's signal
% sent (and the decoy intenisty).
%
% Input parameters:
% * decoys: a cell of the intensities used in decoy analysis. These are the
%   mean photon numbers that Alice can choose from when performing the
%   decoy protocol. The first element in the cell array is treated as the
%   intensity used for key generation.
% * eta (1): the transmissivity of the quantum channel; Must be between 0 
%   and 1 inclusive.
% * detectorEfficiency ([1,1,1,1]): Vector of efficiencies of Bob's detectors. Must be
%   between 0 and 1 inclusive. Must be indexed in order for each detector
%   - (H,V,A,D).
% * misalignmentAngle (0):  Physical angle of misalignment between Alice
%   and Bob's measurements around Y axix. This angle is measured as the
%   physical rotation of the device (period 2pi). Although calculations are
%   done on the Bloch sphere, angles should not be given in that form
%   (period 4pi).
% * darkCountRate ([0,0,0,0]): Vector of probability that a detector that recieves no
%   photons will still randomly click anyway. Must be between 0 and 1. Must
%   be indexed in order for each detector - (H,V,A,D).
% Output parameters:
% * expectationsConditional: The conditional expectations (as a 3D array)
%   from Alice and Bob's measurements. This should be organized as a 4 x 8
%   x n array, where 4 = number of signals Alice sent, 8 = Bob's outcomes
% , and n = the number of intensities used in the decoy
%   protocol. The Table is conditioned on the signal Alice sent, and the
%   intensity she chose. Therefore, each row should sum to 1.
% Options:
% * None.
% DebugInfo:
% * transMat: The linear operator that trasforms the mode operators from
%   what Alice sent, to what Bob recieves. This includes Bob's detector
%   setup, except for the non-linear dark counts. See Coherent for more
%   details.
% * probDetectorClickCon: The probability of each individual detector
%   clicking given the signal choice Alice sent (and intensity in dim 3).
%
% 
% See also QKDChannelModule, Coherent
arguments
    params (1,1) struct
    options (1,1) struct
    debugInfo (1,1) DebugInfo
end

%% options parser
optionsParser = makeGlobalOptionsParser(mfilename);
optionsParser.parse(options);
options = optionsParser.Results;

%% module parser
modParser = moduleParser(mfilename);

%Decoy intensities
modParser.addRequiredParam("decoys", @(x) mustBeCellOf(x, 'numeric'));
modParser.addAdditionalConstraint(@(x) allCells(x,@isscalar),"decoys");
modParser.addAdditionalConstraint(@(x) allCells(x,@(y) y>=0),"decoys");

%Z-basis choice
modParser.addRequiredParam("pz", @(x) mustBeInRange(x, 0, 1));
modParser.addAdditionalConstraint(@isscalar,"pz");

%Channel loss
modParser.addRequiredParam("eta", @(x) mustBeInRange(x, 0, 1));
modParser.addAdditionalConstraint(@isscalar,"eta");

%Detector efficiency
modParser.addOptionalParam("detectorEfficiency", [1,1,1,1], @(x) mustBeInRange(x, 0, 1));

%Misalingment angle
modParser.addOptionalParam("misalignmentAngle",0,@mustBeReal);
modParser.addAdditionalConstraint(@isscalar,"misalignmentAngle");

%Darkcount rate
modParser.addOptionalParam("darkCountRate", [0,0,0,0], @(x) mustBeInRange(x, 0, 1));

modParser.parse(params);

params = modParser.Results;

pz = params.pz;
px = 1-pz;
dc = params.darkCountRate;
eff = params.detectorEfficiency;
theta = params.misalignmentAngle;

DarkCountMatrixZ = DarkCountMatrix(dc(1),dc(2)); % Dark count post-processing in the Z basis
DarkCountMatrixX = DarkCountMatrix(dc(3),dc(4)); % Dark count post-processing in the X basis

DarkCountPostProcessing = blkdiag(DarkCountMatrixZ,DarkCountMatrixX); % Full post-processing in both bases

for index = 1:numel(params.decoys)
    intensity = params.eta*params.decoys{index}; % Intensity entering Bob's detection setup after channel loss
    intensityH = eff(1)*intensity; % Intensity entering the H detector after detector loss
    intensityV = eff(2)*intensity; % Intensity entering the V detector after detector loss
    intensityA = eff(3)*intensity; % Intensity entering the A detector after detector loss
    intensityD = eff(4)*intensity; % Intensity entering the D detector after detector loss

    % If Alice sends H. Represented by the first index of p
    % Accounting for misalignment
    intensityHH = intensityH*cos(theta)^2;
    intensityVH = intensityV*sin(theta)^2;
    intensityDH = intensityD*((cos(theta)+sin(theta))^2)/2;
    intensityAH = intensityA*((cos(theta)-sin(theta))^2)/2;
    [probZ,probX] = ActiveBB84Setup(intensityHH, intensityVH, intensityAH, intensityDH);
    p(1,1:4) = pz*probZ;
    p(1,5:8) = px*probX;

    % If Alice sends V. Represented by the first index of p
    % Accounting for misalignment
    intensityHV = intensityH*sin(theta)^2;
    intensityVV = intensityV*cos(theta)^2;
    intensityDV = intensityD*((cos(theta)-sin(theta))^2)/2;
    intensityAV = intensityA*((cos(theta)+sin(theta))^2)/2;
    [probZ,probX] = ActiveBB84Setup(intensityHV, intensityVV, intensityAV, intensityDV);
    p(2,1:4) = pz*probZ;
    p(2,5:8) = px*probX;

    % If Alice sends A. Represented by the first index of p
    % Accounting for misalignment
    intensityHA = intensityH*((cos(theta)-sin(theta))^2)/2;
    intensityVA = intensityV*((cos(theta)+sin(theta))^2)/2;
    intensityDA = intensityD*sin(theta)^2;
    intensityAA = intensityA*cos(theta)^2;
    [probZ,probX] = ActiveBB84Setup(intensityHA, intensityVA, intensityAA, intensityDA);
    p(3,1:4) = pz*probZ;
    p(3,5:8) = px*probX;

    % If Alice sends D. Represented by the first index of p
    % Accounting for misalignment
    intensityHD = intensityH*((cos(theta)+sin(theta))^2)/2;
    intensityVD = intensityV*((cos(theta)-sin(theta))^2)/2;
    intensityDD = intensityD*cos(theta)^2;
    intensityAD = intensityA*sin(theta)^2;
    [probZ,probX] = ActiveBB84Setup(intensityHD, intensityVD, intensityAD, intensityDD);
    p(4,1:4) = pz*probZ;
    p(4,5:8) = px*probX;

    % Now apply dark counts
    for aliceIndex = 1:4
        p(aliceIndex,:) = (DarkCountPostProcessing*p(aliceIndex,:)')';
    end
    expectationsCon(:,:,index) = p;
end

newParams.expectationsConditional = expectationsCon;


end

% Function to create the dark count post-processing for a setup with 2
% detectors
function dcPP = DarkCountMatrix(dc1, dc2)

   dcPP = [(1-dc1)*(1-dc2),0,0,0;...
                 dc1*(1-dc2),(1-dc2),0,0;...
                 dc2*(1-dc1),0,(1-dc1),0;...
                 dc1*dc2,dc2,dc1,1];
end


function [probZ,probX] = ActiveBB84Setup(intensityH,intensityV,intensityA,intensityD)
% Input parameters:
% * intensityH: Intensity entering Bob's 'H' detector
% * intensityV: Intensity entering Bob's 'V' detector
% * intensityA: Intensity entering Bob's 'A' detector
% * intensityD: Intensity entering Bob's 'D' detector
% Output parameters:
% * probZ : Vector of the probability of the different detection events
% conditioned on Bob's choice of Z-basis
% * probX : Vector of the probability of the different detection events
% conditioned on Bob's choice of X-basis

    % First Z-basis
    probZ(1) = exp(-(intensityH+intensityV)); % no-click in z-basis
    probZ(2) = (1-exp(-intensityH))*exp(-intensityV); % click in H
    probZ(3) = (1-exp(-intensityV))*exp(-intensityH); % click in V
    probZ(4) = (1-exp(-intensityH))*(1-exp(-intensityV)); % double-click in z-basis
    % Now X-basis
    probX(1) = exp(-(intensityA+intensityD)); % no-click in x-basis
    probX(2) = (1-exp(-intensityA))*exp(-intensityD); % click in A
    probX(3) = (1-exp(-intensityD))*exp(-intensityA); % click in D
    probX(4) = (1-exp(-intensityA))*(1-exp(-intensityD)); % double-click in x-basis
end
