function [newParams,modParser] = BB84WCPDecoyEURDescriptionFunc(params, options, debugInfo)
% BB84WCPDecoyEURDescriptionFunc : This is a trivial function. Since we do 
% utilize the numerical framework for key rate computations, this function 
% is not required. It is only included to be compatible with the software. 
%
arguments
    params (1,1) struct
    options (1,1) struct
    debugInfo (1,1) DebugInfo
end
%Parsing technical options for the module
optionsParser = makeGlobalOptionsParser(mfilename);
optionsParser.parse(options);
options = optionsParser.Results;

%% module parser
%Parsing parameters for the module
modParser = moduleParser(mfilename);
modParser.addRequiredParam("pz",@(x) mustBeInRange(x,0,1));
modParser.parse(params)
params = modParser.Results;
%% Do NOTHING

pz = params.pz;

probSignalsA = [pz/2,pz/2,(1-pz)/2,(1-pz)/2]; %probability of Alice sending each signal.

newParams.probSignalsA = probSignalsA;

end

