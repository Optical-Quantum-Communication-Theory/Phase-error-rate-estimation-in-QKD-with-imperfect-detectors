%%  FUNCTION NAME: binaryEntropy
% binary entropy
%%
function entr = binaryEntropyModified(probability)
if probability == 0 
    entr = 0;
elseif probability >= 0.5
    entr = 1;
else
    entr = - probability * log2(probability) - (1-probability) * log2(1 - probability);
end