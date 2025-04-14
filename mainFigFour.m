%Preset file for BB84 protocol using weak coherent pulses with decoy using
%EUR security proof technique

qkdInput = FigFourPreset(1e12);

%run the QKDSolver with this input
results = MainIteration(qkdInput);


folderpath = "./dataFigFour";
filename = "FigFour_" + ".mat";


pathname = folderpath + "/" + filename;
%save the results and preset to a file in the folder "data"


% Check if the folder exists, and create it if it doesn't
if ~exist(folderpath, 'dir')
    mkdir(folderpath);
end

save(pathname,"results","qkdInput");



PlotFigFour();

