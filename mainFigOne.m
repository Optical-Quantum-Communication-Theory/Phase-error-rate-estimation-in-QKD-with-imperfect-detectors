%Preset file for BB84 protocol using weak coherent pulses with decoy using
%EUR security proof technique


InputList = [0.00, 0.00, 0 
            0.02, 0.02, 0; 
            0.02, 0.02, 1;
             0.05, 0.05, 0;
             0.05, 0.05, 1]; 
% Delta for loss and dark counts in detectors and randomSwaps
% Last column specifies random swapping (1 means we swap).







for i = 1:size(InputList,1)

    qkdInput = FigOnePreset(InputList(i,:,:));

    %run the QKDSolver with this input
    results = MainIteration(qkdInput);
  


    folderpath = "./dataFigOne";
    filename = "FigOne_" + num2str(i) + ".mat";
    
    
    pathname = folderpath + "/" + filename;
    %save the results and preset to a file in the folder "data"
    

    % Check if the folder exists, and create it if it doesn't
    if ~exist(folderpath, 'dir')
        mkdir(folderpath);
    end

    save(pathname,"results","qkdInput");

end


PlotFigOne(InputList);

