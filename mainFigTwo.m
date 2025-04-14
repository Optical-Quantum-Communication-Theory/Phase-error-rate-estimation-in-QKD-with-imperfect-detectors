%Preset file for BB84 protocol using weak coherent pulses with decoy using
%EUR security proof technique

FixedDeltaValue = 0.05;
InputList = [1e8,FixedDeltaValue, FixedDeltaValue, 0;
            1e8, FixedDeltaValue, FixedDeltaValue, 1;
            1e12, FixedDeltaValue, FixedDeltaValue, 0;
            1e12, FixedDeltaValue, FixedDeltaValue, 1;
            1e16, FixedDeltaValue, FixedDeltaValue, 0;
            1e16, FixedDeltaValue, FixedDeltaValue, 1];






for i = 1:size(InputList,1)

    qkdInput = FigTwoPreset(InputList(i,:));

    %run the QKDSolver with this input
    results = MainIteration(qkdInput);
  


    folderpath = "./dataFigTwo";
    filename = "FigTwo_" + num2str(i) + ".mat";
    
    
    pathname = folderpath + "/" + filename;
    %save the results and preset to a file in the folder "data"
    

    % Check if the folder exists, and create it if it doesn't
    if ~exist(folderpath, 'dir')
        mkdir(folderpath);
    end

    save(pathname,"results","qkdInput");
end


PlotFigTwo(InputList);

