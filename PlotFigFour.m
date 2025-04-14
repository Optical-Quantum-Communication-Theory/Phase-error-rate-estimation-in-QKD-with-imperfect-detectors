function PlotFigFour()



%% plot the result

% Colour = {"#0072BD"	,"#D95319" ,"#EDB120", "#7E2F8E","#77AC30",	"#4DBEEE", "#A2142F"};
% styles = {'-', '--', '-', '--','-', '--', ':'}; % Line styles
% markers = {'o', 'x', '*', '+','square', 'diamond', "."}; % Marker types
% 
% 



filename = "FigFour_" + ".mat";

pathname = "./dataFigFour"  + "/" + filename;
data = load(pathname);

dataStruct = data.results;


% Preallocate arrays for the fields
numEntries = numel(dataStruct);
keyRate = zeros(numEntries, 1);
deltaLoss = zeros(numEntries, 1);
deltaDarkCounts = zeros(numEntries, 1);

% Extract the data from the struct
for i = 1:numEntries
    keyRate(i) = dataStruct(i).keyRate;
    deltaLoss(i) = dataStruct(i).currentParams.DeltaLoss;
    deltaDarkCounts(i) = dataStruct(i).currentParams.DeltaDarkCounts;
end

% Create a grid for DeltaLoss and DeltaDarkCounts
numGridPoints = 50; % Adjust grid density as needed
[gridX, gridY] = meshgrid(linspace(min(deltaLoss), max(deltaLoss), numGridPoints), ...
                          linspace(min(deltaDarkCounts), max(deltaDarkCounts), numGridPoints));

% Interpolate keyRate onto the grid
gridZ = griddata(deltaLoss, deltaDarkCounts, keyRate, gridX, gridY, 'linear');

% Plot the heatmap
figure;
imagesc(linspace(min(deltaLoss), max(deltaLoss), numGridPoints), ...
        linspace(min(deltaDarkCounts), max(deltaDarkCounts), numGridPoints), ...
        gridZ);
set(gca, 'YDir', 'normal'); % Flip y-axis to display correctly
colorbar;




% Other things

fsize = 50;
lw = 4;
ms = 13;


title('Key Rate','Interpreter','latex');
xlabel("$\Delta_{\eta}$",'Interpreter','latex');
ylabel("$\Delta_{\mathrm{dc}}$", 'Interpreter', 'latex');
% legend('Location','best','Interpreter','latex');
fontsize(gca,fsize,"pixels")

set(findall(gcf,'Type','line'),'LineWidth',lw);
set(findall(gcf,'Type','line'),'MarkerSize',ms);

yl = ylim;        % Get current y-axis limits
ylim([1e-7 yl(2)]); %sets lower limit, keeps upper limit unchanged.
grid on;
hold off;

savefig("./dataFigFour/dataFigFour"  + ".fig");
end



