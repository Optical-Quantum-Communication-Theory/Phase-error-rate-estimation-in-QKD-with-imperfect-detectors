function PlotFigThree(InputList)



%% plot the result

Colour = {"#0072BD"	,"#D95319" ,"#EDB120", "#7E2F8E","#77AC30",	"#4DBEEE", "#A2142F"};
styles = {'-', '--', '-', '--','-', '--', ':'}; % Line styles
markers = {'o', 'x', '*', '+','square', 'diamond', "."}; % Marker types



for i = 1: size(InputList,1)

  
    filename = "FigThree_" + num2str(i) + ".mat";
    
    pathname = "./dataFigThree"  + "/" + filename;
    data = load(pathname);

    
    keyRate = [data.results.keyRate];
    keyRate(keyRate<0) = 0;
   

    DeltaList = cell2mat(data.qkdInput.scanParameters.DeltaLoss);
    
    % Create string for legend
    exponent = num2str(log10(InputList(i,1)) );

    stringLegend = "N_{\mathrm{tot}}=10^{"+exponent+"}";

    if InputList(i,2) == 1
        stringRandomSwap = "\mathrm{Random Swap}";
    else
        stringRandomSwap = "";
    end
    
   
    if ~isempty(stringRandomSwap)
        stringLegend = stringLegend + ", " + stringRandomSwap; % Concatenate only if non-empty
    end

    semilogy(DeltaList , keyRate, [markers{i},styles{i}], "Color",Colour{i}, 'DisplayName', "$"+stringLegend +"$");
    hold on;

end


% Other things

fsize = 50;
lw = 4;
ms = 13;


xlabel("$\Delta_{\eta} = \Delta_{\mathrm{dc}}$",'Interpreter','latex');
ylabel('Key Rate','Interpreter','latex');
legend('Location','best','Interpreter','latex');
fontsize(gca,fsize,"pixels")

set(findall(gcf,'Type','line'),'LineWidth',lw);
set(findall(gcf,'Type','line'),'MarkerSize',ms);

yl = ylim;        % Get current y-axis limits
ylim([1e-7 yl(2)]); %sets lower limit, keeps upper limit unchanged.
grid on;
hold off;

savefig("./dataFigThree/dataFigThree"  + ".fig");
end



%plotting(numel(DeltaList)); % Looks at the key rate data and plots.


%QKDPlot.simple1DPlot(qkdInput,results,"xScaleStyle","dB","yScaleStyle","log")
