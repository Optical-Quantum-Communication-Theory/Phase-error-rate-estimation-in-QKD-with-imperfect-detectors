function PlotFigOne(InputList)



%% plot the result

Colour = {"#0072BD"	,"#D95319" ,"#EDB120", "#7E2F8E","#77AC30",	"#4DBEEE", "#A2142F"};
styles = {'-', '-', '--', '-','--', '--', ':'}; % Line styles
markers = {'o', 'x', '*', '+','square', 'diamond', "."}; % Marker types



for i = 1: size(InputList,1)

  
    filename = "FigOne_" + num2str(i) + ".mat";
    
    pathname = "./dataFigOne"  + "/" + filename;
    data = load(pathname);

    
    keyRate = [data.results.keyRate];
    keyRate(keyRate<0) = 0;
   

    LossList = -10*log10(cell2mat(data.qkdInput.scanParameters.eta));
    
    % Create string for legend
    stringDeltaEta = "\Delta_{\eta} =";
    stringDeltaDarkCounts = "\Delta_{\mathrm{dc}} =";

    if InputList(i,3) == 1
        stringRandomSwap = "\mathrm{Random Swap}";
    else
        stringRandomSwap = "";
    end
    
    stringLegend = stringDeltaEta + num2str(InputList(i,1)) + ", " + stringDeltaDarkCounts + num2str(InputList(i,2));
    if ~isempty(stringRandomSwap)
        stringLegend = stringLegend + ", " + stringRandomSwap; % Concatenate only if non-empty
    end

    semilogy(LossList , keyRate, [markers{i},styles{i}], "Color",Colour{i}, 'DisplayName', "$"+stringLegend +"$");
    hold on;

end


% Other things

fsize = 50;
lw = 4;
ms = 13;


xlabel('Loss (dB)','Interpreter','latex');
ylabel('Key Rate','Interpreter','latex');
legend('Location','best','Interpreter','latex');
fontsize(gca,fsize,"pixels")

set(findall(gcf,'Type','line'),'LineWidth',lw);
set(findall(gcf,'Type','line'),'MarkerSize',ms);
grid on;
hold off;


savefig("./dataFigOne/dataFigOne"  + ".fig");
end



%plotting(numel(DeltaList)); % Looks at the key rate data and plots.


%QKDPlot.simple1DPlot(qkdInput,results,"xScaleStyle","dB","yScaleStyle","log")
