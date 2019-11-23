%   COPYRIGHT
%   This file is part of TSSA: https://github.com/ayrna/tssa
%   Original authors: Antonio M. Duran Rosal, Pedro A. Gutierrez
%   Copyright:
%       This software is released under the The GNU General Public License v3.0 licence
%       available at http://www.gnu.org/licenses/gpl-3.0.html
%   Citation: If you use this code, please cite any of the following papers:
%     [1] A.M. Durán-Rosal, P.A. Gutiérrez, Á. Carmona-Poyato and C. Hervás-Martínez.
%         "A hybrid dynamic exploitation barebones particle swarm optimisation
%         algorithm for time series segmentation", Neurocomputing,
%         Vol. 353, August, 2019, pp. 45-55.
%		  https://doi.org/10.1016/j.neucom.2018.05.129
%     [2] A.M. Durán-Rosal, P.A. Gutiérrez, F.J. Martínez-Estudillo and C. Hervás-Martínez.
%         "Simultaneous optimisation of clustering quality and approximation error
%         for time series segmentation", Information Sciences, Vol. 442-443, May, 2018, pp. 186-201.
%         https://doi.org/10.1016/j.ins.2018.02.041
%     [3] A.M. Durán-Rosal, P.A. Gutiérrez, S. Salcedo-Sanz and C. Hervás-Martínez.
%         "A statistically-driven Coral Reef Optimization algorithm for optimal
%         size reduction of time series", Applied Soft Computing, 
%         Vol. 63. 2018, pp. 139-153.
%         https://doi.org/10.1016/j.asoc.2017.11.037
%
%% plotApproximatedTimeSeries
% Function: Plot the original time series without the cut points and with the cut points
%
% Input:
%     model:      contains the information necessary for the reporter
%     xlb:        x asis label for the graphics
%     ylb:        y asis label for the graphics
%     dataset:    name of the dataset
%     repsuffix:  path of the output file
%     serie:      time series of approximation
%     ori_serie:  to have the same scale with the original
%     
% Output:
%     No output variables. Only two files which contain the graphics of the approximation (without and with cutpoints)
function plotApproximatedTimeSeries(model,xlb,ylb,dataset,repsuffix,serie,ori_serie)
    
    outputPDF = true;
    outputSVG = true;
    outputFile = [repsuffix filesep dataset];
    
    cuts = model.cuts;
    X(:,1)=transpose(1:numel(serie(:,1)));
    X(:,2)=serie(:,1);
    X2(:,1)=X(:,1);
    X2(:,2)=ori_serie(:,1);
    
    
    %% GRAPHICS %%
    % Set lim of all graphs
    xmin=0;
    xmax=max(max(X(:,1)),max(X2(:,1)));
    ymin=min(min(X(:,2)),min(X2(:,2)));
    ymax=max(max(X(:,2)),max(X2(:,2)));


    %% Plot time series without cut points
    f=figure;
    set(f, 'Position', [50 50 1200 400])
    hold on;
    set(gca,'fontsize',14,'LineWidth',1);
    ylabel(ylb,'fontsize',14)
    xlabel(xlb,'fontsize',14)

    xlim([xmin xmax]);
    ylim([ymin ymax]);

    plot(X(:,1),X(:,2),'linewidth',1);

    hold off;

    if outputPDF
        export_fig([outputFile '.pdf'],'-pdf','-transparent');
    end
    if outputSVG
        plot2svg([outputFile '.svg']);
    end
    close all;

    %% Plot time series with cut points
    f2=figure;
    set(f2, 'Position', [50 50 1200 400])
    hold on;
    set(gca,'fontsize',14,'LineWidth',1);
    ylabel(ylb,'fontsize',14)
    xlabel(xlb,'fontsize',14)

    xlim([xmin xmax]);
    ylim([ymin ymax]);

    plot(X(:,1),X(:,2),'linewidth',1);
    
    timecuts = X(cuts,1);
    timecuts = timecuts/1.^3;
    for i=1:numel(timecuts),
        plot([timecuts(i) timecuts(i)], [min((X(:,2))) max((X(:,2)))], 'k--');
    end

    if outputPDF
        export_fig([outputFile 'cuts.pdf'],'-pdf','-transparent');
    end
    if outputSVG
        plot2svg([outputFile 'cuts.svg']);
    end
    hold off;
    
    close all;
end
