%   COPYRIGHT
%   This file is part of TSSA: https://github.com/ayrna/tssa
%   Original authors: Antonio M. Duran Rosal, Pedro A. Gutierrez
%   Copyright:
%       This software is released under the The GNU General Public License v3.0 licence
%       available at http://www.gnu.org/licenses/gpl-3.0.html
%   Citation: If you use this code, please cite any of the following papers:
%     [1] M. Pérez-Ortiz, A.M. Durán-Rosal, P.A. Gutiérrez, et al.
%         "On the use of evolutionary time series analysis for segmenting paleoclimate data"
%         Neurocomputing, Vol. 326-327, January, 2019, pp. 3-14
%         https://doi.org/10.1016/j.neucom.2016.11.101
%     [2] A.M. Durán-Rosal, P.A. Gutiérrez, F.J. Martínez-Estudillo and C. Hervás-Martínez.
%         "Simultaneous optimisation of clustering quality and approximation error
%         for time series segmentation", Information Sciences, Vol. 442-443, May, 2018, pp. 186-201.
%         https://doi.org/10.1016/j.ins.2018.02.041
%     [3] A.M. Durán-Rosal, J.C. Fernández, P.A. Gutiérrez and C. Hervás-Martínez.
%         "Detection and prediction of segments containing extreme significant wave heights"
%         Ocean Engineering, Vol. 142, September, 2017, pp. 268-279.
%         https://doi.org/10.1016/j.oceaneng.2017.07.009
%     [4] A.M. Durán-Rosal, M. de la Paz Marín, P.A. Gutiérrez and C. Hervás-Martínez.
%         "Identifying market behaviours using European Stock Index time series by 
%         a hybrid segmentation algorithm", Neural Processing Letters,
%         Vol. 46, December, 2017, pp. 767–790.
%         https://doi.org/10.1007/s11063-017-9592-8
%
%% plotSegmentedTimeSeries
% Function: Plot the original time series without the cut points and with the cut points
%
% Input:
%     model:      contains the information necessary for the reporter
%     xlb:        x asis label for the graphics
%     ylb:        y asis label for the graphics
%     dataset:    name of the dataset
%     repsuffix:  path of the output file
%     serie:      time series to be plotted
%     serie2:     time series of approximation (if needed)
%     
% Output:
%     No output variables. Only two files which contain the graphics of the segmentation (without and with cutpoints)
function plotSegmentedTimeSeries(model,xlb,ylb,dataset,repsuffix,serie,serie2)
    
    outputPDF = true;
    outputSVG = true;
    outputFile = [repsuffix filesep dataset];
    
    cuts = model.cuts;
    X(:,1)=transpose(1:numel(serie(:,1)));
    X(:,2)=serie(:,1);
    X2(:,1)=X(:,1);
    X2(:,2)=serie2(:,1);
    centroids = model.C;
    nOfClusters = numel(centroids(:,1));
    
    %% GRAPHICS %%
    markers = {'r-','g-','b-','m-','c-','y-','k-'};
    % Set lim of all graphs
    xmin=0;
    xmax=max(max(X(:,1)),max(X2(:,1)));
    ymin=min(min(X(:,2)),min(X2(:,2)));
    ymax=max(max(X(:,2)),max(X2(:,2)));


    %% Plot time series without cut points
    initialPoint = 1;
    f=figure;
    set(f, 'Position', [50 50 1200 400])
    hold on;
    set(gca,'fontsize',14,'LineWidth',1);
    ylabel(ylb,'fontsize',14)
    xlabel(xlb,'fontsize',14)

    xlim([xmin xmax]);
    ylim([ymin ymax]);

    h = zeros((nOfClusters),1);
    cuts = [cuts length(X)];
    nOfCuts = size(cuts,2);

    for c=1:(nOfCuts)
        ht = plot(X(initialPoint:cuts(c),1),...
             X(initialPoint:cuts(c),2),markers{model.L(c)},'linewidth',1);
        initialPoint = cuts(c);

        % To build the correct legend
        h(model.L(c))=ht;
    end

    % Build legend string
    legendStr = cell(nOfClusters,1);
    for j=1:nOfClusters
        legendStr{j,1} = sprintf('Cluster %d',j); 
    end

    legend(h,legendStr,'Location','NorthWest')

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
    initialPoint = 1;
    set(f2, 'Position', [50 50 1200 400])
    hold on;
    set(gca,'fontsize',14,'LineWidth',1) ;
    ylabel(ylb,'fontsize',14)
    xlabel(xlb,'fontsize',14)
    % Set lim of graph
    xlim([xmin xmax]);
    ylim([ymin ymax]);

    h2 = zeros((nOfClusters),1);
    for c=1:(nOfCuts)
        ht2 = plot(X(initialPoint:cuts(c),1),...
             X(initialPoint:cuts(c),2),markers{model.L(c)},'linewidth',1);
        initialPoint = cuts(c);

        % To build the correct legend
        h2(model.L(c))=ht2;
    end
    legend(h2,legendStr,'Location','NorthWest')


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
