%   COPYRIGHT
%   This file is part of TSSA: https://github.com/ayrna/tssa
%   Original authors: Antonio M. Duran Rosal, Pedro A. Gutierrez
%   Copyright:
%       This software is released under the The GNU General Public License v3.0 licence
%       available at http://www.gnu.org/licenses/gpl-3.0.html
%   Citation: If you use this code, please cite the following paper:
%     [1] A.M. Durán-Rosal, P.A. Gutiérrez, F.J. Martínez-Estudillo and C. Hervás-Martínez.
%         "Simultaneous optimisation of clustering quality and approximation error
%         for time series segmentation", Information Sciences, Vol. 442-443, May, 2018, pp. 186-201.
%         https://doi.org/10.1016/j.ins.2018.02.041
%
%% plotParetoFront
% Function: Plot the pareto Front of a GMO for segmentation
% 
% Input:
%     model:      contains the information necessary for the reporter
%     dataset:    name of the dataset
%     repsuffix:  path of the output file
%     
% Output:
%     No output variables. Only a file which contains the graphic of the Pareto Front
function plotParetoFront(model,dataset,repsuffix)

    outputPDF = true;
    outputSVG = true;
    outputFile = [repsuffix filesep dataset];
    f=figure;
    xmin=min(model.fitness(1,:));
    xmax=max(model.fitness(1,:));
    ymin=min(model.fitness(2,:));
    ymax=max(model.fitness(2,:));
    xlabel('Clustering Fitness','fontsize',13);
    ylabel('Error Fitness','fontsize',13);
    hold on;
    xlim([xmin xmax]);
    ylim([ymin ymax]);

    plot(model.fitness(1,1:model.number_of_firstFront),model.fitness(2,1:model.number_of_firstFront),'bo');
    plot(model.fitness(1,model.fbest:model.fbest),model.fitness(2,model.fbest:model.fbest),'b*');
    if model.number_of_firstFront < numel(model.fitness(1,:)),
        plot(model.fitness(1,model.number_of_firstFront+1:end),model.fitness(2,model.number_of_firstFront+1:end),'co');
    end

    if outputPDF
        export_fig([outputFile '_FRONT.pdf'],'-pdf','-transparent');
    end
    if outputSVG
        plot2svg([outputFile '_FRONT.svg']);
    end

    hold off;
    close all;
end
