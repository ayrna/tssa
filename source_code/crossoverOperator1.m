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
%% crossoverOperator1
% Function: Single point cross over operator with size restriction
% Input:
%     individual1:    first parent
%     individual2:    second parent
%     minSeg:         mimimum size of segment
%     maxSeg:         maximum size of segment
% 
% Output:
%     crossedIndividual1: first offspring individual
%     crossedIndividual2: second offspring individual
%     flag:
%         true - the crossover is sucesfully applied
%         false - the crossover is not sucesfully applied
function [crossedIndividual1,crossedIndividual2,flag] = crossoverOperator1(individual1,individual2,minSeg,maxSeg)
    sizeChromosome = numel(individual1);
    crossPoint = randi(sizeChromosome-3,1,1) + 1;
    crossedIndividual1(1,:) = [individual1(1,1:crossPoint) individual2(1,crossPoint+1:end)];
    crossedIndividual2(1,:) = [individual2(1,1:crossPoint) individual1(1,crossPoint+1:end)];
    point1= find(crossedIndividual1(1,1:crossPoint)==1,1,'last');
    if isempty(point1)
        point1=1;
    end
    point2= find(crossedIndividual1(1,crossPoint+1:end)==1,1,'first');
    point2= point2 + crossPoint;
    if isempty(point2)
        point2= sizeChromosome;
    end
    point3= find(crossedIndividual2(1,1:crossPoint)==1,1,'last');
    if isempty(point3)
        point3=1;
    end
    point4= find(crossedIndividual2(1,crossPoint+1:end)==1,1,'first');
    point4= point4 + crossPoint;
    if isempty(point4)
        point4= sizeChromosome;
    end

    if((point2 - point1)< minSeg) || ((point4 - point3)< minSeg),
        flag = false;
        crossedIndividual1(1,:)=individual1(1,:);
        crossedIndividual2(1,:)=individual2(1,:);
    elseif ((point2 - point1)> maxSeg) || ((point4 - point3)> maxSeg),
        flag = false;
        crossedIndividual1(1,:)=individual1(1,:);
        crossedIndividual2(1,:)=individual2(1,:);
    else
        flag = true;
    end

end
