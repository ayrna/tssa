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
%     [2] M. Pérez-Ortiz, A.M. Durán-Rosal, P.A. Gutiérrez, et al.
%         "On the use of evolutionary time series analysis for segmenting paleoclimate data"
%         Neurocomputing, Vol. 326-327, January, 2019, pp. 3-14
%         https://doi.org/10.1016/j.neucom.2016.11.101
%     [3] A.M. Durán-Rosal, P.A. Gutiérrez, F.J. Martínez-Estudillo and C. Hervás-Martínez.
%         "Simultaneous optimisation of clustering quality and approximation error
%         for time series segmentation", Information Sciences, Vol. 442-443, May, 2018, pp. 186-201.
%         https://doi.org/10.1016/j.ins.2018.02.041
%     [4] A.M. Durán-Rosal, P.A. Gutiérrez, S. Salcedo-Sanz and C. Hervás-Martínez.
%         "A statistically-driven Coral Reef Optimization algorithm for optimal
%         size reduction of time series", Applied Soft Computing, 
%         Vol. 63. 2018, pp. 139-153.
%         https://doi.org/10.1016/j.asoc.2017.11.037
%     [5] A.M. Durán-Rosal, J.C. Fernández, P.A. Gutiérrez and C. Hervás-Martínez.
%         "Detection and prediction of segments containing extreme significant wave heights"
%         Ocean Engineering, Vol. 142, September, 2017, pp. 268-279.
%         https://doi.org/10.1016/j.oceaneng.2017.07.009
%     [6] A.M. Durán-Rosal, M. de la Paz Marín, P.A. Gutiérrez and C. Hervás-Martínez.
%         "Identifying market behaviours using European Stock Index time series by 
%         a hybrid segmentation algorithm", Neural Processing Letters,
%         Vol. 46, December, 2017, pp. 767–790.
%         https://doi.org/10.1007/s11063-017-9592-8
%
%% crossoverStr2Op3
% Function: The algorithm determines pairs of parents
%           Operator3: Single point cross over operator which respects the number of segments
% Input:
%     population:     set of segmentations
%     fitness:        fitness value for each segmentation
%     maxAttempts:    maximum number of attempts to re-apply failed crossover
%     
% Output:
%     crossedPopulation:  population after applying crossover
%     fitnessChanged:     fitness of each segmentation (NaN in the case of changes)
function [crossedPopulation,fitnessChanged] = crossoverStr2Op3(population,fitness,pCross,maxAttempts)
    crossedPopulation = population;
    fitnessChanged = fitness;
    [nPop, nCutPoints] = size(population);
    indexes = 1:nPop;
    randIndexes = indexes(randperm(length(indexes)));
    i=1;
    while i<nPop,
        if rand() < pCross,
            ind1 = randIndexes(i);
            ind2 = randIndexes(i+1);
            attempt2=0;
            while attempt2<maxAttempts,
                [crossedPopulation(ind1,:),crossedPopulation(ind2,:),flag] = crossoverOperator3(population(ind1,:),population(ind2,:));
                if flag == false,
                    attempt2=attempt2+1;
                else
                    attempt2=maxAttempts+1;
                end
            end
            population(ind1,:) = crossedPopulation(ind1,:);
            population(ind2,:) = crossedPopulation(ind2,:);
            fitnessChanged(ind1) = NaN;
            fitnessChanged(ind2) = NaN;
        end
        i=i+2;
    end
end
