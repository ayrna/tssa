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
%% crossoverStr1Op2
% Function: The algorithm determines if a parent is selected to be crossed. Another parent is selected randomly
%           crossoverOperator2: Single point cross over operator without size restriction
% Input:
%     population:     set of segmentations
%     fitness:        fitness value for each segmentation
%     pCross:         cross probability
%     maxAttempts:    maximum number of attempts to re-apply failed crossover
%     
% Output:
%     crossedPopulation:  population after applying crossover
%     fitnessChanged:     fitness of each segmentation (NaN in the case of changes)
function [crossedPopulation,fitnessChanged] = crossoverStr1Op2(population,fitness,pCross,maxAttempts)
    crossedPopulation = population;
    fitnessChanged = fitness;
    [nPop, nCutPoints] = size(population);
    for i=1:nPop
        %Crossover
        if rand()<pCross,
            % Find individuals to apply crossover
            ind1 = i;
            attempt = 1;
            while attempt == 1,
                ind2 = randi(nPop,1,1);
                while ind2==i,
                    ind2 = randi(nPop,1,1);
                end
                attempt2=0;
                while attempt2<maxAttempts,
                    [crossedPopulation(ind1,:),crossedPopulation(ind2,:),flag] = crossoverOperator2(population(ind1,:),population(ind2,:));
                    if flag == false,
                        attempt2=attempt2+1;
                    else
                        attempt2=maxAttempts+1;
                        attempt=0;
                    end
                end
            end
            population(ind1,:) = crossedPopulation(ind1,:);
            population(ind2,:) = crossedPopulation(ind2,:);
            fitnessChanged(ind1) = NaN;
            fitnessChanged(ind2) = NaN;                    
        end
    end
end
