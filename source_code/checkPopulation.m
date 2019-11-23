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
%     [2] A.M. Durán-Rosal, P.A. Gutiérrez, S. Salcedo-Sanz and C. Hervás-Martínez.
%         "A statistically-driven Coral Reef Optimization algorithm for optimal
%         size reduction of time series", Applied Soft Computing, 
%         Vol. 63. 2018, pp. 139-153.
%         https://doi.org/10.1016/j.asoc.2017.11.037
%
%% checkPopulation
% Function: Check the chromosomes are in the feasible area of the problem
% 
% Input:
%     newCurrentPopulationInt:  currentPopulation
%     currentPopulationInt:     currentPopulation of the previous iteration
%     oldFitness:               fitness
%     nPobl:                    population size
%     sizeChromosome:           binary chromosome length
%     
% Output:
%     newCurrentPopulationInt:  new integer population
%     newFitness:               new fitness
function [newCurrentPopulationInt, newFitness] = checkPopulation(newCurrentPopulationInt, currentPopulationInt,oldFitness,nPobl,sizeChromosome)
    newFitness = oldFitness*NaN;
    for i=1:nPobl,
        ind = round(newCurrentPopulationInt(i,:));
        if numel(find(ind > sizeChromosome-1)) > 0 || numel(find(ind < 2)) > 0,
            minNew=min(newCurrentPopulationInt(i,:));
            maxNew=max(newCurrentPopulationInt(i,:));
            minCurrent = min(currentPopulationInt(i,:));
            maxCurrent = max(currentPopulationInt(i,:));
            newCurrentPopulationInt(i,:)=(newCurrentPopulationInt(i,:)-minNew)/(maxNew - minNew)*(maxCurrent-minCurrent)+minCurrent;
        end
    end
end
