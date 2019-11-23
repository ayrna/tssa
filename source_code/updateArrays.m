%   COPYRIGHT
%   This file is part of TSSA: https://github.com/ayrna/tssa
%   Original authors: Antonio M. Duran Rosal, Pedro A. Gutierrez
%   Copyright:
%       This software is released under the The GNU General Public License v3.0 licence
%       available at http://www.gnu.org/licenses/gpl-3.0.html
%   Citation: If you use this code, please cite the following paper:
%     [1] A.M. Durán-Rosal, P.A. Gutiérrez, Á. Carmona-Poyato and C. Hervás-Martínez.
%         "A hybrid dynamic exploitation barebones particle swarm optimisation
%         algorithm for time series segmentation", Neurocomputing,
%         Vol. 353, August, 2019, pp. 45-55.
%		  https://doi.org/10.1016/j.neucom.2018.05.129
%
%% updateArrays
% Function: Update the best global position and the best personal positions of each particle
% 
% Input:
%     currentPopulationInt:     population
%     currentFitness:           fitness of the population
%     bestLocalPopulationInt:   best local positions
%     bestLocalFitness:         best local fitness
%     bestIndividual:           best global position
%     bestIndividualFitness:    best global fitness
%     
% Output:
%     bestLocalPopulationInt:   updated best local positions
%     bestLocalFitness:         updated best local fitness
%     bestIndividual:           updated best global position
%     bestIndividualFitness:    updated best global fitness
function [bestLocalPopulationInt, bestLocalFitness, bestIndividual, bestIndividualFitness] = updateArrays(currentPopulationInt, currentFitness, bestLocalPopulationInt, bestLocalFitness, bestIndividual, bestIndividualFitness)
    ind = find(currentFitness > bestLocalFitness);
    bestLocalPopulationInt(ind,:)=currentPopulationInt(ind,:);
    bestLocalFitness(ind)=currentFitness(ind);

    [fbest, indMax] = max(bestLocalFitness);
    if fbest > bestIndividualFitness,
        bestIndividual = bestLocalPopulationInt(indMax,:);
        bestIndividualFitness = fbest;                
    end 
end
