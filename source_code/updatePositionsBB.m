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
%% updatePositionsBB
% Function: Make an iteration of the barebones PSO algorithm
% 
% Input:
%     currentPopulationInt:     population
%     bestLocalPopulationInt:   best local positions
%     bestIndividualInt:        best global position
%     oldFitness:               current fitness
%     nPobl:                    population size
%     sizeChromosome:           binary chromosome length
%     sizeChromosomeInt:        integer chromosome length
%     
% Output:
%     newCurrentPopulationInt:  new integer population
%     newFitness:               new fitness
function [newCurrentPopulationInt,newFitness] = updatePositionsBB(currentPopulationInt,bestLocalPopulationInt,bestIndividualInt,oldFitness,nPobl,sizeChromosome,sizeChromosomeInt)

    newCurrentPopulationInt = currentPopulationInt;
    for i=1:nPobl,
        for j=1:sizeChromosomeInt,
            if rand > 0.5,
                a = (bestLocalPopulationInt(i,j) + bestIndividualInt(j))/2;
                b = abs(bestLocalPopulationInt(i,j) - bestIndividualInt(j));
                newCurrentPopulationInt(i,j) = normrnd(a,b);
            else
                newCurrentPopulationInt(i,j) = bestLocalPopulationInt(i,j);
            end
        end
    end

    newCurrentPopulationInt = sort(newCurrentPopulationInt,2,'ascend');
    [newCurrentPopulationInt, newFitness] = checkPopulation(newCurrentPopulationInt,currentPopulationInt,oldFitness,nPobl,sizeChromosome);
end
