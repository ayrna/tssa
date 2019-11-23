%   COPYRIGHT
%   This file is part of TSSA: https://github.com/ayrna/tssa
%   Original authors: Antonio M. Duran Rosal, Pedro A. Gutierrez
%   Copyright:
%       This software is released under the The GNU General Public License v3.0 licence
%       available at http://www.gnu.org/licenses/gpl-3.0.html
%   Citation: If you use this code, please cite the following paper:
%     [1] A.M. Durán-Rosal, D. Guijo-Rubio, P.A. Gutiérrez and C. Hervás-Martínez.
%		  "Hybrid Weighted Barebones Exploiting Particle Swarm Optimization Algorithm
%		  for Time Series Representation". BIOMA2018. 16th-18th May. 2018. 
%		  Paris (France). LNCS, vol. 10835. pp. 126-137
%		  https://doi.org/10.1007/978-3-319-91641-5_11	
%
%% updatePositionsWBB
% Function: Make an iteration of the weighted barebones PSO algorithm
% 
% Input:
%     currentPopulationInt:     population
%     bestLocalPopulationInt:   best local positions
%     acumulateBestLocal:       weighted acumulated best local positions
%     acumulateBestGlobal:      weighted acumulated best global position
%     oldFitness:               current fitness
%     generation:               denominator of the weigthed acumulation
%     nPobl:                    population size
%     sizeChromosome:           binary chromosome length
%     sizeChromosomeInt:        integer chromosome length
%     
% Output:
%     newCurrentPopulationInt:  new integer population
%     newFitness:               new fitness
function [newCurrentPopulationInt,newFitness] = updatePositionsWBB(currentPopulationInt,bestLocalPopulationInt,acumulateBestLocal,acumulateBestGlobal,oldFitness,generation,nPobl,sizeChromosome,sizeChromosomeInt)
    
    newCurrentPopulationInt = currentPopulationInt;
    for i=1:nPobl,
        for j=1:sizeChromosomeInt,
            if rand > 0.5,
                mediaI = acumulateBestLocal(i,j)/generation;
                mediaG = acumulateBestGlobal(j)/generation;
                a = (mediaI+mediaG)/2;
                b = abs(mediaI-mediaG);
                newCurrentPopulationInt(i,j) = normrnd(a,b);
            else
                newCurrentPopulationInt(i,j) = bestLocalPopulationInt(i,j);
            end
        end
    end
    
    newCurrentPopulationInt = sort(newCurrentPopulationInt,2,'ascend');
    [newCurrentPopulationInt, newFitness] = checkPopulation(newCurrentPopulationInt,currentPopulationInt,oldFitness,nPobl,sizeChromosome);
end
