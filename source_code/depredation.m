%   COPYRIGHT
%   This file is part of TSSA: https://github.com/ayrna/tssa
%   Original authors: Antonio M. Duran Rosal, Pedro A. Gutierrez
%   Copyright:
%       This software is released under the The GNU General Public License v3.0 licence
%       available at http://www.gnu.org/licenses/gpl-3.0.html
%   Citation: If you use this code, please cite any the following paper:
%     [1] A.M. Durán-Rosal, P.A. Gutiérrez, S. Salcedo-Sanz and C. Hervás-Martínez.
%         "A statistically-driven Coral Reef Optimization algorithm for optimal
%         size reduction of time series", Applied Soft Computing, 
%         Vol. 63. 2018, pp. 139-153.
%         https://doi.org/10.1016/j.asoc.2017.11.037
%
%% depredation
% Function: depredate a little percentage of a population
% 
% Input:
%     population: population to be depredated
%     fitness1:   fitness of the population
%     Fd:         percentage of the population to be depredated
%     pDep:       probability of depredation
%     
% Output:
%     newPopulation:  depredated population
%     newFitness1:    updated fitness
function [newPopulation, newFitness1] = depredation(population,fitness1,Fd,pDep)
    newPopulation = population;
    newFitness1 = fitness1;
    [sortedFitness, sortedIndexes] = sort(fitness1,'ascend');
    maxValue = round(Fd*numel(find(fitness1~=-1)));
    indStart = find(sortedFitness~=-1);
    
    for i=indStart(1):indStart(1)+maxValue-1,
        if rand()<pDep,
            newPopulation(sortedIndexes(i),:) = -1;
            newFitness1(sortedIndexes(i)) = -1;
        end
    end
end
