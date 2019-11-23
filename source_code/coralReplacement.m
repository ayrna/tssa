%   COPYRIGHT
%   This file is part of TSSA: https://github.com/ayrna/tssa
%   Original authors: Antonio M. Duran Rosal, Pedro A. Gutierrez
%   Copyright:
%       This software is released under the The GNU General Public License v3.0 licence
%       available at http://www.gnu.org/licenses/gpl-3.0.html
%   Citation: If you use this code, please cite the following paper:
%     [1] A.M. Durán-Rosal, P.A. Gutiérrez, S. Salcedo-Sanz and C. Hervás-Martínez.
%         "A statistically-driven Coral Reef Optimization algorithm for optimal
%         size reduction of time series", Applied Soft Computing, 
%         Vol. 63. 2018, pp. 139-153.
%         https://doi.org/10.1016/j.asoc.2017.11.037
%
%% Coral Reef replacement
% Function: Performance a selection by a coral reef replacement
% 
% Input:
%     population: set of chromosomes
%     fitness:    fitness of each individual
%     nPobl:      population size
%     poolPopulation: solutions of a pool
%     poolFitness:    fitness of the pool population
%     Natt:           max attempts to replacement
%     
% Output:
%     newPopulation:  selected population
%     newFitness:     fitness of the new population
function [newPopulation,newFitness1] = coralReplacement(population,fitness,nPobl,poolPopulation,poolFitness,Natt)
    newPopulation = population;
    newFitness1 = fitness;
    
    for i=1:numel(poolPopulation(:,1)),
        nAttempts = 0;
        while nAttempts < Natt,
            randPosition = randi([1 nPobl],1);
            if poolFitness(i) > newFitness1(randPosition),
                newPopulation(randPosition,:) = poolPopulation(i,:);
                newFitness1(randPosition) = poolFitness(i);
                nAttempts = Natt+1;
            else
                nAttempts=nAttempts+1;
            end
        end
    end

end
