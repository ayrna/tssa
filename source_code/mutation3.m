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
%% mutation3
% Function: Cyclic mutation operator with number of segments restriction
% 
% Input:
%     population: set of segmentations
%     fitness:    fitness value of each segmentation
%     pMut:       mutation probability
%     
% Output:
%     mutatedPopulation:  population after applying mutation
%     changedFitness:     fitness of each segmentation (NaN in the case of changes)
function [mutatedPopulation,changedFitness] = mutation3(population,fitness,pMut)
    mutatedPopulation = population;
    changedFitness = fitness;
    sizeChromosome = numel(population(1,:));
    %t= cputime
    for i=1:size(population,1),
        %Mutate?
        if rand()<pMut,
            type = randi(2,1,1);
            if type == 1,  %Displacement to the right
                attempts=0;
                while attempts<=2;
                    displacement = randi(sizeChromosome-1,1,1);
                    mutatedPopulation(i,1+displacement:end)=population(i,1:end-displacement);
                    mutatedPopulation(i,1:displacement)=population(i,end-displacement+1:end);

                    if mutatedPopulation(i,1)==1 || mutatedPopulation(i,end)==1,
                        attempts=attempts+1;
                        mutatedPopulation(i,:)=population(i,:);
                    else
                        attempts=10;
                        changedFitness(i) = NaN;
                    end
                end
            else  %Displacement to the left
                attempts=0;
                while attempts<=2;
                    displacement = randi(sizeChromosome-1,1,1);
                    mutatedPopulation(i,1:end-displacement)=population(i,displacement+1:end);
                    mutatedPopulation(i,end-displacement+1:end)=population(i,1:displacement);

                    if mutatedPopulation(i,1)==1 || mutatedPopulation(i,end)==1,
                        attempts=attempts+1;
                        mutatedPopulation(i,:)=population(i,:);
                    else
                        attempts=10;
                        changedFitness(i) = NaN;
                    end
                end
            end
        end
    end
end
