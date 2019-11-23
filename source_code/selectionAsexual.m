%   COPYRIGHT
%   This file is part of TSSA: https://github.com/ayrna/tssa
%   Original authors: Antonio M. Duran Rosal, Pedro A. Gutierrez
%   Copyright:
%       This software is released under the The GNU General Public License v3.0 licence
%       available at http://www.gnu.org/licenses/gpl-3.0.html
%   Citation: If you use this code, please cite any of the following papers:
%     [1] A.M. Durán-Rosal, P.A. Gutiérrez, S. Salcedo-Sanz and C. Hervás-Martínez.
%         "A statistically-driven Coral Reef Optimization algorithm for optimal
%         size reduction of time series", Applied Soft Computing, 
%         Vol. 63. 2018, pp. 139-153.
%         https://doi.org/10.1016/j.asoc.2017.11.037
%     [2] A.M. Durán-Rosal, P.A. Gutiérrez, S. Salcedo-Sanz and C. Hervás-Martínez.
%         "Dynamical Memetization in Coral Reef Optimization Algorithms
%         for Optimal Time Series Approximation",
%         Progress in Artificial Intelligence, Vol. 8, June, 2019, pp. 253-262.
%         https://doi.org/10.1007/s13748-019-00176-0
%     [3] A.M. Durán-Rosal, D. Guijo-Rubio, P.A. Gutiérrez, S. Salcedo-Sanz and C. Hervás-Martínez.
%         "A coral reef optimization algorithm for wave height time series segmentation problems".
%         International Work-Conference on Artificial and Natural Neural Networks (IWANN2017).
%         14th-16th June. 2017. Cádiz (Spain). LNCS, vol. 10305. pp. 673-684
%         https://doi.org/10.1007/978-3-319-59153-7_58
%
%% Asexual selection
% Function: Performance a selection by an asexual procedure
% 
% Input:
%     population: set of chromosomes
%     fitness:    fitness of each individual
%     nPobl:      population size
%     Fa:         percentage of asexual reproduction (selection)
%     
% Output:
%     newPopulation:  selected population
%     newFitness:     fitness of the new population
function [newPopulation, newFitness1] = selectionAsexual(population,fitness1,Fa)
    [sortedFitness, sortedIndexes] = sort(fitness1,'descend');
    minValue = 1;
    maxValue = round(Fa*numel(find(fitness1~=-1)));
    if maxValue ~= 0,
        indx = randi([minValue maxValue],1);
    else
        indx = 1;
    end
    newPopulation = population(sortedIndexes(indx),:);
    newFitness1 = fitness1(sortedIndexes(indx));
        
end
