%   COPYRIGHT
%   This file is part of TSSA: https://github.com/ayrna/tssa
%   Original authors: Antonio M. Duran Rosal, Pedro A. Gutierrez
%   Copyright:
%       This software is released under the The GNU General Public License v3.0 licence
%       available at http://www.gnu.org/licenses/gpl-3.0.html
%   Citation: If you use this code, please cite the following paper:
%     [1] A.M. Durán-Rosal, P.A. Gutiérrez, F.J. Martínez-Estudillo and C. Hervás-Martínez.
%         "Simultaneous optimisation of clustering quality and approximation error
%         for time series segmentation", Information Sciences, Vol. 442-443, May, 2018, pp. 186-201.
%         https://doi.org/10.1016/j.ins.2018.02.041
%
%% selection2GMO
% Function: Selection of the best individuals situated in the Pareto Front. Only for GMO.
% 
% Input:
%     population:       set of chromosomes
%     fitness1:         fitness of the first objective of each individual
%     fitness2:         fitness of the second objective of each individual
%     nPobl:            population size
%     
% Output:
%     newPopulation:  selected population
%     newFitness1:    fitness of the first objective of the new population
%     newFintess2:    fitness of the second objective of the new population (NaN in case of moobjective algorithmm)
function [newPopulation,newFitness1,newFitness2] = selection2GMO(population,fitness1,fitness2,nPobl)
    sizeChromosome = numel(population(1,:));
    newFitness1 = zeros(1,nPobl);
    newFitness2 = zeros(1,nPobl);
    newPopulation = zeros(nPobl,sizeChromosome);

    newPopulation(1:nPobl,:)=population(1:nPobl,:);
    newFitness1(1,1:nPobl)=fitness1(1,1:nPobl);
    newFitness2(1,1:nPobl)=fitness2(1,1:nPobl);
end 
