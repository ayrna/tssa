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
%% selectionbinaryTournament
% Function: Selection of individuals performing a binary tournament. Only for GMO.
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
function [newPopulation, newFitness1, newFitness2] = selectionBinaryTournament(population,fitness1,fitness2,nPobl)
    nRows = nPobl;
    nCols = numel(population(1,:));
    newPopulation = zeros(nRows,nCols);
    newFitness1 = zeros(1, nRows);
    newFitness2 = zeros(1, nRows);

    for i=1:numel(population(:,1)),
        ind1 = randi(nRows);
        ind2 = randi(nRows);
        while ind1 == ind2,
            ind2 = randi(nRows);
        end
        if ind1 < ind2,
            newPopulation(i,:)=population(ind1,:);
            newFitness1(1,i)=fitness1(1,ind1);
            newFitness2(1,i)=fitness2(1,ind1);
        else
            newPopulation(i,:)=population(ind2,:);
            newFitness1(1,i)=fitness1(1,ind2);
            newFitness2(1,i)=fitness2(1,ind2);
        end
    end
end
