%   COPYRIGHT
%   This file is part of TSSA: https://github.com/ayrna/tssa
%   Original authors: Antonio M. Duran Rosal, Pedro A. Gutierrez
%   Copyright:
%       This software is released under the The GNU General Public License v3.0 licence
%       available at http://www.gnu.org/licenses/gpl-3.0.html
%   Citation: If you use this code, please cite any of the following papers:
%     [1] M. Pérez-Ortiz, A.M. Durán-Rosal, P.A. Gutiérrez, et al.
%         "On the use of evolutionary time series analysis for segmenting paleoclimate data"
%         Neurocomputing, Vol. 326-327, January, 2019, pp. 3-14
%         https://doi.org/10.1016/j.neucom.2016.11.101
%     [2] A.M. Durán-Rosal, J.C. Fernández, P.A. Gutiérrez and C. Hervás-Martínez.
%         "Detection and prediction of segments containing extreme significant wave heights"
%         Ocean Engineering, Vol. 142, September, 2017, pp. 268-279.
%         https://doi.org/10.1016/j.oceaneng.2017.07.009
%     [3] A.M. Durán-Rosal, M. de la Paz Marín, P.A. Gutiérrez and C. Hervás-Martínez.
%         "Identifying market behaviours using European Stock Index time series by 
%         a hybrid segmentation algorithm", Neural Processing Letters,
%         Vol. 46, December, 2017, pp. 767–790.
%         https://doi.org/10.1007/s11063-017-9592-8
%
%% selection1Roulette
% Function: Performance a selection by a probabilistic roulette wheel
% 
% Input:
%     population: set of chromosomes
%     fitness:    fitness of each individual
%     nPobl:      population size
%     
% Output:
%     newPopulation:  selected population
%     newFitness:     fitness of the new population
function [newPopulation,newFitness] = selection1Roulette(population,fitness,nPobl)
    [fbest,indBestSegmentation] = max(fitness);
    newPopulation(1,:) = population(indBestSegmentation,:);
    newFitness = zeros(1,nPobl);
    newFitness(1) = fbest;
    cumFitness = cumsum(fitness);

    randCums = cumFitness(end).*rand(1,nPobl);
    for i=2:nPobl,
        % Roulette
        index1 = find((cumFitness > randCums(i))==1);
        newPopulation(i,:) = population(index1(1),:);
        newFitness(i) = fitness(index1(1));
    end
end
