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
%% transformPopulationToBi
% Function: transform the population representation from integer to binary
% 
% Input:
%     populationInt:     integer population
%     nPobl:             population size
%     sizeChromosome:    binary chromosome length
% 
% Output:
%     population:  binary population
function [population] = transformPopulationToBi(populationInt,nPobl,sizeChromosome)
    population = false(nPobl,sizeChromosome);
    for i=1:nPobl,
        ind = populationInt(i,:);
        ind=round(ind);
        population(i,ind)=1;
    end
end
