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
%     [2] A.M. Durán-Rosal, P.A. Gutiérrez, S. Salcedo-Sanz and C. Hervás-Martínez.
%         "Dynamical Memetization in Coral Reef Optimization Algorithms
%         for Optimal Time Series Approximation",
%         Progress in Artificial Intelligence, Vol. 8, June, 2019, pp. 253-262.
%         https://doi.org/10.1007/s13748-019-00176-0
%
%% initialisePopulation4
% Function: Initialisation with free positions
% 
% Input:
%     nPobl:          population size
%     sizeChromosome: chromosome size
%     numSeg:         number of segments
%     freePositions:  number of free positions
% 
% Output:
%     newPopulation:  initialised population
function [newPopulation] = initialisePopulation4(nPobl,sizeChromosome,numSeg,freePositions)
    newPopulation = ones(nPobl,sizeChromosome) * -1;
    ind = noRepetitionRand(1,nPobl,nPobl-freePositions);
    
    for i=1:numel(ind),
        ind2=noRepetitionRand(2,sizeChromosome-1,numSeg-1);
        newPopulation(ind(i),:)=0;
        newPopulation(ind(i),ind2)=1;
    end
end
