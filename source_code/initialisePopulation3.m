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
%     [2] A.M. Durán-Rosal, D. Guijo-Rubio, P.A. Gutiérrez, S. Salcedo-Sanz and C. Hervás-Martínez.
%         "A coral reef optimization algorithm for wave height time series segmentation problems".
%         International Work-Conference on Artificial and Natural Neural Networks (IWANN2017).
%         14th-16th June. 2017. Cádiz (Spain). LNCS, vol. 10305. pp. 673-684
%         https://doi.org/10.1007/978-3-319-59153-7_58
%
%% initialisePopulation3
% Function: Initialisation with free positions
% 
% Input:
%     nPobl:          population size
%     sizeChromosome: chromosome size
%     minSeg:         mimimum segment size
%     maxSeg:         maximum segment size
%     freePositions:  number of free positions
% 
% Output:
%     newPopulation:  initialised population
function [newPopulation] = initialisePopulation3(nPobl,sizeChromosome,minSeg,maxSeg,freePositions)
    newPopulation = ones(nPobl,sizeChromosome) * -1;
    ind = noRepetitionRand(1,nPobl,nPobl-freePositions);
    
    for i=1:numel(ind),
        cont = 0;
        newPopulation(ind(i),:) = 0;
        while (sizeChromosome - cont) > maxSeg
            cont = cont + randi([minSeg,maxSeg]);
            newPopulation(ind(i),cont)=1;
        end
    end
end
