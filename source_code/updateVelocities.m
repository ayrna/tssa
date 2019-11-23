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
%
%% updateVelocities
% Function: update the velocities of a particle swarm
% 
% Input:
%     velocities:               current velocities of the particles
%     currentPopulationInt:     population
%     bestLocalPopulationInt:   best local positions
%     bestIndividualInt:        best global position
%     nPobl:                    population size
%     sizeChromosomeInt:        integer chromosome length
%     W:                        weight param
%     C1,C2:                    acceleration constants
%     
% Output:
%     newVelocities:            new velocities of the particles
function [newVelocities] = updateVelocities(velocities,currentPopulationInt,bestLocalPopulationInt,bestIndividualInt,nPobl,sizeChromosomeInt,W,C1,C2)
    newVelocities = velocities;
    ro1 = unifrnd(0,1,[1 sizeChromosomeInt]);
    ro2 = unifrnd(0,1,[1 sizeChromosomeInt]);

    for i=1:nPobl,
        newVelocities(i,:) = W * velocities(i,:) + C1 * ro1 .* (bestLocalPopulationInt(i,:) - currentPopulationInt(i,:)) + C2 * ro2 .* (bestIndividualInt(1,:) - currentPopulationInt(i,:));
    end
end 
