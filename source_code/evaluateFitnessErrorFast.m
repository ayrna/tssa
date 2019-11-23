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
%     [2] A.M. Durán-Rosal, P.A. Gutiérrez, F.J. Martínez-Estudillo and C. Hervás-Martínez.
%         "Simultaneous optimisation of clustering quality and approximation error
%         for time series segmentation", Information Sciences, Vol. 442-443, May, 2018, pp. 186-201.
%         https://doi.org/10.1016/j.ins.2018.02.041
%     [3] A.M. Durán-Rosal, P.A. Gutiérrez, S. Salcedo-Sanz and C. Hervás-Martínez.
%         "A statistically-driven Coral Reef Optimization algorithm for optimal
%         size reduction of time series", Applied Soft Computing, 
%         Vol. 63. 2018, pp. 139-153.
%         https://doi.org/10.1016/j.asoc.2017.11.037
%
%% evaluateFitnessErrorFast
% Function: Evaluation method for approximation error
% 
% Input:
%     population:         set of chromosomes
%     oldFitness:         fitness of population
%     nPobl:              population size
%     sizeChromosome:     chromosome length
%     x:                  time indexes (horizontal)
%     y:                  time indexes (vertical)
%     matrix:             precomputed matrix distances
%     
% Output:
%     fitness:            fitness of current population
function [fitness] = evaluateFitnessErrorFast(population,oldFitness,nPobl,sizeChromosome,x,y,matrix)
    fitness = zeros(1,size(population,1));
    errors = zeros(1,1+sum(population(:,1)));
    %fitness = oldFitness;
    for i=1:nPobl,
        if isnan(oldFitness(i)),
            ind=find(population(i,:)==1);
            lefts = [1 ind];
            rights = [ind sizeChromosome];
            for j=1:numel(lefts),
                left = lefts(j);
                right = rights(j);
                %if (right - left) > 1,
                    a = (y(right) - y(left))/(x(right)-x(left));
                    b= y(right) - (a * x(right));
                    % errors j stores the quatratic sum of the errors of the j-segment
                    errors(j) = ( (a^2 * matrix.suma_xx(left,right)) + ((right-left+1) * b^2) + matrix.suma_yy(left,right) + (2 * a * matrix.suma_x(left,right) * b) - (2 * a * matrix.suma_xy(left,right))  - (2 * b * matrix.suma_y(left,right)) ); 
                %end
            end
            fitness(i) = 1/(1+sqrt(sum(errors)/sizeChromosome));
        elseif oldFitness(i) == -1,
            fitness(i) = -1;
        end 
    end
end
