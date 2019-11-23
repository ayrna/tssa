%   COPYRIGHT
%   This file is part of TSSA: https://github.com/ayrna/tssa
%   Original authors: Antonio M. Duran Rosal, Pedro A. Gutierrez
%   Copyright:
%       This software is released under the The GNU General Public License v3.0 licence
%       available at http://www.gnu.org/licenses/gpl-3.0.html
%   Citation: If you use this code, please cite any of the following papers:
%     [1] A.M. Durán-Rosal, P.A. Gutiérrez, F.J. Martínez-Estudillo and C. Hervás-Martínez.
%         "Simultaneous optimisation of clustering quality and approximation error
%         for time series segmentation", Information Sciences, Vol. 442-443, May, 2018, pp. 186-201.
%         https://doi.org/10.1016/j.ins.2018.02.041
%
%% evaluateBothFitness
% Function: Evaluation method for both fitness
% 
% Input:
%     typeError:          type of error (RMSE, RMSEp, MAXe) (see function computeErrors)
%     typeClustMeasure:   type of clustering fitness measure (see function fitnessF)
%     population:         set of chromosomes
%     oldFitness1:        clustering fitness of population
%     oldFitness2:        error fitness of population
%     k:                  number of clusters
%     iterClust:          maximum number of iteration for k-means (see function clustering)
%     serie:              time series
%     characActivation:   flag array to decide which characteristics are used
%     degree:             degree of approximation
%     
% Output:
%     fitness1:            clustering fitness of current population
%     fitness2:            error fitness of current population
function [fitness1,fitness2] = evaluateBothFitness(typeError,typeClustMeasure,population,oldFitness1,oldFitness2,k,iterClust,serie,characActivation,degree)
    fitness1 = zeros(1,size(population,1));
    fitness2 = fitness1;
    %fitness = oldFitness;
    for i=1:size(population,1),
        if isnan(oldFitness1(i)),
            [charac] = computeMetrics(population(i,:),serie,characActivation,degree,typeError);
            [normCharac] = normalizeFunction(charac);
            %Clustering
            [assignation,centroids] = clusteringKmeans(normCharac,k,iterClust);
            fitness1(i) = fitnessF(typeClustMeasure,normCharac,assignation,centroids,k);
            %Errors
            [errors] = computeErrors(population(i,:),serie,degree);
            fitness2(i) = fitnessError(errors,typeError);
        elseif oldFitness1(i) == -1,
            fitness1(i) = -1;
            fitness2(i) = -1;
        end   
    end
end
