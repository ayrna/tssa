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
%     [2] A.M. Durán-Rosal, P.A. Gutiérrez, F.J. Martínez-Estudillo and C. Hervás-Martínez.
%         "Simultaneous optimisation of clustering quality and approximation error
%         for time series segmentation", Information Sciences, Vol. 442-443, May, 2018, pp. 186-201.
%         https://doi.org/10.1016/j.ins.2018.02.041
%     [3] A.M. Durán-Rosal, J.C. Fernández, P.A. Gutiérrez and C. Hervás-Martínez.
%         "Detection and prediction of segments containing extreme significant wave heights"
%         Ocean Engineering, Vol. 142, September, 2017, pp. 268-279.
%         https://doi.org/10.1016/j.oceaneng.2017.07.009
%     [4] A.M. Durán-Rosal, M. de la Paz Marín, P.A. Gutiérrez and C. Hervás-Martínez.
%         "Identifying market behaviours using European Stock Index time series by 
%         a hybrid segmentation algorithm", Neural Processing Letters,
%         Vol. 46, December, 2017, pp. 767–790.
%         https://doi.org/10.1007/s11063-017-9592-8
%
%% evaluateFitnessClusteringKmeans
% Function: Evaluation method for clustering Kmeans
% 
% Input:
%     typeError:          type of error (RMSE, RMSEp, MAXe) (see function computeErrors)
%     typeClustMeasure:   type of clustering fitness measure (see function fitnessF)
%     population:         set of chromosomes
%     oldFitness:         fitness of population
%     k:                  number of clusters
%     iterClust:          maximum number of iteration for k-means (see function clustering)
%     serie:              time series
%     characActivation:   flag array to decide which characteristics are used
%     degree:             degree of approximation
%     
% Output:
%     fitness:            fitness of current population
function [fitness] = evaluateFitnessClusteringKmeans(typeError,typeClustMeasure,population,oldFitness,k,iterClust,serie,characActivation,degree)
    fitness = zeros(1,size(population,1));
    %fitness = oldFitness;
    for i=1:size(population,1),
        if isnan(oldFitness(i)),
            [charac] = computeMetrics(population(i,:),serie,characActivation,degree,typeError);
            [normCharac] = normalizeFunction(charac);
            %Clustering
            [assignation,centroids] = clusteringKmeans(normCharac,k,iterClust);
            fitness(i) = fitnessF(typeClustMeasure,normCharac,assignation,centroids,k);
        elseif oldFitness(i) == -1,
            fitness(i) = -1;
        end   
    end
end
