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
%% saveCentroids
% Function: Save the centroids of a segmentation in two files (normalized and unnormalized)
% 
% Input:
%     model:      contains the information necessary for the reporter
%     dataset:    name of the dataset
%     repsuffix:  path of the output file
%     
% Output:
%     No output variables. Only two files which contain the saved centroids (normalized and unnormalized)
function saveCentroids(model, dataset, repsuffix)

    outputFile = [repsuffix filesep dataset];
    centroids = model.C;
    centroidsNorm = centroids;
    nOfClusters = numel(centroids(:,1));
    nOfFeatures = numel(centroids(1,:));
    
    %Unnormalize centroids
    for j=1:nOfFeatures,
       maximo=max(model.features(:,j));
       minimo=min(model.features(:,j));
       for i=1:nOfClusters,
           centroids(i,j)=(centroids(i,j)*(maximo-minimo))+minimo;
       end
    end
    
    % Save centroids
    f = fopen([outputFile '_centroids.csv'], 'wt');
    for i=1:nOfClusters
        fprintf(f, '%d;', i);
        for j=1:nOfFeatures-1,
            fprintf(f, '%f;', centroids(i,j));
        end
        fprintf(f, '%f\n', centroids(i,end));
    end
    fclose(f);
    
    % Save Normalised centroids
    f = fopen([outputFile '_centroidsNorm.csv'], 'wt');
    for i=1:nOfClusters
        fprintf(f, '%d;', i);
        for j=1:nOfFeatures-1,
            fprintf(f, '%f;', centroidsNorm(i,j));
        end
        fprintf(f, '%f\n', centroidsNorm(i,end));
    end
    fclose(f);
    
end
