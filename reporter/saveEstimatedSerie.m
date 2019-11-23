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
%% saveEstimatedSerie
% Function: Save the approximation time series produced by a segmentation
% 
% Input:
%     model:      contains the information necessary for the reporter
%     dataset:    name of the dataset
%     repsuffix:  path of the output file
%     
% Output:
%     No output variables. Only a file which contains the approximated time series
function saveEstimatedSerie(model, dataset, repsuffix)   
    X(:,1)=transpose(1:numel(model.estimatedSerie(:,1)));
    X(:,2)=model.estimatedSerie(:,1);

    outputFile = [repsuffix filesep dataset];
    % Estimated Time Series
    f = fopen([outputFile '_estimated.txt'], 'wt');
    for i=1:numel(X(:,1)),
        fprintf(f, '%d\t%f\n', X(i,1), X(i,2)); 
    end
    fclose(f);
end
