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
%     [2] M. Pérez-Ortiz, A.M. Durán-Rosal, P.A. Gutiérrez, et al.
%         "On the use of evolutionary time series analysis for segmenting paleoclimate data"
%         Neurocomputing, Vol. 326-327, January, 2019, pp. 3-14
%         https://doi.org/10.1016/j.neucom.2016.11.101
%     [3] A.M. Durán-Rosal, P.A. Gutiérrez, F.J. Martínez-Estudillo and C. Hervás-Martínez.
%         "Simultaneous optimisation of clustering quality and approximation error
%         for time series segmentation", Information Sciences, Vol. 442-443, May, 2018, pp. 186-201.
%         https://doi.org/10.1016/j.ins.2018.02.041
%     [4] A.M. Durán-Rosal, P.A. Gutiérrez, S. Salcedo-Sanz and C. Hervás-Martínez.
%         "A statistically-driven Coral Reef Optimization algorithm for optimal
%         size reduction of time series", Applied Soft Computing, 
%         Vol. 63. 2018, pp. 139-153.
%         https://doi.org/10.1016/j.asoc.2017.11.037
%     [5] A.M. Durán-Rosal, J.C. Fernández, P.A. Gutiérrez and C. Hervás-Martínez.
%         "Detection and prediction of segments containing extreme significant wave heights"
%         Ocean Engineering, Vol. 142, September, 2017, pp. 268-279.
%         https://doi.org/10.1016/j.oceaneng.2017.07.009
%     [6] A.M. Durán-Rosal, M. de la Paz Marín, P.A. Gutiérrez and C. Hervás-Martínez.
%         "Identifying market behaviours using European Stock Index time series by 
%         a hybrid segmentation algorithm", Neural Processing Letters,
%         Vol. 46, December, 2017, pp. 767–790.
%         https://doi.org/10.1007/s11063-017-9592-8
%
%% computeErrors
% Function: Error values for a segmentation
% 
% Input:
%     individual: chromosome (segmentation)
%     serie:      time series
%     degree:     degree of approximation
%     
% Output:
%     errors:     array formed by [RMSE, RMSEp, MAXe]
function [errors] = computeErrors(individual,serie,degree)
    ind = find(individual==1);
    nOfSegments = size(ind,2);
    errors = zeros(1,3);
    errorsRMSE = zeros(nOfSegments+1,1);
    errorsRMSEp = zeros(nOfSegments+1,1);
    errorsMAXe = zeros(nOfSegments+1,1);
    [errorsRMSE(1,1), errorsRMSEp(1,1), errorsMAXe(1,1)] = errorsSegment(serie(1:ind(1)),degree);

    % Number of Segments - 1
    for j=1:nOfSegments-1,
        [errorsRMSE(j+1,1), errorsRMSEp(j+1,1), errorsMAXe(j+1,1) ] = errorsSegment(serie(ind(j):ind(j+1)),degree);
    end
    [errorsRMSE(end,1), errorsRMSEp(end,1), errorsMAXe(end,1)] = errorsSegment(serie(ind(end):end),degree);
    
    MSE=mean(errorsRMSE);
    RMSE=sqrt(MSE);

    MSEp=sum(errorsRMSEp)/numel(individual);
    RMSEp=sqrt(MSEp);

    MAXe=max(errorsMAXe);
    MAXe=sqrt(MAXe);
    
    errors(1)=RMSE;
    errors(2)=RMSEp;
    errors(3)=MAXe;
end
