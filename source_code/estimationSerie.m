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
%% estimationSerie
% Function: Estimating time series
% 
% Input:
%     individual: chromosome (segmentation)
%     serie:      time series
%     degree:     degree of approximation for regression (0 for interpolation)
%     
% Output:
%     yEstimated: approximated time series
function [yEstimated] = estimationSerie(individual,serie,degree)
    ind = find(individual==1);
    nOfSegments = size(ind,2);
    yEstimated = zeros(numel(serie),1);
    yEstimated(1:ind(1))= estimationSegment(serie(1:ind(1)),degree);

    % Number of Segments - 1
    for j=1:nOfSegments-1,
        yEstimated(ind(j):ind(j+1)) = estimationSegment(serie(ind(j):ind(j+1)),degree);
    end
    yEstimated(ind(end):end) = estimationSegment(serie(ind(end):end),degree);
end
