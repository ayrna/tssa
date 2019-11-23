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
%% estimationSegment
% Function: Estimating segment
% 
% Input:
%     segment:    time series values of the segment
%     degree:     degree of approximation for regression (0 for interpolation)
%     
% Output:
%     yEstimated: approximated segment
function [yEstimated] = estimationSegment(segment,degree)
    yEstimated = zeros(numel(segment),1);
    X=1:numel(segment);
    X=transpose(X);
    if degree == 0,
        yEstimated(:,1)=interp1([X(1,1) X(end,1)], [segment(1,1) segment(end,1)], X(:,1));
    else
        p = polyfit(X,segment,degree);
        yEstimated(:,1) = polyval(p,X(:,1));
        % yEstimated(:,1) = p(1)*X(:,1).*X(:,1) + p(2)*X(:,1) + p(3);
    end
end
