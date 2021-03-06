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
%% errorsSegment
% Function: Error values for a given segment
% 
% Input:
%     segment:    time series values of the segment
%     degree:     degree of approximation
%     
% Output:
%     errorMSE:   Mean Squared Error of the segment
%     errorSSE:   Sum of Squared Error of the segment
%     errorMAXe:  Maximum Error of the segment
function [errorMSE, errorSSE, errorMAXe] = errorsSegment(segment,degree)
    if numel(segment)>1,
        X=1:numel(segment);
        X=transpose(X);
        if degree == 0,
            estimated(:,1)=interp1([X(1,1) X(end,1)], [segment(1,1) segment(end,1)], X(:,1));
        else
            p = polyfit(X,segment,degree);
            estimated(:,1) = polyval(p,X(:,1));
        end
        error = estimated(:,1) - segment(:,1);
        error = error.*error;
        %RMSE
        N=numel(error);
        errorMSE=sum(error);
        errorMSE=errorMSE/N;
        %RMSEp
        errorSSE=sum(error);
        %MAXe
        errorMAXe=max(error);
    else
        errorMSE=0;
        errorSSE=0;
        errorMAXe=0;
    end
end
