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
%% metrics
% Function: Values of the metrics for a given segment
% 
% Input:
%     segment:            time series values of the segment
%     characActivation:   flag array to decide which characteristics are used
%     degree:             degree of approximation
%     typeError:          type of error (RMSE, RMSEp, MAXe) (see function computeErrors)
%     
% Output:
%     values:             array of characteristics of the segment
function [values] = metrics(segment,characActivation,degree,typeError)
    % To extract statistic metrics
    values = zeros(1,(sum(characActivation(1:end-2))+degree*characActivation(end-1)+characActivation(end)));
    X = 1:numel(segment);
    X = transpose(X);
    c = (1/numel(segment));
    m = sum(segment)*c;
    a = (segment-m);
    varsegment = c * sum(a.^2);
    s = sqrt(varsegment);

    if s==0,
        counterValues = 1;
        if characActivation(1)==1,
            values(counterValues) = 0; %Variance
            counterValues = counterValues + 1;
        end
        if characActivation(2)==1,
            values(counterValues) = 0; %Skewness
            counterValues = counterValues + 1;
        end
        if characActivation(3)==1,
            values(counterValues) = -3; %Kurtosis
            counterValues = counterValues + 1;
        end
        if characActivation(4)==1,
            values(counterValues) = 0; %Autocorrelation
            counterValues = counterValues + 1;
        end
        if characActivation(5)==1,
            for i=1:degree,
                values(counterValues) = 0;
                counterValues = counterValues + 1;
            end
        end
        if characActivation(6)==1,
            values(counterValues) = 0; %Error
        end
            
    else
        counterValues = 1;
        if characActivation(1)==1,
            values(counterValues) = varsegment; %Variance
            counterValues = counterValues + 1;
        end
        if characActivation(2)==1,
            values(counterValues) = c*(sum(a.^3)/(s.^3)); %Skewness
            counterValues = counterValues + 1;
        end
        if characActivation(3)==1,
            values(counterValues) = c*(sum(a.^4)/varsegment.^2) - 3; %Kurtosis
            counterValues = counterValues + 1;
        end
        if characActivation(4)==1,
            values(counterValues) = sum((segment(1:end-1) - m) .* (segment(2:end) - m))/varsegment; %Autocorrelation
            counterValues = counterValues + 1;
        end
        p = polyfit(X,segment,degree);
        if characActivation(5)==1,
            for i=1:degree,
                values(counterValues) = p(i);
                counterValues = counterValues + 1;
            end
        end
        if characActivation(6)==1,
            values(counterValues) = 0; %Error
            % Error
            % estimated(:,1)= p(1)*X(:,1).*X(:,1) + p(2)*X(:,1) + p(3);
            estimated(:,1) = polyval(p,X(:,1));
            error = estimated(:,1) - segment(:,1);
            error = error.*error;
            if typeError == 1,
                N=numel(error);
                error=sum(error);
                values(counterValues)=error/N; %MSE
            elseif typeError == 2,
                values(counterValues)=sum(error); %SSE
            else
                values(counterValues)=max(error); %MAXe
            end
        end 
    end
end
