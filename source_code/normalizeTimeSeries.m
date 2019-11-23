%   COPYRIGHT
%   This file is part of TSSA: https://github.com/ayrna/tssa
%   Original authors: Antonio M. Duran Rosal, Pedro A. Gutierrez
%   Copyright:
%       This software is released under the The GNU General Public License v3.0 licence
%       available at http://www.gnu.org/licenses/gpl-3.0.html
%   Citation: If you use this code, please cite any of the following papers:
%     [1] A.M. Durán-Rosal, J.C. Fernández, P.A. Gutiérrez and C. Hervás-Martínez.
%         "Detection and prediction of segments containing extreme significant wave heights"
%         Ocean Engineering, Vol. 142, September, 2017, pp. 268-279.
%         https://doi.org/10.1016/j.oceaneng.2017.07.009
%     [2] A.M. Durán-Rosal, M. de la Paz Marín, P.A. Gutiérrez and C. Hervás-Martínez.
%         "Identifying market behaviours using European Stock Index time series by 
%         a hybrid segmentation algorithm", Neural Processing Letters,
%         Vol. 46, December, 2017, pp. 767–790.
%         https://doi.org/10.1007/s11063-017-9592-8
%
%% normalizeTimeSeries
% Function: Scales time series between intervalLeft and intervalRight
% 
% Input:
%     serie:          time series
%     intervalLeft:   left interval of the scaled
%     intervalRight:  right interval of the scaled
%     
% Output:
%     scaled_serie:   scaled time series
%     minimum:        mimimum value of the input time series
%     maximum:        maximum value of the input time series
function [scaled_serie, minimum, maximum] = normalizeTimeSeries(serie,intervalLeft,intervalRight)
   minimum = min(serie);
   maximum = max(serie);
   scaled_serie=(serie - minimo)/(maximo-minimo);
   scaled_serie=(scaled_serie*(intervalRight-intervalLeft))+intervalLeft;
end
