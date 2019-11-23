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
%% unNormalizeTimeSeries
% Function: Re-scales time series between intervalLeft and intervalRight to the original before being scaled
% 
% Input:
%     scaled_serie:   scaled time series
%     minimum:        mimimum value of the no scales time series
%     maximum:        maximum value of the no scales time series
%     intervalLeft:   left interval of the scaled
%     intervalRight:  right interval of the scaled
%     
% Output:
%     serie:   no scaled time series
function [serie] = unNormalizeTimeSeries(scaled_serie,minimum,maximum,intervalLeft,intervalRight)
    serie=(((scaled_serie-intervalLeft)/(intervalRight-intervalLeft))*(maximum-minimum))+maximum;
end
