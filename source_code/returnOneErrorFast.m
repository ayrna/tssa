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
%     [4] A.M. Durán-Rosal, P.A. Gutiérrez, S. Salcedo-Sanz and C. Hervás-Martínez.
%         "A statistically-driven Coral Reef Optimization algorithm for optimal
%         size reduction of time series", Applied Soft Computing, 
%         Vol. 63. 2018, pp. 139-153.
%         https://doi.org/10.1016/j.asoc.2017.11.037
%
%% returnOneErrorFast
% Function: A specific error value for a given segment
% 
% Input:
%     left:       left index of the segment
%     right:      right index of the segment
%     x:          time indexes (horizontal)
%     y:          time indexes (vertical)
%     matrix:     precomputed matrix distances
% Output:
%     error:      error of the segment
function [error] = returnOneErrorFast(left,right,x,y,matrix)
    %if (right - left) > 1,
        a = (y(right) - y(left))/(x(right)-x(left));
        b= y(right) - (a * x(right));
        % error stores the quatratic sum of the errors of the segment delimited
        % by left and right
        error = ( (a^2 * matrix.suma_xx(left,right)) + ((right-left+1) * b^2) + matrix.suma_yy(left,right) + (2 * a * matrix.suma_x(left,right) * b) - (2 * a * matrix.suma_xy(left,right))  - (2 * b * matrix.suma_y(left,right)) ); 
    %
end
