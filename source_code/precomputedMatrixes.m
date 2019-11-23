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
%     [2] A.M. Durán-Rosal, P.A. Gutiérrez, S. Salcedo-Sanz and C. Hervás-Martínez.
%         "A statistically-driven Coral Reef Optimization algorithm for optimal
%         size reduction of time series", Applied Soft Computing, 
%         Vol. 63. 2018, pp. 139-153.
%         https://doi.org/10.1016/j.asoc.2017.11.037
%
%% precomputedMatrixes
% Function: Precomputes the matrix of distances
% 
% Input:
%     x:                  time indexes (horizontal)
%     y:                  time indexes (vertical)
%     sizeChromosome:     chromosome length
%     
% Output:
%     matrix:             precomputed matrix distances
function [matrix] = precomputedMatrixes(x,y,sizeChromosome)  
    matrix.suma_x = zeros(sizeChromosome,sizeChromosome);
    matrix.suma_y = zeros(sizeChromosome,sizeChromosome);
    matrix.suma_xx = zeros(sizeChromosome,sizeChromosome);
    matrix.suma_yy = zeros(sizeChromosome,sizeChromosome);
    matrix.suma_xy = zeros(sizeChromosome,sizeChromosome);

    % Sumas acumuladas desde el primer punto a cualquiera
    matrix.suma_x(1,:) = cumsum(x);
    matrix.suma_y(1,:) = cumsum(y);
    matrix.suma_xx(1,:) = cumsum(x.^2);
    matrix.suma_yy(1,:) = cumsum(y.^2);
    matrix.suma_xy(1,:) = cumsum(x.*y);

    % Para el segundo punto, se le resta la del primero
    % Para el tercero, la del segundo
    % etc...
    for i=2:(sizeChromosome-1)
        matrix.suma_x(i,(i+1):end) = matrix.suma_x(i-1,(i+1):end) - x(i-1);
        matrix.suma_y(i,(i+1):end) = matrix.suma_y(i-1,(i+1):end) - y(i-1);
        matrix.suma_xx(i,(i+1):end) = matrix.suma_xx(i-1,(i+1):end) - (x(i-1)).^2;
        matrix.suma_yy(i,(i+1):end) = matrix.suma_yy(i-1,(i+1):end)- (y(i-1)).^2;
        matrix.suma_xy(i,(i+1):end) = matrix.suma_xy(i-1,(i+1):end)- (x(i-1)*y(i-1));
    end
end
