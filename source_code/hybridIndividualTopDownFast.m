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
%% hybridIndividualTopDownFast
% Function: Add some cut points in a chromosome with Top Down algorithm
% 
% Input:
%     chromosome: segmentation
%     max_iter:   number of points to be added
%     x:          time indexes (horizontal)
%     y:          time indexes (vertical)
%     matrix:     precomputed matrix distances
%     
% Output:
%     chromosome: hybridized chromosome
function [chromosome] = hybridIndividualTopDownFast(chromosome,max_iter,x,y,matrix)
    chromosome(1,1)=0;
    indexes = find(chromosome==1);
    gain_cost = zeros(1,numel(chromosome));

    base = returnOneErrorFast(1,indexes(1),x,y,matrix);
    for j=2:indexes(1),
        gain_cost(1,j) = base - (returnOneErrorFast(1,j,x,y,matrix) + returnOneErrorFast(j,indexes(1),x,y,matrix));
    end
    for i=1:numel(indexes)-1,
        base = returnOneErrorFast(indexes(i),indexes(i+1),x,y,matrix);
        for j=indexes(i):indexes(i+1),
            gain_cost(1,j) = base - (returnOneErrorFast(indexes(i),j,x,y,matrix) + returnOneErrorFast(j,indexes(i+1),x,y,matrix));
        end
    end
    base = returnOneErrorFast(indexes(end),numel(chromosome),x,y,matrix);
    for j=indexes(end):numel(chromosome),
        gain_cost(1,j) = base - (returnOneErrorFast(indexes(end),j,x,y,matrix) + returnOneErrorFast(j,numel(chromosome),x,y,matrix));
    end
    gain_cost(1,indexes)=-Inf;
    gain_cost(1,1)=-Inf;
    gain_cost(1,end)=-Inf;

    iterations = 0;

    while iterations < max_iter,
        ind = find(gain_cost==max(gain_cost),1);
        chromosome(1,ind)=1;
        left = find(chromosome(1:ind-1)==1,1,'last');
        if isempty(left),
            left=1;
        end
        right = find(chromosome(ind+1:end)==1,1,'first');
        if isempty(right),
            right=numel(chromosome);
        else
            right=right+ind;
        end
        base = returnOneErrorFast(left,ind,x,y,matrix);
        for j=left:ind,
            gain_cost(1,j)=base - (returnOneErrorFast(left,j,x,y,matrix)+returnOneErrorFast(j,ind,x,y,matrix));
        end
        base = returnOneErrorFast(ind,right,x,y,matrix);
        for j=ind:right,
            gain_cost(1,j)=base - (returnOneErrorFast(ind,j,x,y,matrix)+returnOneErrorFast(j,right,x,y,matrix));
        end
        gain_cost(1,ind)=-Inf;
        gain_cost(1,left)=-Inf;
        gain_cost(1,right)=-Inf;
        iterations=iterations+1;
    end          
end
