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
%% hybridIndividualBottomUp
% Function: Remove some cut points in a chromosome with Bottom-Up algorithm
% 
% Input:
%     chromosome: segmentation
%     serie:      time series
%     max_iter:   number of points to be deleted
%     degree:     degree of approximation
%     typeError:  type of error (see returnOneError)
%     
% Output:
%     chromosome: hybridized chromosome
function [chromosome] = hybridIndividualBottomUp(chromosome,serie,max_iter,degree,typeError)
    chromosome(1,1)=0;
    chromosome(1,end)=0;
    indexes = find(chromosome == 1);
    if numel(indexes) > 10,                
        merge_cost = zeros(1,numel(indexes));
        merge_cost(1,1)=returnOneError(serie(1:indexes(2)),degree,typeError);
        for i=2:numel(indexes)-1,
            merge_cost(i)=returnOneError(serie(indexes(i-1):indexes(i+1)),degree,typeError); 
        end
        merge_cost(1,end)=returnOneError(serie(indexes(i-1):end),degree,typeError);

        iterations=0;
        while iterations < max_iter,
            ind = find(merge_cost==min(merge_cost),1);
            if ind == 2,
                merge_cost(ind-1)=returnOneError(serie(1:indexes(ind+1)),degree,typeError);
            elseif ind > 2,
                if ind == numel(merge_cost),
                    merge_cost(ind-1)=returnOneError(serie(indexes(ind-2):end),degree,typeError);
                else
                    merge_cost(ind-1)=returnOneError(serie(indexes(ind-2):indexes(ind+1)),degree,typeError);
                end
            end

            if ind == (numel(merge_cost)-1),
                merge_cost(ind+1)=returnOneError(serie(indexes(ind-1):end),degree,typeError);
            elseif ind < (numel(merge_cost)-1),
                if ind == 1,
                    merge_cost(ind+1)=returnOneError(serie(1:indexes(ind+2)),degree,typeError);
                else
                    merge_cost(ind+1)=returnOneError(serie(indexes(ind-1):indexes(ind+2)),degree,typeError);
                end
            end
            chromosome(1,indexes(ind))=0;
            merge_cost(ind)=[];
            indexes(ind)=[];
            iterations=iterations+1;
        end
        chromosome(1,1)=0;
        chromosome(1,end)=0;
    end
end
