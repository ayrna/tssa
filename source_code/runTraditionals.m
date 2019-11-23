%   COPYRIGHT
%   This file is part of TSSA: https://github.com/ayrna/tssa
%   Original authors: Antonio M. Duran Rosal, Pedro A. Gutierrez
%   Copyright:
%       This software is released under the The GNU General Public License v3.0 licence
%       available at http://www.gnu.org/licenses/gpl-3.0.html
%   Citation: If you use this code, please cite any of the following papers:
%     [1] E. Keogh, S. Chu, D. Hart and M. Pazzani.
%         "Segmenting time series: A survey and novel approach",
%		  In Data mining in time series databases, 2004, pp.1-21.
%         https://doi.org/10.1142/9789812565402_0001
%     [2] A.M. Durán-Rosal, P.A. Gutiérrez, Á. Carmona-Poyato and C. Hervás-Martínez.
%         "A hybrid dynamic exploitation barebones particle swarm optimisation
%         algorithm for time series segmentation", Neurocomputing,
%         Vol. 353, August, 2019, pp. 45-55.
%		  https://doi.org/10.1016/j.neucom.2018.05.129
%     [3] A.M. Durán-Rosal, P.A. Gutiérrez, S. Salcedo-Sanz and C. Hervás-Martínez.
%         "A statistically-driven Coral Reef Optimization algorithm for optimal
%         size reduction of time series", Applied Soft Computing, 
%         Vol. 63. 2018, pp. 139-153.
%         https://doi.org/10.1016/j.asoc.2017.11.037
%
%% runTraditionals
% Main Traditional function: It is a selector of the type of traditional algorithm
% 
% Input:
%     typeAlgorithm:  specifiy the type of algorithm
%         1 - Sliding_Window
%         2 - Top_Down
%         3 - Bottom_Up
%         4 - SWAB
%         5 - Top_Down (given numSeg)
%         6 - Bottom_Up (given numSeg)
%     serie:          time series
%     maxError:       maximum error allowed in the segment
%     numSeg:         number of segment for the segmentation
%     degree:         degree of approximation
%     typeError:      type of error (see returnOneError)
%     buffer_size:    size of the initial window (only for SWAB)
%     
% Output:
%     chromosome:     segmentation by the traditional algorithm
function [chromosome] = runTraditionals(typeAlgorithm,serie,maxError,numSeg,degree,typeError,buffer_size)
    switch typeAlgorithm
        case 1
            chromosome = Sliding_Window(serie,maxError,degree,typeError);
        case 2
            chromosome = Top_Down(serie,maxError,degree,typeError);
        case 3
            chromosome = Bottom_Up(serie,maxError,degree,typeError);
        case 4
            chromosome = SWAB(serie,maxError,degree,typeError,buffer_size);
        case 5
            chromosome = Top_Down_numSeg(serie,numSeg,degree,typeError);
        case 6
            chromosome = Bottom_Up_numSeg(serie,numSeg,degree,typeError);
        otherwise
            chromosome = Sliding_Window(serie,maxError,degree,typeError);
    end
end

%% SWAB
% Function: Performance of the SWAB algorithm
% 
% Input:
%     segment:        values of the time series
%     maxError:       maximum error allowed in the segment
%     degree:         degree of approximation
%     typeError:      type of error (see returnOneError)
%     buffer_size:    size of the initial window
%     
% Output:
%     chromosome:     segmentation by the SWAB algorithm
function [chromosome] = SWAB(segment,max_error,degree,typeError,buffer_size)
    chromosome = zeros(1,numel(segment));
    left = 1;
    right = buffer_size;
    salida = 0;
    while salida==0 && left~=right,
        initial_left = left;
        initial_right = right;
        %chromosome(left:right) = 0;
        chromosome(left:right) = Bottom_Up(segment(left:right),max_error,degree,typeError);
        chromosome(left)=1;
        chromosome(right)=1;
        ind = find(chromosome(left+1:right)==1,1);
        chromosome(left+ind)=1;
        left=left+ind;
        %chromosome(right:end) = 0;
        chromosome(right:end) = Sliding_Window(segment(right:end),max_error,degree,typeError);
        chromosome(end) = 1;
        ind = find(chromosome(right:end)==1,1);
        right = right + ind - 1; 
        if(initial_left == left) && (initial_right == right),
            salida=1;
        end
    end
    chromosome(1,1)=0;
    chromosome(1,end)=0;
end

%% Sliding_Window
% Function: Performance of the Sliding Window algorithm
% 
% Input:
%     segment:        values of the time series
%     maxError:       maximum error allowed in the segment
%     degree:         degree of approximation
%     typeError:      type of error (see returnOneError)
%     
% Output:
%     chromosome:     segmentation by the Sliding Window algorithm
function [chromosome] = Sliding_Window(segment,max_error,degree,typeError)
    chromosome = zeros(1,numel(segment));
    left=1;
    right=2;
    while right < numel(segment),
        while (right <= numel(segment)) && (returnOneError(segment(left:right),degree,typeError) < max_error),
            right=right+1;
        end
        left=right-1;
        chromosome(left)=1;
    end
    chromosome(numel(segment))=0;
end

%% Top_Down
% Function: Performance of the Top Down algorithm
% 
% Input:
%     segment:        values of the time series
%     maxError:       maximum error allowed in the segment
%     degree:         degree of approximation
%     typeError:      type of error (see returnOneError)
%     
% Output:
%     chromosome:     segmentation by the Top Down algorithm
function [chromosome] = Top_Down(segment,max_error,degree,typeError)
    chromosome = zeros(1,numel(segment));
    best_so_far = inf;
    for i=2:numel(segment)-1,
        error_partition = returnOneError(segment(1:i),degree,typeError) + returnOneError(segment(i:end),degree,typeError);
        if error_partition < best_so_far,
            breakpoint = i;
            best_so_far = error_partition;
        end
    end

    if returnOneError(segment(1:breakpoint),degree,typeError) > max_error,
        chromosome(1:breakpoint) = Top_Down(segment(1:breakpoint),max_error);
    end

    if returnOneError(segment(breakpoint:end),degree,typeError) > max_error,
        chromosome(breakpoint:end) = Top_Down(segment(breakpoint:end),max_error);
    end
    chromosome(1,breakpoint)=1;            
end

%% Bottom_Up
% Function: Performance of the Bottom_Up algorithm
% 
% Input:
%     segment:        values of the time series
%     maxError:       maximum error allowed in the segment
%     degree:         degree of approximation
%     typeError:      type of error (see returnOneError)
%     
% Output:
%     chromosome:     segmentation by the Bottom_Up algorithm
function [chromosome] = Bottom_Up(segment,max_error,degree,typeError)
    chromosome = ones(1,numel(segment));
    merge_cost = zeros(1,numel(segment));
    merge_cost(1,1)=inf;
    merge_cost(1,end)=inf;
    indexes = find(chromosome==1);
    for i=2:numel(segment)-1,
        merge_cost(i)=returnOneError(segment(indexes(i-1):indexes(i+1)),degree,typeError); 
    end

    while min(merge_cost) < max_error,
        ind = find(merge_cost==min(merge_cost),1);
        if ind > 2,
            merge_cost(ind-1)=returnOneError(segment(indexes(ind-2):indexes(ind+1)),degree,typeError);
        end
        if ind < (numel(merge_cost)-1),
            merge_cost(ind+1)=returnOneError(segment(indexes(ind-1):indexes(ind+2)),degree,typeError);
        end
        chromosome(1,indexes(ind))=0;
        merge_cost(ind)=[];
        indexes(ind)=[];                
    end
    chromosome(1,1)=0;
    chromosome(1,end)=0;           
end


%% Top_Down_numSeg
% Function: Performance of the Top Down algorithm given a number of segments
% 
% Input:
%     segment:        values of the time series
%     numSeg:         number of segments
%     degree:         degree of approximation
%     typeError:      type of error (see returnOneError)
%     
% Output:
%     chromosome:     segmentation by the Top Down algorithm
function [chromosome] = Top_Down_numSeg(segment,numSeg,degree,typeError)
    chromosome=zeros(1,numel(segment));
    gain_cost = zeros(1,numel(chromosome));
    max_iter=numSeg-1;

    base = returnOneError(segment(1:end),degree,typeError);
    for j=2:numel(chromosome),
        gain_cost(1,j) = base - (returnOneError(segment(1:j),degree,typeError) + returnOneError(segment(j:end),degree,typeError));
    end

    gain_cost(1,1)=-Inf;
    gain_cost(1,end)=-Inf;

    iterations = 0;

    while iterations < max_iter,
        ind = find(gain_cost==max(gain_cost),1);
        if ind ==1,
            fprintf('Error');
        end
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
        base = returnOneError(segment(left:ind),degree,typeError);
        for j=left:ind,
            gain_cost(1,j)=base - (returnOneError(segment(left:j),degree,typeError)+returnOneError(segment(j:ind),degree,typeError));
        end
        base = returnOneError(segment(ind:right),degree,typeError);
        for j=ind:right,
            gain_cost(1,j)=base - (returnOneError(segment(ind:j),degree,typeError)+returnOneError(segment(j:right),degree,typeError));
        end
        gain_cost(1,ind)=-Inf;
        gain_cost(1,left)=-Inf;
        gain_cost(1,right)=-Inf;
        iterations=iterations+1;
    end          
end

%% Bottom_Up
% Function: Performance of the Bottom_Up algorithm given a number of segments
% 
% Input:
%     segment:        values of the time series
%     numSeg:         number of segments
%     degree:         degree of approximation
%     typeError:      type of error (see returnOneError)
%     
% Output:
%     chromosome:     segmentation by the Bottom_Up algorithm
function [chromosome] = Bottom_Up_numSeg(segment,numSeg,degree,typeError)
    chromosome=ones(1,numel(segment));
    chromosome(1,1)=0;
    chromosome(1,end)=0;
    indexes = find(chromosome == 1);
    max_iter=numel(indexes)-numSeg+1;
    if numel(indexes) > 10,                
        merge_cost = zeros(1,numel(indexes));
        merge_cost(1,1)=returnOneError(segment(1:indexes(2)),degree,typeError);
        for i=2:numel(indexes)-1,
            merge_cost(i)=returnOneError(segment(indexes(i-1):indexes(i+1)),degree,typeError); 
        end
        merge_cost(1,end)=returnOneError(segment(indexes(i-1):end),degree,typeError);

        iterations=0;
        while iterations < max_iter,
            ind = find(merge_cost==min(merge_cost),1);
            if ind == 2,
                merge_cost(ind-1)=returnOneError(segment(1:indexes(ind+1)),degree,typeError);
            elseif ind > 2,
                if ind == numel(merge_cost),
                    merge_cost(ind-1)=returnOneError(segment(indexes(ind-2):end),degree,typeError);
                else
                    merge_cost(ind-1)=returnOneError(segment(indexes(ind-2):indexes(ind+1)),degree,typeError);
                end
            end

            if ind == (numel(merge_cost)-1),
                merge_cost(ind+1)=returnOneError(segment(indexes(ind-1):end),degree,typeError);
            elseif ind < (numel(merge_cost)-1),
                if ind == 1,
                    merge_cost(ind+1)=returnOneError(segment(1:indexes(ind+2)),degree,typeError);
                else
                    merge_cost(ind+1)=returnOneError(segment(indexes(ind-1):indexes(ind+2)),degree,typeError);
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
