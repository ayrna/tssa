%   COPYRIGHT
%   This file is part of TSSA: https://github.com/ayrna/tssa
%   Original authors: Antonio M. Duran Rosal, Pedro A. Gutierrez
%   Copyright:
%       This software is released under the The GNU General Public License v3.0 licence
%       available at http://www.gnu.org/licenses/gpl-3.0.html
%   Citation: If you use this code, please cite any of the following papers:
%     [1] A.M. Durán-Rosal, P.A. Gutiérrez, F.J. Martínez-Estudillo and C. Hervás-Martínez.
%         "Simultaneous optimisation of clustering quality and approximation error
%         for time series segmentation", Information Sciences, Vol. 442-443, May, 2018, pp. 186-201.
%         https://doi.org/10.1016/j.ins.2018.02.041
%
%% evaluateFronts
% Function: Metrics for the evaluation of the Pareto front
%     
% Input:
%     fitness:    Consists in two rows which saves the fitness for each objective
%     
% Output:
%     spacing:    Metric Spacing (minimize)
%     m3:         Metric M3 (maximize)
%     hv:         Metric Hyperarea (maximize)
%     HRS:        Metric HRS (minimize)
function [spacing,m3,hv,HRS] = evaluateFronts(fitness)
    spacing = evaluateSpacing(fitness);
    m3 = evaluateM3(fitness);
    hv = evaluateHyperarea(fitness);
    HRS = evaluateHRS(fitness);
end

%% evaluateSpacing
% Function: Calculates the Spacing metric. Lower values indicates that the Pareto front are perfect distributed, on the contrary case, the distribution is irregular
% 
% Input:
%     fitness:    Consists in two rows which saves the fitness for each objective
%     
% Output:
%     spacing:    Value of the Spacing metric
function [spacing] = evaluateSpacing(fitness)
    num = numel(fitness(1,:));
    auxiliar_distances = zeros(1,num);
    d = zeros(1,num);

    for i=1:num,
        auxiliar_distances(1,:)=0;
        for j=1:num,
            auxiliar_distances(1,j)=abs(fitness(1,i)-fitness(1,j))+abs(fitness(2,i)-fitness(2,j));
            if i==j,
               auxiliar_distances(1,j)=inf;
            end
        end
        d(i)=min(auxiliar_distances);
    end

    spacing = sum((d - mean(d)).*(d - mean(d)));
    spacing = sqrt(spacing/(num-1));    
end

%% evaluateM3
% Function: Calculates the M3 metric. To higher value, better hypervolume
% 
% Input:
%     fitness:    Consists in two rows which saves the fitness for each objective
%     
% Output:
%     m3:    Value of the M3 metric
function [m3] = evaluateM3(fitness)
    maximo1 = max(fitness(1,:));
    minimo1 = min(fitness(1,:));

    maximo2 = max(fitness(2,:));
    minimo2 = min(fitness(2,:));

    m3 = sqrt((maximo1-minimo1)^2 + (maximo2-minimo2)^2);
end

%% evaluateHyperarea
% Function: Calculates the Hyperarea metric. To higher value, better hypervolume
% 
% Input:
%     fitness:    Consists in two rows which saves the fitness for each objective
%     
% Output:
%     H:    Value of the Hyperarea metric
function [H] = evaluateHyperarea(fitness)
    [fit1 ind]=sort(fitness(1,:));

    fit2 = fitness(2,ind);

    H = fit1(1)*fit2(1);
    for i=2:numel(fit1),
        H = H + ((fit1(i) - fit1(i-1))*fit2(i)); 
    end

end

%% evaluateHRS
% Function: Calculates the HRS metric. Lower values indicates that the Pareto front are perfect distributed, on the contrary case, the distribution is irregular. It is a generalization of the Spacing.
% 
% Input:
%     fitness:    Consists in two rows which saves the fitness for each objective
%     
% Output:
%     HRS:    Value of the HRS metric
function [HRS] = evaluateHRS(fitness)
    % Minimizar (como el Spacing, es una generalizacion)
    num = numel(fitness(1,:));
    auxiliar_distances = zeros(1,num);
    d = zeros(1,num);

    for i=1:num,
        auxiliar_distances(1,:)=0;
        for j=1:num,
            auxiliar_distances(1,j)=abs(fitness(1,i)-fitness(1,j))+abs(fitness(2,i)-fitness(2,j));
            if i==j,
               auxiliar_distances(1,j)=inf;
            end
        end
        d(i)=min(auxiliar_distances);
    end

    HRS = max(d)/mean(d);   
end
