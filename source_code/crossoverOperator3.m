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
%% crossoverOperator3
% Function: Single point cross over operator which respects the number of segments
% Input:
%     individual1:    first parent
%     individual2:    second parent
% 
% Output:
%     crossedIndividual1: first offspring individual
%     crossedIndividual2: second offspring individual
%     flag:
%         true - the crossover is sucesfully applied
%         false - the crossover is not sucesfully applied
function [crossedIndividual1,crossedIndividual2,flag] = crossoverOperator3(individual1,individual2)
    sizeChromosome = numel(individual1);
    crossPoint = randi(sizeChromosome-3,1,1) + 1;
    ind_1=find(individual1(1,crossPoint:end)==1)+crossPoint-1;
    ind_0=find(individual2(1,crossPoint:end)==0)+crossPoint-1;
    inter10 = intersect(ind_1,ind_0);
    count_10=numel(inter10);

    ind_00=find(individual1(1,crossPoint:end)==0)+crossPoint-1;
    ind_11=find(individual2(1,crossPoint:end)==1)+crossPoint-1;
    inter01 = intersect(ind_00,ind_11);
    count_01=numel(inter01);

    if(count_01~=0 && count_10~=0),
        minimo = min(count_10,count_01);
        crossedIndividual1(1,:) = individual1(1,:);
        crossedIndividual2(1,:) = individual2(1,:);

        inter10=inter10(1:minimo);
        inter01=inter01(1:minimo);

        crossedIndividual1(1,inter10)=0;
        crossedIndividual2(1,inter10)=1;
        crossedIndividual1(1,inter01)=1;
        crossedIndividual2(1,inter01)=0;
        flag = true;
    else
        flag = false;
        crossedIndividual1(1,:)=individual1(1,:);
        crossedIndividual2(1,:)=individual2(1,:);
    end
end
