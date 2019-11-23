%   COPYRIGHT
%   This file is part of TSSA: https://github.com/ayrna/tssa
%   Original authors: Antonio M. Duran Rosal, Pedro A. Gutierrez
%   Copyright:
%       This software is released under the The GNU General Public License v3.0 licence
%       available at http://www.gnu.org/licenses/gpl-3.0.html
%   Citation: If you use this code, please cite the following paper:
%     [1] A.M. Durán-Rosal, P.A. Gutiérrez, F.J. Martínez-Estudillo and C. Hervás-Martínez.
%         "Simultaneous optimisation of clustering quality and approximation error
%         for time series segmentation", Information Sciences, Vol. 442-443, May, 2018, pp. 186-201.
%         https://doi.org/10.1016/j.ins.2018.02.041
%
%% Non domination sort
% Function: Sort the population by the domination of the solution, and then by the crowding distance in each front
% 
% Input:
%     population:     set of solutions (segmentations)
%     fitness1:       fitness in the first objective
%     fitness2:       fitness in the second objective
%     
% Output:
%     f:                  indicates the number of front of each solution
%     sorted_population:  sorted population
%     changedFitness1:    sorted fitness in the first objective
%     changedFitness2:    sorted fitness in the second objective
%     crowding_distances: crowding_distances of each solution
function [f, sorted_population, changedFitness1, changedFitness2, crowding_distances] = non_domination_sort(population, fitness1, fitness2)
    sizeChromosome = numel(population(1,:));
    N = numel(fitness1(1,:));
    crowding_distances = zeros(N,1);
    f = zeros(N,1);
    front = 1;
    F(front).f = [];
    individual = [];
    for i = 1 : N,
        individual(i).n = 0;
        individual(i).p = [];
        for j = 1 : N,
            dom_less = 0;
            dom_equal = 0;
            dom_more = 0;
            %We compare fitness 1 (clustering)
            if fitness1(1,i) < fitness1(1,j),
                dom_less = dom_less + 1;
            elseif fitness1(1,i) == fitness1(1,j),
                dom_equal = dom_equal+1;
            else
                dom_more = dom_more +1;
            end
            % We compare fitness 2 (error approximation)
            if fitness2(1,i) < fitness2(1,j),
                dom_less = dom_less + 1;
            elseif fitness2(1,i) == fitness2(1,j),
                dom_equal = dom_equal+1;
            else
                dom_more = dom_more +1;
            end

            if dom_more == 0 && dom_equal ~= 2,
                individual(i).n = individual(i).n + 1;
            elseif dom_less == 0 && dom_equal ~=2,
                individual(i).p = [individual(i).p j];
            end
        end
        if individual(i).n == 0,
            f(i,1)=front;
            F(front).f = [F(front).f i];
        end
    end

    while ~isempty(F(front).f)
        Q = [];
        for i = 1 : length(F(front).f)
            if ~isempty(individual(F(front).f(i)).p)
                for j = 1 : length(individual(F(front).f(i)).p)
                    individual(individual(F(front).f(i)).p(j)).n = individual(individual(F(front).f(i)).p(j)).n - 1;
                    if individual(individual(F(front).f(i)).p(j)).n == 0
                        f(individual(F(front).f(i)).p(j),1) = front + 1;
                        Q = [Q individual(F(front).f(i)).p(j)];
                    end
                end
            end
        end
        front = front + 1;
        F(front).f = Q;
    end
    sorted_population = zeros(numel(population(:,1)),numel(population(1,:)));
    changedFitness1 = zeros(1,numel(fitness1(1,:)));
    changedFitness2 = zeros(1,numel(fitness2(1,:)));
    [temp, index_of_fronts] = sort(f(:,1));
    for i = 1 : length(index_of_fronts),
        sorted_population(i,:)=population(index_of_fronts(i),:);
        changedFitness1(1,i)=fitness1(1,index_of_fronts(i));
        changedFitness2(1,i)=fitness2(1,index_of_fronts(i));
        f(i,1)=temp(i,1);
    end

    % Apply crowding
    current_index = 0;

    for front = 1 : (length(F) - 1),
        distance = 0;
        y = [];
        yFitness1 = [];
        yFitness2 = [];
        yf = [];
        previous_index = current_index + 1;
        for i = 1: length(F(front).f)
            y(i,:) = sorted_population(current_index + i, :);
            yFitness1(1,i) = changedFitness1(1,current_index + i);
            yFitness2(1,i) = changedFitness2(1,current_index + i);
            yf(i,1) = f(current_index + i,1);                    
        end
        current_index = current_index + i;

        %%FITNESS 1
        [sorted_based_on_objective, index_of_objectives] = sort(yFitness1(1,:));
        sorted_based_on_objective = [];
        sorted_based_on_objectiveFitness1 = [];
        sorted_based_on_objectiveFitness2 = [];
        for j = 1 : length(index_of_objectives)
            sorted_based_on_objective(j,:) = y(index_of_objectives(j),:);
            sorted_based_on_objectiveFitness1(1,j)=yFitness1(1,index_of_objectives(j));
            sorted_based_on_objectiveFitness2(1,j)=yFitness2(1,index_of_objectives(j));
        end
        f_max = sorted_based_on_objectiveFitness1(1,length(index_of_objectives));
        f_min = sorted_based_on_objectiveFitness1(1,1);
        y(index_of_objectives(length(index_of_objectives)),sizeChromosome+1) = Inf;
        y(index_of_objectives(1), sizeChromosome + 1) = Inf;
        for j = 2 : length(index_of_objectives) - 1
           next_obj  = sorted_based_on_objectiveFitness1(1,j + 1);
           previous_obj  = sorted_based_on_objectiveFitness1(1,j - 1);
           if (f_max - f_min == 0)
               y(index_of_objectives(j),sizeChromosome+1) = Inf;
           else
               y(index_of_objectives(j),sizeChromosome+1) = (next_obj - previous_obj)/(f_max - f_min);
           end
        end
        %%FITNESS 2
        [sorted_based_on_objective, index_of_objectives] = sort(yFitness2(1,:));
        sorted_based_on_objective = [];
        sorted_based_on_objectiveFitness1 = [];
        sorted_based_on_objectiveFitness2 = [];
        for j = 1 : length(index_of_objectives)
            sorted_based_on_objective(j,:) = y(index_of_objectives(j),:);
            sorted_based_on_objectiveFitness1(1,j)=yFitness1(1,index_of_objectives(j));
            sorted_based_on_objectiveFitness2(1,j)=yFitness2(1,index_of_objectives(j));
        end
        f_max = sorted_based_on_objectiveFitness2(1,length(index_of_objectives));
        f_min = sorted_based_on_objectiveFitness2(1,1);
        y(index_of_objectives(length(index_of_objectives)),sizeChromosome+2) = Inf;
        y(index_of_objectives(1), sizeChromosome + 2) = Inf;
        for j = 2 : length(index_of_objectives) - 1
           next_obj  = sorted_based_on_objectiveFitness2(1,j + 1);
           previous_obj  = sorted_based_on_objectiveFitness2(1,j - 1);
           if (f_max - f_min == 0)
               y(index_of_objectives(j),sizeChromosome+2) = Inf;
           else
               y(index_of_objectives(j),sizeChromosome+2) = (next_obj - previous_obj)/(f_max - f_min);
           end
        end
        distance = [];
        distance(:,1)= zeros(length(F(front).f),1);
        distance(:,1) = distance(:,1) + y(:,sizeChromosome+1) + y(:,sizeChromosome+1);

        y(:,sizeChromosome+1)=[];
        y(:,sizeChromosome+1)=[];

        sorted_populationAux=[];
        changedFitness1Aux=[];
        changedFitness2Aux=[];
        fAux=[];

        [sorted_based_on_distances, index_of_distances] = sort(distance(:,1),'descend');
        for j = 1: length(index_of_distances),
            sorted_populationAux(j,:) = y(index_of_distances(j),:);
            changedFitness1Aux(1,j) = yFitness1(1,index_of_distances(j));
            changedFitness2Aux(1,j) = yFitness2(1,index_of_distances(j));
            fAux(j,1) = yf(index_of_distances(j),1);
        end
        sorted_population(previous_index:current_index,:)=sorted_populationAux(:,:);
        changedFitness1(1,previous_index:current_index)=changedFitness1Aux(1,:);
        changedFitness2(1,previous_index:current_index)=changedFitness2Aux(1,:);
        f(previous_index:current_index,:)=fAux(:,1);
        crowding_distances(previous_index:current_index,:)=sorted_based_on_distances(:,1);
    end
    sorted_population; changedFitness1; changedFitness2; f; crowding_distances;
    % When the function is finished we have everything in:
    % sorted_population, changedFitness1, changedFitness2, f, crowding_distances
end
