%   COPYRIGHT
%   This file is part of TSSA: https://github.com/ayrna/tssa
%   Original authors: Antonio M. Duran Rosal, Pedro A. Gutierrez
%   Copyright:
%       This software is released under the The GNU General Public License v3.0 licence
%       available at http://www.gnu.org/licenses/gpl-3.0.html
%   Citation: If you use this code, please cite any of the following papers:
%     [1] M. Pérez-Ortiz, A.M. Durán-Rosal, P.A. Gutiérrez, et al.
%         "On the use of evolutionary time series analysis for segmenting paleoclimate data"
%         Neurocomputing, Vol. 326-327, January, 2019, pp. 3-14
%         https://doi.org/10.1016/j.neucom.2016.11.101
%     [2] A.M. Durán-Rosal, P.A. Gutiérrez, F.J. Martínez-Estudillo and C. Hervás-Martínez.
%         "Simultaneous optimisation of clustering quality and approximation error
%         for time series segmentation", Information Sciences, Vol. 442-443, May, 2018, pp. 186-201.
%         https://doi.org/10.1016/j.ins.2018.02.041
%     [3] A.M. Durán-Rosal, J.C. Fernández, P.A. Gutiérrez and C. Hervás-Martínez.
%         "Detection and prediction of segments containing extreme significant wave heights"
%         Ocean Engineering, Vol. 142, September, 2017, pp. 268-279.
%         https://doi.org/10.1016/j.oceaneng.2017.07.009
%     [4] A.M. Durán-Rosal, M. de la Paz Marín, P.A. Gutiérrez and C. Hervás-Martínez.
%         "Identifying market behaviours using European Stock Index time series by 
%         a hybrid segmentation algorithm", Neural Processing Letters,
%         Vol. 46, December, 2017, pp. 767–790.
%         https://doi.org/10.1007/s11063-017-9592-8
%
%% mutation1
% Function: Four mutation types Add/Remove/MoveLeft/MoveRight segment with size restriction
% 
% Input:
%     population: set of segmentations
%     fitness:    fitness value of each segmentation
%     pMut:       mutation probability
%     mutedPoints:    percentage of mutation
%     minSeg:         mimimum segment size
%     maxSeg:         maximum segment size
%     
% Output:
%     mutatedPopulation:  population after applying mutation
%     changedFitness:     fitness of each segmentation (NaN in the case of changes)
function [mutatedPopulation,changedFitness] = mutation1(population,fitness,pMut,mutedPoints,minSeg,maxSeg)
    mutatedPopulation = population;
    changedFitness = fitness;
    sizeChromosome = numel(population(1,:));
    %t= cputime
    for i=1:size(population,1),
        %Mutate?
        if rand()<pMut,
            type = randi(2,1,1);
            cLength=sum(population(i,:));
            if type == 1,
                for j=1:(mutedPoints*cLength),
                    % add point
                    if rand()>0.5,
                        [ind] = find(mutatedPopulation(i,:)==0);
                        attempts=0;
                        flag = 1;
                        while attempts<3 && flag==1,
                            point = ind(randi(numel(ind),1,1));
                            mutatedPopulation(i,point) = 1;
                            previousPoint = find(mutatedPopulation(i,1:point-1)==1,1,'last');
                            nextPoint = find(mutatedPopulation(i,point+1:end)==1,1,'first');
                            nextPoint = nextPoint + point;
                            if isempty(previousPoint),
                                previousPoint=1;
                            end
                            if isempty(nextPoint),
                                nextPoint = sizeChromosome;
                            end
                            if ((nextPoint - point)< minSeg || (point - previousPoint)< minSeg),
                                mutatedPopulation(i,point)=0;
                                attempts=attempts+1;
                            else
                                flag=0;
                            end
                        end
                        if attempts < 3,
                            changedFitness(i) = NaN;
                        end
                    % remove point
                    else
                        [ind] = find(mutatedPopulation(i,:)==1);
                        attempts=0;
                        flag=1;

                        while attempts < 3 && flag==1,
                            choice = randi(numel(ind),1,1);
                            point = ind(choice);
                            mutatedPopulation(i,point) = 0;
                            if choice == 1,
                                previousPoint = 1;
                            else
                                previousPoint = ind(choice-1);
                            end
                            if choice == numel(ind),
                                nextPoint = sizeChromosome;
                            else
                                nextPoint = ind(choice+1);
                            end
                            if nextPoint - previousPoint > maxSeg,
                                mutatedPopulation(i,point)=1;
                                attempts=attempts+1;
                            else
                                flag=0;
                            end

                        end
                        if attempts < 3,
                            changedFitness(i) = NaN;
                        end
                    end
                end
            else
                for j=1:(mutedPoints*cLength),
                    [ind] = find(mutatedPopulation(i,:)==1);
                    [choice] = randi(numel(ind),1);
                    % Desplacement to the right
                    if rand()>0.5,
                        if choice == numel(ind),
                            difference = numel(mutatedPopulation(i,:)) - ind(choice);
                        else
                            difference = ind(choice+1) - ind(choice);
                        end
                        if (difference > minSeg),
                            attempts=0;
                            flag=1;
                            while flag==1 && attempts<3,
                                mutatedPopulation(i,ind(choice)) = 0;
                                displacement = randi(difference-minSeg,1);
                                mutatedPopulation(i,ind(choice)+displacement) = 1;
                                if choice==1,
                                    previousPoint = 1;
                                else
                                    previousPoint = ind(choice-1);
                                end
                                if ind(choice) - previousPoint + displacement > maxSeg,
                                    mutatedPopulation(i,ind(choice)) = 1;
                                    mutatedPopulation(i,ind(choice)+displacement) = 0;
                                    attempts=attempts+1;
                                else
                                    flag=0;
                                end
                            end
                            if attempts < 3,
                            changedFitness(i) = NaN;
                            end
                        end 
                    else %Desplacement to the left
                        if choice == 1,
                            difference = ind(choice);
                        else
                            difference = ind(choice) - ind(choice-1);
                        end
                        if (difference > minSeg),
                            attempts=0;
                            flag=1;
                            while flag==1 && attempts<3,
                                mutatedPopulation(i,ind(choice)) = 0;
                                displacement = randi(difference-minSeg,1);
                                mutatedPopulation(i,ind(choice)-displacement) = 1;
                                if choice == numel(ind),
                                    nextPoint = sizeChromosome;
                                else
                                    nextPoint = ind(choice+1);
                                end
                                if nextPoint - ind(choice) - displacement > maxSeg,
                                    mutatedPopulation(i,ind(choice)) = 1;
                                    mutatedPopulation(i,ind(choice)-displacement) = 0;
                                    attempts=attempts+1;
                                else
                                    flag=0;
                                end
                            end
                            if attempts < 3,
                            changedFitness(i) = NaN;
                            end
                        end
                    end

                end


            end
        end

    end
end
