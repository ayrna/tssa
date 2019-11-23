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
%% fitnessF
% Function: Calculate the fitness function for a clustering assignment
% 
% Input:
%     typeClustMeasure: type of measure (quality goodness)
%         1 - Davies-Bouldin
%         2 - Dunn
%         3 - Calinski and Harabasz
%         4 - SSE
%         5 - NSSE
%         6 - Silhouette
%         7 - GD33
%         8 - GD43
%         9 - GD53
%         10 - COP
%     charac:         set of patterns (mapped segments, characteristics)
%     assignation:    array of clustering assignation
%     centroids:      centroids of the resulting clustering
%     k:              number of clusters
%     
% Output:
%     fitness:    value of the chosen fitness
function [fitness] = fitnessF(typeClustMeasure,charac, assignation, centroids,k)

    realTargets = unique(assignation);
    newk = numel(realTargets);
    distanceToCentroid = zeros(1,newk);
    fitness = 0;
    for i=1:newk
        distanceToCentroid(i) = (mean(pdist2(charac(assignation==realTargets(i),:),centroids(realTargets(i),:))));
    end

    if newk < k,
        fitness = 0;

    elseif typeClustMeasure == 1,
        %%  typeFitness == 1 => Davies-Bouldin Criterion => Minimize
        for i=1:newk
            selected = 0;
            for j=(i+1):newk
                aux = (distanceToCentroid(i) + distanceToCentroid(j))/pdist2(centroids(realTargets(i),:),centroids(realTargets(j),:));
                if aux>selected,
                    selected = aux;
                end
            end
            fitness = fitness + selected;
        end
        fitness = fitness / newk;
        % Change the fitness so as we have to maximise it
        fitness = 1/(1+fitness);

    elseif typeClustMeasure == 2,
        %%  typeFitness == 2 => Dunn index => Maximize
        bad_solution=0;
        for i=1:newk,
            num_elem = numel(find(assignation==i));
            if num_elem < 2,
                bad_solution=1;
            end
        end
        if bad_solution==0,
            selected = inf;
            diam = zeros(1,newk);
            indexesI = false(newk,size(assignation,1));
            allDistances = squareform(pdist(charac));
            for i=1:newk
                indexesI(i,:) = assignation==realTargets(i);
                numElI = sum(indexesI(i,:));
                diam(i) = (1/(2*numElI*(numElI-1)))*sum(sum(allDistances(indexesI(i,:),indexesI(i,:))));
            end
            maxintraclusterdistance = max(diam);
            for i=1:newk
                for j=(i+1):newk
                    aux = min(min(allDistances(indexesI(i,:),indexesI(j,:))))/maxintraclusterdistance;
                    if aux < selected
                        selected = aux;
                    end

                end
            end
            fitness = selected;
        else
            fitness = 0;
        end

    elseif typeClustMeasure == 3,
        %% typeFitness == 3 => CALINSKI AND HARABASZ INDEX 1974 => Maximize
        Sw=zeros(numel(charac(1,:)),numel(charac(1,:))); %Degree
        Sb=zeros(numel(charac(1,:)),numel(charac(1,:)));
        m = mean(charac);
        totalNum = numel(assignation);
        for currentClass = 1:newk,
            indexesI = assignation == realTargets(currentClass);
            mk = mean(charac(indexesI,:));
            Sw = Sw + (1/totalNum)*cov( (charac(indexesI,:)),1);
            Sb = Sb + sum(indexesI)/totalNum*(mk-m)'*(mk-m);
        end


        fitness = (trace(Sb)/(k-1))/(trace(Sw)/(totalNum-k));

    elseif typeClustMeasure == 4,
        %% typeFitness == 4 => SSE => Minimize
        for i=1:newk
            fitness = fitness + sum(pdist2(charac(assignation==realTargets(i),:),centroids(realTargets(i),:)));
        end
        fitness = fitness/numel(assignation);
        fitness = 1/(1+fitness);

    elseif typeClustMeasure == 5,
        %% typeFitness == 5 => SSE/DentreK => Minimize
        for i=1:newk
            fitness = fitness + sum(pdist2(charac(assignation==realTargets(i),:),centroids(realTargets(i),:)));
        end
        dist_entre=0;
        for i=1:newk,
            for j=(i+1):newk,
                dist_entre = dist_entre + sum(pdist2(centroids(i,:),centroids(j,:)));  
            end
        end
        fitness = fitness/numel(assignation);
        dist_entre = dist_entre / factorial(newk-1);

        fitness = fitness / dist_entre;

        fitness = 1/(1+fitness);

    elseif typeClustMeasure == 6,
        %% typeFitness == 6 => Silhouette index => Maximize
        fitness = 0;
        allDistances = squareform(pdist(charac));
        for k=1:newk,
           indK = find(assignation==realTargets(k));
           for i=1:numel(indK),
               b=zeros(1,newk);
               for l=1:newk,
                   b(1,l)=sum(allDistances(indK(i),assignation==realTargets(l)));
                   b(1,l)=b(1,l)/(sum(assignation==realTargets(l)));
               end
               a=b(1,k);
               b(k)=[];
               b=min(b);
               fitness=fitness + ( (b-a)/(max(a,b)));
           end                                 
        end
        fitness = fitness / numel(assignation);

    elseif typeClustMeasure == 7,
        %% typeFitness == 7 => gD33 index => Maximize
        fitness=generalizedDunn(charac, assignation, centroids, 3, 3);

    elseif typeClustMeasure == 8,
        %% typeFitness == 8 => gD43 index => Maximize
        fitness=generalizedDunn(charac, assignation, centroids, 4, 3);

    elseif typeClustMeasure == 9,
        %% typeFitness == 9 => gD53 index => Maximize
        fitness=generalizedDunn(charac, assignation, centroids, 5, 3);

    else
        %% typeFitness == 10 => COP index => Minimize
        allDistances = squareform(pdist(charac));
        for k=1:newk,
           num = sum(pdist2(charac(assignation==k,:),centroids(k,:)));
           num = num / sum(assignation==k);

           maximos=zeros(1,newk);

           for i=1:newk,
               maximos(1,i)=max(max(allDistances(assignation==k,assignation==i)));
           end
           maximos(k)=[];
           den=min(maximos);
           fitness=fitness + num/den;
        end
        fitness=fitness/numel(assignation);
        fitness = 1/(1+fitness);

    end

end

%% generalizedDunn
% Function: Calculate the generalized Dunn
% 
% Input:
%     charac:         set of patterns (mapped segments, characteristics)
%     assignation:    array of clustering assignation
%     centroids:      centroids of the resulting clustering
%     type_delta:     lower delta
%     type_Delta:     capital delta
%     
% Output:
%     gDunn:      value of the generalized Dunn
function [gDunn] = generalizedDunn(charac, assignation, centroids, type_delta, type_Delta)
   newk=numel(unique(assignation));
   min = deltaValue(charac,assignation,centroids,type_delta,1,2);
   max = deltaValueG(charac,assignation,centroids,type_Delta,1);
   for i=1:newk,
       for j=1:newk,
          if i~=j,
              newMin=deltaValue(charac,assignation,centroids,type_delta,i,j);
              if newMin < min,
                  min=newMin;
              end
          end
       end
       newMax=deltaValueG(charac,assignation,centroids,type_Delta,i);
       if newMax > max,
           max=newMax;
       end
   end
   gDunn=min/max;
end

%% deltaValue
% Function: Calculate the value for \delta
% 
% Input:
%     charac:         set of patterns (mapped segments, characteristics)
%     assignation:    array of clustering assignation
%     centroids:      centroids of the resulting clustering
%     type:           type of calculation
%     k:              index of first centroid
%     l:              index of second centroid
%     
% Output:
%     value:          value for \delta
function [value] = deltaValue(charac,assignation,centroids,type,k,l)
    value=0;
    if type == 3,
        value=pdist2(charac(assignation==k,:),charac(assignation==l,:));
        value=sum(sum(value))/(sum(assignation==k)*sum(assignation==l));
    elseif type == 4,
        value=pdist2(centroids(k,:),centroids(l,:));    
    else
        value=sum(pdist2(charac(assignation==k,:),centroids(k,:)));
        value2=sum(pdist2(charac(assignation==l,:),centroids(l,:)));
        value=value+value2;
        value=value/(sum(assignation==k)+sum(assignation==l));
    end


end

%% deltaValueG
% Function: Calculate the value for \Delta
% 
% Input:
%     charac:         set of patterns (mapped segments, characteristics)
%     assignation:    array of clustering assignation
%     centroids:      centroids of the resulting clustering
%     type:           type of calculation
%     k:              index of first centroid
%     
% Output:
%     value:          value for \Delta
function [value] = deltaValueG(charac, assignation, centroids, type, k)
    value=0;
    % Types 1 and 3, is correct for generalizedDunn
    if type==1,
        allDistances=squareform(pdist(charac(assignation==k)));
        value=max(max(allDistances));
    else
        value=sum(pdist2(charac(assignation==k,:),centroids(k,:)));
        value=2*value/(sum(assignation==k));
    end

end
