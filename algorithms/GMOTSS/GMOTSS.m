%GMOTSS Genetic Multiobjective Time series segmentation [1]
%
%   GMOTSS methods:
%      runAlgorithm               - runs the corresponding algorithm (GMOTSS in [1])
%      saveInformation            - specific information of the algorithm
%      saveAll                    - save all information of the algorithm
%
%   References:
%     [1] A.M. Durán-Rosal, P.A. Gutiérrez, F.J. Martínez-Estudillo and C. Hervás-Martínez.
%         "Simultaneous optimisation of clustering quality and approximation error
%         for time series segmentation", Information Sciences, Vol. 442-443, May, 2018, pp. 186-201.
%         https://doi.org/10.1016/j.ins.2018.02.041
%
%   This file is part of TSSA: https://github.com/ayrna/tssa
%   Original authors: Antonio M. Duran Rosal, Pedro A. Gutierrez Peña
%   Citation: If you use this code, please cite the associated paper [1]
%   Copyright:
%       This software is released under the The GNU General Public License v3.0 licence
%       available at http://www.gnu.org/licenses/gpl-3.0.html 
classdef GMOTSS < handle
    
    properties
        name_parameters = {'numIt','nPobl','k','pCross','pMut','seed','minSeg','maxSeg','sizeChromosome','iterClust','polyDegree','characActivation'}
        dataFile
        data
        idealSegFile
        idealseg
        parameters
    end
    
    methods
        
        %% Constructor
        function obj = GMOTSS()
            obj.defaultParameters();
        end

        %% Default parameters of the class
        function obj = defaultParameters(obj)
            % Number of generations
            obj.parameters.numIt = 200;
            % Population size
            obj.parameters.nPobl = 80;
            % crossover probability
            obj.parameters.pCross = 0.8;
            % mutation probability
            obj.parameters.pMut = 0.2;
            % percentage cut points to be mutated
            obj.parameters.mutedPoints = 0.1;
            % Random number generation seed
            obj.parameters.seed = 1;
            % Minimum length of each segment
            obj.parameters.minSeg = 2;
            % Maximum length of each segment
            obj.parameters.maxSeg = 100;
            % number of cluster for k-means
            obj.parameters.k = 5;
            % max number of iterations for k-means
            obj.parameters.iterClust = 20;
            % fitness function
            obj.parameters.typeFitness = 1;
            % degree for approximations (0 - Interpolation, >=1 degree)
            obj.parameters.polyDegree = 1;
            % bool array for characteristics:
            obj.parameters.characActivation = [1 1 1 1 1 0];
            % Fitness for error: 1-RMSE 2-RMSEp 3-MAXe
            obj.parameters.typeError = 1;
        end
        
        %% Parameters of the algorithm
        function [parameters_as_str] = getParameters(obj)
            parameters = obj.parameters;
            
            fields = fieldnames(parameters);
            parameters_as_str = '';
            
            for i = 1:numel(fields)
                if ~strcmp(fields{i},'characActivation'),
                    parameters_as_str = [parameters_as_str sprintf('%s;%f\n', fields{i}, parameters.(fields{i}))];
                else
                    parameters_as_str = [parameters_as_str sprintf('%s;%s\n', fields{i}, num2str(parameters.(fields{i})))];
                end
            end
        end
        
        %% Main algorithm
        function [information,informationBestClustering,informationBestError] = runAlgorithm(obj, serie)
            addpath(['..' filesep '..' filesep 'source_code' filesep]);
            addpath(['..' filesep '..' filesep 'source_code' filesep 'kmeans' filesep]);
            
            obj.data = serie;
            nOfData = length(serie);
            obj.parameters.sizeChromosome = nOfData;
            
            if strcmp(version('-release'),'2013a')
                s = RandStream('mt19937ar','Seed',obj.parameters.seed);
                RandStream.setGlobalStream(s);
            else
                s = RandStream.create('mt19937ar','seed',obj.parameters.seed);
                RandStream.setDefaultStream(s);
            end
            
            %'Initialisation'
            initialPopulation = initialisePopulation1(obj.parameters.nPobl,obj.parameters.sizeChromosome,obj.parameters.minSeg,obj.parameters.maxSeg);
            
            %'Evaluation'
            oldFitness = zeros(1,obj.parameters.nPobl)*NaN;
            [initialFitness1,initialFitness2] = evaluateBothFitness(obj.parameters.typeError,obj.parameters.typeFitness,initialPopulation,...
                                                                                oldFitness,oldFitness,obj.parameters.k,obj.parameters.iterClust,obj.data,...
                                                                                obj.parameters.characActivation,obj.parameters.polyDegree);
            [f,currentPopulation,currentFitness1,currentFitness2,crowding_distances] = non_domination_sort(initialPopulation, initialFitness1, initialFitness2);
            
            
            for i=1:obj.parameters.numIt,
                
                parentPopulation = currentPopulation;
                parentFitness1 = currentFitness1;
                parentFitness2 = currentFitness2;
                
                %'Crossover'
                [newPopulation, newFitness] = crossoverStr1Op1(parentPopulation,parentFitness1,obj.parameters.pCross,...
                                                               obj.parameters.minSeg,obj.parameters.maxSeg,3);
                
                %'Mutation'
                [newPopulation, newFitness] = mutation1(newPopulation,newFitness,obj.parameters.pMut,obj.parameters.mutedPoints,...
                                                        obj.parameters.minSeg,obj.parameters.maxSeg);

                %'Evaluation'
                copyParentFitness2 = parentFitness2;
                copyParentFitness2(isnan(newFitness)==1)=NaN;
                [newFitness1,newFitness2] = evaluateBothFitness(obj.parameters.typeError,obj.parameters.typeFitness,newPopulation,newFitness,...
                                                              copyParentFitness2,obj.parameters.k,obj.parameters.iterClust,obj.data,...
                                                              obj.parameters.characActivation,obj.parameters.polyDegree);
                [f,resultantPopulation,resultantFitness1,resultantFitness2,crowding_distances] = non_domination_sort([parentPopulation; newPopulation], [parentFitness1 newFitness1], [parentFitness2 newFitness2]);
                
                 %'Selection'
                [currentPopulation, currentFitness1, currentFitness2] = selection2GMO(resultantPopulation,resultantFitness1,resultantFitness2,obj.parameters.nPobl);
            end
            
            % Best global solution: Fitness closer to 1,1
            minimo = min(currentFitness1(1,:));
            maximo = max(currentFitness1(1,:));
            normFit1(1,:) = (currentFitness1(1,:)-minimo) / (maximo -minimo);
            minimo = min(currentFitness2(1,:));
            maximo = max(currentFitness2(1,:));
            normFit2(1,:) = (currentFitness2(1,:)-minimo) / (maximo -minimo);
            normFit1(1,:)= 1 - normFit1(1,:);
            normFit1(1,:)= normFit1(1,:).*normFit1(1,:);
            normFit2(1,:)= 1 - normFit2(1,:);
            normFit2(1,:)= normFit2(1,:).*normFit2(1,:);
            distances_to_ideal(1,:) = sqrt(normFit1(1,:)+normFit2(1,:));
            [mini, fbestidx] = min(distances_to_ideal(1,:));
            
            fprintf('*******BEST SEGMENTATION IN GLOBAL TERMS*******\n');
            fprintf('Indx: %d\n', fbestidx);
            fprintf('Fitness Clustering: %.15f\n', currentFitness1(1,fbestidx));
            fprintf('Fitness Error: %.15f\n', currentFitness2(1,fbestidx));
            % Information for the reporter
            information.fbestClustering = currentFitness1(1,fbestidx);
            information.fbestError = currentFitness2(1,fbestidx);
            information.segmentation = currentPopulation(fbestidx,:);
            information.errors = computeErrors(information.segmentation,obj.data,obj.parameters.polyDegree);
            information.features = computeMetrics(information.segmentation,obj.data,obj.parameters.characActivation,obj.parameters.polyDegree,obj.parameters.typeError);
            information.estimatedSerie = estimationSerie(information.segmentation,obj.data,obj.parameters.polyDegree);
            information.cuts = find(information.segmentation==1);
            [normCharac] = normalizeFunction(information.features);
            [information.L, information.C] = clusteringKmeans(normCharac,obj.parameters.k,obj.parameters.iterClust);
            information.parameters = obj.parameters;
            information.fitness = [currentFitness1; currentFitness2];
            information.number_of_firstFront = numel(find(f==1));
            if information.number_of_firstFront > obj.parameters.nPobl,
                information.number_of_firstFront = obj.parameters.nPobl;
            end
            [spacing,m3,hv,HRS]=evaluateFronts(information.fitness(:,1:information.number_of_firstFront));
            information.spacing = spacing;
            information.m3 = m3;
            information.hv = hv;
            information.HRS = HRS;
            information.fbest = fbestidx;
            information.degree = obj.parameters.polyDegree;
            
            % Best solution in Clustering terms
            [maxClustering, fbestidx] = max(currentFitness1(1,:));
            fprintf('*******BEST SEGMENTATION IN CLUSTERING TERMS*******\n');
            fprintf('Indx: %d\n', fbestidx);
            fprintf('Fitness Clustering: %.15f\n', currentFitness1(1,fbestidx));
            fprintf('Fitness Error: %.15f\n', currentFitness2(1,fbestidx));
            % Information for the reporter
            informationBestClustering.fbestClustering = currentFitness1(1,fbestidx);
            informationBestClustering.fbestError = currentFitness2(1,fbestidx);
            informationBestClustering.segmentation = currentPopulation(fbestidx,:);
            informationBestClustering.errors = computeErrors(informationBestClustering.segmentation,obj.data,obj.parameters.polyDegree);
            informationBestClustering.features = computeMetrics(informationBestClustering.segmentation,obj.data,obj.parameters.characActivation,obj.parameters.polyDegree,obj.parameters.typeError);
            informationBestClustering.estimatedSerie = estimationSerie(informationBestClustering.segmentation,obj.data,obj.parameters.polyDegree);
            informationBestClustering.cuts = find(informationBestClustering.segmentation==1);
            [normCharac] = normalizeFunction(informationBestClustering.features);
            [informationBestClustering.L, informationBestClustering.C] = clusteringKmeans(normCharac,obj.parameters.k,obj.parameters.iterClust);
            informationBestClustering.parameters = obj.parameters;
            informationBestClustering.fitness = [currentFitness1; currentFitness2];
            informationBestClustering.number_of_firstFront = numel(find(f==1));
            if informationBestClustering.number_of_firstFront > obj.parameters.nPobl,
                informationBestClustering.number_of_firstFront = obj.parameters.nPobl;
            end
            informationBestClustering.spacing = spacing;
            informationBestClustering.m3 = m3;
            informationBestClustering.hv = hv;
            informationBestClustering.HRS = HRS;
            informationBestClustering.fbest = fbestidx;
            informationBestClustering.degree = obj.parameters.polyDegree;
            
           % Best solution in Error terms
            [maxError, fbestidx] = max(currentFitness2(1,:));
            fprintf('*******BEST SEGMENTATION IN ERROR TERMS*******\n');
            fprintf('Indx: %d\n', fbestidx);
            fprintf('Fitness Clustering: %.15f\n', currentFitness1(1,fbestidx));
            fprintf('Fitness Error: %.15f\n', currentFitness2(1,fbestidx));
            % Information for the reporter
            informationBestError.fbestClustering = currentFitness1(1,fbestidx);
            informationBestError.fbestError = currentFitness2(1,fbestidx);
            informationBestError.segmentation = currentPopulation(fbestidx,:);
            informationBestError.errors = computeErrors(informationBestError.segmentation,obj.data,obj.parameters.polyDegree);
            informationBestError.features = computeMetrics(informationBestError.segmentation,obj.data,obj.parameters.characActivation,obj.parameters.polyDegree,obj.parameters.typeError);
            informationBestError.estimatedSerie = estimationSerie(informationBestError.segmentation,obj.data,obj.parameters.polyDegree);
            informationBestError.cuts = find(informationBestError.segmentation==1);
            [normCharac] = normalizeFunction(informationBestError.features);
            [informationBestError.L, informationBestError.C] = clusteringKmeans(normCharac,obj.parameters.k,obj.parameters.iterClust);
            informationBestError.parameters = obj.parameters;
            informationBestError.fitness = [currentFitness1; currentFitness2];
            informationBestError.number_of_firstFront = numel(find(f==1));
            if informationBestError.number_of_firstFront > obj.parameters.nPobl,
                informationBestError.number_of_firstFront = obj.parameters.nPobl;
            end
            informationBestError.spacing = spacing;
            informationBestError.m3 = m3;
            informationBestError.hv = hv;
            informationBestError.HRS = HRS;
            informationBestError.fbest = fbestidx;
            informationBestError.degree = obj.parameters.polyDegree;
        end
        
        %% Specific information of the algorithm
        function saveInformation(obj,model,dataset,repsuffix)
            outputFile = [repsuffix filesep dataset];
            f = fopen([outputFile '_info.csv'], 'wt');
            fprintf(f, 'Number of Cuts;%d\n', numel(model.cuts));
            fprintf(f, 'Number of Segments;%d\n',numel(model.cuts)+1);
            fprintf(f, 'Clustering Fitness Value;%f\n',model.fbestClustering);
            fprintf(f, 'Error Fitness Value;%f\n',model.fbestError);
            fprintf(f, 'RMSE;%f\n',model.errors(1));
            fprintf(f, 'RMSE;%f\n',model.errors(2));
            fprintf(f, 'RMSE;%f\n',model.errors(3));
            fprintf(f, 'Spacing;%f\n', model.spacing);
            fprintf(f, 'M3;%f\n', model.m3);
            fprintf(f, 'HV;%f\n', model.hv);
            fprintf(f, 'HRS;%f\n', model.HRS);
            fprintf(f, 'GMO parameters\n');
            fprintf(f, '%s\n', obj.getParameters());
            fclose(f);
        end
        
        %% Save all information of the algorithm
        function saveAll(obj,model,dataset,repsuffix)
            addpath(['..' filesep '..' filesep 'reporter' filesep]);
            addpath(['..' filesep '..' filesep 'reporter' filesep 'external_tools' filesep 'export_fig' filesep]);
            addpath(['..' filesep '..' filesep 'reporter' filesep 'external_tools' filesep 'plot2svg' filesep]);
            saveSegments(model,dataset,repsuffix);
            saveCentroids(model,dataset,repsuffix);
            saveEstimatedSerie(model,dataset,repsuffix);
            obj.saveInformation(model,dataset,repsuffix);
            plotSegmentedTimeSeries(model,'xlabel','ylabel',dataset,repsuffix,obj.data,model.estimatedSerie);
            plotSegmentedTimeSeries(model,'xlabel','ylabel',[dataset '_estimated'],repsuffix,model.estimatedSerie,obj.data);
            plotParetoFront(model,dataset,repsuffix);
        end
        
    end
    
    
end
