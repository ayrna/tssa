classdef CROTSS < handle
%CROTSS Clustering CRO Time series segmentation [1]
%
%   CRO methods:
%      runAlgorithm               - runs the corresponding algorithm  (CRO in [1])
%      saveInformation            - specific information of the algorithm
%      saveAll                    - save all information of the algorithm
%
%   References:
%     [1] A.M. Durán-Rosal, D. Guijo-Rubio, P.A. Gutiérrez, S. Salcedo-Sanz and C. Hervás-Martínez.
%         "A coral reef optimization algorithm for wave height time series segmentation problems".
%         International Work-Conference on Artificial and Natural Neural Networks (IWANN2017).
%         14th-16th June. 2017. Cádiz (Spain). LNCS, vol. 10305. pp. 673-684
%         https://doi.org/10.1007/978-3-319-59153-7_58
%
%
%   This file is part of TSSA: https://github.com/ayrna/tssa
%   Original authors: Antonio M. Duran Rosal, Pedro A. Gutierrez Peña
%   Citation: If you use this code, please cite the associated paper [1]
%   Copyright:
%       This software is released under the The GNU General Public License v3.0 licence
%       available at http://www.gnu.org/licenses/gpl-3.0.html
    properties
        name_parameters = {'numIt','nPobl','k','seed','minSeg','maxSeg','sizeChromosome','iterClust','polyDegree','characActivation','umbralEntropy','intervalLeft','intervalRight','freePositions','Fa','Fb','Fd','pDep','Natt'}
        dataFile
        data
        parameters
    end
    
    methods
        %% Constructor
        function obj = CROTSS()
            obj.defaultParameters();
        end
        
        %% Default parameters
        function obj = defaultParameters(obj)
            % Number of generations
            obj.parameters.numIt = 200;
            % Population size
            obj.parameters.nPobl = 80;
            % crossover probability
            % obj.parameters.pCross = 0.8;
            % mutation probability
            % obj.parameters.pMut = 0.2;
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
            % umbral for beta hybridization
            obj.parameters.umbralEntropy = 0.1;
            % left interval for the normalization
            obj.parameters.intervalLeft = 0.05;
            % right interval for the normalization
            obj.parameters.intervalRight = 0.95;
            % Specific parameters for CRO
            % Number of free positions
            obj.parameters.freePositions = 20;
            % Percentage of asexual reproduction
            obj.parameters.Fa = 0.2;
            % Percentage of sexual reproduction (ext)
            obj.parameters.Fb = 0.5;
            % Percentage of depredation
            obj.parameters.Fd = 0.1;
            % Probability of depredetation
            obj.parameters.pDep = 0.01;
            % Maximum number of attempts to replacement corals
            obj.parameters.Natt = 3;
            
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
        function [information] = runAlgorithm(obj, serie)
            addpath(['..' filesep '..' filesep 'source_code' filesep]);
            addpath(['..' filesep '..' filesep 'source_code' filesep 'kmeans' filesep]);

            obj.data = serie;
            nOfData = length(serie);
            obj.parameters.sizeChromosome = nOfData;
            
            % Seed
            if strcmp(version('-release'),'2013a')
                s = RandStream('mt19937ar','Seed',obj.parameters.seed);
                RandStream.setGlobalStream(s);
            else
                s = RandStream.create('mt19937ar','seed',obj.parameters.seed);
                RandStream.setDefaultStream(s);
            end
                        
            %'Initialisation'
            currentPopulation = initialisePopulation3(obj.parameters.nPobl,obj.parameters.sizeChromosome,obj.parameters.minSeg,...
                                                      obj.parameters.maxSeg,obj.parameters.freePositions);
            freeIndexes = currentPopulation(:,1)==-1;
            
            %'Evaluation'
            oldFitness = zeros(1,obj.parameters.nPobl)*NaN;
            oldFitness(freeIndexes)=-1;
            currentFitness = evaluateFitnessClusteringKmeans(1,obj.parameters.typeFitness,currentPopulation,oldFitness,obj.parameters.k,...
                                                             obj.parameters.iterClust,obj.data,obj.parameters.characActivation,obj.parameters.polyDegree);
            
            occupiedIndexes = currentFitness~=-1;
            information.meanFitness(1) = mean(currentFitness(occupiedIndexes));
            information.bestFitness(1) = max(currentFitness(occupiedIndexes));

            for i=1:obj.parameters.numIt,
                
                % Asexual reproduction
                [asexualIndividual,asexualFitness]=selectionAsexual(currentPopulation,currentFitness,obj.parameters.Fa);
                
                % Sexual reproduction
                occupiedIndexes = find(currentFitness~=-1);
                randIndexes = occupiedIndexes(randperm(length(occupiedIndexes)));
                occupiedIndexes = randIndexes;
                nOccupied = numel(occupiedIndexes);
                numberOfCrossed = round(nOccupied*obj.parameters.Fb);
                if mod(numberOfCrossed,2) == 1,
                    numberOfCrossed = numberOfCrossed - 1;
                end
                
                %'Crossover'
                [poolCrossPopulation, poolCrossFitness] = crossoverStr2Op1(currentPopulation(occupiedIndexes(1:numberOfCrossed),:),...
                                                                          currentFitness(occupiedIndexes(1:numberOfCrossed)),1.1,...
                                                                          obj.parameters.minSeg,obj.parameters.maxSeg,3);
                
                %'Mutation'
                [poolMutPopulation, poolMutFitness] = mutation1(currentPopulation(occupiedIndexes(numberOfCrossed+1:end),:),...
                                                                currentFitness(occupiedIndexes(numberOfCrossed+1:end)),1.1,...
                                                                obj.parameters.mutedPoints,obj.parameters.minSeg,obj.parameters.maxSeg);

                %'Evaluation'
                poolPopulation = [asexualIndividual; poolCrossPopulation; poolMutPopulation];
                poolFitness = [asexualFitness poolCrossFitness poolMutFitness];
                poolFitness = evaluateFitnessClusteringKmeans(1,obj.parameters.typeFitness,poolPopulation,poolFitness,obj.parameters.k,...
                                                             obj.parameters.iterClust,obj.data,obj.parameters.characActivation,obj.parameters.polyDegree);
                
                %'Selection'
                [resultantPopulation, resultantFitness] = coralReplacement(currentPopulation,currentFitness,obj.parameters.nPobl,...
                                                                                         poolPopulation,poolFitness,obj.parameters.Natt);
                
                %'Depredation'
                [currentPopulation,currentFitness]=depredation(resultantPopulation,resultantFitness,obj.parameters.Fd,obj.parameters.pDep);
                
                occupiedIndexes = find(currentFitness~=-1);
                information.meanFitness(i+1) = mean(currentFitness(occupiedIndexes));
                information.bestFitness(i+1) = max(currentFitness(occupiedIndexes));
            end
            
            validIndexes = find(currentFitness~=-1);
            finalPopulation = currentPopulation(validIndexes,:);
            finalFitness = currentFitness(validIndexes);
            currentPopulation = finalPopulation;
            currentFitness = finalFitness;
            
            % Hybridization assuming normal distribution
            [fbest,fbestidxGA] = max(currentFitness);
            chromosome = currentPopulation(fbestidxGA,:);

            
            % Information for the reporter
            information.fbest = fbest;
            information.segmentation = chromosome;
            information.features = computeMetrics(information.segmentation,obj.data,obj.parameters.characActivation,obj.parameters.polyDegree,1);
            information.cuts = find(information.segmentation==1);
            [normCharac] = normalizeFunction(information.features);
            [information.L, information.C] = clusteringKmeans(normCharac,obj.parameters.k,obj.parameters.iterClust);
            information.parameters = obj.parameters;
            information.degree = obj.parameters.polyDegree;
        end
        
        %% Specific information of the algorithm
        function saveInformation(obj,model,dataset,repsuffix)
            outputFile = [repsuffix filesep dataset];
            f = fopen([outputFile '_info.csv'], 'wt');
            fprintf(f, 'Number of Cuts;%d\n', numel(model.cuts));
            fprintf(f, 'Number of Segments;%d\n',numel(model.cuts)+1);
            fprintf(f, 'Initial Fitness Value;%f\n',model.bestFitness(1));
            fprintf(f, 'CRO Fitness Value;%f\n',model.fbest);
            fprintf(f, 'CRO parameters\n');
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
            obj.saveInformation(model,dataset,repsuffix);
            plotSegmentedTimeSeries(model,'xlabel','ylabel',dataset,repsuffix,obj.data,obj.data);
        end
    end
    
end
