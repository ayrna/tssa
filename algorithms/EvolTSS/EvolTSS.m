%EvolTSS Clustering Genetic Time series segmentation [1]
%
%   EvolTSS methods:
%      runAlgorithm               - runs the corresponding algorithm and (GA in [1])
%      saveInformation            - specific information of the algorithm
%      saveAll                    - save all information of the algorithm
%
%   References:
%     [1] M. Pérez-Ortiz, A.M. Durán-Rosal, P.A. Gutiérrez, et al.
%         "On the use of evolutionary time series analysis for segmenting paleoclimate data"
%         Neurocomputing, Vol. 326-327, January, 2019, pp. 3-14
%         https://doi.org/10.1016/j.neucom.2016.11.101
%
%   This file is part of TSSA: https://github.com/ayrna/tssa
%   Original authors: Antonio M. Duran Rosal, Pedro A. Gutierrez Peña
%   Citation: If you use this code, please cite the associated paper [1]
%   Copyright:
%       This software is released under the The GNU General Public License v3.0 licence
%       available at http://www.gnu.org/licenses/gpl-3.0.html
classdef EvolTSS < handle
    
    properties
        name_parameters = {'numIt','nPobl','k','pCross','pMut','seed','minSeg','maxSeg','sizeChromosome','iterClust','polyDegree','characActivation'}
        dataFile
        data
        parameters
    end
    
    methods
        %% Constructor
        function obj = EvolTSS()
            obj.defaultParameters();
        end
        
        %% Default parameters
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
            currentPopulation = initialisePopulation1(obj.parameters.nPobl,obj.parameters.sizeChromosome,obj.parameters.minSeg,obj.parameters.maxSeg);
            
            %'Evaluation'
            oldFitness = zeros(1,obj.parameters.nPobl)*NaN;
            currentFitness = evaluateFitnessClusteringKmeans(1,obj.parameters.typeFitness,currentPopulation,oldFitness,obj.parameters.k,...
                                                             obj.parameters.iterClust,obj.data,obj.parameters.characActivation,obj.parameters.polyDegree);
            
            information.meanFitness(1) = mean(currentFitness);
            information.bestFitness(1) = max(currentFitness);

            for i=1:obj.parameters.numIt,
                
                %'Crossover'
                [newPopulation, newFitness] = crossoverStr1Op1(currentPopulation,currentFitness,obj.parameters.pCross,...
                                                               obj.parameters.minSeg,obj.parameters.maxSeg,3);
                
                %'Mutation'
                [newPopulation, newFitness] = mutation1(newPopulation,newFitness,obj.parameters.pMut,obj.parameters.mutedPoints,...
                                                        obj.parameters.minSeg,obj.parameters.maxSeg);

                %'Evaluation'
                newFitness = evaluateFitnessClusteringKmeans(1,obj.parameters.typeFitness,newPopulation,newFitness,obj.parameters.k,...
                                                             obj.parameters.iterClust,obj.data,obj.parameters.characActivation,obj.parameters.polyDegree);
                
                %'Selection'
                [currentPopulation, currentFitness]= selection1Roulette([currentPopulation; newPopulation],[currentFitness newFitness],obj.parameters.nPobl);
                
                information.meanFitness(i+1) = mean(currentFitness);
                information.bestFitness(i+1) = max(currentFitness);
            end

            % Information for the reporter
            [information.fbest,fbestidx] = max(currentFitness);
            information.segmentation = currentPopulation(fbestidx,:);
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
            fprintf(f, 'GA Fitness Value;%f\n',model.fbest);
            fprintf(f, 'GA parameters\n');
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