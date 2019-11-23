 classdef AMCROTSS < handle
    %AMCROTSS Approximation MCRO Time series segmentation [1]
    %
    %   AMCROTSS methods:
    %      runAlgorithm               - runs the corresponding algorithm (DMCRO in [1])
    %      saveInformation            - specific information of the algorithm
    %      saveAll                    - save all information of the algorithm
    %
    %   References:
    %     [1] A.M. Durán-Rosal, P.A. Gutiérrez, S. Salcedo-Sanz and C. Hervás-Martínez.
    %         "Dynamical Memetization in Coral Reef Optimization Algorithms
    %         for Optimal Time Series Approximation",
    %         Progress in Artificial Intelligence, Vol. 8, June, 2019, pp. 253-262.
    %         https://doi.org/10.1007/s13748-019-00176-0
    %
    %   This file is part of TSSA: https://github.com/ayrna/tssa
    %   Original authors: Antonio M. Duran Rosal, Pedro A. Gutierrez Peña
    %   Citation: If you use this code, please cite the associated paper [1]
    %   Copyright:
    %       This software is released under the The GNU General Public License v3.0 licence
    %       available at http://www.gnu.org/licenses/gpl-3.0.html
    properties
        name_parameters = {'numIt','nPobl','numSeg','pCross','pMut','seed','sizeChromosome','polyDegree','percentage_hybridation','freePositions','Fa','Fb','Fd','pDep','Natt'}
        dataFile
        data
        parameters
    end
    
    methods
        %% Constructor
        function obj = AMCROTSS()
            obj.defaultParameters();
        end
        
        %% Default parameters
        function obj = defaultParameters(obj)
            % Number of generations
            obj.parameters.numIt = 200;
            % Population size
            obj.parameters.nPobl = 80;
            % Number of segments
            obj.parameters.numSeg = 80;
            % Crossover probability
            obj.parameters.pCross = 0.8;
            % Mutation probability
            obj.parameters.pMut = 0.2;
            % Random number generation seed
            obj.parameters.seed = 1;
            % degree for approximations (0 - Interpolation, >=1 degree)
            obj.parameters.polyDegree = 0;
            % Percentage hybridation
            obj.parameters.percentage_hybridation = 0.40;
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
                parameters_as_str = [parameters_as_str sprintf('%s;%f\n', fields{i}, parameters.(fields{i}))];
            end
            
        end
        
        %% Main algorithm
        function [information] = runAlgorithm(obj, serie)
            addpath(['..' filesep '..' filesep 'source_code' filesep]);
            
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
            x=1:obj.parameters.sizeChromosome;
            y=transpose(serie);
            tic;
            [matrix] = precomputedMatrixes(x,y,obj.parameters.sizeChromosome);
            currentPopulation = initialisePopulation4(obj.parameters.nPobl,obj.parameters.sizeChromosome,obj.parameters.numSeg,obj.parameters.freePositions);
            freeIndexes = currentPopulation(:,1)==-1;
            occupiedIndexes = find(currentPopulation(:,1)~=-1);
            for iter=1:(numel(occupiedIndexes)/2),
                randNumber=randi([1 numel(occupiedIndexes)]);
                i=occupiedIndexes(randNumber);
                % Bottom-Up solution
                max_iters = round(obj.parameters.percentage_hybridation*(numel(find(currentPopulation(i,:)==1))));
                [chromosomeBU] = hybridIndividualBottomUpFast(currentPopulation(i,:),max_iters,x,y,matrix);
                % Top-Down solution (HA solution)
                currentPopulation(i,:) = hybridIndividualTopDownFast(chromosomeBU,max_iters,x,y,matrix);
            end
            
            %'Evaluation'
            oldFitness = zeros(1,obj.parameters.nPobl)*NaN;
            oldFitness(freeIndexes)=-1;
            numberEvaluations = 0;
            numberEvaluations = numberEvaluations + numel(find(isnan(oldFitness)));
            currentFitness = evaluateFitnessErrorFast(currentPopulation,oldFitness,obj.parameters.nPobl,obj.parameters.sizeChromosome,x,y,matrix);
            
            m=1;
            i=1;
            while numberEvaluations < obj.parameters.numIt,
                % Asexual reproduction
                [asexualIndividual,asexualFitness] = selectionAsexual(currentPopulation,currentFitness,obj.parameters.Fa);
                [asexualIndividual,asexualFitness] = mutation3(asexualIndividual,NaN,obj.parameters.pMut);
                
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
                [poolCrossPopulation, poolCrossFitness] = crossoverStr2Op3(currentPopulation(occupiedIndexes(1:numberOfCrossed),:),...
                                                                           currentFitness(occupiedIndexes(1:numberOfCrossed)),...
                                                                           obj.parameters.pCross,3);
                
                %'Mutation'
                [poolMutPopulation, poolMutFitness] = mutation3(currentPopulation(occupiedIndexes(numberOfCrossed+1:end),:),...
                                                                currentFitness(occupiedIndexes(numberOfCrossed+1:end)),obj.parameters.pMut);
                
                %'Evaluation'
                poolPopulation = [asexualIndividual; poolCrossPopulation; poolMutPopulation];
                poolFitness = [asexualFitness poolCrossFitness poolMutFitness];
                % Memetic
                if ((numberEvaluations >= 0.25 * obj.parameters.numIt && m == 1) ||...
                    (numberEvaluations >= 0.50 * obj.parameters.numIt && m == 2) ||...
                    (numberEvaluations >= 0.75 * obj.parameters.numIt && m == 3)),
                      [trush,indBestSegmentationGA] = max(currentFitness);
                      [chromosome] = hybridIndividualBottomUpFast(currentPopulation(indBestSegmentationGA,:),max_iters,x,y,matrix);
                      chromosome = hybridIndividualTopDownFast(chromosome,max_iters,x,y,matrix);
                      [poolPopulation] = [poolPopulation; chromosome];
                      [poolFitness] = [poolFitness NaN];
                      m=m+1;
                end
                numberEvaluations = numberEvaluations + numel(find(isnan(poolFitness)));              
                poolFitness = evaluateFitnessErrorFast(poolPopulation,poolFitness,numel(poolFitness),obj.parameters.sizeChromosome,x,y,matrix);
                
                %'Coral replacement'
                [resultantPopulation, resultantFitness] = coralReplacement(currentPopulation,currentFitness,obj.parameters.nPobl,...
                                                                                         poolPopulation,poolFitness,2);
                
                %'Depredation'
                [currentPopulation,currentFitness]=depredation(resultantPopulation,resultantFitness,obj.parameters.Fd,obj.parameters.pDep);

                i=i+1;
            end
            
            % GA solution
            [fbestGA,indBestSegmentationGA] = max(currentFitness);
            chromosomeGA = currentPopulation(indBestSegmentationGA,:);
            timeGA=toc;
            
            % Bottom-Up solution
            max_iters = round(obj.parameters.percentage_hybridation*(numel(find(chromosomeGA==1))));
            [chromosomeBU] = hybridIndividualBottomUpFast(chromosomeGA,max_iters,x,y,matrix);
            %fbestBU = evaluateFitnessErrorFast(chromosomeBU,NaN,obj.parameters.nPobl,obj.parameters.sizeChromosome,x,y,matrix);
            
            % Top-Down solution (HA solution)
            chromosomeHA = hybridIndividualTopDownFast(chromosomeBU,max_iters,x,y,matrix);
            fbestHA = evaluateFitnessErrorFast(chromosomeHA,NaN,1,obj.parameters.sizeChromosome,x,y,matrix);
            timeHA=toc;
            
            % Information for the reporter
            information.errorsGA = (1/fbestGA)-1;
            information.errorsHA = (1/fbestHA)-1;
            information.fitnessGA = fbestGA;
            information.fitnessHA = fbestHA;
            information.segmentation = chromosomeHA;
            information.estimatedSerie = estimationSerie(information.segmentation,obj.data,0);
            information.cuts = find(information.segmentation==1);
            information.parameters = obj.parameters;
            
            information.timeGA=timeGA;
            information.timeHA=timeHA;
            information.numberEvaluations = numberEvaluations;
            information.generations = i-1;
        end
        
        %% Specific information of the algorithm
        function saveInformation(obj,model,dataset,repsuffix)
            outputFile = [repsuffix filesep dataset];
            f = fopen([outputFile '_info.csv'], 'wt');
            fprintf(f, 'Number of Cuts;%d\n', numel(model.cuts));
            fprintf(f, 'Number of Segments;%d\n',numel(model.cuts)+1);
            fprintf(f, 'Number of Evaluations;%d\n',model.numberEvaluations);
            fprintf(f, 'Number of Generations;%d\n',model.generations);
            fprintf(f, 'Solution;RMSEp;fitness;time\n');
            fprintf(f, 'GA solution;%f;%f;%f\n',model.errorsGA,model.fitnessGA,model.timeGA);
            fprintf(f, 'HA solution;%f;%f;%f;%f\n',model.errorsHA,model.fitnessHA,model.timeHA);
            fprintf(f, 'HA parameters\n');
            fprintf(f, '%s\n', obj.getParameters());
            fclose(f);
        end
        
         %% Save all information of the algorithm
        function saveAll(obj,model,dataset,repsuffix)
            addpath(['..' filesep '..' filesep 'reporter' filesep]);
            addpath(['..' filesep '..' filesep 'reporter' filesep 'external_tools' filesep 'export_fig' filesep]);
            addpath(['..' filesep '..' filesep 'reporter' filesep 'external_tools' filesep 'plot2svg' filesep]);
            obj.saveInformation(model,dataset,repsuffix);
            saveEstimatedSerie(model,dataset,repsuffix);
            plotApproximatedTimeSeries(model,'xlabel','ylabel',dataset,repsuffix,model.estimatedSerie,obj.data);
        end
    end
    
end
