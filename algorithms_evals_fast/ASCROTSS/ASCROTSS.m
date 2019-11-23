 classdef ASCROTSS < handle
    %ASCROTSS Approximation SCRO Time series segmentation [1]
    %
    %   ASCROTSS methods:
    %      runAlgorithm               - runs the corresponding algorithm and its hybrid version (SCRO and HSCRO in [1])
    %      saveInformation            - specific information of the algorithm
    %      saveAll                    - save all information of the algorithm
    %
    %   References:
    %     [1] A.M. Durán-Rosal, P.A. Gutiérrez, S. Salcedo-Sanz and C. Hervás-Martínez.
    %         "A statistically-driven Coral Reef Optimization algorithm for optimal
    %         size reduction of time series", Applied Soft Computing, 
    %         Vol. 63. 2018, pp. 139-153.
    %         https://doi.org/10.1016/j.asoc.2017.11.037
    %
    %   This file is part of TSSA: https://github.com/ayrna/tssa
    %   Original authors: Antonio M. Duran Rosal, Pedro A. Gutierrez Peña
    %   Citation: If you use this code, please cite the associated paper [1]
    %   Copyright:
    %       This software is released under the The GNU General Public License v3.0 licence
    %       available at http://www.gnu.org/licenses/gpl-3.0.html   
    properties
        name_parameters = {'numIt','nPobl','numSeg','pCross','pMut','seed','sizeChromosome','polyDegree','percentage_hybridation'}
        dataFile
        data
        parameters
    end
    
    methods
        %% Constructor
        function obj = ASCROTSS()
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
            % Percentage hybridation
            obj.parameters.percentage_hybridation = 0.40;
            % Type of error (MSE, SSE, MAXe)
            obj.parameters.typeError = 1;
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
            x=1:obj.parameters.sizeChromosome;
            y=transpose(serie);
            tic;
            [matrix] = precomputedMatrixes(x,y,obj.parameters.sizeChromosome);
            currentPopulation = initialisePopulation2(obj.parameters.nPobl,obj.parameters.sizeChromosome,obj.parameters.numSeg);
            
            %'Evaluation'
            oldFitness = zeros(1,obj.parameters.nPobl)*NaN;
            numberEvaluations = 0;
            numberEvaluations = numberEvaluations + numel(find(isnan(oldFitness)));
            currentFitness = evaluateFitnessErrorFast(currentPopulation,oldFitness,obj.parameters.nPobl,obj.parameters.sizeChromosome,x,y,matrix);
            
            med = mean(currentFitness);
            dev = std(currentFitness);
            
            ind = find(currentFitness < (med - dev));
            currentFitness(ind) = -1;
            currentPopulation(transpose(ind),:)=-1;
            
            i=1;
            while numberEvaluations < obj.parameters.numIt,
                %saveFitness(i,:) = currentFitness;
                % Asexual reproduction
                occupiedIndexes = find(currentFitness~=-1);
                med = mean(currentFitness(occupiedIndexes));
                dev = std(currentFitness(occupiedIndexes));
                Fa = numel(find(currentFitness > (med + 1*dev)));
                [asexualIndividual,asexualFitness]=selectionAsexual(currentPopulation,currentFitness,Fa/numel(occupiedIndexes));
                [asexualIndividual,asexualFitness] = mutation3(asexualIndividual,NaN,obj.parameters.pMut);
                
                % Sexual reproduction
                indCross = find(currentFitness > (med - dev));
                indMut = find(currentFitness <= (med - dev) & currentFitness ~= -1);
                if mod(numel(indCross),2) == 1,
                    indCross = [indCross indCross(1)];
                end
                
                %'Crossover'
                [poolCrossPopulation, poolCrossFitness] = crossoverStr2Op3(currentPopulation(indCross,:),currentFitness(indCross),...
                                                                           obj.parameters.pCross,3);
                
                %'Mutation'
                if (numel(indMut)==0),
                    poolMutPopulation = [];
                    poolMutFitness = [];
                else
                    [poolMutPopulation, poolMutFitness] = mutation3(currentPopulation(indMut,:),currentFitness(indMut),obj.parameters.pMut);
                end
                
                %'Evaluation'
                poolPopulation = [asexualIndividual; poolCrossPopulation; poolMutPopulation];
                poolFitness = [asexualFitness poolCrossFitness poolMutFitness];
                numberEvaluations = numberEvaluations + numel(find(isnan(poolFitness)));              
                poolFitness = evaluateFitnessErrorFast(poolPopulation,poolFitness,numel(poolFitness),obj.parameters.sizeChromosome,x,y,matrix);
                
                %'Coral replacement'
                [resultantPopulation, resultantFitness] = coralReplacement(currentPopulation,currentFitness,obj.parameters.nPobl,...
                                                                           poolPopulation,poolFitness,2);
                
                %'Depredation'
                occupiedIndexes = find(resultantFitness~=-1);
                med = mean(resultantFitness(occupiedIndexes));
                dev = std(resultantFitness(occupiedIndexes));
                Fd = numel(find(resultantFitness < (med - 2*dev) & resultantFitness ~= -1));
                Fd = Fd/numel(occupiedIndexes);
                [currentPopulation,currentFitness]=depredation(resultantPopulation,resultantFitness,Fd,1);

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
