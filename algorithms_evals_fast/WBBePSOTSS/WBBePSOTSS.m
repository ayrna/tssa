 classdef WBBePSOTSS < handle
    %WBBePSOTSS Approximation WBBePSO Time series segmentation [1]
    %
    %   WBBePSOTSS methods:
    %      runAlgorithm               - runs the corresponding algorithm and its hybrid versions (WBBePSO and HWBBePSO in [1])
    %      saveInformation            - specific information of the algorithm
    %      saveAll                    - save all information of the algorithm
    %
    %   References:
    %     [1] A.M. Durán-Rosal, D. Guijo-Rubio, P.A. Gutiérrez and C. Hervás-Martínez.
    %		  "Hybrid Weighted Barebones Exploiting Particle Swarm Optimization Algorithm
    %		  for Time Series Representation". BIOMA2018. 16th-18th May. 2018. 
    %		  Paris (France). LNCS, vol. 10835. pp. 126-137
    %		  https://doi.org/10.1007/978-3-319-91641-5_11	
    %
    %   This file is part of TSSA: https://github.com/ayrna/tssa
    %   Original authors: Antonio M. Duran Rosal, Pedro A. Gutierrez Peña
    %   Citation: If you use this code, please cite the associated paper [1]
    %   Copyright:
    %       This software is released under the The GNU General Public License v3.0 licence
    %       available at http://www.gnu.org/licenses/gpl-3.0.html   
    properties
        name_parameters = {'numIt','nPobl','numSeg','seed','sizeChromosome','sizeChromosomeInt','polyDegree','percentage_hybridation','typeError'}
        dataFile
        data
        parameters
    end
    
    methods
        %% Constructor
        function obj = WBBePSOTSS()
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
            % Random number generation seed
            obj.parameters.seed = 1;
            % degree for approximations (0 - Interpolation, >=1 degree)
            obj.parameters.polyDegree = 0;
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
            obj.parameters.sizeChromosomeInt = obj.parameters.numSeg-1;
            
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
            matrix=precomputedMatrixes(x,y,obj.parameters.sizeChromosome);
            currentPopulation = initialisePopulation2(obj.parameters.nPobl,obj.parameters.sizeChromosome,obj.parameters.numSeg);
            currentPopulationInt = transformPopulationToInt(currentPopulation,obj.parameters.nPobl,obj.parameters.sizeChromosomeInt);
            bestLocalPopulationInt = currentPopulationInt;
            acumulateBestLocal = bestLocalPopulationInt;
            
            %'Evaluation'
            oldFitness = zeros(1,obj.parameters.nPobl)*NaN;
            numberEvaluations = 0;
            numberEvaluations = numberEvaluations + numel(find(isnan(oldFitness)));
            currentFitness = evaluateFitnessErrorFast(currentPopulation,oldFitness,obj.parameters.nPobl,obj.parameters.sizeChromosome,x,y,matrix);
            bestLocalFitness = currentFitness;
            [bestIndividualFitness, bidx]= max(currentFitness);
            bestIndividualInt = currentPopulationInt(bidx,:);
            acumulateBestGlobal = bestIndividualInt;
            
            %'Initial velocities'
            i=1;
            generations = 1;
            while numberEvaluations < obj.parameters.numIt,
                
                %'PSO'                
                [currentPopulationInt,newFitness] = updatePositionsWBB(currentPopulationInt,bestLocalPopulationInt,acumulateBestLocal,acumulateBestGlobal,...
                                                                       currentFitness,generations,obj.parameters.nPobl,obj.parameters.sizeChromosome,...
                                                                       obj.parameters.sizeChromosomeInt);
    
                %'Evaluation'
                currentPopulation = transformPopulationToBi(currentPopulationInt,obj.parameters.nPobl,obj.parameters.sizeChromosome);
                numberEvaluations = numberEvaluations + numel(find(isnan(newFitness)));
                currentFitness = evaluateFitnessErrorFast(currentPopulation,newFitness,obj.parameters.nPobl,obj.parameters.sizeChromosome,x,y,matrix);
                
                %'Update arrays'
                [bestLocalPopulationInt, bestLocalFitness, bestIndividualInt, bestIndividualFitness] = updateArrays(currentPopulationInt, currentFitness,...
                                                                                                                bestLocalPopulationInt, bestLocalFitness,...
                                                                                                                bestIndividualInt, bestIndividualFitness);
                acumulateBestLocal = acumulateBestLocal + (i+1) *bestLocalPopulationInt;
                acumulateBestGlobal = acumulateBestGlobal + (i+1) *bestIndividualInt;
                
                i=i+1;
                generations = generations + i;
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
