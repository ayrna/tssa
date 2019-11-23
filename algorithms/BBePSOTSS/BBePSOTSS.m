 classdef BBePSOTSS < handle
    %BBePSOTSS Approximation BBePSO Time series segmentation [1]
    %
    %   BBePSOTSS methods:
    %      runAlgorithm               - runs the corresponding algorithm and its hybrid versions (BBePSO and HBBePSO in [1])
    %      saveInformation            - specific information of the algorithm
    %      saveAll                    - save all information of the algorithm
    %
    %   References:
    %     [1] A.M. Durán-Rosal, P.A. Gutiérrez, Á. Carmona-Poyato and C. Hervás-Martínez.
    %         "A hybrid dynamic exploitation barebones particle swarm optimisation
    %         algorithm for time series segmentation", Neurocomputing,
    %         Vol. 353, August, 2019, pp. 45-55.
    %		  https://doi.org/10.1016/j.neucom.2018.05.129
    %
    %   This file is part of TSSA: https://github.com/ayrna/tssa
    %   Original authors: Antonio M. Duran Rosal, Pedro A. Gutierrez Peña
    %   Citation: If you use this code, please cite the associated paper [1]
    %   Copyright:
    %       This software is released under the The GNU General Public License v3.0 licence
    %       available at http://www.gnu.org/licenses/gpl-3.0.html
    properties
        name_parameters = {'numIt','nPobl','numSeg','C1','C2','W','seed','sizeChromosome','sizeChromosomeInt','polyDegree','percentage_hybridation','typeError'}
        dataFile
        data
        parameters
    end
    
    methods
        %% Constructor
        function obj = BBePSOTSS()
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
            % PSO constants
            obj.parameters.C1 = 0.5;
            obj.parameters.C2 = 0.5;
            obj.parameters.W = 0.72;
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
            currentPopulation = initialisePopulation2(obj.parameters.nPobl,obj.parameters.sizeChromosome,obj.parameters.numSeg);
            currentPopulationInt = transformPopulationToInt(currentPopulation,obj.parameters.nPobl,obj.parameters.sizeChromosomeInt);
            bestLocalPopulationInt = currentPopulationInt;
            
            
            %'Evaluation'
            oldFitness = zeros(1,obj.parameters.nPobl)*NaN;
            currentFitness = evaluateFitnessError(obj.parameters.typeError,currentPopulation,oldFitness,obj.data,obj.parameters.polyDegree);
            bestLocalFitness = currentFitness;
            [bestIndividualFitness, bidx]= max(currentFitness);
            bestIndividualInt = currentPopulationInt(bidx,:);
            
            
            information.meanFitness(1) = mean(currentFitness);
            information.stdFitness(1) = std(currentFitness);
            [information.bestFitness(1), idx] = max(currentFitness);
            information.bestAll(1) = information.bestFitness(1);
            chromosomeInit = currentPopulation(idx,:);
            
            %'Initial velocities'
            tic;
            for i=1:obj.parameters.numIt,
                
                %'PSO'
                [currentPopulationInt, newFitness] = updatePositionsBB(currentPopulationInt,bestLocalPopulationInt,bestIndividualInt,currentFitness,...
                                                                       obj.parameters.nPobl,obj.parameters.sizeChromosome,obj.parameters.sizeChromosomeInt);
                
                %'Evaluation'
                currentPopulation = transformPopulationToBi(currentPopulationInt,obj.parameters.nPobl,obj.parameters.sizeChromosome);
                currentFitness = evaluateFitnessError(obj.parameters.typeError,currentPopulation,newFitness,obj.data,obj.parameters.polyDegree);
                
                %'Update arrays'
                [bestLocalPopulationInt, bestLocalFitness, bestIndividualInt, bestIndividualFitness] = updateArrays(currentPopulationInt, currentFitness,...
                                                                                                                bestLocalPopulationInt, bestLocalFitness,...
                                                                                                                bestIndividualInt, bestIndividualFitness);
                
                
                information.meanFitness(i+1) = mean(currentFitness);
                information.stdFitness(i+1) = std(currentFitness);
                information.bestFitness(i+1) = max(currentFitness);
                information.bestAll(i+1) = bestIndividualFitness;
            end
            
            % Initial solution
            [errorsInit] = computeErrors(chromosomeInit,obj.data,obj.parameters.polyDegree);
            
            % GA solution
            fbestGA = bestIndividualFitness;
            chromosomeGA = transformPopulationToBi(bestIndividualInt,1,obj.parameters.sizeChromosome);
            [errorsGA] = computeErrors(chromosomeGA,obj.data,obj.parameters.polyDegree);
            timeGA=toc;
            
            % Bottom-Up solution
            max_iters = round(obj.parameters.percentage_hybridation*(numel(find(chromosomeGA==1))));
            [chromosomeBU] = hybridIndividualBottomUp(chromosomeGA,obj.data,max_iters,obj.parameters.polyDegree,obj.parameters.typeError);
            [errorsBU] = computeErrors(chromosomeBU,obj.data,obj.parameters.polyDegree);
            fbestBU = evaluateFitnessError(obj.parameters.typeError,chromosomeBU,NaN,obj.data,obj.parameters.polyDegree);
            
            % Top-Down solution (HA solution)
            chromosomeHA = hybridIndividualTopDown(chromosomeBU,obj.data,max_iters,obj.parameters.polyDegree,obj.parameters.typeError);
            [errorsHA] = computeErrors(chromosomeHA,obj.data,obj.parameters.polyDegree);
            fbestHA = evaluateFitnessError(obj.parameters.typeError,chromosomeHA,NaN,obj.data,obj.parameters.polyDegree);
            timeHA=toc;
            
            % Information for the reporter
            information.errorsInit = errorsInit;
            information.errorsGA = errorsGA;
            information.errorsBU = errorsBU;
            information.errorsHA = errorsHA;
            information.fitnessGA = fbestGA;
            information.fitnessBU = fbestBU;
            information.fitnessHA = fbestHA;
            information.segmentation = chromosomeHA;
            information.estimatedSerie = estimationSerie(information.segmentation,obj.data,obj.parameters.polyDegree);
            information.cuts = find(information.segmentation==1);
            information.parameters = obj.parameters;
            information.degree = obj.parameters.polyDegree;
            
            information.timeGA=timeGA;
            information.timeHA=timeHA;
            information.generations = i-1;
        end        

        
        %% Specific information of the algorithm
        function saveInformation(obj,model,dataset,repsuffix)
            outputFile = [repsuffix filesep dataset];
            f = fopen([outputFile '_info.csv'], 'wt');
            fprintf(f, 'Number of Cuts;%d\n', numel(model.cuts));
            fprintf(f, 'Number of Segments;%d\n',numel(model.cuts)+1);
            fprintf(f, 'Number of Generations;%d\n', model.generations);
            fprintf(f, 'Solution;RMSE;RMSEp;MAXe;fitness\n');
            fprintf(f, 'Initial solution;%f;%f;%f;%f\n',model.errorsInit(1),model.errorsInit(2),model.errorsInit(3),model.bestFitness(1));
            fprintf(f, 'GA solution;%f;%f;%f;%f\n',model.errorsGA(1),model.errorsGA(2),model.errorsGA(3),model.fitnessGA);
            fprintf(f, 'Bottom-Up solution;%f;%f;%f;%f\n',model.errorsBU(1),model.errorsBU(2),model.errorsBU(3),model.fitnessBU);
            fprintf(f, 'HA solution;%f;%f;%f;%f\n',model.errorsHA(1),model.errorsHA(2),model.errorsHA(3),model.fitnessHA);
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
