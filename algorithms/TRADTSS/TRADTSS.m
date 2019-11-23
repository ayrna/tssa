 classdef TRADTSS < handle
%TRADTSS Approximation algorithms extracted from [1] and changed in [2,3]
%
%   TRADSS methods:
%      runAlgorithm               - runs the corresponding algorithm (Bottom-Up, Top-Down and SWAB in [1,2,3])
%      saveInformation            - specific information of the algorithm
%      saveAll                    - save all information of the algorithm
%
%   References:
%     [1] E. Keogh, S. Chu, D. Hart and M. Pazzani.
%         "Segmenting time series: A survey and novel approach",
%		  In Data mining in time series databases, 2004, pp.1-21.
%         https://doi.org/10.1142/9789812565402_0001
%     [2] A.M. Durán-Rosal, P.A. Gutiérrez, S. Salcedo-Sanz and C. Hervás-Martínez.
%         "A statistically-driven Coral Reef Optimization algorithm for optimal
%         size reduction of time series", Applied Soft Computing, 
%         Vol. 63. 2018, pp. 139-153.
%         https://doi.org/10.1016/j.asoc.2017.11.037
%     [3] A.M. Durán-Rosal, P.A. Gutiérrez, Á. Carmona-Poyato and C. Hervás-Martínez.
%         "A hybrid dynamic exploitation barebones particle swarm optimisation
%         algorithm for time series segmentation", Neurocomputing,
%         Vol. 353, August, 2019, pp. 45-55.
%		  https://doi.org/10.1016/j.neucom.2018.05.129
%
%   This file is part of TSSA: https://github.com/ayrna/tssa
%   Original authors: Antonio M. Duran Rosal, Pedro A. Gutierrez Peña
%   Citation: If you use this code, please cite the associated papers [1,2,3]
%   Copyright:
%       This software is released under the The GNU General Public License v3.0 licence
%       available at http://www.gnu.org/licenses/gpl-3.0.html  
    properties
        name_parameters = {'typeAlgorithm','maxError','numSeg','polyDegree','typeError','buffer_size','sizeChromosome'}
        dataFile
        data
        parameters
    end
    
    methods
        %% Constructor
        function obj = TRADTSS()
            obj.defaultParameters();
        end
        
        %% Default parameters
        function obj = defaultParameters(obj)
            % Type Algorithm
            obj.parameters.typeAlgorithm = 1;
            % Maximum error allowed
            obj.parameters.maxError = 80;
            % Number of segments
            obj.parameters.numSeg = 20;
            % Degree for approximations (0 - Interpolation, >=1 degree)
            obj.parameters.polyDegree = 0;
            % Type of error (MSE, SSE, MAXe)
            obj.parameters.typeError = 1;
            % Size of buffer for SWAB
            obj.parameters.buffer_size = 20;
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
                       
            chromosome = runTraditionals(obj.parameters.typeAlgorithm,obj.data,obj.parameters.maxError,obj.parameters.numSeg,obj.parameters.polyDegree,obj.parameters.typeError,obj.parameters.buffer_size);
            [errors] = computeErrors(chromosome,obj.data,obj.parameters.polyDegree);
            
            %'Evaluation'
            paramsEval.typeError = obj.parameters.typeError;
            paramsEval.degree = obj.parameters.polyDegree;
            paramsEval.serie = obj.data;
            paramsEval.population = chromosome;
            paramsEval.oldFitness = NaN;
            fbest = evaluate(2,paramsEval);
                        
            % Information for the reporter
            information.errors = errors;
            information.fbest = fbest;
            information.segmentation = chromosome;
            information.estimatedSerie = estimationSerie(information.segmentation,obj.data,obj.parameters.polyDegree);
            information.cuts = find(information.segmentation==1);
            information.parameters = obj.parameters;
            information.degree = obj.parameters.polyDegree;
        end
        
        %% Specific information of the algorithm
        function saveInformation(obj,model,dataset,repsuffix)
            outputFile = [repsuffix filesep dataset];
            f = fopen([outputFile '_info.csv'], 'wt');
            fprintf(f, 'Number of Cuts;%d\n', numel(model.cuts));
            fprintf(f, 'Number of Segments;%d\n',numel(model.cuts)+1);
            fprintf(f, 'RMSE;RMSEp;MAXe;fitness\n');
            fprintf(f, '%f;%f;%f;%f\n',model.errors(1),model.errors(2),model.errors(3),model.fbest);
            fprintf(f, 'TRAD parameters\n');
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
