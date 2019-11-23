%   COPYRIGHT
%   This file is part of TSSA: https://github.com/ayrna/tssa
%   Original authors: Antonio M. Duran Rosal, Pedro A. Gutierrez
%   Citation: If you use this code, please cite the associated paper [1]
%   Copyright:
%       This software is released under the The GNU General Public License v3.0 licence
%       available at http://www.gnu.org/licenses/gpl-3.0.html
%
%   References:
%     [1] M. Pérez-Ortiz, A.M. Durán-Rosal, P.A. Gutiérrez, et al.
%         "On the use of evolutionary time series analysis for segmenting paleoclimate data"
%         Neurocomputing, Vol. 326-327, January, 2019, pp. 3-14
%         https://doi.org/10.1016/j.neucom.2016.11.101
%
% MASTER SAVE ALL
% This function saves the information of multiple runnings given an
% experimental design for EvolTSS algorithm.
% Summary information are saved in resultsMultipleRunnings.csv file
% Summary models are saved in informationMultipleRunnings.mat file
function masterSaveAll(folder,model)
    fid = fopen([folder filesep 'resultsMultipleRunnings.csv'],'wt');
    fprintf(fid,'#Run;NumberSegments;InitialFitness;GAFitness\n');
    if (numel(model) > 1),
        for i=1:numel(model),
            numSegments(i) = numel(model(i).cuts)+1;
            fitnessI(i) = model(i).bestFitness(1);
            fitness(i) = model(i).fbest;
            fprintf(fid,'%d;%d;%f;%f\n',i,numSegments(i),fitnessI(i),fitness(i));
        end
        fprintf(fid,'Mean;%f;%f;%f\n',mean(numSegments),mean(fitnessI),mean(fitness));
        fprintf(fid,'Std;%f;%f;%f\n',std(numSegments),std(fitnessI),std(fitness));
    else
        fprintf(fid,'1;%d;%f;%f\n',numel(model.cuts)+1,model.bestFitness(1),model.fbest);
    end
    fclose(fid);
    save([folder filesep 'informationMultipleRunnings.mat'], 'model');
end