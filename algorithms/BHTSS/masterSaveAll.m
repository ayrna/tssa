%   COPYRIGHT
%   This file is part of TSSA: https://github.com/ayrna/tssa
%   Original authors: Antonio M. Duran Rosal, Pedro A. Gutierrez
%   Citation: If you use this code, please cite the associated paper [1]
%   Copyright:
%       This software is released under the The GNU General Public License v3.0 licence
%       available at http://www.gnu.org/licenses/gpl-3.0.html
%
%   References:
%     [1] A.M. Durán-Rosal, J.C. Fernández, P.A. Gutiérrez and C. Hervás-Martínez.
%         "Detection and prediction of segments containing extreme significant wave heights"
%         Ocean Engineering, Vol. 142, September, 2017, pp. 268-279.
%         https://doi.org/10.1016/j.oceaneng.2017.07.009
%
% MASTER SAVE ALL
% This function saves the information of multiple runnings given an
% experimental design for BHTSS algorithm.
% Summary information are saved in resultsMultipleRunnings.csv file
% Summary models are saved in informationMultipleRunnings.mat file
function masterSaveAll(folder,model)
    fid = fopen([folder filesep 'resultsMultipleRunnings.csv'],'wt');
    fprintf(fid,'#Run;NumberSegments;InitialFitness;GAFitness;GALSFitness\n');
    if (numel(model)>1),
        for i=1:numel(model),
            numSegments(i) = numel(model(i).cuts)+1;
            fitnessI(i) = model(i).bestFitness(1);
            fitnessGA(i) = model(i).fbestGA;
            fitness(i) = model(i).fbest;
            fprintf(fid,'%d;%d;%f;%f;%f\n',i,numSegments(i),fitnessI(i),fitnessGA(i),fitness(i));
        end
        fprintf(fid,'Mean;%f;%f;%f;%f\n',mean(numSegments),mean(fitnessI),mean(fitnessGA),mean(fitness));
        fprintf(fid,'Std;%f;%f;%f;%f\n',std(numSegments),std(fitnessI),std(fitnessGA),std(fitness));
    else
        fprintf(fid,'1;%d;%f;%f;%f',numel(model.cuts)+1,model.bestFitness(1),model.fbestGA,model.fbest);
    end
    fclose(fid);
    save([folder filesep 'informationMultipleRunnings.mat'], 'model');
end
