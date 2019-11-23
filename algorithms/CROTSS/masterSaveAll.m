%   COPYRIGHT
%   This file is part of TSSA: https://github.com/ayrna/tssa
%   Original authors: Antonio M. Duran Rosal, Pedro A. Gutierrez
%   Citation: If you use this code, please cite the associated paper [1]
%   Copyright:
%       This software is released under the The GNU General Public License v3.0 licence
%       available at http://www.gnu.org/licenses/gpl-3.0.html
%
%   References:
%     [1] A.M. Durán-Rosal, D. Guijo-Rubio, P.A. Gutiérrez, S. Salcedo-Sanz and C. Hervás-Martínez.
%         "A coral reef optimization algorithm for wave height time series segmentation problems".
%         International Work-Conference on Artificial and Natural Neural Networks (IWANN2017).
%         14th-16th June. 2017. Cádiz (Spain). LNCS, vol. 10305. pp. 673-684
%         https://doi.org/10.1007/978-3-319-59153-7_58
%
% MASTER SAVE ALL
% This function saves the information of multiple runnings given an
% experimental design for CROTSS algorithm.
% Summary information are saved in resultsMultipleRunnings.csv file
% Summary models are saved in informationMultipleRunnings.mat file
function masterSaveAll(folder,model)
    fid = fopen([folder filesep 'resultsMultipleRunnings.csv'],'wt');
    fprintf(fid,'#Run;NumberSegments;InitialFitness;CROFitness\n');
    if (numel(model)>1),
        for i=1:numel(model),
            numSegments(i) = numel(model(i).cuts)+1;
            fitnessI(i) = model(i).bestFitness(1);
            fitness(i) = model(i).fbest;
            fprintf(fid,'%d;%d;%f;%f\n',i,numSegments(i),fitnessI(i),fitness(i));
        end
        fprintf(fid,'Mean;%f;%f;%f;%f\n',mean(numSegments),mean(fitnessI),mean(fitness));
        fprintf(fid,'Std;%f;%f;%f\n',std(numSegments),std(fitnessI),std(fitness));
    else
        fprintf(fid,'1;%d;%f;%f',numel(model.cuts)+1,model.bestFitness(1),model.fbest);
    end
    fclose(fid);
    save([folder filesep 'informationMultipleRunnings.mat'], 'model');
end
