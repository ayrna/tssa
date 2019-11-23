%   COPYRIGHT
%   This file is part of TSSA: https://github.com/ayrna/tssa
%   Original authors: Antonio M. Duran Rosal, Pedro A. Gutierrez
%   Citation: If you use this code, please cite the associated paper [1,2,3]
%   Copyright:
%       This software is released under the The GNU General Public License v3.0 licence
%       available at http://www.gnu.org/licenses/gpl-3.0.html
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
% MASTER SAVE ALL
% This function saves the information of multiple runnings given an
% experimental design for TRADTSS algorithm.
% Summary information are saved in resultsMultipleRunnings.csv file
% Summary models are saved in informationMultipleRunnings.mat file
function masterSaveAll(folder,model)
    fid = fopen([folder filesep 'resultsMultipleRunnings.csv'],'wt');
    fprintf(fid,'#Run;NumberSegments;RMSE;RMSEp;MAXe;Fitnnes\n');
    if (numel(model)>1),
        for i=1:numel(model),
            numSegments(i) = numel(model(i).cuts)+1;

            RMSE(i) = model(i).errors(1);
            RMSEp(i) = model(i).errors(2);
            MAXe(i) = model(i).errors(3);
            Fitness(i) = model(i).fbest;

            fprintf(fid,'%d;%d;%f;%f;%f;%f\n',i,numSegments(i),RMSE(i),RMSEp(i),MAXe(i),Fitness(i));
        end
        fprintf(fid,'Mean;%f;',mean(numSegments));
        fprintf(fid,'%f;%f;%f;%f\n',mean(RMSE),mean(RMSEp),mean(MAXe),mean(Fitness));
        fprintf(fid,'Std;%f;',std(numSegments));
        fprintf(fid,'%f;%f;%f;%f\n',std(RMSE),std(RMSEp),std(MAXe),std(Fitness));
    else
        fprintf(fid,'1;%d;%f;%f;%f;%f\n',numel(model.cuts)+1,model.errors(1),model.errors(2),model.errors(3),model.fbest);
    end
    
    fclose(fid);
    save([folder filesep 'informationMultipleRunnings.mat'], 'model');
end
