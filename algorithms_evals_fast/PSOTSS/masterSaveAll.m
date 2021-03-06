%   COPYRIGHT
%   This file is part of TSSA: https://github.com/ayrna/tssa
%   Original authors: Antonio M. Duran Rosal, Pedro A. Gutierrez
%   Citation: If you use this code, please cite the associated paper [1]
%   Copyright:
%       This software is released under the The GNU General Public License v3.0 licence
%       available at http://www.gnu.org/licenses/gpl-3.0.html
%
%   References:
%     [1] A.M. Durán-Rosal, P.A. Gutiérrez, Á. Carmona-Poyato and C. Hervás-Martínez.
%         "A hybrid dynamic exploitation barebones particle swarm optimisation
%         algorithm for time series segmentation", Neurocomputing,
%         Vol. 353, August, 2019, pp. 45-55.
%		  https://doi.org/10.1016/j.neucom.2018.05.129
%
% MASTER SAVE ALL
% This function saves the information of multiple runnings given an
% experimental design for PSOTSS algorithm.
% Summary information are saved in resultsMultipleRunnings.csv file
% Summary models are saved in informationMultipleRunnings.mat file
function masterSaveAll(folder,model)
    fid = fopen([folder filesep 'resultsMultipleRunnings.csv'],'wt');
    fprintf(fid,'#Run;NumberSegments;');
    fprintf(fid,'GA_RMSEp;GA_Fitnnes;');
    fprintf(fid,'HA_RMSEp;HA_Fitnnes\n');
    if (numel(model)>1),
        for i=1:numel(model),
            numSegments(i) = numel(model(i).cuts)+1;

            GA_RMSEp(i) = model(i).errorsGA;
            GA_Fitness(i) = model(i).fitnessGA;

            HA_RMSEp(i) = model(i).errorsHA;
            HA_Fitness(i) = model(i).fitnessHA;

            fprintf(fid,'%d;%d;',i,numSegments(i));
            fprintf(fid,'%f;%f;',GA_RMSEp(i),GA_Fitness(i));
            fprintf(fid,'%f;%f\n',HA_RMSEp(i),HA_Fitness(i));
        end
        fprintf(fid,'Mean;%f;',mean(numSegments));
        fprintf(fid,'%f;%f;',mean(GA_RMSEp),mean(GA_Fitness));
        fprintf(fid,'%f;%f\n',mean(HA_RMSEp),mean(HA_Fitness));

        fprintf(fid,'Std;%f;',std(numSegments));
        fprintf(fid,'%f;%f;',std(GA_RMSEp),std(GA_Fitness));
        fprintf(fid,'%f;%f;\n',std(HA_RMSEp),std(HA_Fitness));
    else
        fprintf(fid,'1;%d;',numel(model.cuts)+1);
        fprintf(fid,'%f;%f;%f;%f;',model.errorsGA,model.fitnessGA);
        fprintf(fid,'%f;%f;%f;%f\n',model.errorsHA,model.fitnessHA);
    end
    fclose(fid);
    
    fid = fopen([folder filesep 'times.csv'],'wt');
    fprintf(fid,'#Run;timeGA;timeHA\n');
    if (numel(model)>1),
        for i=1:numel(model),
            timeGA(i)=model(i).timeGA;
            timeHA(i)=model(i).timeHA;
            fprintf(fid,'%d;%f;%f\n',i,timeGA(i),timeHA(i));
        end
        fprintf(fid,'Mean;%f;%f\n',mean(timeGA),mean(timeHA));
        fprintf(fid,'Std;%f;%f\n',std(timeGA),std(timeHA));

    else
        fprintf(fid,'1;%f;%f\n',model.timeGA,model.timeHA);
    end
    fclose(fid);
    
    save([folder filesep 'informationMultipleRunnings.mat'], 'model');
end
