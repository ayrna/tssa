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
    fprintf(fid,'#Run;NumberSegments;Init_RMSE;Init_RMSEp;Init_MAXe;Init_Fitnnes;');
    fprintf(fid,'GA_RMSE;GA_RMSEp;GA_MAXe;GA_Fitnnes;');
    fprintf(fid,'BU_RMSE;BU_RMSEp;BU_MAXe;BU_Fitnnes;');
    fprintf(fid,'HA_RMSE;HA_RMSEp;HA_MAXe;HA_Fitnnes\n');
    if (numel(model)>1),
        for i=1:numel(model),
            numSegments(i) = numel(model(i).cuts)+1;

            Init_RMSE(i) = model(i).errorsInit(1);
            Init_RMSEp(i) = model(i).errorsInit(2);
            Init_MAXe(i) = model(i).errorsInit(3);
            Init_Fitness(i) = model(i).bestFitness(1);

            GA_RMSE(i) = model(i).errorsGA(1);
            GA_RMSEp(i) = model(i).errorsGA(2);
            GA_MAXe(i) = model(i).errorsGA(3);
            GA_Fitness(i) = model(i).fitnessGA;

            BU_RMSE(i) = model(i).errorsBU(1);
            BU_RMSEp(i) = model(i).errorsBU(2);
            BU_MAXe(i) = model(i).errorsBU(3);
            BU_Fitness(i) = model(i).fitnessBU;

            HA_RMSE(i) = model(i).errorsHA(1);
            HA_RMSEp(i) = model(i).errorsHA(2);
            HA_MAXe(i) = model(i).errorsHA(3);
            HA_Fitness(i) = model(i).fitnessHA;

            fprintf(fid,'%d;%d;%f;%f;%f;%f;',i,numSegments(i),Init_RMSE(i),Init_RMSEp(i),Init_MAXe(i),Init_Fitness(i));
            fprintf(fid,'%f;%f;%f;%f;',GA_RMSE(i),GA_RMSEp(i),GA_MAXe(i),GA_Fitness(i));
            fprintf(fid,'%f;%f;%f;%f;',BU_RMSE(i),BU_RMSEp(i),BU_MAXe(i),BU_Fitness(i));
            fprintf(fid,'%f;%f;%f;%f\n',HA_RMSE(i),HA_RMSEp(i),HA_MAXe(i),HA_Fitness(i));
        end
        fprintf(fid,'Mean;%f;',mean(numSegments));
        fprintf(fid,'%f;%f;%f;%f;',mean(Init_RMSE),mean(Init_RMSEp),mean(Init_MAXe),mean(Init_Fitness));
        fprintf(fid,'%f;%f;%f;%f;',mean(GA_RMSE),mean(GA_RMSEp),mean(GA_MAXe),mean(GA_Fitness));
        fprintf(fid,'%f;%f;%f;%f;',mean(BU_RMSE),mean(BU_RMSEp),mean(BU_MAXe),mean(BU_Fitness));
        fprintf(fid,'%f;%f;%f;%f\n',mean(HA_RMSE),mean(HA_RMSEp),mean(HA_MAXe),mean(HA_Fitness));

        fprintf(fid,'Std;%f;',std(numSegments));
        fprintf(fid,'%f;%f;%f;%f;',std(Init_RMSE),std(Init_RMSEp),std(Init_MAXe),std(Init_Fitness));
        fprintf(fid,'%f;%f;%f;%f;',std(GA_RMSE),std(GA_RMSEp),std(GA_MAXe),std(GA_Fitness));
        fprintf(fid,'%f;%f;%f;%f;',std(BU_RMSE),std(BU_RMSEp),std(BU_MAXe),std(BU_Fitness));
        fprintf(fid,'%f;%f;%f;%f\n',std(HA_RMSE),std(HA_RMSEp),std(HA_MAXe),std(HA_Fitness));
    else
        fprintf(fid,'1;%d;',numel(model.cuts)+1);
        fprintf(fid,'%f;%f;%f;%f;',model.errorsInit(1),model.errorsInit(2),model.errorsInit(3),model.bestFitness(1));
        fprintf(fid,'%f;%f;%f;%f;',model.errorsGA(1),model.errorsGA(2),model.errorsGA(3),model.fitnessGA);
        fprintf(fid,'%f;%f;%f;%f;',model.errorsBU(1),model.errorsBU(2),model.errorsBU(3),model.fitnessBU);
        fprintf(fid,'%f;%f;%f;%f\n',model.errorsHA(1),model.errorsHA(2),model.errorsHA(3),model.fitnessHA);
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
