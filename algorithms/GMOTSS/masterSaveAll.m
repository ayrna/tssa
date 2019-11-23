%   COPYRIGHT
%   This file is part of TSSA: https://github.com/ayrna/tssa
%   Original authors: Antonio M. Duran Rosal, Pedro A. Gutierrez
%   Citation: If you use this code, please cite the associated paper [1]
%   Copyright:
%       This software is released under the The GNU General Public License v3.0 licence
%       available at http://www.gnu.org/licenses/gpl-3.0.html
%
%     [1] A.M. Durán-Rosal, P.A. Gutiérrez, F.J. Martínez-Estudillo and C. Hervás-Martínez.
%         "Simultaneous optimisation of clustering quality and approximation error
%         for time series segmentation", Information Sciences, Vol. 442-443, May, 2018, pp. 186-201.
%         https://doi.org/10.1016/j.ins.2018.02.041
%
% MASTER SAVE ALL
% This function saves the information of multiple runnings given an
% experimental design for GMOTSS algorithm.
% Summary information are saved in resultsMultipleRunnings.csv file
% Summary models are saved in informationMultipleRunnings.mat file
function masterSaveAll(folder,filename,model)
    fid = fopen([folder filesep filename '.csv'],'wt');
    fprintf(fid,'#Run;NumberSegments;FitnessClustering;FitnessError;RMSE;RMSEp;MAXe;Spacing;M3;HV;HRS\n');
    if numel(model)>1,
        for i=1:numel(model),
            resultsRMSE(i) = model(i).errors(1);
            resultsRMSEp(i) = model(i).errors(2);
            resultsMAXe(i) = model(i).errors(3);
            resultsClustering(i) = model(i).fbestClustering;
            resultsError(i) = model(i).fbestError;
            numSegments(i) = size(model(i).cuts,2) + 1;
            spacing(i) = model(i).spacing;
            m3(i) = model(i).m3;
            hv(i) = model(i).hv;
            HRS(i) = model(i).HRS;
            fprintf(fid,'%d;%d;%f;%f;%f;%f;%f;%f;%f;%f;%f\n',i,numSegments(i),resultsClustering(i),resultsError(i),resultsRMSE(i),resultsRMSEp(i),resultsMAXe(i),spacing(i),m3(i),hv(i),HRS(i));
        end

        fprintf(fid,'Mean;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f\n',mean(numSegments),mean(resultsClustering),mean(resultsError),mean(resultsRMSE),mean(resultsRMSEp),mean(resultsMAXe),mean(spacing),mean(m3),mean(hv),mean(HRS));
        fprintf(fid,'Std;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f\n',std(numSegments),std(resultsClustering),std(resultsError),std(resultsRMSE),std(resultsRMSEp),std(resultsMAXe),std(spacing),std(m3),std(hv),std(HRS));
    else
        fprintf(fid,'1;%d;%f;%f;%f;%f;%f;%f;%f;%f;%f\n',size(model.cuts,2)+1,model.fbestClustering,model.fbestError,model.errors(1),model.errors(2),model.errors(3),model.spacing,model.m3,model.hv,model.HRS);
    end
    save([folder filesep filename '.mat'], 'model');
    fclose(fid);
end
