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
%   MASTER EXPERIMENTER
%   This script runs the TRADTSS algorithm given an experimental design
%   Results are saved in reports/ folder

% You should fill the files array with your time series
files ={'MALLAT_.txt','Donoho-Johnstone.txt'};
warning('off')

c = clock;
folder = ['reports' filesep num2str(c(1)) '-' num2str(c(2)) '-'  num2str(c(3)) '-' num2str(c(4)) '-' num2str(c(5)) '-' num2str(uint8(c(6)))];
mkdir('.',folder)
f = fopen([folder filesep 'info.csv'], 'wt');
fprintf(f,'Type;BBDD;# Segments; RMSE; Time\n');
for type=1:2,
    for file = 1:numel(files),
        input = ['..' filesep '..' filesep 'time_series' filesep char(files(file))];
        serie = load(input);
        alg = TRADTSS;
        alg.dataFile = char(files(file));
        alg.parameters.typeAlgorithm = type;
        alg.parameters.numSeg = round(0.025*numel(serie(:,2)));
        information = alg.runAlgorithm(serie(:,2));
        fprintf(f,'%d;%s,%d;%.3f;%.2f\n',type,alg.dataFile,numel(information.cuts)+1,information.errorsHA,information.timeHA);
        clear information;
    end
    fprintf(f,'\n');
end
fclose(f);



