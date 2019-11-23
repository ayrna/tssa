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
nOfruns = 1;
typeE = [1, 2];
for file = 1:numel(files),
    input = ['..' filesep '..' filesep 'time_series' filesep char(files(file))];
    serie = load(input);
    c = clock;
    folder = ['reports' filesep num2str(c(1)) '-' num2str(c(2)) '-'  num2str(c(3)) '-' num2str(c(4)) '-' num2str(c(5)) '-' num2str(uint8(c(6)))];
    mkdir('.',folder)
    for i=1:nOfruns,
        alg(i) = TRADTSS;
        alg(i).dataFile = char(files(file));
        alg(i).parameters.typeAlgorithm = 6;
        alg(i).parameters.maxError = 24000;
        alg(i).parameters.numSeg = round(0.025*numel(serie(:,2)));
        alg(i).parameters.polyDegree = 0;
        alg(i).parameters.typeError = 2;
        alg(i).parameters.buffer_size = 188;
        information(i) = alg(i).runAlgorithm(serie(:,2));
        mkdir('.',[folder filesep num2str(i)])
        alg(i).saveAll(information(i),char(files(file)),[folder filesep num2str(i)]);
    end
    masterSaveAll(folder,information)
    clear information;
end



