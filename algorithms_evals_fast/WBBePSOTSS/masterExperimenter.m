%   COPYRIGHT
%   This file is part of TSSA: https://github.com/ayrna/tssa
%   Original authors: Antonio M. Duran Rosal, Pedro A. Gutierrez
%   Citation: If you use this code, please cite the associated paper [1]
%   Copyright:
%       This software is released under the The GNU General Public License v3.0 licence
%       available at http://www.gnu.org/licenses/gpl-3.0.html
%
%   References:
%     [1] A.M. Durán-Rosal, D. Guijo-Rubio, P.A. Gutiérrez and C. Hervás-Martínez.
%		  "Hybrid Weighted Barebones Exploiting Particle Swarm Optimization Algorithm
%		  for Time Series Representation". BIOMA2018. 16th-18th May. 2018. 
%		  Paris (France). LNCS, vol. 10835. pp. 126-137
%		  https://doi.org/10.1007/978-3-319-91641-5_11
%
%   MASTER EXPERIMENTER
%   This script runs the WBBePSOTSS algorithm given an experimental design
%   Results are saved in reports/ folder

% You should fill the files array with your time series
files ={'MALLAT_.txt','Donoho-Johnstone.txt'};
warning('off')
nOfruns = 30;

for file = 1:numel(files),
    for percentage=4:4,
        input = ['..' filesep '..' filesep 'time_series' filesep char(files(file))];
        serie = load(input);
        c = clock;
        folder = ['reports' filesep num2str(c(1)) '-' num2str(c(2)) '-'  num2str(c(3)) '-' num2str(c(4)) '-' num2str(c(5)) '-' num2str(uint8(c(6)))];
        mkdir('.',folder)
        for i=1:nOfruns,
            alg(i) = WBBePSOTSS;
            alg(i).dataFile = char(files(file));
            alg(i).parameters.numIt = numel(serie(:,2))*3.5;
            alg(i).parameters.nPobl = 100;
            alg(i).parameters.numSeg = round(0.025*numel(serie(:,2)));
            alg(i).parameters.seed = i*10;
            alg(i).parameters.percentage_hybridation=0.20+percentage*0.05;
            information(i) = alg(i).runAlgorithm(serie(:,2));
            mkdir('.',[folder filesep num2str(i)])
            alg(i).saveInformation(information(i),char(files(file)),[folder filesep num2str(i)]);
        end
        masterSaveAll(folder,information)
        clear information;
    end
end
