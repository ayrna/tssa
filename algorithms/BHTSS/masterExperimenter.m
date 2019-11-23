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
%   MASTER EXPERIMENTER
%   This script runs the BHTSS algorithm given an experimental design
%   Results are saved in reports/ folder

% You should fill the files array with your time series
files ={'MALLAT_.txt','Donoho-Johnstone.txt'};
clusterArray = [4, 4];
warning('off')
nOfruns = 30;

for file = 1:numel(files),
    for fit = 3:3,
        input = ['..' filesep '..' filesep 'time_series' filesep char(files(file))];
        serie = load(input);
        c = clock;
        folder = ['reports' filesep num2str(c(1)) '-' num2str(c(2)) '-'  num2str(c(3)) '-' num2str(c(4)) '-' num2str(c(5)) '-' num2str(uint8(c(6)))];
        mkdir('.',folder)
        for i=1:nOfruns,
            alg(i) = BHTSS;
            alg(i).dataFile = char(files(file));
            alg(i).parameters.numIt = 100;
            alg(i).parameters.nPobl = 200;
            alg(i).parameters.pCross = 0.8;
            alg(i).parameters.pMut = 0.2;
            alg(i).parameters.mutedPoints = 0.2;
            alg(i).parameters.seed = i*10;
            alg(i).parameters.minSeg = 20;
            alg(i).parameters.maxSeg = 120;
            alg(i).parameters.k = clusterArray(file);
            alg(i).parameters.iterClust = 20;
            alg(i).parameters.typeFitness = fit;
            alg(i).parameters.polyDegree = 1;
            alg(i).parameters.characActivation = [1 1 1 1 1 0];
            alg(i).parameters.umbralEntropy = 0.05;
            alg(i).parameters.intervalLeft = 0.05;
            alg(i).parameters.intervalRight = 0.95;
            information(i) = alg(i).runAlgorithm(serie(:,2));
            mkdir('.',[folder filesep num2str(i)])
            alg(i).saveAll(information(i),char(files(file)),[folder filesep num2str(i)]);
        end
        masterSaveAll(folder,information)
        clear information;
    end
end



