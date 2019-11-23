%   COPYRIGHT
%   This file is part of TSSA: https://github.com/ayrna/tssa
%   Original authors: Antonio M. Duran Rosal, Pedro A. Gutierrez
%   Citation: If you use this code, please cite the associated paper [1]
%   Copyright:
%       This software is released under the The GNU General Public License v3.0 licence
%       available at http://www.gnu.org/licenses/gpl-3.0.html
%
%   References:
%     [1] A.M. Durán-Rosal, P.A. Gutiérrez, F.J. Martínez-Estudillo and C. Hervás-Martínez.
%         "Simultaneous optimisation of clustering quality and approximation error
%         for time series segmentation", Information Sciences, Vol. 442-443, May, 2018, pp. 186-201.
%         https://doi.org/10.1016/j.ins.2018.02.041
%
%   MASTER EXPERIMENTER
%   This script runs the GMOTSS algorithm given an experimental design
%   Results are saved in reports/ folder

% You should fill the files array with your time series
files ={'MALLAT_.txt','Donoho-Johnstone.txt'};
warning('off')
nOfruns = 30;

for file=1:numel(files),    
    input = ['..' filesep '..' filesep 'time_series' filesep char(files(file))];
    serie = load(input);
    c = clock;
    folder = ['reports' filesep num2str(c(1)) '-' num2str(c(2)) '-'  num2str(c(3)) '-' num2str(c(4)) '-' num2str(c(5)) '-' num2str(uint8(c(6)))];
    mkdir('.',folder)
    for i=1:nOfruns,
        alg(i) = GMOTSS;
        alg(i).dataFile = char(files(file));
        alg(i).parameters.numIt = 200;
        alg(i).parameters.nPobl = 200;
        alg(i).parameters.pCross = 0.8;
        alg(i).parameters.pMut = 0.2;
        alg(i).parameters.mutedPoints = 0.2;
        alg(i).parameters.seed = i*10;
        alg(i).parameters.maxSeg = 150;
        alg(i).parameters.k = 5;
        alg(i).parameters.iterClust = 20;
        alg(i).parameters.typeFitness = 6;
        alg(i).parameters.polyDegree = 2;
        alg(i).parameters.minSeg = alg(i).parameters.polyDegree+2;
        alg(i).parameters.characActivation = [1 1 1 1 1 0];
        alg(i).parameters.typeError = 2;
        [information(i),informationClustering(i),informationError(i)] = alg(i).runAlgorithm(serie(:,2));
        mkdir('.',[folder filesep num2str(i)])
        mkdir('.',[folder filesep num2str(i) filesep 'bestGlobal'])
        mkdir('.',[folder filesep num2str(i) filesep 'bestClustering'])
        mkdir('.',[folder filesep num2str(i) filesep 'bestError'])
        alg(i).saveAll(information(i),char(files(file)),[folder filesep num2str(i) filesep 'bestGlobal']);
        alg(i).saveAll(informationClustering(i),char(files(file)),[folder filesep num2str(i) filesep 'bestClustering']);
        alg(i).saveAll(informationError(i),char(files(file)),[folder filesep num2str(i) filesep 'bestError']);
    end

    masterSaveAll(folder,'resultsBestGlobal',information);
    masterSaveAll(folder,'resultsBestClustering',informationClustering);
    masterSaveAll(folder,'resultsBestError',informationError);
    clear information;
    clear informationClustering;
    clear informationError;
end
