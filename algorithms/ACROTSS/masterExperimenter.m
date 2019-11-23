%   COPYRIGHT
%   This file is part of TSSA: https://github.com/ayrna/tssa
%   Original authors: Antonio M. Duran Rosal, Pedro A. Gutierrez
%   Citation: If you use this code, please cite the associated paper [1]
%   Copyright:
%       This software is released under the The GNU General Public License v3.0 licence
%       available at http://www.gnu.org/licenses/gpl-3.0.html
%
%   References:
%     [1] A.M. Durán-Rosal, P.A. Gutiérrez, S. Salcedo-Sanz and C. Hervás-Martínez.
%         "A statistically-driven Coral Reef Optimization algorithm for optimal
%         size reduction of time series", Applied Soft Computing, 
%         Vol. 63. 2018, pp. 139-153.
%         https://doi.org/10.1016/j.asoc.2017.11.037
%
%   MASTER EXPERIMENTER
%   This script runs the ACROTSS algorithm given an experimental design
%   Results are saved in reports/ folder

% You should fill the files array with your time series
files ={'MALLAT_.txt','Donoho-Johnstone.txt'};
warning('off')
nOfruns = 30;

% Experimental setting
for file = 1:numel(files),
    for percentage=4:4,
        input = ['..' filesep '..' filesep 'time_series' filesep char(files(file))];
        serie = load(input);
        c = clock;
        folder = ['reports' filesep num2str(c(1)) '-' num2str(c(2)) '-'  num2str(c(3)) '-' num2str(c(4)) '-' num2str(c(5)) '-' num2str(uint8(c(6)))];
        mkdir('.',folder)
        for i=1:nOfruns,
            alg(i) = ACROTSS;
            alg(i).dataFile = char(files(file));
            alg(i).parameters.numIt = 100;
            alg(i).parameters.nPobl = 200;
            alg(i).parameters.numSeg = round(0.025*numel(serie(:,2)));
            alg(i).parameters.pCross = 1.1;
            alg(i).parameters.pMut = 1.1;
            alg(i).parameters.seed = i*10;
            alg(i).parameters.polyDegree = 0;
            alg(i).parameters.percentage_hybridation=0.20+percentage*0.05;
            alg(i).parameters.typeError = 2;
            % parametros http://onlinelibrary.wiley.com/doi/10.1002/ett.2759/full
            alg(i).parameters.freePositions = 0.2*alg(i).parameters.nPobl;
            alg(i).parameters.Fa = 0.05;
            alg(i).parameters.Fb = 0.98;
            alg(i).parameters.Fd = 0.05;
            alg(i).parameters.pDep = 0.01;
            alg(i).parameters.Natt = 2;
            information(i) = alg(i).runAlgorithm(serie(:,2));
            mkdir('.',[folder filesep num2str(i)])
            alg(i).saveAll(information(i),char(files(file)),[folder filesep num2str(i)]);
        end
        masterSaveAll(folder,information)
        clear information;
    end
end



