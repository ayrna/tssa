%   COPYRIGHT
%   This file is part of TSSA: https://github.com/ayrna/tssa
%   Original authors: Antonio M. Duran Rosal, Pedro A. Gutierrez
%   Copyright:
%       This software is released under the The GNU General Public License v3.0 licence
%       available at http://www.gnu.org/licenses/gpl-3.0.html
%   Citation: If you use this code, please cite any of the following papers:
%     [1] A.M. Durán-Rosal, M. de la Paz Marín, P.A. Gutiérrez and C. Hervás-Martínez.
%         "Identifying market behaviours using European Stock Index time series by 
%         a hybrid segmentation algorithm", Neural Processing Letters,
%         Vol. 46, December, 2017, pp. 767–790.
%         https://doi.org/10.1007/s11063-017-9592-8
%
%% hybridIndividualNormal
% Function: Hybridization using entropy functions with Normal distribution
% 
% Input:
%     bestIndividual: chromosome to be hybridized
%     serie:          time series
%     umbralEntropy:  determine if the division is satisfactory
%     minSeg:         minimum segment size
%     
% Output:
%     individual:     hybridized segmentation
function [individual] = hybridIndividualNormal(bestIndividual,serie,umbralEntropy,minSeg)
    % Parche provisional puede funcionar (comprobar)
    bestIndividual(1,1)=1;
    bestIndividual(1,end)=1;
    cutsGenetics = find(bestIndividual==1);
    for i=1:numel(cutsGenetics)-1,
        segment = serie(cutsGenetics(i):cutsGenetics(i+1));
        bestIndividual2 = divideSegm(segment, numel(segment), bestIndividual(cutsGenetics(i):cutsGenetics(i+1)), umbralEntropy, minSeg);
        bestIndividual(cutsGenetics(i):cutsGenetics(i+1))=bestIndividual2;  
    end
    bestIndividual(1,1)=0;
    bestIndividual(1,end)=0;
    individual = bestIndividual;
end

%% divideSegm
% Function: Recursive division of a given segment
% 
% Input:
%     segment:        time series values of the segment
%     ind:            number of points of the segment
%     individual:     piece of the chrosomome corresponding to the segment
%     umbralEntropy:  determine if the division is satisfactory
%     minSeg:         minimum segment size
%     
% Output:
%     individual:     hybridized segment chromosome
function [individual] = divideSegm(segment, ind, individual,umbralEntropy,minSeg)
    umbral = -2*log(1-umbralEntropy);  %Determine if the division is satisfactory. alfa=0.05
    flag = 0;

    tamMinSeg = minSeg;

    %si el segmento puede ser dividido
    if(numel(segment) > (tamMinSeg*2)),
        %entropia del segmento
        e_segment = calculateEntropy(segment);
        %fprintf('e segmento-------------------------> %f\n',e_segment);
        %entropia de cada una de las particiciones posibles
        for j=tamMinSeg+1:numel(segment)-tamMinSeg,
           e_izq = calculateEntropy(segment(1:j));
           e_der = calculateEntropy(segment(j:end));
           if (~(isnan(e_izq)) && ~(isnan(e_der))),
           %se almacena la mejor division
               if(flag == 0),
                   flag = 1;
                   e_division = e_izq + e_der;
                   punto = j;
               elseif((e_izq + e_der) < e_division),
                       %fprintf('Actualizo\n');
                       e_division = e_izq + e_der;
                       punto = j;
               end
           end
        end


        %fprintf('e segmento-------------------------> %f\n',e_segment);
        %fprintf('e division-------------------------> %f\n\n',e_division);

        %Si la division es satisfactoria, se guarda
        if((e_segment - e_division) > umbral),
           individual(ind-numel(segment)+punto) = 1;
           %recursivo parte izq
           individual2 = divideSegm(segment(1:punto), (ind-numel(segment)+punto), individual, umbralEntropy, minSeg);
           %recursivo parte der
           individual = divideSegm(segment(punto:end), ind, individual2, umbralEntropy, minSeg);
        end
    end
end

%% calculateEntropy
% Function: Calculate the entropy of a given segment
% 
% Input:
%     segment:          time series values of the segment
% Output:
%     entropy:     entropy value of the segment
function [entropy] = calculateEntropy(segment)
        x = sum(segment);
        x2 = sum(segment.^2);
        entropy = numel(segment) * log(sqrt((x2/numel(segment)) - (x/numel(segment))^2));
end
