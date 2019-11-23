%   COPYRIGHT
%   This file is part of TSSA: https://github.com/ayrna/tssa
%   Original authors: Antonio M. Duran Rosal, Pedro A. Gutierrez
%   Copyright:
%       This software is released under the The GNU General Public License v3.0 licence
%       available at http://www.gnu.org/licenses/gpl-3.0.html
%   Citation: If you use this code, please cite the following paper:
%     [1] A.M. Durán-Rosal, J.C. Fernández, P.A. Gutiérrez and C. Hervás-Martínez.
%         "Detection and prediction of segments containing extreme significant wave heights"
%         Ocean Engineering, Vol. 142, September, 2017, pp. 268-279.
%         https://doi.org/10.1016/j.oceaneng.2017.07.009
%
%% hybridIndividual
% Function: Hybridization using entropy functions
% 
% Input:
%     bestIndividual: chromosome to be hybridized
%     serie:          time series
%     intervalLeft:   left interval of the scaled
%     intervalRight:  right interval of the scaled
%     umbralEntropy:  determine if the division is satisfactory
%     minSeg:         minimum segment size
%     
% Output:
%     individual:     hybridized segmentation
function [individual] = hybridIndividualBeta(bestIndividual,serie,intervalLeft,intervalRight,umbralEntropy,minSeg)
    % Parche provisional puede funcionar (comprobar)
    bestIndividual(1,1)=1;
    bestIndividual(1,end)=1;
    cutsGenetics = find(bestIndividual==1);
    for i=1:numel(cutsGenetics)-1,
        segment = serie(cutsGenetics(i):cutsGenetics(i+1));
        bestIndividual2 = divideSegm(segment, numel(segment), bestIndividual(cutsGenetics(i):cutsGenetics(i+1)), intervalLeft, intervalRight, umbralEntropy, minSeg);
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
%     intervalLeft:   left interval of the scaled
%     intervalRight:  right interval of the scaled
%     umbralEntropy:  determine if the division is satisfactory
%     minSeg:         minimum segment size
%     
% Output:
%     individual:     hybridized segment chromosome
function [individual] = divideSegm(segment, ind, individual,intervalLeft,intervalRight,umbralEntropy,minSeg)
    umbral = -2*log(1-umbralEntropy);  %Determine if the division is satisfactory. alfa=0.05
    flag = 0;

    tamMinSeg = minSeg;

    %si el segmento puede ser dividido
    if(numel(segment) > (tamMinSeg*2)),
        %entropia del segmento
        e_segment = calculateEntropy(segment,intervalLeft,intervalRight);
        %fprintf('e segmento-------------------------> %f\n',e_segment);
        %entropia de cada una de las particiciones posibles
        for j=tamMinSeg+1:numel(segment)-tamMinSeg,
           e_izq = calculateEntropy(segment(1:j),intervalLeft,intervalRight);
           e_der = calculateEntropy(segment(j:end),intervalLeft,intervalRight);
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
           individual2 = divideSegm(segment(1:punto), (ind-numel(segment)+punto), individual, intervalLeft, intervalRight, umbralEntropy, minSeg);
           %recursivo parte der
           individual = divideSegm(segment(punto:end), ind, individual2, intervalLeft, intervalRight, umbralEntropy, minSeg);
        end
    end
end

%% calculateEntropy
% Function: Calculate the entropy of a given segment
% 
% Input:
%     segment:          time series values of the segment
%     intervalLeft:     left interval of the scaled
%     intervalRight:    right interval of the scaled
%     
% Output:
%     entropy:     entropy value of the segment
function [entropy] = calculateEntropy(segment,intervalLeft,intervalRight)
   %Transformamos el vector entre 0.05 y 0.05 
    minimo = min(segment);
    maximo = max(segment);
    segment=(segment - minimo)/(maximo-minimo);
    segment=(segment*(intervalRight-intervalLeft))+intervalLeft;

    %Se calcula entropia asumiendo distribucion gamma
    [alpha,beta]=calculateFinalsAlfaBeta(segment);
    resultado=-(alpha-1)*(gammaEuler(alpha)*(log(alpha)-(1/(2*alpha)))-gammaEuler(alpha+beta)*(log(alpha+beta)-(1/(2*(alpha+beta)))));
    resultado=resultado+ (  -(beta-1)* (gammaEuler(beta)*(log(beta)-(1/(2*beta)))-gammaEuler(alpha+beta)*(log(alpha+beta)-(1/(2*(alpha+beta)))))  );
    resultado=resultado+ ( log(gammaEuler(alpha)) + log(gammaEuler(beta)) - log(gammaEuler(alpha+beta)) );
    entropy=numel(segment)*resultado;
%             if(entropia<0),
%                 fprintf('ENTROPIA ES NEGATIVA\n');
%                 fprintf('Alpha %f\n', alpha);
%                 fprintf('Beta %f\n', beta);
%             end

    %Desnormalizamos a sus valores reales
    segment=(((segment-intervalLeft)/(intervalRight-intervalLeft))*(maximo-minimo))+minimo;

end

%% calculateFinalsAlfaBeta
% Function: Calculate alpha and beta
% 
% Input:
%     segment:    time series values of the segment
%     
% Output:
%     alpha:      value of estimated alpha
%     beta:       value of estimated beta
function [alpha,beta] = calculateFinalsAlfaBeta(segment)
    [alpha,beta]=calculateInitialsAlfaBeta(segment);
    s=10;
    for k=1:s,
        alpha_a=alpha;
        beta_a=beta;

        % Para alfa
        num1=sum(log(segment))/numel(segment);
        num2=log((s+alpha_a+beta_a-0.5)/(s+alpha_a-0.5));
        num3=0;
        den=0;
        for j=1:s,
           num3=num3+((beta_a*(j+alpha_a))/(j*(j+alpha_a-1)*(j+alpha_a+beta_a-1)));
           den=den+((beta_a)/(j*(j+alpha_a-1)*(j+alpha_a+beta_a-1))); 
        end
        alpha=(num1+num2+num3)/(den);

        % Para beta
        num1=sum(log(1-segment))/numel(segment);
        num2=log((s+alpha+beta_a-0.5)/(s+beta_a-0.5));
        num3=0;
        den=0;
        for j=1:s,
            num3=num3+((alpha*(j+beta_a))/(j*(j+beta_a-1)*(j+alpha+beta_a-1)));
            den=den+((alpha)/(j*(j+beta_a-1)*(j+alpha+beta_a-1)));
        end
        beta=(num1+num2+num3)/(den);
    end
%             if(isnan(alpha)),
%                 fprintf('ALPHA ES NAN\n');
%             end
%             if(isnan(beta)),
%                 fprintf('BETA ES NAN\n');
%             end
end

%% calculateInitialsAlfaBeta
% Function: Calculate initial alpha and beta
% 
% Input:
%     segment:    time series values of the segment
%     
% Output:
%     alpha_0:      value of estimated initial alpha
%     beta_0:       value of estimated initial beta
function [alpha_0,beta_0] = calculateInitialsAlfaBeta(segment)
    tau = (sum(segment)/numel(segment));
    gamma = ((sum(segment)/numel(segment))-(sum(segment.*segment)/numel(segment)))/((sum(segment.*segment)/numel(segment))-((sum(segment)/numel(segment))*(sum(segment)/numel(segment))));
    num1=sum(segment)/numel(segment);
    num2=sum(segment.*segment)/numel(segment);
    num=num1-num2;
    den1=sum(segment.*segment)/numel(segment);
    den2=(sum(segment)/numel(segment))*(sum(segment)/numel(segment));
    den=den1-den2;
    gamma2=num/den;
    alpha_0=tau*gamma;
    beta_0=gamma*(1-tau);
%             if(isnan(alpha_0)),
%                 fprintf('ALPHA 0 ES NAN\n');
%             end
%             if(isnan(beta_0)),
%                 fprintf('BETA 0 ES NAN\n');
%             end

end

%% gammaEuler
% Function: Calculate the Gamma de Euler of a value
% 
% Input:
%     value:  input value
%     
% Output:
%     result: gamma Euler value
function [result] = gammaEuler(value)
   value=value-1;
   result=sqrt(2*pi*value)*((value^value)*(exp(-value)));            
end
