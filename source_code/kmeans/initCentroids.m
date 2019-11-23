function C=initCentroids(X,k, method)
%INITCENTROIDS Initialice centroids for k-means in a deterministic way. 
%   C = initCentroids(X,k) produces a 1-by-k matrix containing several centroids
if nargin < 2
    error(message('initCentroids:TooFewInputs'));
end

if nargin < 3
    method = 'extremevalues';
end

switch(method)
    case {'maxvariance'}
        % TODO: Finish this method
        
        % unnecesary
        % var=sqrt(std(X)

        % sort by variance/std
        [y,idx]=sort(std(X));

        % get the highest variance
        [y,idx] = sort(X(:,idx(1)),'descend');

        clear y;

        C = X(idx(1:k),:);
    case {'extremevalues'}
        % get the highest range (No funciona si normalizamos)
        %[y,idx]=max(max(X)-min(X));
        % get the highest standard deviation variable
        [y,idx]=max(std(X));
        [y,idx] = max(X(:,idx(1)));
        clear y;

        % First centroid
        C(1,:) = X(idx,:);
        
%         figure;hold on;
%         plot(X(:,1),X(:,2),'.','markersize',5)
%         plot(C(:,1),C(:,2),'r*','markersize',10)
        
        for i=2:k         
            % Original (Javi)
            %[y,idx] = max(sum(pdist2(X,C),2));
            
            distances = pdist2(X,C);            
            [y,idx] = max(sum(distances,2).*all(distances,2));
            
            C(i,:) = X(idx,:);
%             plot(C(:,1),C(:,2),'r*','markersize',10)
        end
       
%         hold off;

end