function w = WCSS(X,y,C)
% Returns the within-cluster sum of squares (WCSS)

k = length(C);
w = zeros(size(1,k));

for i=1:k
    idx = (y==i);
    w(i) = sum( sum( ( X(idx,:)-repmat(C(i,:),size(X(idx),1),1) ).^2 )) ;
end

