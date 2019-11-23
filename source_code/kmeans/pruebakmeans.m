% Pequeño ejemplo de prueba de varias implementaciones de k-means con datos sintéticos. 
% No he sacado el SSE de las otras implementaciones TODO ;)

repeat = 30;
k=5;

%Let's make some fake data with two groups
RAW =load('dat/Synthetic-Gaussian-K-2-sig-0.250-modes3.dat');

x1=RAW(:,1);
x2=RAW(:,2);
X=RAW(:,1:2);

y=RAW(:,end);

n=size(RAW,1);

% scatter(x1,x2,[],y,'filled');
% 
% SSE = zeros(repeat, 1);
% times = zeros(repeat, 1);
% for i=1:repeat
%     tic;
%     C=initCentroids(X,k);
%     [L,C,S] = kmeans(X,k,'Start',C);
%     times(i,1)=toc;
%     SSE(i,1) = mean(sum(S));
% end
% 
% disp(['kmeans:   Mean SSE ' num2str(mean(SSE))])
% disp(['kmeans:   STD SSE ' num2str(std(SSE))])
% disp(['kmeans:   Mean Time ' num2str(mean(times))])


SSE = zeros(repeat, 1);
times = zeros(repeat, 1);
for i=1:repeat
    tic;
    [L,C] = kmeanspp(X,k);
    S=WCSS(X,L,C);
    times(i,1)=toc;
    SSE(i,1) = mean(sum(S));
end

disp(['kmeanspp: Mean SSE ' num2str(mean(SSE))])
disp(['kmeanspp: STD SSE ' num2str(std(SSE))])
disp(['kmeanspp: Mean Time ' num2str(mean(times))])

SSE = zeros(repeat, 1);
times = zeros(repeat, 1);
for i=1:repeat
    tic;
    C=initCentroids(X,k);
    [L,C] = dcKMeans(X,k,C);
    S=WCSS(X,L,C);
    times(i,1)=toc;
    SSE(i,1) = mean(sum(S));
end

disp(['dcKMeans: Mean SSE ' num2str(mean(SSE))])
disp(['dcKMeans: STD SSE ' num2str(std(SSE))])
disp(['dcKMeans: Mean Time ' num2str(mean(times))])
