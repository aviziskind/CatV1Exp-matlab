
N = 15;
% X = rand(N,2); 
for i = 1:10;
X(1,1) = X(1,1) + .05; 


Y = pdist(X); 
Z = linkage(Y,'average'); 

figure(1); 
[H, T, perm] = dendrogram(Z, 'colorthreshold','default');

figure(2); 
D = calcDendrogram(Z, 'colorthreshold','default');
[h, cols] = plotDendrogram(h, D);    

figure(3); clf; hold on;
for i = 1:N
%     plot(X(i,1), X(i,2), 'o', 'color', cols(i,:) );
    text(X(i,1), X(i,2), num2str(i), 'color', cols(i,:), 'horiz', 'center', 'vert', 'mid' );
end



%     pause(.05);
