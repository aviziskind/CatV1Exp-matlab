

N = 15;
rand('state', 0);
X = rand(N,2); 

Y = pdist(X); 
Z = linkage(Y,'average'); 

figure(1);  clf;
[H, T, perm, col] = myDendrogram(Z, 'colorthreshold','default');

figure(2);  clf;
D = calcDendrogram(Z, 'colorthreshold','default');
[h, cols] = plotDendrogram([], D);

figure(3); clf;
for i = 1:N
    h_txt(i) = text(X(i,1), X(i,2), num2str(i), 'color', cols(i,:), 'horiz', 'center', 'vert', 'mid' ); %#ok<SAGROW>
end

% return;
for i = 1:20;
    
%     X(1,:) = X(1,:) + randn(1,2)/10; 
%     X = X + randn(size(X))/100;
%     X(X<0) = 0;
%     X(X>1) = 1;
    n_cur = size(X,1);
    if i <= 3
        idx = setdiff(1:n_cur, randi(n_cur));
        X = X(idx,:);
    elseif i > 3 && i <= 12        
        X(end+1,1:2) = rand(1,2);
    elseif i > 12    
        idx = setdiff(1:n_cur, randi(n_cur));
        X = X(idx,:);
    end

    Y = pdist(X); 
    Z = linkage(Y,'average');     
    D = calcDendrogram(Z, 'colorthreshold','default');
    [h, cols] = plotDendrogram(h, D);    
    
    figure(3); clf;
    N = size(X,1);

    for j= 1:N
        h_txt(j) = text(X(j,1), X(j,2), num2str(j), 'color', cols(j,:), 'horiz', 'center', 'vert', 'mid' ); %#ok<SAGROW>
    end
    
%     for j = 1:N
%         set(h_txt(j), 'position', X(j,:), 'color', cols(j,:) );
%     end
    3;
%     pause(.1);
end

