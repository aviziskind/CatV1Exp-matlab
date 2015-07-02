function testWgtKmeans


    rand('state', 0);
    K = 4;
    Cnt_pole = [-1, 2, 1, -2;
                2, 1, -2, -1];
    Cnt_cent = Cnt_pole;
            
%     dist = 'correlation';
    dist = 'sqEuclidean';
    n = 10;
%     Ntot = n*K;
    pts = cell(1,K);    
    s = 2;
    for i = 1:K
        pts{i} = [Cnt_pole(:,i), bsxfun(@plus, Cnt_cent(:,i), rand(2,n)*s)];
    end
    pts_v = cat(2, pts{:});

    sizeData = repmat([6^2, (2^2)*ones(1,n)], 1, K);
    colorData = ones(n+1,1)*[1:K];
    wgts = repmat([50, ones(1,n)], 1, K);
    
    figure(1); clf;
    subplot(1,3,1); 
    scatter(pts_v(1,:), pts_v(2,:), sizeData, colorData(:), 'o');

    [ids, Cnt] = kmeans(pts_v', K, 'start', Cnt_cent', 'distance', dist, 'emptyaction', 'singleton');    
    subplot(1,3,2); 
    scatter(pts_v(1,:), pts_v(2,:), sizeData, ids(:), 'o'); hold on;
    for i = 1:K
        plot(Cnt(i,1), Cnt(i,2), [color_s(i), 's'], 'markersize', 10)
    end
    
    [ids, Cnt_wgt] = wgtKmeans(pts_v', K, wgts, 'start', Cnt_cent', 'distance', dist, 'emptyaction', 'singleton');    
    subplot(1,3,3); 
    scatter(pts_v(1,:), pts_v(2,:), sizeData, ids(:), 'o'); hold on;
    for i = 1:K
        plot(Cnt_wgt(i,1), Cnt_wgt(i,2), [color_s(i), 's'], 'markersize', 10, 'markerfacecolor',color_s(i))
    end
    axis tight;
    
3;


end


