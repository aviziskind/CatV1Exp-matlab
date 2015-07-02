function slides2

doSlides = 3;


if any(doSlides == 1)   % show effect of not removing spikes.

    % get stuff ready
    fn = getFileName('properties', 4470);
    S = load(fn);

    amps = S.negAmps';
    amps_norm = bsxfun(@rdivide, amps, normV(amps, 1));
    
    [M,C] = getChannelMeansAndCovariance(4470);
    [V,D] = eig(C); 
    EvectMax = V(:,end);
    N = size(amps,2);    
    combos = nchoosek([1:4], 2);
    dots = zeros(6,N);
    for j = 1:6
        Ev = EvectMax([combos(j,:)]);
        V = amps_norm(combos(j,:), :);
        dots(j,:) = sum( bsxfun(@times, Ev, V), 1);
    end
    3;
    idx_spks_perp = find(abs(dots(1,:)) < .001);
    s1 = round(S.position(idx_spks_perp)); %#ok<FNDSB>
    
    t1 = s1(5) - 30;
    
    
    Gid = 4470;
    siteData = siteDataFor('Gid', Gid);
    dfInfo = siteData.dataFileInfo;
    nChannelsTot = dfInfo.nChannels;
    nSamplesTot = dfInfo.filesize/(2*nChannelsTot);            

    channelIds = dfInfo.channelIds;
    nChannelsToUse = length(channelIds);
    ticksPerMs = dfInfo.samplingRateHz/1000;
    [tmp, actualCov, whitenMtx] = getChannelMeansAndCovariance(Gid, 0);

    ii = 5;
    id1 = combos(ii, 1);
    id2 = combos(ii, 2);
    ids = [id1 id2];

    dfname = [dfInfo.dataFileName 'F'];
    pathname = rawDataDir(dfname);    
    dfname = [pathname dfname];

    removeAroundSpk_ms = [0 0];
    idx_aroundSpk = removeAroundSpk_ms(1)*ticksPerMs:removeAroundSpk_ms(2)*ticksPerMs;

%     t_starts = 5100+150;  % blank
%     chunkSize = 400;
    t_starts = 4398253;
    chunkSize = 500;
     
    
%     t_starts = t1+400;  % sideways spike
%     chunkSize = 500;

%     t_starts = t1+5900;  % downwards spike
%     chunkSize = 450;
% 
%     t_starts = t1+2900;  % 3 downwards spikes
%     chunkSize = 500;
% 
%     t_starts = t1+10000;
%     chunkSize = 10000;


    nStdTh = 3;
    sclMtx = diag(1./sqrt(diag(C)));
    chunk_i = double(readWindowsFromRawDataFile(dfname, t_starts, chunkSize, nChannelsTot, channelIds));     
%     chunk_i = bsxfun(@minus, chunk_i, mean(chunk_i, 1));
    chunkCov = cov(chunk_i); 
    [chunkEvect,chunkEval] = eig(chunkCov);        
    chunkWhitenMtx = sqrt(chunkEval)\chunkEvect';
    whitenMtx = chunkWhitenMtx;
    
    chunk_i_raw = double(chunk_i);
    chunk_i_scl = (sclMtx * chunk_i_raw')';
    chunk_i_ccw = (whitenMtx * chunk_i_raw')';
    
    
    chunk_i_toUse = chunk_i_ccw;
    figure(33); plot([1:chunkSize]/20, chunk_i_toUse);
    3;
    
    [~, chunk_i_itp] = fourierInterp(1:chunkSize, chunk_i_toUse, 30);
%     fs = 20000;
    plotRawData = 1;
    
    if plotRawData
        nChPlot = 2; spc = 6;        
        chunk_i_plot = bsxfun(@plus, chunk_i_toUse(:,1:nChPlot), -[0:nChPlot-1]*spc);
        figure(10+nChPlot); clf; plot([1:chunkSize]/20, chunk_i_plot);
        drawHorizontalLine(-[0:nChPlot-1]*spc)
        set(gca, 'ytick', [-nChPlot+1:0]*spc, 'yticklabel', arrayfun(@(i) sprintf('Ch %d', i), nChPlot:-1:1, 'un', 0), ...
            'units', 'pixels', 'position', [55, 50, 330, 80*(nChPlot+1)])
        ylim([-nChPlot*spc, spc])
        xlabel('ms')
    end

    mean_tmp = mean(chunk_i_toUse,1)';
    chunk_i_meansub = bsxfun(@minus, chunk_i_toUse, mean_tmp');
    
    C_samp = cov(chunk_i_meansub);
    [V,D] = eig(C_samp); %#ok<NASGU>
    EvectMax_samp = V(:,end);
    invCovMtx_samp = inv(C_samp);

    [V,D] = eig(actualCov); %#ok<NASGU>
    EvectMax_ac = V(:,end);
    invCovMtx_ac = inv(actualCov);
    
    nStdTh_ext = 5;
    [aboveTh_tmp, mndist] = aboveEllipsoidalThreshold(chunk_i_meansub', nStdTh_ext, EvectMax_ac, invCovMtx_ac); %#ok<ASGLU>
    aboveTh = mndist > nStdTh;   %aboveTh_tmp also takes into account the plane perp to the eigenvector, which we don't care about here.
    idx_inspks = find(aboveTh);

    if length(idx_aroundSpk) > 1
        idxToDelete = [];
        while ~isempty(idx_inspks)            
            idxToDelete = [idxToDelete, idx_inspks(1)+idx_aroundSpk]; %#ok<AGROW>
            remain = find( idx_inspks > idxToDelete(end), 1);
            idx_inspks = idx_inspks(remain:end);
        end
        idxToDelete(idxToDelete < 1 | idxToDelete > chunkSize) = []; 
    else
        idxToDelete = idx_inspks;
    end    
    
    chunk_i_all = chunk_i_toUse;
    chunk_i_aboveTh = chunk_i_toUse(idxToDelete,:);        
    chunk_i(idxToDelete,:) = [];        

    chunk_i_itp_meansub = bsxfun(@minus, chunk_i_itp, mean_tmp');
    [aboveTh_tmp, mndist_itp] = aboveEllipsoidalThreshold(chunk_i_itp_meansub', nStdTh_ext, EvectMax_samp, invCovMtx_samp);
    aboveTh_itp = mndist_itp > nStdTh;
    
    
%     id1 = 1;
%     id2 = 2;

%     figure(1); clf;
% %     plot(chunk_i_tmp(:,id1), chunk_i_tmp(:,id2), 'o'); hold on;
%     plot(chunk_i_belowTh(:,id1), chunk_i_belowTh(:,id2), 'g.', 'markersize', 1);
%     plot(chunk_i_aboveTh(:,id1), chunk_i_aboveTh(:,id2), 'r.', 'markersize', 1);
        markSpikeTracesRed = 0;
        drawEvects = 1;

        figure(2); clf; hold on;
        hax = gca;
        box on; axis equal;
%         set(hax, 'units', 'pixels', 'position', [50 50 300 300]); 
        if markSpikeTracesRed
            plot(chunk_i_aboveTh(:,id1), chunk_i_aboveTh(:,id2), 'ro', 'markersize', 3, 'markerfacecolor', 'r');
            plot(chunk_i_belowTh(:,id1), chunk_i_belowTh(:,id2), 'bo', 'markersize', 3, 'markerfacecolor', 'b'); 
            3;
    %         plot(chunk_i_all(:,id1), chunk_i_all(:,id2), 'k:'); 
            [idx_above, idx_below] = idxAboveBelow(aboveTh_itp);
            chunk_i_itp(end+1,:) = nan;

            plot(chunk_i_itp(idx_above,id1), chunk_i_itp(idx_above,id2), 'r-'); 
            plot(chunk_i_itp(idx_below,id1), chunk_i_itp(idx_below,id2), 'b-'); 
        else
            plot(chunk_i_all(:,id1), chunk_i_all(:,id2), 'bo', 'markersize', 3, 'markerfacecolor', 'b');
            plot(chunk_i_itp(:,id1), chunk_i_itp(:,id2), 'b-');             
        end
        3;
        
%         axis equal;

        ids = [id1, id2];
        
        cov_withSpikes = cov(chunk_i_toUse(:, ids));
        cov_withoutSpikes = cov(chunk_i_toUse(:, ids));
        
        actualCov = actualCov(ids, ids);
%         [cov_withSpikes, cov_withoutSpikes] = deal(actualCov);
        
        [x_ws, y_ws] = ellipsoidFromCov(mean_tmp(ids), cov_withSpikes, 3, 100);
        [x_ns, y_ns] = ellipsoidFromCov(mean_tmp(ids), cov_withoutSpikes, 5, 100);
        [x_ac, y_ac] = ellipsoidFromCov(mean_tmp(ids), actualCov, 5, 100);

        xlabel('Channel 1'); ylabel('Channel 2');
%         axis([-60, 60, -60, 60])
        3;
        if drawEvects
            [Evcts, D] = eig(cov_withoutSpikes);
            Evct1 = Evcts(:,1)*3*.95 * sqrt(D(1,1));
            Evct2 = Evcts(:,2)*3*.95 * sqrt(D(2,2));
            quiver(0,0,Evct1(1), Evct1(2), 'r', 'linewidth', 4, 'MaxHeadSize', 1)
            quiver(0,0,Evct2(1), Evct2(2), 'r', 'linewidth', 4, 'MaxHeadSize', 1)
            3;
            
        end
        
        plot(x_ws, y_ws, 'color', [.0 .85 0], 'linewidth', 2); hold on;
% %         plot(x_ns, y_ns, 'g-', 'linewidth', 2);
%         plot(x_ac, y_ac, 'k', 'linewidth', 3);
        
    
%         ax = [-580, 529, -580, 529];
%         axis(ax);
    3;
end

if any(doSlides == 2)  % circular mean & var colors
    
   
    N = 401;
    C = zeros(N,N,3);
    
    nOri = 360; 
    allCols = hsv(nOri);
    allR = zeros(N,N);
    allTh= zeros(N,N);

    ps = 0:.01:1;
    th = -pi:.01:pi;
    nP = length(ps);
    nTh = length(th);
    C1 = zeros(nP, nTh, 3);
    for i = 1:nP
        for j = 1:nTh
            col = meanVarToCol(th(j), ps(i), allCols);
            C1(i,j,:) = col;
        end
    end
    figure(4); clf;
    image(C1);
    3;
    
    
    mid = ceil(N/2);
    for i = 1:N
        for j = 1:N
            p = 2*([i,j]-mid)/(N);
            r = norm(p);
            th = atan2(p(1), p(2));
         
            col = meanVarToCol(th, r, allCols);
            
            if r > 1
                C(i,j,:) = [1 1 1];
            else
                C(i,j,:) = col;
            end
        end
    end
        
    
    mid = ceil(N/2);
    for i = 1:N
        for j = 1:N
            p = 2*([i,j]-mid)/(N);
            r = norm(p);
            th = atan2(p(1), p(2));
         
            col = meanVarToCol(th, 1-r, allCols);
            
            if r > 1
                C(i,j,:) = [1, 1, 1];                           
            else
                C(i,j,:) = col;
            end
            allR(i,j) = r;
            allTh(i,j) = th;

        end
    end
    figure(2); clf;
    image(C);
    set(gca, 'xtick', [], 'ytick', [], 'xcolor', [1 1 1], 'ycolor', [1 1 1]);
    box off;
    axis equal tight
    
%     figure(10); surf(C(:,:,1));
%     figure(11); surf(C(:,:,2));
%     figure(12); surf(C(:,:,3));
    3;
end


if any(doSlides == 3) % cluster pruning
    
%     Gid = 
    
    randn('state', 0);
    rand('state', 34);

    N1 = 1000;
    M1 = [ -270 -405]; 
    C1 = [4468 4360; 4360 5485]; 
    
    N2 = 200;
    M2 = [-230 -300];  %M2 = [-230 -245]; 
    C2 = [ 2237 1950; 1950 2273]; 
    
    X1 = mvnrnd(M1, C1, N1);
    X2 = mvnrnd(M2, C2, N2);
    figure(1); clf; hold on; box on;
    plot(X1(:,1), X1(:,2), 'bo', 'markersize', 3) 
    plot(X2(:,1), X2(:,2), 'go', 'markersize', 3);
    xlabel('Channel 1'); ylabel('Channel 2');
    set(gca, 'units', 'pixels', 'position', [60 60 300 300])
    
    title('Spike Amplitudes from Cell 1 and Cell 2');
    legend('Cell 1', 'Cell 2', 'location', 'NW');
    axis equal;
    
    [Evect,Eval] = eig(C1);  whitenMtx = sqrt(Eval)\Evect'; % channel whitening matrix                            
    whiten = @(X, M, W) [whitenMtx * bsxfun(@minus, X', M(:))]';
    X1w = whiten(X1, M1, whitenMtx);
    X2w = whiten(X2, M1, whitenMtx);
    figure(2); clf; hold on; box on;
    plot(X1w(:,1), X1w(:,2), 'bo', 'markersize', 3);
    plot(X2w(:,1), X2w(:,2), 'go', 'markersize', 3);
    xlabel('Eigenvector 1'); ylabel('Eigenvector 2');
    set(gca, 'units', 'pixels', 'position', [60 60 300 300])
    axis equal;
    legend('Cell 1', 'Cell 2', 'location', 'NW');
    title('Spike Amplitudes in eigenbasis of Cell 1');
    
    nRefr = 5;
    idx1 = randi(N1, 1, nRefr);
    idx2 = randi(N2, 1, nRefr);
    h_refr = zeros(1,nRefr);
    for i = 1:nRefr
        Xr1 = [X1w(idx1(i), 1), X1w(idx1(i), 2)];
        Xr2 = [X2w(idx2(i), 1), X2w(idx2(i), 2)];
        
        h_refr(i) = plot([Xr1(1), Xr2(1)], [Xr1(2), Xr2(2)], 'k-', 'linewidth', 2);        
        h_refr(i) = plot([Xr1(1)], [Xr1(2)], 'ko', 'markerfacecolor', 'b', 'linewidth', 2);        
        h_refr(i) = plot([Xr2(1)], [Xr2(2)], 'ko', 'markerfacecolor', 'g', 'linewidth', 2);        
    end
    3;
    
    
   
    
    
end


if any(doSlides == 33) % 3D hyperplane picture
   
    randn('state', 0);
    rand('state', 34);

    N1 = 1000;
    M1 = [0;0;0]; 
    C1 = [eye(3)]; 
    
    N2 = 200;
    M2 = [2; 1.2; 0.9];  %M2 = [-230 -245]; 
    C2 = [eye(3)]/3; 
    
    X1 = mvnrnd(M1, C1, N1);
    X2 = mvnrnd(M2, C2, N2);
    figure(1); clf; hold on; box on;
    plot3(X1(:,1), X1(:,2), X1(:,3), 'bo', 'markersize', 3); 
    plot3(X2(:,1), X2(:,2), X2(:,3), 'go', 'markersize', 3);
    
    view(3);
    3;
    Z = @(x,y) 6-3*x-1*y;
    xlims = xlim;
    ylims = ylim;
    h = fmesh(Z, xlims, ylims);
    set(h, 'edgecolor', [1 .2 .2], 'edgealpha', .9, 'facealpha', .2);
    
    args = {{'a', [-10:.1:10], 6}, {'b', [-10:.1:10], 3}, {'c', [-10:.1:10], 1}};
    manipulate(@updatePlot, args, 'FigId', 5);
    3;
end

function updatePlot(a,b,c)
    Z2 = @(x,y) a-b*x-c*y;    
    [xs, ys, zs] = fmesh(Z2, xlims, ylims);
    set(h, 'xdata', xs, 'ydata', ys, 'zdata', zs);
    set(h, 'edgecolor', [1 .2 .2], 'edgealpha', .9, 'facealpha', .2);
end



if any(doSlides == 4) % wgt kmeans
    
   % regular k-means   
   randn('state', 3);
   
   ns = [3, 4, 5, 4];
   Ms = {[-2; 0], [0; 4], [2; 2], [5; -2]};
   N = length(ns);   
   X = cellfun(@(M,n) bsxfun(@plus, randn(2,n)/1.7, M), Ms, num2cell(ns), 'un', 0);
   idxs = arrayfun(@(i,n) i*ones(1,n), 1:N, ns, 'un', 0);
   idxs = [idxs{:}];
   
   npts = sum(ns);
   
   figure(43); clf; hold on; box on; axis equal;
   
   cnt = 1; h = zeros(1,npts);
   for i = 1:N
       for j = 1:ns(i)
            h(cnt) = plot(X{i}(1,j), X{i}(2,j), 'ko', 'markersize', 4, 'markerfacecolor', 'auto');
            cnt = cnt +1;
       end
   end
   3;
   X_all = [X{:}]';
   
   allColors = lines(7);
   
   K = 3;
   [k_idx, C] = kmeans(X_all, K, 'replicates', 100);
   
   for i = 1:npts
        set(h(i), 'color', allColors(k_idx(i),:), 'markerfacecolor', allColors(k_idx(i),:));
   end
   
   showPtToCent = 0;
   showWithinLines = 0;
   showBetweenLines = 1;
   for i = 1:K
%        col_idx_i = mode( idxs( k_idx == i) );
        plot(C(i,1), C(i,2), 'o', 'color', 'k', 'markersize', 9);               
        if showPtToCent 
            for j = find(k_idx(:)' == i)
                plot([X_all(j,1) C(i,1)], [X_all(j,2) C(i,2)], 'k-')
            end            
        end
        
        [idx_sameCells, idx_diffCells] = getIdxSameDiffCells(k_idx, npts);
        if showWithinLines
            [i1, i2] = find(idx_sameCells);
            for jj = 1:length(i1)                
                plot([X_all(i1(jj),1) X_all(i2(jj),1)], [X_all(i1(jj),2) X_all(i2(jj),2)], 'k-')                
            end
        end
        if showBetweenLines
            [i1, i2] = find(idx_diffCells);
            for jj = 1:length(i1)
                plot([X_all(i1(jj),1) X_all(i2(jj),1)], [X_all(i1(jj),2) X_all(i2(jj),2)], 'k-')
            end
        end
        
   end
    3;    
    
end

if any(doSlides == 5)  % kmeans with only 4 points (a little clearer)
    
% regular k-means   
   randn('state', 3);
   
   X = [1   .5,  2  6;
        1, 1.5,  3, 1.5];

    figure(44); clf; hold on; box on;
    plot(X(1,:), X(2,:), 'o'); axis equal
    axis([0 6.5 0 4]); 
    
   Ms = {[-2; 0], [0; 4], [2; 2], [5; -2]};
   N = length(ns);   
   X = cellfun(@(M,n) bsxfun(@plus, randn(2,n)/1.7, M), Ms, num2cell(ns), 'un', 0);
   idxs = arrayfun(@(i,n) i*ones(1,n), 1:N, ns, 'un', 0);
   idxs = [idxs{:}];
   
   npts = sum(ns);
   
   figure(43); clf; hold on; box on; axis equal;
   
   cnt = 1; h = zeros(1,npts);
   for i = 1:N
       for j = 1:ns(i)
            h(cnt) = plot(X{i}(1,j), X{i}(2,j), 'ko', 'markersize', 4, 'markerfacecolor', 'auto');
            cnt = cnt +1;
       end
   end
   3;
   X_all = [X{:}]';
   
   allColors = lines(7);
   
   K = 3;
   [k_idx, C] = kmeans(X_all, K, 'replicates', 100);
   
   for i = 1:npts
        set(h(i), 'color', allColors(k_idx(i),:), 'markerfacecolor', allColors(k_idx(i),:));
   end
   
   showPtToCent = 0;
   showWithinLines = 0;
   showBetweenLines = 1;
   for i = 1:K
%        col_idx_i = mode( idxs( k_idx == i) );
        plot(C(i,1), C(i,2), 'o', 'color', 'k', 'markersize', 9);               
        if showPtToCent 
            for j = find(k_idx(:)' == i)
                plot([X_all(j,1) C(i,1)], [X_all(j,2) C(i,2)], 'k-')
            end            
        end
        
        [idx_sameCells, idx_diffCells] = getIdxSameDiffCells(k_idx, npts);
        if showWithinLines
            [i1, i2] = find(idx_sameCells);
            for jj = 1:length(i1)                
                plot([X_all(i1(jj),1) X_all(i2(jj),1)], [X_all(i1(jj),2) X_all(i2(jj),2)], 'k-')                
            end
        end
        if showBetweenLines
            [i1, i2] = find(idx_diffCells);
            for jj = 1:length(i1)
                plot([X_all(i1(jj),1) X_all(i2(jj),1)], [X_all(i1(jj),2) X_all(i2(jj),2)], 'k-')
            end
        end
        
   end
    3;   
    
    
    
end


if any(doSlides == 6)  % circular kmeans 
    
% regular k-means   
   randn('state', 3);

   th = 2*pi*[.1 .12  .2, .5];
   X = [cos(th); sin(th)];
   
    figure(44); clf; hold on; box on;
    plot(X(1,:), X(2,:), 'o'); axis equal
    axis([0 6.5 0 4]); 
    
   Ms = {[-2; 0], [0; 4], [2; 2], [5; -2]};
   N = length(ns);   
   X = cellfun(@(M,n) bsxfun(@plus, randn(2,n)/1.7, M), Ms, num2cell(ns), 'un', 0);
   idxs = arrayfun(@(i,n) i*ones(1,n), 1:N, ns, 'un', 0);
   idxs = [idxs{:}];
   
   npts = sum(ns);
   
   figure(43); clf; hold on; box on; axis equal;
   
   cnt = 1; h = zeros(1,npts);
   for i = 1:N
       for j = 1:ns(i)
            h(cnt) = plot(X{i}(1,j), X{i}(2,j), 'ko', 'markersize', 4, 'markerfacecolor', 'auto');
            cnt = cnt +1;
       end
   end
   3;
   X_all = [X{:}]';
   
   allColors = lines(7);
   
   K = 3;
   [k_idx, C] = kmeans(X_all, K, 'replicates', 100);
   
   for i = 1:npts
        set(h(i), 'color', allColors(k_idx(i),:), 'markerfacecolor', allColors(k_idx(i),:));
   end
   
   showPtToCent = 0;
   showWithinLines = 0;
   showBetweenLines = 1;
   for i = 1:K
%        col_idx_i = mode( idxs( k_idx == i) );
        plot(C(i,1), C(i,2), 'o', 'color', 'k', 'markersize', 9);               
        if showPtToCent 
            for j = find(k_idx(:)' == i)
                plot([X_all(j,1) C(i,1)], [X_all(j,2) C(i,2)], 'k-')
            end            
        end
        
        [idx_sameCells, idx_diffCells] = getIdxSameDiffCells(k_idx, npts);
        if showWithinLines
            [i1, i2] = find(idx_sameCells);
            for jj = 1:length(i1)                
                plot([X_all(i1(jj),1) X_all(i2(jj),1)], [X_all(i1(jj),2) X_all(i2(jj),2)], 'k-')                
            end
        end
        if showBetweenLines
            [i1, i2] = find(idx_diffCells);
            for jj = 1:length(i1)
                plot([X_all(i1(jj),1) X_all(i2(jj),1)], [X_all(i1(jj),2) X_all(i2(jj),2)], 'k-')
            end
        end
        
   end
    3;   
    
    
    
end


if any(doSlides == 7)  % mean waveforms (illustrating negative amplitudes
    Gid = 2342;
    S_wvfm = load(getFileName('mwvfm_cell', Gid));
    mwvfms = S_wvfm.meanWaveforms(2:end);
    t_ms = S_wvfm.t_ms;
    nClust = length(mwvfms);
    [nT, nChannels] = deal(43,4);
    figure(57); clf;
    for i = 1:nClust
        subplot(1,nClust,i);
        offset_ms = 100;
        wvfm_ch = reshape(mwvfms(i).wvfm_raw, [nT, nChannels]);
        plot(t_ms, wvfm_ch)
    end
    
    3;
    
    
    
end





end


function col = meanVarToCol(circMn, circVr, allCols)
    
%     if fd_flag == 1 % flashed gratings
%         nOri = 180;
%     else            % drifting gratings.
%         nOri = 360; 
%     end
    nOri = size(allCols, 1);
    
    col_id = round( rad2deg( circMn ) );
    col_id = mod(col_id, nOri)+1;

    col = allCols(col_id,:);
   
%     circVr = max((circVr-0.5),0);
    col = bound( col*(1- circVr), 0, 1);
    
end

function [idx_above, idx_below] = idxAboveBelow(tf)
    n_end = length(tf);    

    idx_above_raw = find(tf);
    grp_above = continuousGroupings(idx_above_raw)';           
    grp_above(2,:) = {n_end+1};
    idx_above = cat(1, grp_above{:});
    
    idx_below_raw = find(~tf);
    grp_below = continuousGroupings(idx_below_raw)';           
    grp_below(2,:) = {n_end+1};
    idx_below = cat(1, grp_below{:});
        
    
end


function [idx_sameCells, idx_diffCells] = getIdxSameDiffCells(cell_clustIds, nClust)

    if isnumeric(cell_clustIds)
        [~, cell_clustIds] = uniqueList(cell_clustIds);
    end            
%     nClustsPerCell = cellfun(@length, cell_clustIds);
%     nClusts = sum(nClustsPerCell);
    nCells = length(cell_clustIds);

    tf_sameCell = false(nClust, nClust);
    for cell_i = 1:nCells
        tf_sameCell(cell_clustIds{cell_i}, cell_clustIds{cell_i}) = true;
    end    
    idxBelowDiag_clust = tril(true(nClust), -1);
    idx_sameCells = ( tf_sameCell & idxBelowDiag_clust);
    idx_diffCells = (~tf_sameCell & idxBelowDiag_clust);
end

