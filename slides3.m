function slides3



doSlides = [1];


if any(doSlides == 1) % various cell-cell distance measures
    
   S = load('driftingGratingCells_GLFcuw8_phase.mat');
   allGids = [S.allCells.Gid];
   allCellIds = [S.allCells.cellId];

   figure(1); % ori overlaps
   Gid1 = 422;
   Gid2 = 434;
   cellId1 = 1;
   cellId2 = 3;
   
   
   id1 = 2;
   id2 = 3;
   ids = [id1, id2];
   C1_idx1 = 53; %find(allGids == Gid & allCellIds == cellId1, 1);
   C1_idx2 = 59; %find(allGids == Gid & allCellIds == cellId2, 1);

   fet1 = S.allCells(C1_idx1).spkFeatures;
   fet2 = S.allCells(C1_idx2).spkFeatures;
   
   doOverlap = 1;
   if doOverlap
   
    M1 = double(fet1.negAmps_ccw_mean);  C1 = double(fet1.negAmps_ccw_cov);
    M2 = double(fet2.negAmps_ccw_mean);  C2 = double(fet2.negAmps_ccw_cov);                
    
    xa = min( M1(id1)-C1(id1,id1)*3, M2(id1)-C2(id1,id1)*3 );
    xb = max( M1(id1)+C1(id1,id1)*3, M2(id1)+C2(id1,id1)*3 );

    ya = min( M2(id2)-C2(id2,id2)*3, M2(id2)-C2(id2,id2)*3 );
    yb = max( M2(id2)+C2(id2,id2)*3, M2(id2)+C2(id2,id2)*3 );
    
    xs = linspace(xa, xb, 40);
    ys = linspace(ya, yb, 40);
    [xs_grid, ys_grid] = meshgrid(xs, ys);
    XY = [xs_grid(:), ys_grid(:)]';
    Zs1 = reshape( gaussianN(XY,  M1(ids), C1(ids, ids) ), size(xs_grid) );
    Zs2 = reshape( gaussianN(XY,  M2(ids), C2(ids, ids) ), size(xs_grid) );
%     surf(xs, ys, Zs1, 'facecolor', 'b'); hold on;
%     surf(xs, ys, Zs2, 'facecolor', 'g'); hold on;
    nCont = 8;
    figure(17); clf;
    [~, h1] = contour(xs, ys, Zs1, nCont, 'color', 'b'); hold on;
    [~, h2] = contour(xs, ys, Zs2, nCont, 'color', 'r');
    axis equal tight;
    3;
    
    figure(18); clf;
    h = surf(xs, ys, Zs1); hold on;
    surf(xs, ys, Zs2, 'color', 'r');
    3;
    
    
   end
   return;
    3;
    %}
    figure(88); clf;
    subplotGap(1,1,1,1, [.01 0 .01], [.01 0 .01])
    plot(fet1.wvfm_scl_mean, 'b', 'linewidth', 2); hold on;
    plot(fet2.wvfm_scl_mean, 'r', 'linewidth', 2); 
    xlim([1 172]);
    set(gca, 'xtick', [], 'ytick', [])
    set(gca, 'color', 'none', 'position', [.01 .01 .98, .98])
    
    
    figure(89); clf;
    n = 43; wvfm_idxs = arrayfun(@(i) i*n+1:(i+1)*n, 0:3, 'un', 0);
    for i = 1:4
        h(i) = subplotGap(1,4,1, i, [.01, .003, .01], [.01, .003, .01]);
        plot(fet1.wvfm_scl_mean(wvfm_idxs{i}), 'b', 'linewidth', 2); hold on;
        plot(fet2.wvfm_scl_mean(wvfm_idxs{i}), 'r', 'linewidth', 2); hold on;
        axis tight;
        set(gca, 'xtick', [], 'ytick', [])
        set(gca, 'color', 'none', 'position', get(gca, 'outerposition'))
    end
    matchAxes('Y', h);
    3;
    
    idx0 = 19+[0:3]*43;
        
    figure(90); clf;
    fwhm_idxs = [3, 170];
    cols = 'br';
    for i = 1:2
        h = subplotGap(1,2,1, i, [.01 .003, .01], [.01, .00, .01]);
        idx_ch1 = indmax(abs(fet1.wvfm_scl_mean(idx0)));
        f = S.allCells(fwhm_idxs(i)).spkFeatures;
        plot(f.wvfm_scl_mean(wvfm_idxs{idx_ch1}), cols(i), 'linewidth', 2);
        set(gca, 'xtick', [], 'ytick', []);
        set(gca, 'color', 'none');
        set(gca, 'position', get(gca, 'outerPosition'))
        axis tight;
%         subplotGap(1,2,1, i, [.003 .00, .003], [.002, .00, .002], h);
        3;
    end
        
    3;
        
    
    3;
%    Gid = 
    for i = 40:80;
        figure(100+i); plot(S.allCells(i).spkFeatures.wvfm_scl_mean);
        8
    end

%     32 - wide, peak at 1
%     25 - narrow, peak at 3
%     53 = narrow-ish, peak at 2
%     59 - wide-ish, peak at 1 & 4
%     
%     72 - narrow, peak at 3
%     71 - narrow peak at 2
    
    
    figure(2); 
   
    
    
end

