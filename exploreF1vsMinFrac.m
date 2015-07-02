function exploreF1vsMinFrac
    
    [gratingType, gratingType_s] = curGratingType;  % FLASHED_GRATING = 1; DRIFTING_GRATING = 2;
    [cmpType, cmpType_s] = curCmpType;    
    [pt_ids, pairTypes] = curPairTypes;
    pairTypes_str = [pairTypes{:}];
    
    ospDatafile = [CatV1Path gratingType_s 'GratingCells_DB.mat'];    
    cmpDatafile = [CatV1Path gratingType_s 'GratingComparisonData_DB_' cmpType_s '_' pairTypes_str '.mat'];     

    S1 = load(ospDatafile);        
    allCells = S1.allCells;
    allGids = [allCells.Gid];
    allCellIds = [allCells.cellId];
    nUnits = length(allCells);    
    
    S3 = load(cmpDatafile); 
    [pairData, S, pairTypes, measures, locations] = deal(S3.pairData, S3.allStatsC, S3.pairTypes, S3.measures, S3.locations);
    nCmp = length(pairData);
%     loc_names = {'maxMinFracR', 'maxR1xR2'};
    loc_names = {'maxMinFracR'};
%     ms_names = {'cc', 'rho'};
    ms_names = {'cc'};
    nLoc = length(loc_names);
    nM = length(ms_names);
    vals = cell(nLoc,nM); 
    hV_ax = zeros(nLoc, nM); hV_col = cell(nLoc, nM); 
    hH_ax = zeros(nLoc, nM); hH_im = zeros(nLoc, nM);
    allLocs_minFracs = cat(1, pairData.loc_minFracOfMaxes);        
    allF1oDCs_cmp  = cat(3, pairData.loc_F1oDCs_cmp);
    F1oDCs_pref = cat(1, pairData.F1oDCs_pref);    
    F1oDC1s_pref = F1oDCs_pref(:,1);
    F1oDC2s_pref = F1oDCs_pref(:,2);    
    %gather data
    pairGids = cat(1, pairData.Gids);
    pairCellIds = cat(1, pairData.cellIds);
    [MNdists1, MNdists2] = deal( zeros(nCmp, 1) );
    for pair_i = 1:nCmp 
        cell1_idx = find(allGids == pairGids(pair_i,1) & allCellIds == pairCellIds(pair_i,1) , 1);
        cell2_idx = find(allGids == pairGids(pair_i,2) & allCellIds == pairCellIds(pair_i,2) , 1);
        
        MNdists1(pair_i) = allCells(cell1_idx).spkFeatures.meanMNdist;
        MNdists2(pair_i) = allCells(cell2_idx).spkFeatures.meanMNdist;
    end
    allMN = [MNdists1(:); MNdists2(:)];
    f1odc_range = [0 2];
    mn_range = [0, round(max(allMN))*1.05];
    
    nbins0 = 10;
    binEdges0 = linspace(f1odc_range(1), f1odc_range(2), nbins0+1);
    binCent0 = binEdge2cent(binEdges0);
    
%     Gid1s_idx = arrayfun(@(x) find(allGids == pairGids(:    
%     MNdists1 = allCells(i).spkFeatures.meanMNdist;    

    [F1oDC1s_cmp, F1oDC2s_cmp, minFracs] = deal(zeros(nCmp, nLoc));

    for loc_i = 1:nLoc
        loc(loc_i) = find(strcmp(loc_names{loc_i}, locations)); %#ok<AGROW>

        F1oDC1s_cmp(:,loc_i) = squeeze(allF1oDCs_cmp(1,loc(loc_i),:));
        F1oDC2s_cmp(:,loc_i) = squeeze(allF1oDCs_cmp(2,loc(loc_i),:));
        
        minFracs(:,loc_i) = allLocs_minFracs(:,loc(loc_i));
        for m_i = 1:nM
            m(m_i) = find(strcmp(ms_names{m_i}, measures));           %#ok<AGROW>
            vals{loc_i, m_i} = S{loc(loc_i),m(m_i)}.val;
                        
        end
    end
        
%     allF1oDCs_pref = cat(3, pairData.F1oDCs_cmp);


    
    % plot
    idx = minFracs > 0;
    nCols = 30;
%     xlabel('F1/DC_1'); ylabel('F1/DC_2')
    hCbar_tmp = colorbar;
    pos_hc = get(hCbar_tmp, 'position'); set(hCbar_tmp, 'visible', 'off');

    figure(1); clf; 
    figure(2); clf;
    for loc_i = 1:nLoc
        for m_i = 1:nM
                      
            figure(1);
            hV_ax(loc_i, m_i) = subplot(nM, nLoc, loc_i+nLoc*(m_i-1));
            hV_col{loc_i,m_i} = colorPlot([], {F1oDC1s_cmp, F1oDC2s_cmp}, vals{loc_i, m_i}, {nCols, @jet, [-1, 1]}, '.');    
            axis equal tight;
            axis ([0 2 0 2]);
            
            jet1 = [1 1 1; jet(100)];
            figure(2);
            hH_ax(loc_i, m_i) = subplot(nM, nLoc, loc_i+nLoc*(m_i-1));            
%             p = get(hH_ax(loc_i, m_i), 'position');
            Z = wgtHist2D(F1oDC1s_cmp,F1oDC2s_cmp, vals{loc_i, m_i}, binEdges0, binEdges0);
            hH_im(loc_i,m_i) = imagesc(binCent0, binCent0, Z); colorbar; axis xy; 
            tk = [0, .5, 1, 1.5, 2];
            axis equal tight;
            axis ([0 2 0 2]);
            caxis([-1 1])
            set(gca, 'xtick', tk, 'ytick', tk);
            colormap(jet1);
%             set(hH_ax(loc_i, m_i), 'position', p);
%             xlabel('F1/DC_1'); ylabel('F1/DC_2');
            
        end
    end

    3;
    
%     p = get(gca, 'position');
%     dum_ax = axes('position', [p(1)+p(3), p(2)+p(4), .01, .01]);
%     dum_im = imagesc([-1; 1]);
%     hCbar = colorbar; 
%     colormap('jet');
%     set(hCbar, 'position', pos_hc);
    
%     figure(2); clf;
%     Z = wgtHist2D(F1oDC1s_cmp,F1oDC2s_cmp,vals, xBinEdges, yBinEdges);
%     hIm = imagesc(xBinCent, yBinCent, Z); colorbar; axis xy;
%     tk = [0, .5, 1, 1.5, 2];
%     set(gca, 'xtick', tk, 'ytick', tk);
%     xlabel('F1/DC_1'); ylabel('F1/DC_2');
    
    function plotData(xyaxes, minFrac_th, maxFrac_th, minVal_th, maxVal_th, nbins)        
                        
        for loc_j = 1:nLoc
            idx_frac = ibetween(  minFracs(:,loc_j), minFrac_th, maxFrac_th);            
            
            for m_j = 1:nM
                
                idx_vals = ibetween( vals{loc_j,m_j}, minVal_th, maxVal_th);
                idx = idx_frac & idx_vals;
                
                switch xyaxes
                    case 'Mahalanobis', Xs = MNdists1(idx,loc_j);
                                        Ys = MNdists2(idx,loc_j);
                                        xy_range = mn_range;
                    
                    case 'F1/DC @cmp', Xs = F1oDC1s_cmp(idx,loc_j);
                                       Ys = F1oDC2s_cmp(idx,loc_j);
                                       xy_range = f1odc_range;
                        
                    case 'F1/DC @pref', Xs = F1oDC1s_pref(idx);
                                        Ys = F1oDC2s_pref(idx);
                                        xy_range = f1odc_range;

                end
                binEdges = linspace( xy_range(1), xy_range(2), nbins+1);
                binCents = binEdge2cent(binEdges);
                
                hV_col{loc_j, m_j} = colorPlot(hV_col{loc_j, m_j}, {Xs, Ys}, ...
                                               vals{loc_j,m_j}(idx), {nCols, @jet, [-1, 1]}, '.');
                Z = wgtHist2D(Xs, Ys, vals{loc_j,m_j}(idx), binEdges, binEdges);
                set(hH_im(loc_j,m_j), 'xdata', binCents, 'ydata', binCents, 'cdata', Z');
                set([hV_ax(:); hH_ax(:)], 'xlim', xy_range, 'ylim', xy_range, 'xtickmode', 'auto', 'ytickmode', 'auto');
            end
        end
        
    end
    
    
    args = { {'axes', {'Mahalanobis', 'F1/DC @cmp', 'F1/DC @pref'}}, ...
             {'minFrac', [0:.01:1], 0}, {'maxFrac', [0:.01:1], 1}, ...
             {'minVal', [-1:.01:1], -1}, {'maxVal', [-1:.01:1], 1}, ...
             {'nbins', [5:40], 10} };
        
    manipulate(@plotData, args, 'FigId', 10);

3;


end


%------------------------------------
function Z = wgtHist2D(x,y,v, xBinEdges, yBinEdges)
    nx = length(xBinEdges)-1;
    ny = length(yBinEdges)-1;
    Z = zeros(nx, ny);
    cnt = zeros(nx, ny);
    [binnedX, xBinIds] = histcnt(x, xBinEdges);
    [binnedY, yBinIds] = histcnt(y, yBinEdges);        
    for i = 1:length(x)
        if (xBinIds(i) > 0) && (yBinIds(i) > 0)
            Z(xBinIds(i), yBinIds(i)) = Z(xBinIds(i), yBinIds(i)) + v(i);
    %         Z(xBinIds(i), yBinIds(i)) = Z(xBinIds(i), yBinIds(i)) + 1;
            cnt(xBinIds(i), yBinIds(i)) = cnt(xBinIds(i), yBinIds(i)) +1;
        end
    end  
    idx = (cnt > 0);
    Z(idx) = Z(idx)./cnt(idx);
    Z(~idx) = nan;
end

