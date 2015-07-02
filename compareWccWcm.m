function compareWccWcm
    
    gratingType = curGratingType;  % FLASHED_GRATING = 1; DRIFTING_GRATING = 2;

    ospDatafile = [CatV1Path gratingType 'GratingCells_DB.mat'];    
    cmpDatafile = [CatV1Path gratingType 'GratingComparisonData_DB.mat']; 

    S1 = load(ospDatafile);        
    allOSPs = S1.allOSPs;
    nUnits = length(allOSPs);    
    
    S3 = load(cmpDatafile); 
    [pairData, S, pairTypes, measures, locations] = deal(S3.pairData, S3.allStatsC, S3.pairTypes, S3.measures, S3.locations);
    nCmp = length(pairData);
    loc_name = 'maxR1xR2';
    ms_name = 'cc';
    
    loc = find(strcmp(loc_name, locations));
    m   = find(strcmp(ms_name, measures));
    
    vals = S{loc,m}.val;
    
    allGids = unique([allOSPs.GroupId]);
    nGroups = length(allGids);
    Wcc_corrs = nan(1,nGroups);
    Wcm_corrs = nan(1,nGroups);
    
    siteGids = cat(1, pairData.Gids);
    siteCellIds = cat(1, pairData.cellIds);
    Wcc_inds = cell(1,nGroups);
    Wcm_inds = cell(1,nGroups);
    
    for gi = 1:nGroups
        Wsite_inds = find(siteGids(:,1) == allGids(gi) & siteGids(:,2) == allGids(gi));
        Wcc_idx = (siteCellIds(Wsite_inds,1) > 0 & siteCellIds(Wsite_inds,2) > 0);
        Wcm_idx = (siteCellIds(Wsite_inds,1) == 0 | siteCellIds(Wsite_inds,2) == 0);
        
        Wcc_inds{gi} = Wsite_inds(Wcc_idx); 
        Wcm_inds{gi} = Wsite_inds(Wcm_idx);
    end
    
    ntotal_cc = 0; ntotal_cm = 0;
    rep_pval_th = 2;
    minF1oDC_cmp = 0;
    loc_minF1oDC = cat(1, pairData.loc_minF1oDCs_cmp);
    loc_minF1oDC = loc_minF1oDC(:,loc);
    min_respPvalData = cat(1, pairData.min_respPval);
    minN = 3;
    for gi = 1:nGroups
%             (min_respPvalData(idx) > rep_pval_th) & (loc_minF1oDC(idx) > minF1oDC_cmp);
        
        if ~isempty(Wcc_inds{gi})
            wcc_inds = (min_respPvalData(Wcc_inds{gi}) > rep_pval_th) & (loc_minF1oDC(Wcc_inds{gi}) > minF1oDC_cmp);
            Wcc_vals = vals(Wcc_inds{gi}(wcc_inds)); 
            if length(Wcc_vals) >= minN
                Wcc_corrs(gi) = nanmean( Wcc_vals );
                ntotal_cc = ntotal_cc + length(Wcc_vals);
            end
        end
        if ~isempty(Wcm_inds{gi})
            wcm_inds = (min_respPvalData(Wcm_inds{gi}) > rep_pval_th) & (loc_minF1oDC(Wcm_inds{gi}) > minF1oDC_cmp);
            Wcm_vals = vals(Wcm_inds{gi}(wcm_inds));
            if length(Wcm_vals) >= minN
                Wcm_corrs(gi) = nanmean( Wcm_vals );        
                ntotal_cm = ntotal_cm + length(Wcm_vals);
            end
        end
    end
    fprintf('ntotal cc = %d, cm = %d\n', ntotal_cc, ntotal_cm);
        
    figure(15); clf;
    plot(Wcm_corrs, Wcc_corrs, 'bo', 'markersize', 2); hold on;
    id = ~isnan(Wcm_corrs) & ~isnan(Wcc_corrs);
    p = polyfit(Wcm_corrs(id), Wcc_corrs(id), 2);
%     xs = [min(Wcm_corrs):.01:max(Wcm_corrs)];
    xs = [-1:.01:1];
    h2 = plot(xs, polyval(p, xs), 'r:');
    xlabel('Wcm'); ylabel('Wcc')
    title([titleCase(gratingType) ' gratings']);
    axis([-1 1, -1 1]);
    zeroaxes;
    
    

    3;
end