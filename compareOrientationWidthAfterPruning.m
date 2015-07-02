function compareOrientationWidthAfterPruning
% function analyzePruningInfo

refrRange_ms = [.8:.05:1.2];

allGids = getAllGids('o');

% gids = gids(1:100);

nGroups = length(allGids);

allData_c = cell(1,nGroups);
matchDB = 0;

    minNSpkRemoved = 1;

    wgtCCs_flag = 1;
    nRand = 1;
%     whichPlots = 'remaining';
    whichPlots = 'removed';

    curGroupingType('clusters')
    clustGrp_fn = getFileName('Groups', 'grating_dOr');
    clustPGrp_fn = strrep(clustGrp_fn, 'clusters', 'clustersPruned');
    
    S_grp_c = load(clustGrp_fn);  grps_c = S_grp_c.gratingGroups_dOr;
    S_grp_cp = load(clustPGrp_fn); grps_cp = S_grp_cp.gratingGroups_dOr;

%     clustsToUse{gi} = arrayfun(@(s1, s2) intersect(s1.cellIds, s2.cellIds), 0), S_grp_c, S_grp_cp, 'un', 0);    
    assert(length(grps_c) == nGroups);
    clustsToUse = cell(1,nGroups);
    grp_nSpikesRemoved_stim = cell(1,nGroups);
    grp_nSpikesRemoved_tot = cell(1,nGroups);
    grp_fracSpikesRemoved = cell(1,nGroups);
    nSpikeOrig = cell(1,nGroups);
%     clustOK = cell(1,nGroups);
    
    %%    
    
    for gi = 1:nGroups
        [clustIds_common, ia, ib] = intersect(grps_c(gi).cellIds, grps_cp(gi).cellIds);
        idx_clusts = grps_c(gi).cellIds(ia) > 0;
        ns_c = grps_c(gi).nSpikes(ia);
        ns_cp = grps_cp(gi).nSpikes(ib);
        nSpikesRemoved = ns_c - ns_cp;
        idx_pruned = nSpikesRemoved >= minNSpkRemoved;                
        idx_use = idx_clusts & idx_pruned;
        clustsToUse{gi} = clustIds_common(idx_use);                
        grp_nSpikesRemoved_tot{gi} = nSpikesRemoved(idx_use);
        assert(all(nSpikesRemoved(idx_use) > 0));
        grp_fracSpikesRemoved{gi} = nSpikesRemoved(idx_use) ./ ns_c(idx_use);
        
        nSpikeOrig{gi} = ns_c(idx_use);
        grp_fracSpikesRemoved{gi} = nSpikesRemoved(idx_use) ./ ns_c(idx_use);
    end
       
    %%
%         , S_grp_c, S_grp_cp, 'un', 0);    
                
    nClustEachGrp = cellfun(@length, clustsToUse);    
    nClustTot = sum(nClustEachGrp);
%     clustSpikes = 
%     nCellsTot = grps_c.cellIds
    
    frFactor = round(1000/(8+1/3)); % for drifting gratings (1 bin);
    
    clust_Gids_C = cell(1,nGroups);    
    clust_oriStats_C = cell(1,nGroups);
    clustP_oriStats_C = cell(1,nGroups);
    clustR_oriStats_C = cell(1,nGroups);
    clustP_diff_oriStats_C = cell(1,nGroups);
    clustR_diff_oriStats_C = cell(1,nGroups);
    clustFs_oriStats_C = cell(1,nGroups);

    clustFs_phaseCCs_C = cell(1,nGroups);
    clustR_phaseCCs_C = cell(1,nGroups);
    clustP_phaseCCs_C = cell(1,nGroups);
    
    otc_clust_C = cell(1,nGroups);
    otc_clustP_C = cell(1,nGroups);
    otc_clustR_C = cell(1,nGroups);
    
    
    clusts_ok_C = cell(1,nGroups);
%     jj = 0;
    progressBar('init-', nClustTot, 40);
    for gi = 1:nGroups
        Gid = allGids(gi);
        
        if nClustEachGrp(gi) == 0
            continue;
        end

%         clust_oriStats_C = cell(1,nGroups);
%         clustP_oriStats_C = cell(1,nGroups);
        
        [uClusterIds, clusterIdxs] = getCellSorting(Gid, 'clusters', matchDB);
        spikeStimIds = getSpikesStimIds(Gid);

        pruningFile = getFileName('clusterPrunings', Gid, matchDB, 1);
        S_prune = load(pruningFile);        
        clustSpikesRemoved_idx = S_prune.clustSpikesRemoved_idx;
        assert(isequal(S_prune.clustIds, uClusterIds));
        
        [nOri, nSpf, nPh, nTrials, ~] = getGratingStimType(Gid);
        R_full_dims = [nOri, nSpf, nPh, nTrials];
        
        clust_Gids_C{gi} = Gid(ones(1, length(clustsToUse{gi})));        
        
        for ci = 1:length(clustsToUse{gi})
            progressBar;
            clustId = clustsToUse{gi}(ci);
            c_idx = find(clustId == uClusterIds);
%             curGroupingType('clusters')
%             [R1_c, R_full1_c] = getOriSpfPhaseProfile_simple(Gid, uClusterIds(ci), [], 1);

%             curGroupingType('clustersPruned')
%             [R1_cp, R_full1_cp] = getOriSpfPhaseProfile_simple(Gid, uClusterIds(ci), []);            

%             figure(3); imageOSP(R1_c); figure(4); imageOSP(R2_c);
%             figure(13); imageOSP(R1_cp); figure(14); imageOSP(R2_cp);
                        
            idxs_clust = clusterIdxs{c_idx};
            cellSpikeStims_all = spikeStimIds(idxs_clust);            
            stimulus_idxs = find(cellSpikeStims_all);
            cellSpikeStims = cellSpikeStims_all(stimulus_idxs);
            R_full_clust = countIdxs(R_full_dims, cellSpikeStims) * frFactor;                        
            
            idxs_clustP = idxs_clust;
            idxs_clustP(clustSpikesRemoved_idx{c_idx}) = [];
            cellSpikeStimsP = nonzeros( spikeStimIds(idxs_clustP) );            
            R_full_clustP = countIdxs(R_full_dims, cellSpikeStimsP) * frFactor;
            R_full_clustP_diff = R_full_clust - R_full_clustP;                        
            assert(length(idxs_clust)-length(idxs_clustP) == grp_nSpikesRemoved_tot{gi}(ci));
                        
            
            %%
%             nSpksInClust_tot = length(idxs_clust);
            nSpksInClust_stim = nnz( cellSpikeStims );
%             nSpksRemoved = length(clustSpikesRemoved_idx{c_idx});
            nSpksAfterPruning = length(cellSpikeStimsP);
            nSpksRemoved = nSpksInClust_stim - nSpksAfterPruning;
            
            grp_nSpikesRemoved_stim{gi}(ci) = nSpksRemoved;
            
            
            clust_i_ok = nSpksRemoved > 1;
            clusts_ok_C{gi}(ci) = clust_i_ok;
            if ~clust_i_ok
                continue;
            end
            
            %%
            % correlation coefficients:
            clustP_ptc = reshape(mean(R_full_clustP,4), [nOri, nPh])';
            clustP_diff_ptc = reshape(mean(R_full_clustP_diff,4), [nOri, nPh])';          
            cc_wgt = getWeightedAverageCC(clustP_ptc, clustP_diff_ptc, wgtCCs_flag);
            clustP_phaseCCs_C{gi}{ci} = cc_wgt;

            
            otc_clustR_C{gi}{ci} = zeros(nOri, nRand);
            for ri = 1:nRand
                idxs_clustR = idxs_clust;
                idx_perm = randperm(nSpksInClust_stim); 
                idx_rm_rand = stimulus_idxs(idx_perm(1:nSpksRemoved));
                idxs_clustR(idx_rm_rand) = [];
                cellSpikeStimsR = nonzeros( spikeStimIds(idxs_clustR) );            
                R_full_clustR = countIdxs(R_full_dims, cellSpikeStimsR) * frFactor;
                
                otc_clustR_C{gi}{ci}(:,ri) = mean(mean(R_full_clustR, 3), 4);

                R_full_clustR_diff = R_full_clust - R_full_clustR;
                cc_wgt_rnd = getWeightedAverageCC(R_full_clustR, R_full_clustR_diff, wgtCCs_flag);                
                clustR_phaseCCs_C{gi}{ci}(ri) = cc_wgt_rnd;                
            end
            R_full_clustR_diff = R_full_clust - R_full_clustR;

            
            
            for ri = 1:nRand                                
                idx_choose_rand = randperm(length(cellSpikeStimsP)); 
                idx_choose_rand = idx_choose_rand(1:nSpksRemoved);                
                cellSpikeStimsFs = cellSpikeStimsP(idx_choose_rand);
                R_full_clustFs = countIdxs(R_full_dims, cellSpikeStimsFs) * frFactor;
                stats_clustFs = calcDegreeOfTuningStats(R_full_clustFs, [], Gid);                
                clustFs_oriStats_C{gi}{ci}(ri) = stats_clustFs.oriStats_si;
                
                clustFs_ptc = reshape(mean(R_full_clustFs,4), [nOri, nPh])';
                cc_wgt_rnd = getWeightedAverageCC(clustP_ptc, clustFs_ptc, wgtCCs_flag);                
                clustFs_phaseCCs_C{gi}{ci}(ri) = cc_wgt_rnd;                
                
            end
            
            
             
            assert( sum(R_full_clustP(:)) == sum(R_full_clustR(:)) );
                               
            stats_clust = calcDegreeOfTuningStats(R_full_clust, [], Gid);
            clust_oriStats_C{gi}(ci) = stats_clust.oriStats_si;
            clust_i_ok = stats_clust.oriStats_si.cellOK;
            clusts_ok_C{gi}(ci) = clust_i_ok;
                     
            if ~clust_i_ok
                continue;
            end

            otc_clust_C{gi}{ci} = mean(mean(R_full_clust, 3), 4);
            otc_clustP_C{gi}{ci} = mean(mean(R_full_clustP, 3), 4);
%             if nRand
%             otc_clustR_C{gi}{ci} = mean(mean(R_full_clustR, 3), 4);
            
            stats_clustP = calcDegreeOfTuningStats(R_full_clustP, [], Gid);
            clustP_oriStats_C{gi}(ci) = stats_clustP.oriStats_si;

            stats_clustR = calcDegreeOfTuningStats(R_full_clustR, [], Gid);
            clustR_oriStats_C{gi}(ci) = stats_clustR.oriStats_si;                

            stats_clustP_diff = calcDegreeOfTuningStats(R_full_clustP_diff, [], Gid);
            clustP_diff_oriStats_C{gi}(ci) = stats_clustP_diff.oriStats_si;

            stats_clustR_diff = calcDegreeOfTuningStats(R_full_clustR_diff, [], Gid);                                
            clustR_diff_oriStats_C{gi}(ci) = stats_clustR_diff.oriStats_si;                                

            stats_clustP = calcDegreeOfTuningStats(R_full_clustP, [], Gid);
            clustP_oriStats_C{gi}(ci) = stats_clustP.oriStats_si;


            
            3;
        end                
                      
        
%         [uClusterPrunedIds, clusterPrunedIdxs] = getCellSorting(Gid, 'clusters', matchDB);
        
%         calculateOriStatsFromClustIdxs(
        
%         [bins, allHistVals, meanRate, bckgSamples] = dbGetCellSpkStimHists(Gid, cellId, opt)
        
%         
%         fn = getFileName('prunedClustersStats', Gid);
%         S = load(fn);
%         allData_c{i} = S.pruningAutoData;

    end
    
    
    clust_Gids = selectOK(clust_Gids_C, clusts_ok_C);
    clust_clustIds = selectOK(clustsToUse, clusts_ok_C);
    
    otc_clust = selectOK(otc_clust_C, clusts_ok_C);
    otc_clust_P = selectOK(otc_clustP_C, clusts_ok_C);
    otc_clust_R = selectOK(otc_clustR_C, clusts_ok_C);
    %%
    fracSpikesRemoved = selectOK(grp_fracSpikesRemoved, clusts_ok_C);
    idx = ord(fracSpikesRemoved, 'descend');    
    
    clust_oriStats = selectOK(clust_oriStats_C, clusts_ok_C);    
    
    for j = 1:30; %:length(otc_clust)
       i = idx(j);
       figure(100+j); clf;
       n = length(otc_clust{i}); dOri = 360/n;
       oris = linspace(0, 360, n+1); oris = oris(1:n);
       
       dir_pref_i = clust_oriStats(i).dir_pref_deg;
       
       subplot(4,1,1:3);  
       
       [oris_shifted, otc_clust_shft, idx_shift] = shiftOriTuningCurve(oris, otc_clust{i}, dir_pref_i, -90);
%        shift
%        dir_pref_i       
       plot(oris_shifted, otc_clust{i}(idx_shift), 'ks-', 'markersize', 3); hold on;
       plot(oris_shifted, otc_clust_P{i}(idx_shift), 'bo-', 'markersize', 3);
       otc_clust_R_m = mean(otc_clust_R{i}, 2);
       otc_clust_R_s = std(otc_clust_R{i}, [], 2);
       errorbar(oris_shifted, otc_clust_R_m(idx_shift), otc_clust_R_s(idx_shift), 'r.-');       
       title(sprintf('Gid = %d. clustId = %d. fracRemoved = %.2f', ...
           clust_Gids(i), clust_clustIds(i), fracSpikesRemoved(i)));
       axis tight;
       set(gca, 'xtick', -90:45:270); xlim([-90, 270]);
       
       subplot(4,1,4);              
       clustP_diff = otc_clust{i}-otc_clust_P{i};
       clustR_diffs = bsxfun(@minus, otc_clust{i}, otc_clust_R{i});
       plot(oris_shifted, clustP_diff(idx_shift), 'bo-', 'markersize', 3); hold on;       
       otc_clust_diff_R_m = mean(clustR_diffs, 2);
       otc_clust_diff_R_s = std(clustR_diffs, [], 2);       
       
       errorbar(oris_shifted, otc_clust_diff_R_m(idx_shift), otc_clust_diff_R_s(idx_shift), 'r.-');
       axis tight;
       set(gca, 'xtick', -90:45:270); xlim([-90, 270]);
%        plot(oris, otc_clust_R{i}, 'r.-');

        3;
    end
    %%
    3;
    
    clusts_ok = [clusts_ok_C{:}];
    
    if strcmp(whichPlots, 'removed')

        %%
        
        clust_oriStats = selectOK(clust_oriStats_C, clusts_ok_C);
        clustP_oriStats = selectOK(clustP_oriStats_C, clusts_ok_C);
        clustR_oriStats = selectOK(clustR_oriStats_C, clusts_ok_C);        
        clustP_diff_oriStats = selectOK(clustP_diff_oriStats_C, clusts_ok_C);
        clustR_diff_oriStats = selectOK(clustR_diff_oriStats_C, clusts_ok_C);
        clustFs_oriStats = selectOK(clustFs_oriStats_C, clusts_ok_C);    

        
        clustFs_phaseCCs = selectOK(clustFs_phaseCCs_C, clusts_ok_C);    
        clustP_phaseCCs = selectOK(clustP_phaseCCs_C, clusts_ok_C);    
        clustR_phaseCCs = selectOK(clustR_phaseCCs_C, clusts_ok_C);    
        
        CFs_ori_width = cellfun(@(s) mean([s.w_ori_global]), clustFs_oriStats);
        
        
%         ori_fld = 'dir_pref_deg';
        ori_fld = 'ori_pref_deg';
                
        C_ori_pref = [clust_oriStats.(ori_fld)];
        CP_ori_pref = [clustP_oriStats.(ori_fld)];
        CPd_ori_pref = [clustP_diff_oriStats.(ori_fld)];
        CRd_ori_pref = [clustR_diff_oriStats.(ori_fld)];
        
        C_ori_width = [clust_oriStats.w_ori_global];
        CP_ori_width = [clustP_oriStats.w_ori_global];
        CR_ori_width = [clustR_oriStats.w_ori_global];
        CPd_ori_width = [clustP_diff_oriStats.w_ori_global];        
        CRd_ori_width = [clustR_diff_oriStats.w_ori_global];        
        
        
        %%
        figure(61); clf;
        ori_diffs_p = circDist(C_ori_pref, CPd_ori_pref, 180);
        ori_diffs_r = circDist(C_ori_pref, CRd_ori_pref, 180);
        h = hist2({ori_diffs_p, ori_diffs_r}, 30, 'line');
        set(h(1), 'color', 'b')
        set(h(2), 'color', 'r')
        legend('Pruned', 'Random');
        mn_p = mean(ori_diffs_p);
        mn_r = mean(ori_diffs_r);
        md_p = median(ori_diffs_p);
        md_r = median(ori_diffs_r);
        drawVerticalLine(md_p, 'color', 'b', 'linestyle', ':')
        drawVerticalLine(md_r, 'color', 'r', 'linestyle', ':')
        
        title({'\bfDifference in preferred orientation\rm', ...
                sprintf('Mean pruned : %.2f. Mean control diff: %.2f', mn_p, mn_r), ...
                sprintf('Median pruned: %.2f. Median control: %.2f', md_p, md_r)})

        %% widths: pruned vs random
        figure(62); clf;
        h = hist2({CPd_ori_width, CRd_ori_width}, 20, 'line');       
        set(h(1), 'color', 'b')
        set(h(2), 'color', 'r')
        legend('Pruned', 'Random');
        md_p = median(CPd_ori_width);
        md_r = median(CRd_ori_width);
            drawVerticalLine(md_p, 'color', 'b', 'linestyle', ':')
            drawVerticalLine(md_r, 'color', 'r', 'linestyle', ':')
        
%         title( sprintf('Median pruned width: %.2f. Median control width: %.2f', md_p, md_r) );
        title({'\bfOrientation tuning widths of removed spikes\rm', ...                
                sprintf('Median pruned: %.2f. Median control: %.2f', md_p, md_r)})


        %%  difference in widths: final-pruned    vs     final - (sampling n from final)
        figure(63); clf;
        diff_actual = abs(CP_ori_width-CPd_ori_width);
        diff_control = abs(CP_ori_width-CFs_ori_width);
        h = hist2({diff_actual, diff_control}, 30, 'line');      
        set(h(1), 'color', 'b')
        set(h(2), 'color', 'r')
        legend('Final Spikes Width - Pruned Spikes Width', 'Final Spikes Width - Random Samples Width');
        md_p = median(diff_actual);
        md_r = median(diff_control);        
        drawVerticalLine(md_p, 'color', 'b', 'linestyle', ':')
        drawVerticalLine(md_r, 'color', 'r', 'linestyle', ':')

        title({'\bfDifference in Orientation tuning Widths\rm', ...                
                sprintf('Actual difference in width: %.2f. Control difference in width: %.2f', md_p, md_r) })

            3;
        
        %%  difference in widths: final-pruned    vs     random final vs random-pruned
        figure(64); clf;
        diff_actual = abs(CP_ori_width-CPd_ori_width);
        diff_control = abs(CR_ori_width-CRd_ori_width);
        h = hist2({diff_actual, diff_control}, 30, 'line');      
        set(h(1), 'color', 'b')
        set(h(2), 'color', 'r')
        legend('Final Spikes Width - Pruned Spikes Width', 'Random Final Spikes Width - Random Pruned Spikes Width');
        md_p = median(diff_actual);
        md_r = median(diff_control);
        drawVerticalLine(md_p, 'color', 'b', 'linestyle', ':')
        drawVerticalLine(md_r, 'color', 'r', 'linestyle', ':')
        title({'\bfDifference in Orientation tuning Widths\rm', ...                
                sprintf('Actual difference in width: %.2f. Control difference in width: %.2f', md_p, md_r) })
        
        3;
        
        
        %% correlation coefficients between phase tuning curve: final,pruned vs final-sampled,pruned
        clustFs_phaseCCs_m = cellfun(@anyNonNanIn, clustFs_phaseCCs);
        clustR_phaseCCs_m = cellfun(@anyNonNanIn, clustR_phaseCCs);                        
        clustP_phaseCCs1 = [clustP_phaseCCs{:}];

        figure(65); clf;
        h = hist2({clustP_phaseCCs1, clustFs_phaseCCs_m}, 35, 'stairs');      
        set(h(1), 'color', 'b')
        set(h(2), 'color', 'r')
        legend('ptcCC(Final Spikes, Pruned Spikes)', 'ptcCC(Final Spikes, sampled from Final Spikes)');
        md_p = median(clustP_phaseCCs1);
        md_r = median(clustFs_phaseCCs_m);
        drawVerticalLine(md_p, 'color', 'b', 'linestyle', ':')
        drawVerticalLine(md_r, 'color', 'r', 'linestyle', ':')
        
        title({'\bfPhase Tuning Curve (weighted) mean correlation coefficients\rm', ...                
               sprintf('Median of Pruned: %.2f. Median of Control : %.2f', md_p, md_r)})
        
        
                
        %% correlation coefficients between phase tuning curve: final,pruned vs random final vs random pruned
        figure(66); clf;                
        
        h = hist2({clustP_phaseCCs1, clustR_phaseCCs_m}, 40, 'stairs');      
        set(h(1), 'color', 'b')
        set(h(2), 'color', 'r')
        legend('ptcCC(Final Spikes, Pruned Spikes)', 'ptcCC(Random Final Spikes, Random Pruned Spikes)');
        md_p = median(clustP_phaseCCs1);
        md_r = median(clustR_phaseCCs_m);
        drawVerticalLine(md_p, 'color', 'b', 'linestyle', ':')
        drawVerticalLine(md_r, 'color', 'r', 'linestyle', ':')
        
        title({'\bfPhase Tuning Curve (weighted) mean correlation coefficients\rm', ...                
               sprintf('Median of Pruned: %.2f. Median of Control : %.2f', md_p, md_r)})
                
        %%
        
        
        
        %%
        
        figure(64);
        
        3;
        
        
        
        
%         clust_OSI  = [clust_oriStats.OSI];
%         clustP_OSI = [clustP_oriStats.OSI];        
%         clustR_OSI = [clustR_oriStats.OSI];        
    %%

    %     ori_c = clust_OSI;
    %     ori_cp = clustP_OSI;
    %     ori_cr = clustR_OSI;

        ori_c = clust_w_glob;
        ori_cp = clustP_w_glob;
        ori_cr = clustR_w_glob;

        
       
        
    elseif strcmp(whichPlots, 'remaining')    

        clust_oriStats = selectOK(clust_oriStats_C, clusts_ok_C);
        clustP_oriStats = selectOK(clustP_oriStats_C, clusts_ok_C);
        clustR_oriStats = selectOK(clustR_oriStats_C, clusts_ok_C);
                
        clust_OSI  = [clust_oriStats.OSI];
        clustP_OSI = [clustP_oriStats.OSI];        
        clustR_OSI = [clustR_oriStats.OSI];        
    %%
        clust_w_glob  = [clust_oriStats.w_ori_global];
        clustP_w_glob = [clustP_oriStats.w_ori_global];        
        clustR_w_glob = [clustR_oriStats.w_ori_global];        

    %     ori_c = clust_OSI;
    %     ori_cp = clustP_OSI;
    %     ori_cr = clustR_OSI;

        ori_c = clust_w_glob;
        ori_cp = clustP_w_glob;
        ori_cr = clustR_w_glob;

        fracSpikesRemoved = [grp_fracSpikesRemoved{:}];

        %% plot dOSI vs fracRemoved
        dOSI_real = ori_c-ori_cp;
        dOSI_control = ori_c-ori_cr;

        figure(51); clf;
        h_control = plot(fracSpikesRemoved, dOSI_control, 'ro', 'markersize', 2); hold on;
        [cc_control, cc_p_control] = corr(fracSpikesRemoved(:), dOSI_control(:));
        p_fit = polyfit(fracSpikesRemoved(:), dOSI_control(:), 1);
        fplot(@(x) polyval(p_fit, x), xlim, 'r:');
        xlabel('fraction of cluster removed'); ylabel('Change in Global Ori width');
        str_control = sprintf('\\DeltaW_{Global}^{ORI} vs f_{removed} (control) : cc = %.2f, p = %.3g', cc_control, cc_p_control);

    %     figure(52); clf;    
        h_data = plot(fracSpikesRemoved, dOSI_real, 'o', 'markersize', 2); hold on;
        [cc_real, cc_p_real] = corr(fracSpikesRemoved(:), dOSI_real(:));
        p_fit = polyfit(fracSpikesRemoved(:), dOSI_real(:), 1);
        fplot(@(x) polyval(p_fit, x), xlim, 'b:');
        xlabel('fraction of cluster removed'); ylabel('Change in Global Ori width');
        str_real = sprintf('\\DeltaW_{Global}^{ORI} vs f_{removed} (data) : cc = %.2f, p = %.3g', cc_real, cc_p_real);
        title({str_real, str_control})
        legend([h_data, h_control], 'Pruning', 'Control')


        3;
        %%

        figure(52); clf;
        binE = linspace(-3, 3, 51); binC = binEdge2cent(binE);
        binV = histcnt(ori_cp-ori_c, binE);    
        bar(binC, binV, 1, 'facecolor', 'b');
        xlabel('Ori-Width_{after pruning}-Ori-Width_{before pruning}');
        m = mean(ori_cp-ori_c);
        drawVerticalLine(0, 'color', 'k');
        drawVerticalLine(m, 'color', 'r');


        %%
        figure(50); clf;
        plot(ori_c, ori_cr, 'ro', 'markersize', 2); hold on;
        plot(ori_c, ori_cp, 'bo', 'markersize', 2); hold on;
        mx = roundToNearest( max([ori_c(:); ori_cr(:)]), 5, 'up');
        axis equal;
        axis([0 mx 0 mx]);
        fplot(@(x) x, [0 mx], 'r');
        xlabel('W_{Global}^{ORI} before pruning'); ylabel('W_{Global}^{ORI} after pruning');
        legend('Random', 'Pruning', 'location', 'NW')



        3;

        %%


        figure(52); clf;
        mean_OSI_c = mean(ori_c);     std_OSI_c = std(ori_c);
        mean_OSI_cp = mean(ori_cp);   std_OSI_cp = std(ori_cp);

    %     bar([1, 2], [mean_OSI_c, mean_OSI_cp]); hold on;
        errorbar([1, 2], [mean_OSI_c, mean_OSI_cp], [std_OSI_c, std_OSI_cp], 'ro')
        %%
    3;    
    
        
        
        
    end
    
%     
%     
%     3;
% 
%     clustIds = [data.clustId];
%     data = data(clustIds > 0);
%     
%     clusterFracRemoved = cat(1, data.clusterFracRemoved);
%     
%     fracRemoved = 1 - clusterFracRemain;    
%     
%     idxSomeRefrac = fracRmv(:,end) > 0;
%     
%     fracRmv = fracRmv(idxSomeRefrac,:);
%     
%     Nclust_u = length(fracRmv);
%     
%     rs = zeros(Nclust_u,9);
%     for i = 1:9
%         rs(i,:) = fracRmv(:,9)./fracRmv(:,i);
%     end        
%     
%     rel_dec = (fracRmv(:,9)-fracRmv(:,1))./fracRmv(:,9);
% 
%     3;
%     
end 




function data_ok = selectOK(data_in, idx_ok_in)
data_ok = cellfun(@(data, idx) data(idx), data_in, idx_ok_in, 'un', 0);
data_ok = [data_ok{:}];
end


function wgt_cc = getWeightedAverageCC(R1, R2, wgt_flag)
    
    if size(R1, 4) > 1
        [nOri, nSp, nPh, nTrials] = size(R1); %#ok<NASGU>
        R1 = reshape(mean(R1,4), [nOri, nPh])';
    end
    if size(R2, 4) > 1
        [nOri, nSp, nPh, nTrials] = size(R2); %#ok<NASGU>
        R2 = reshape(mean(R2,4), [nOri, nPh])';
    end
    
    wgt_byFiringRate = exist('wgt_flag', 'var') && ~isempty(wgt_flag) && (wgt_flag == 1);
    
    sumR1 = sum(R1,1);
    sumR2 = sum(R2,1);
    idx_use = (sumR1 > 0) & (sumR2 > 0);
    if nnz(idx_use) == 0
        wgt_cc = nan;
        return;
    end
    
    ccs = pearsonR_v(R1(:,idx_use), R2(:,idx_use));
    
    if wgt_byFiringRate
        wgts = sumR1(idx_use) .* sumR2(idx_use);
        wgt_cc = wgt_nanmean(ccs, wgts);
    else
        wgt_cc = mean(ccs);
    end
    
end



% function X = countIdxs(X, idx)
%     for i = 1:numel(idx)
%         X(idx(i)) = X(idx(i))+1;
%     end
% end

%     
%     
%     clusterFracRemoved = cat(1, data.refrFracRemoved);
%     
%     allBestRefr_ms = [data.bestRefrPeriod_ms];
% allBestRefr_idx = arrayfun(@(x) find(refrRange_ms == x, 1), [data.bestRefrPeriod_ms]);
% allBestCuts = arrayfun(@(s,i) s.refrFracRemoved(i), data, allBestRefr_idx);
% 
% % end
% 3;
% figure(10); i = 9; hist(clusterFracRemoved(:,i), 40); title(sprintf('Refractory period = %.2f ms', refrRange_ms(i))); xlabel('fraction of cluster removed'); ylim([0 500])
% figure(10); i = 5; hist(clusterFracRemoved(:,i), 40); title(sprintf('Refractory period = %.2f ms', refrRange_ms(i))); xlabel('fraction of cluster removed'); ylim([0 500])
% figure(10); i = 1; hist(clusterFracRemoved(:,i), 40); title(sprintf('Refractory period = %.2f ms', refrRange_ms(i))); xlabel('fraction of cluster removed'); ylim([0 500])
% figure(10);  hist(allBestCuts, 40);  title('Best refractory period'); xlabel('fraction of cluster removed'); ylim([0 500])
% 


