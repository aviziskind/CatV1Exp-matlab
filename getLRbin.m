function [l_bin, r_bin, stat_bin_x, stat_val_y, sample_m, sample_s, stat_y_wind_mean, comment, updatedStats] = ...
    getLRbin(Gid, cellId, allBins, statsM_S, intervalSize, statName, stimType, nStdTh, objTh, meanRate)  % psthVals

    minWindowWidth_bins = 3;    
    useFiringRateThreshold = false;
        minFiringRate_Hz = .18;

    useObjThreshold_indivBins = true;
        objTh_indivBins = objTh;%-log10(.01);
    useObjThreshold_entireWindow = true;    
        objTh_entireWindow = objTh; %-log10(.01);
    useNUniqueThreshold = strcmp(curResponseType(''), 'raw'); 
        nUniqueTh_minN = 7;
    useNoisinessThreshold = false;  %&& isfield(statsM_S, 'n_unique');    
        noisinessTh_maxN = .12;
    ifCantFindWindow_SearchAllValidWindows = true;
    
    minWindowAreaFracOfMax = .05;

%     useNzNzAsWeakThreshold = true;
    
    
%     statNameWindowCmp = statName;
%     statNameWindowCmp  = 'rho_p_nznz';
        
    psthOkLims_nonCph = [20 150];
    psthOkLims_Cph    = [20 200];  %4470:2 goes up to 180. 4482:5 goes up to 190
    sampleWindowLims  = [-300, 0]; % i fixed the cph problem, so can always go up to 0;
    l_bin = 0; r_bin = 0;
%     comment = '';
%     stat_y_wind_mean = 0;
    if isfield(statsM_S, 'clustIds')
        statsM_S = rmfield(statsM_S, 'clustIds');
    end
    updatedStats = statsM_S;
    
    binWidth = diff(allBins(1:2));
    psthWindowStart = allBins(1) - binWidth/2;

    binStatOffset = ceil((intervalSize-1)/2); % [0, 1, 1, 2, 2, ...]    
    psthOkLims = iff(stimType.isCphFlashed, psthOkLims_Cph, psthOkLims_nonCph); 
    
    stat_bin_x = allBins(1+binStatOffset:end-binStatOffset);
    nStatBins = length(stat_bin_x);    
    statSampleBin_idx =  stat_bin_x >= sampleWindowLims(1) &  stat_bin_x <= sampleWindowLims(2) ;                       

    useObjThreshold_indivBins_here = useObjThreshold_indivBins && ~strcmp(statName, 'r_entropy');
    useObjThreshold_entireWindow_here = useObjThreshold_entireWindow && ~strcmp(statName, 'r_entropy');        
    idx_rel = find( (stat_bin_x > psthOkLims(1)) & (stat_bin_x < psthOkLims(2)) );    
            
    nStats = 1;
    statsM = statsM_S.(statName);

    %%
    statSize = [120, 120];
    
    lbins_default = [77; 78; 77; 78];
    rbins_default = [100; 100; 101; 101]; 
    idx_default = sub2indV(statSize, [rbins_default'; lbins_default']);
    
    idx_mtx = reshape(1:numel(statsM), size(statsM));
    interval_baseline_idxs = diag( idx_mtx, -(intervalSize-1));
    [rbins_baseline, lbins_baseline] = ind2sub(statSize, interval_baseline_idxs);
    
    all_lbins = [lbins_default; lbins_baseline];
    all_rbins = [rbins_default; rbins_baseline];
    all_idxs = [idx_default; interval_baseline_idxs];
    
    minVals = structfun( @(s) full(min(s(all_idxs))), statsM_S);
    %%
    if any(minVals == 0)        
        all_vals = getPSTHwindowData(Gid, cellId, fieldnames(statsM_S), [], [], all_lbins, all_rbins);
        statsM(all_idxs) = all_vals.(statName)(all_idxs);
    end
    3;
    %%
    
    sigTestMethod = 'stdev'; % 'stdev' or 'prob'    
    
    if nStats == 1
        
        stat_val_y = full( diag(statsM, -(intervalSize-1)) );
        stat_val_y = stat_val_y(1:nStatBins); % for even intervalSize, remove the nan at end.
        stat_val_y_rel = stat_val_y(idx_rel);
        if isnan(stat_val_y(end))
%             keyboard;
            stat_val_y(end) = 0;
        end
        
        sample_1 = stat_val_y(statSampleBin_idx);
%             if ~strcmp(statName, 'r_entropy')
%                 sample_1 = sample_1(sample_1 < 2.5);
%             end
        if any(isnan(sample_1))
            sample_m = nanmean(sample_1);
            sample_s = nanstd(sample_1);
        else
            sample_m = mean(sample_1);
            sample_s = std(sample_1);
        end

        sgn = iff(strcmp(statName, 'r_entropy'), -1, 1);
        
%         if strcmp(sigTestMethod, 'stdev')
            
            idx_significant = find ( sgn*(stat_val_y_rel - sample_m)/sample_s > nStdTh );
            
%         elseif strcmp(sigTestMethod, 'prob')
            erfc1 = @(x) erfc(x/sqrt(2)); % integral of gaussian with variance 1 (instead of 1/2)        
            
            val_y_norm = sgn*(stat_val_y_rel - sample_m)/sample_s;
            stat_pval_y = erfc1(val_y_norm);            
            p_th = erfc1(nStdTh);            
            idx_significant2 = find( stat_pval_y < p_th) ;
%         end
            assert( isequal(idx_significant, idx_significant2) )
            
    end
                
    stat_y_wind_mean = nanmean(stat_val_y_rel);
    
    

    
    
%     meanRate = otherStats.meanRate;    
    if useFiringRateThreshold && (meanRate < minFiringRate_Hz)
        comment = sprintf('Below FR threshold %.2f < %.2f', meanRate, minFiringRate_Hz);
        return;
    end    
    
    if isempty(idx_significant)
        comment = sprintf('No significant bins.');
        return;
    end
    
    
%     findLocalMinima
    idx_bin_offset = idx_rel(1)-1;
    if ~strcmp(statName, 'r_entropy');
        minMinimaWidth = 1;
        stat_val_y_rel_smoothed = gaussSmooth(stat_val_y_rel, 1);
        idx_minima = findLocalMinima( stat_val_y_rel_smoothed, minMinimaWidth );
    else
        idx_minima = [];
    end
    windowIdxs = continuousGroupings(idx_significant + idx_bin_offset,  idx_minima + idx_bin_offset);
    
    areas = cellfun(@(i) sum(abs(stat_val_y(i)-sample_m)), windowIdxs);
    if minWindowAreaFracOfMax > 0
        maxArea = max(areas);
        idx_remove = areas < minWindowAreaFracOfMax*maxArea ; 
        windowIdxs ( idx_remove ) = [];
    end        
    nWindows = length(windowIdxs);

    % consider all possible concatenations of the groups (if more than 1)
    % and adding 1 or 2 bins before / after.
    addNbinsAround = 2;
    subTractNbinsAround = binStatOffset;
    binOptions = [-subTractNbinsAround : 1 : addNbinsAround];
    [xx_arr, yy_arr] = meshgrid(binOptions, binOptions);
    binsAround = [xx_arr(:), yy_arr(:)];
    n_around = size(binsAround,1);

    [xx_cat, yy_cat] = meshgrid([1:nWindows], [1:nWindows]);
    concats = [xx_cat(xx_cat <= yy_cat), yy_cat(xx_cat <= yy_cat)];
    n_concats = size(concats, 1);

    statPerms    = -ones(n_concats, n_around);
    nBinsForPerm = zeros(n_concats, n_around);
    lbins_toGet = zeros(n_concats, n_around);
    rbins_toGet = zeros(n_concats, n_around);
    for cat_i = 1:n_concats
        binL0 = windowIdxs{concats(cat_i,1)}(1);
        binR0 = windowIdxs{concats(cat_i,2)}(end);
        for arr_i = 1:n_around
            binL = binL0 - binsAround(arr_i,1);
            binR = binR0 + binsAround(arr_i,2);
            nBinsForPerm(cat_i, arr_i) = binR-binL+1;
            if (binL >= idx_rel(1)) && (binR <= idx_rel(end)) &&  (binR-binL+1 >= minWindowWidth_bins)
                val = statsM(binR+binStatOffset, binL+binStatOffset);
                if val == 0
                    lbins_toGet(cat_i, arr_i) = binL+binStatOffset; 
                    rbins_toGet(cat_i, arr_i) = binR+binStatOffset;
                else
                    statPerms(cat_i, arr_i) = val;
                end
            else
                3;
            end
        end
    end
        
%     lbins_toGet = [lbins_toGet(:); 77; 78; 77; 78];
%     rbins_toGet = [rbins_toGet(:); 100; 100; 101; 101];
    
    idx_toGet = find(lbins_toGet(:));
    if ~isempty(idx_toGet)
        updatedStats = getPSTHwindowData(Gid, cellId, fieldnames(statsM_S), [], [], lbins_toGet(idx_toGet), rbins_toGet(idx_toGet));
        updatedStatM = updatedStats.(statName);
        sub_x = rbins_toGet(idx_toGet);
        sub_y = lbins_toGet(idx_toGet);
        idx_vals = sub2indV(size(updatedStatM), [sub_x(:), sub_y(:)] );    
        statPerms(idx_toGet) = updatedStatM(idx_vals);
    end

    msk = zeros(size(statPerms));
    msk(statPerms == -1) = nan;    
    
    % add slight bias to statPerms (in case are maxed out at p =0) - in favor of larger window sizes.    
    smallestInc = min(diff(unique(statPerms(statPerms>0))))/2;
    if ~isempty(smallestInc)
        statPerms = statPerms +  smallestInc * (nBinsForPerm/sum(nBinsForPerm(:)));
    end
    
    [rep_p_max, inds] = maxElement((sgn * statPerms + msk) );  % for entropy - find minimum
    if ~any(statPerms > 0)
        comment = {sprintf('No windows were appropriate')};
        return;
%         error('fix algorithm')
    end
    [cat_j, arr_j] = dealV(inds);

%     l_bin_stat = windowIdxs{concats(cat_j,1)}(1)   - binsAround(arr_j,1);
%     r_bin_stat = windowIdxs{concats(cat_j,2)}(end) + binsAround(arr_j,2);        
    bestWindowConcat = concats(cat_j,:);
    bestWindow_binEdges = [windowIdxs{bestWindowConcat(1)}(1), windowIdxs{bestWindowConcat(2)}(end)];
    bestWindow_adjust = binsAround(arr_j,:);
    l_bin_stat = bestWindow_binEdges(1) - bestWindow_adjust(1);
    r_bin_stat = bestWindow_binEdges(2) + bestWindow_adjust(2); 
    l_bin = l_bin_stat + binStatOffset;
    r_bin = r_bin_stat + binStatOffset;

    l_bin_ms = psthWindowStart + (binWidth * [l_bin-1]);
    r_bin_ms = psthWindowStart + (binWidth * [r_bin  ]);
    
    stat_y_selectWind =  stat_val_y(l_bin_stat:r_bin_stat);
    stat_y_wind_mean = nanmean(stat_y_selectWind);
    
    
    if useNUniqueThreshold || useNoisinessThreshold  
        osp_selected = getOspDataForPsthWindow(Gid, cellId, [], [], l_bin, r_bin, [], 'osp');
        
        if useNUniqueThreshold 
            nUniqueHere = length(unique(osp_selected(:)));
        end
        if useNoisinessThreshold
            noisinessHere = ospNoisiness(osp_selected);
        end
        
    end
    
%     stat_y_wind_mean = max(stat_y_wind);
   
        
    
    indivBinsBelowTh =  useObjThreshold_indivBins_here && (nnz(stat_val_y_rel_smoothed >= objTh_indivBins) < 1);
    entireWindowBelowTh =  useObjThreshold_entireWindow_here && ~(rep_p_max >= objTh_entireWindow);
    windowTooSmall = (r_bin-l_bin+1 < minWindowWidth_bins);    
    belowNUniqueThreshold = useNUniqueThreshold && nUniqueHere < nUniqueTh_minN;    
    aboveNoisinessThreshold = useNoisinessThreshold && noisinessHere > noisinessTh_maxN;
    if belowNUniqueThreshold
        3;
    end
    if belowNUniqueThreshold && ~(indivBinsBelowTh || entireWindowBelowTh || windowTooSmall || aboveNoisinessThreshold)
        3;
    end
        
    
    if indivBinsBelowTh || entireWindowBelowTh || windowTooSmall || belowNUniqueThreshold || aboveNoisinessThreshold
        comment = {};
        if indivBinsBelowTh
            comment = {sprintf('No individual bins above min threshold of %.2f', objTh_indivBins)};
        end
        if entireWindowBelowTh
            comment = [comment, {sprintf('Entire window (bins %d-%d; %d-%d ms) was below threshold of %.2f', l_bin, r_bin, round(l_bin_ms), round(r_bin_ms), objTh_entireWindow) }];
        end
        if windowTooSmall
            comment = [comment, {sprintf('Window was too small (only %d bins)', r_bin-l_bin+1) }];
        end        
        if belowNUniqueThreshold
            comment = [comment, {sprintf('Too few unique entries (%d < %d)', nUniqueHere, nUniqueTh_minN) }];
        end                
        if aboveNoisinessThreshold
            comment = [comment, {sprintf('Too noisy (%.3f > %.3f)', noisinessHere, noisinessTh_maxN) }];
        end                
        [l_bin, r_bin] = deal(0);
        return;
    end

    comment = sprintf('Found window: [bins: %d-%d; %d-%d ms]', l_bin, r_bin, round(l_bin_ms), round(r_bin_ms));


end





                % allow for 2 almost continuous group separated by 1 or 2
                % dips below threshold:
                
%                 minDiff = @(i1,i2) min( abs(max(i1)-min(i2)), abs(max(i2)-min(i1)) );
%                 i_max = indmax(areas);   
%                 idx_LR = windowIdxs{i_max}([1, end]);
                
                %               
%                 if length(areas) > 1
%                     [areas_sorted, sort_idx] = sort(areas, 'descend');
%                     maxDistOf2Peaks = iff(isCphFG, 8, 5);
%                     minRatioOf2Areas = iff(isCphFG, .1, .25);
%                     area_ratios = areas_sorted(2) / areas_sorted(1);
%                     area_dist = minDiff(windowIdxs{sort_idx(1)}, windowIdxs{sort_idx(2)});
%                     
%                     if (area_ratios > minRatioOf2Areas) && (area_dist <= maxDistOf2Peaks)
%                         idxs = [windowIdxs{sort_idx(1)}; windowIdxs{sort_idx(2)}];
%                         idx_LR = [min(idxs), max(idxs)];
%                     end                    
%                 end                                     



%         case 'fracArea'           
% %             h_area = varargin{:};
%             areaBinsIdx = find(stat_bin_x >= psthOkLims(1) &  stat_bin_x <= psthOkLims(2)) ;
%             nBinsArea = length(areaBinsIdx);
% 
%             areaRatios = nan(nBinsArea, nBinsArea);
% 
%             areaStats = stat_val_y(areaBinsIdx);
%             
%             if ~strcmp(statName, 'r_entropy')
%                 areaStats = rectified(areaStats - sample_m);
%             else
%                 areaStats = rectified( - (areaStats - sample_m));
%             end
%             
%             for l_bin = 1 : nBinsArea
%                 for r_bin = l_bin : nBinsArea;
%                     areaRatios(r_bin, l_bin) = getPsthFracArea(areaStats, l_bin, r_bin);
%                 end
%             end
% 
%             [tmp_mx, maxInd] = maxElement(areaRatios);
%             [R_auto, L_auto] = dealV( maxInd );
%             l_bin = L_auto + binStatOffset + areaBinsIdx(1);
%             r_bin = R_auto + binStatOffset + areaBinsIdx(1);
%             
%             if isnan(tmp_mx),
%                 [l_bin, r_bin] = deal([]);
%             end
%             
% %             if ~isempty(h_area)
% %                 set(h_area, 'xdata', stat_bin_x(areaBinsIdx), 'ydata', stat_bin_x(areaBinsIdx), 'cdata', areaRatios);
% %                 xylims = [stat_bin_x(areaBinsIdx([1, end]))] + binWidth/2*[-1, 1];
% %                 set(h_area_ax, 'xlim', xylims, 'ylim', xylims);
% %             end
% %             L_bin_auto_ms = stat_bin_x(areaBinsIdx(L_auto));
% %             R_bin_auto_ms = stat_bin_x(areaBinsIdx(R_auto));                                        
% %             set(h_area_mx, 'xdata', L_bin_auto_ms, 'ydata', R_bin_auto_ms);
%             

%{
    if iscell(statName)
        nStats = length(statName);
        error(nchk(1, 2, nStats));        
        statsM = statsM_S.(statName{1});
        statsM2 = statsM_S.(statName{2});
    else
        nStats = 1;
        statsM = statsM_S.(statName);
    end


elseif nStats == 2
            % Let r', h' be zero-mean versions of r, h;  
            % compute covariance matrix  C=( <r'^2>,  <r' h'>; <r' h'> <h'^2> ); 
            % then in a Gaussian approximation, letting v be the vector (r;h) 
            % and m the mean of this vector, we have p(r,h)=( 1/(2\pi sqrt(Det(C)) ) e^{ -.5 (v-m)' C^{-1} (v - m) )
            
            stat_val_y1 = diag(statsM1, -(intervalSize-1));
            stat_val_y1 = stat_val_y1(1:nStatBins); % for even intervalSize, remove the nan at end.
            sample_1 = stat_val_y1(statSampleBin_idx);
            sgn1 = iff(strcmp(statName{1}, 'r_entropy'), -1, 1);
            m1 = nanmean(sample_1); s1 = nanstd(sample_1);            
            sample_1_norm = sgn1*(sample_1 - m1)/s1;
            val_y1_norm = sgn1*(stat_val_y1 - m1)/s1;
            
            stat_val_y2 = diag(statsM2, -(intervalSize-1));
            stat_val_y2 = stat_val_y2(1:nStatBins); % for even intervalSize, remove the nan at end.
            sample_2 = stat_val_y2(statSampleBin_idx);                        
            sgn2 = iff(strcmp(statName{2}, 'r_entropy'), -1, 1);
            m2 = nanmean(sample_2); s2 = nanstd(sample_2);            
            sample_2_norm = sgn1*(sample_1 - m1)/s1;
            val_y2_norm = sgn2*(stat_val_y2 - m2)/s2;

%             if ~strcmp(statName, 'r_entropy')
%                 statSampleVals = statSampleVals(statSampleVals < 2.5);
%             end
            
            C = cov(sample_1_norm(:), sample_1_norm(:));
            
                var1 = mean((sample_1_norm).^2);
                var2 = mean((sample_2_norm).^2);
                var12 = mean((sample_1_norm).*(sample_2_norm));
                C2 = [var1, var12; var12; var2];                
            
            assert(all(C(:) == C2(:)));
            V = [sample_1-m1; sample_2-m2];
            for i = 1:size(V,2);
                stat_p(i) = (1/(2*pi*sqrt(Det(C)))) * exp(-.5 * V(:,i)' * inv(C) * V(:,i) );
            end
            
            end
%}