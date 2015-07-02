function [bestRefrPeriod_ms, fracRemoved, isMU, allFracsRemoved] = getBestRefrPeriodFromDepMtx(...
    depMtx, isis_ms, refrSpikeMarkedIdxs, params, opt, clustIdsCut)

    refr_leeway_ms = opt.refr_leeway_ms;
    refrRange_ms = opt.refrRange_ms;
    
    if isfield(params, 'refrMaxCalc') && ~isempty(params.refrMaxCalc)       
        % for multiunits with lots of refractory period violations, cut off calculations after a
        % certain number. So don't calculate after this point, since dep-matrix doesn't cover it.
        refrRange_ms = refrRange_ms(refrRange_ms <= params.refrMaxCalc);
    end

%     if isfield(params, 'true_refrPeriod_ms') && ~isempty(params.true_refrPeriod_ms);
%         bestRefrPeriod_ms = true_refrPeriod_ms;
%     else
        bestRefrPeriod_ms = refrRange_ms(end);
%     end
    isMU = false;
    nRefrPeriods = length(refrRange_ms);    
    
    isCrossRefr = exist('clustIdsCut', 'var') && ~isempty(clustIdsCut);
        
    
    isMultiUnit = isscalar(depMtx) && isnan(depMtx);
    if isMultiUnit % too many refractory spikes.
        isMU = true;
        fracRemoved = nan;
        allFracsRemoved = nan;
        return;
    end            
    
    nRefrSpikes = length(isis_ms);
    if (nRefrSpikes == 0) % no refractory spikes.        
        fracRemoved = 0;
        allFracsRemoved = 0;
        return;
    end        
    
    allFracsRemoved = zeros(1,nRefrPeriods);
    idx_spksToCut = [];
    sameForRefrPeriodUpTo = 0;
    for ri = 1:nRefrPeriods
        refrRange_ms_i = refrRange_ms(ri);
        if (refrRange_ms_i > sameForRefrPeriodUpTo)
            [idx_spksToCut, sameForRefrPeriodUpTo] = getAutoCutFromDepMtx(depMtx, isis_ms, refrRange_ms_i ) ;
        end           
        
        if isCrossRefr             
            clustI_spkCut_idx = idx_spksToCut & (clustIdsCut == 1);
            clustJ_spkCut_idx = idx_spksToCut & (clustIdsCut == 2);                
            clustI_spksRemoved = unique( vertcat( refrSpikeMarkedIdxs{clustI_spkCut_idx})  );
            clustJ_spksRemoved = unique( vertcat( refrSpikeMarkedIdxs{clustJ_spkCut_idx})  );
            nSpksRemoved = length(clustI_spksRemoved) + length(clustJ_spksRemoved);

            allFracsRemoved(ri) = (nSpksRemoved / sum(params.nSpkInBothClusts));            
        else
            nSpksRemoved = nUniqueVals(cat(1, refrSpikeMarkedIdxs{idx_spksToCut}));            
            
            allFracsRemoved(ri) = (nSpksRemoved / params.nSpksInClust);
        end        
                
        
    end
    
    
        
    allFracsRemoved_null = getFNull_t(refrRange_ms, allFracsRemoved);
    
    
    refrDecreases = (allFracsRemoved_null - allFracsRemoved);

    idx_1ms = indmin( abs(refrRange_ms - 1));
    
%     idx_maxAllowedRefrPeriod = indmin(abs(refrRange_ms - opt.maxRefrPeriod_ms));
%     refrDecreases_allowed = refrDecreases(1:idx_maxAllowedRefrPeriod);

    idx_allowedRefrPeriod = find( ibetween(refrRange_ms, opt.minRefrPeriod_ms, opt.maxRefrPeriod_ms) );
    refrDecreases_allowed = refrDecreases(idx_allowedRefrPeriod);

    haveTrueRefrPeriod = isfield(params, 'true_refrPeriod_ms') && ~isempty(params.true_refrPeriod_ms);

    if haveTrueRefrPeriod
        if isCrossRefr
            3;
        end
        idx_best_tmp = find(refrRange_ms < params.true_refrPeriod_ms, 1, 'last');
        refr_leeway_ms = 0; 
        
    else
    
        if ~isCrossRefr % procedure for single cluster refractory periods

            if ~any(refrDecreases_allowed > 0) % f_removed never goes below f_null
                isMU = true;
                fracRemoved = nan;
                allFracsRemoved = nan;
                return;
            end                  
            
            refrDecMax = max(refrDecreases_allowed );
            idx_best_tmp = find(refrDecreases_allowed == refrDecMax, 1, 'last') + idx_allowedRefrPeriod(1)-1;
            % this often results in a refractory period that is slightly too restrictive - allow for a
            % little leeway (but not less than 1 ms).
            
                      
            
        else  % procedure for cross-cluster refractory periods:
            % instead of probing all possible cross-refractory periods, we will instead use the
            % minimum of the two individual refractory periods.
            % but once again, we can use a bit of leeway (0.15 ms)

            minAutoRefrPeriods = min(params.autoRefrPeriod_ms);
            idx_best_tmp = find( abs(refrRange_ms - minAutoRefrPeriods) < 1e-5, 1);

        end

    end
    
    
    nLeewayPts = round( refr_leeway_ms/diff(refrRange_ms(1:2)) );
    nLeewayPts = min(nLeewayPts, rectified(idx_best_tmp-idx_1ms) );

    allFracsRemoved_nearBest = allFracsRemoved(idx_best_tmp + [-nLeewayPts:0]);
    bestFracRemoved = min ( allFracsRemoved_nearBest );
    idx_best = find(bestFracRemoved == allFracsRemoved_nearBest, 1, 'last') + idx_best_tmp - (nLeewayPts+1);

    bestRefrPeriod_ms = refrRange_ms(idx_best);
    fracRemoved = allFracsRemoved(idx_best);
    
    
    assert(ibetween(bestRefrPeriod_ms, opt.minRefrPeriod_ms, opt.maxRefrPeriod_ms));
    
    haveGoodRefrac = isfield(opt, 'true_refrPeriod_ms') && ~isempty(opt.true_refrPeriod_ms) && ...
        (bestRefrPeriod_ms - opt.true_refrPeriod_ms) >= 1; % opt.true_refrPeriod_ms < refrRange_ms(end);
    
    dbug = 0;% && exist('autoRefrPeriod_ms', 'var');
%     dbug = haveGoodRefrac; %isfield(opt, 'clustId') && opt.clustId == 1;
    
    
    if dbug
        %%
        ori = 'vert';
        showLabels = 1;
        
        nSub = 2;
        
         mergeToOnePlot = 1;
        xlims = [.3, max(refrRange_ms)];
        mn = switchh(ori, {'vert', 'horz'}, {[nSub, 1], [1 nSub]});
        [m,n] = dealV(mn);
        
        figure(5); clf; 
        if ~mergeToOnePlot 
            h_ax(1) = subplot(m,n,1);
            hy1 = plot(refrRange_ms, allFracsRemoved, '.-'); 
            h_ax(2) = subplot(m,n,2);
            hy2 = plot(refrRange_ms, refrDecreases, 'b.:');                  
        else
            [h_ax, hy1, hy2] = plotyy(refrRange_ms, allFracsRemoved, refrRange_ms, refrDecreases); 
        end
        set(h_ax, 'xlim', xlims, 'nextplot', 'add'); 
        set(hy1, 'linewidth', 2);
        set(hy2, 'linewidth', 2, 'color', [0 .7 0], 'linestyle', ':');
        hy2_tmp = plot(h_ax(1), 0, 0, 'linewidth', 2, 'color', [0 .7 0], 'linestyle', ':');
        
        
        
        h_null = plot(h_ax(1), refrRange_ms, allFracsRemoved_null, 'k-'); 
%         drawVerticalLine(bestRefrPeriod_ms);
        ylims1 = [0 allFracsRemoved(end)*1.1];
        ylim(ylims1);
               
        ylims2 = get(h_ax(2), 'ylim');
        
%         drawVerticalLine(opt.maxRefrPeriod_ms, 'linestyle', ':', 'color', 'r', 'linewidth', 2)        
        hB = line(bestRefrPeriod_ms*[1 1], ylims1(2) - [0, diff(ylims1)/2], 'linestyle', '-', 'color', 'r', 'linewidth', 2, 'parent', h_ax(1));
        if ~mergeToOnePlot
            line(bestRefrPeriod_ms*[1 1], ylims2(2) - [0, diff(ylims2)/2], 'linestyle', '-', 'color', 'r', 'linewidth', 2, 'parent', h_ax(2));
        end
        if showLabels
            axis(h_ax(1));  xlabel('t'); 
            if ~mergeToOnePlot
                axis(h_ax(2));  xlabel('t'); ylabel('f_{null}(t) - f_{removed}(t)')
            end
        end
%         drawnow;
%         pause(1)
        if isfield(opt, 'true_refrPeriod_ms') && ~isempty(opt.true_refrPeriod_ms) %&& opt.true_refrPeriod_ms < refrRange_ms(end)
            hT = line(opt.true_refrPeriod_ms*[1 1], ylims1(1) + [0, diff(ylims1)/2], 'linestyle', '-', 'color', [0 .7 0], 'linewidth', 2, 'parent', h_ax(1));
            if ~mergeToOnePlot
                line(opt.true_refrPeriod_ms*[1 1], ylims2(1) + [0, diff(ylims2)/2], 'linestyle', '-', 'color', [0 .7 0], 'linewidth', 2, 'parent', h_ax(2))
            end
        end
        
        3;
        if showLabels
            if ~mergeToOnePlot
                legend(h_ax(1), {'f_{removed}(t)', 'f_{null}(t)', 'Chosen t_{refrac}', 't_{true}'}, 'location', 'best', 'fontsize', 8)                        
            else
                legend([hy1, h_null, hy2_tmp, hB, hT], {'f_{removed}(t)', 'f_{null}(t)', 'f_{removed}(t) - f_{null}(t)' 'Chosen t_{refrac}', 't_{true}'}, 'location', 'SE', 'fontsize', 8)                        
            end
        end
        title(h_ax(1), sprintf('Group %d, cluster %d', opt.Gid, opt.clustId));
       3; 
    end
    
end

%         if exist('autoRefrPeriod_ms', 'var') && ~isempty(autoRefrPeriod_ms)
%             x1 = autoRefrPeriod_ms(1);
%             x2 = autoRefrPeriod_ms(2);
%             line([x1 x1], [0, ylims(2)/2], 'linestyle', ':', 'color', 'r', 'linewidth', 2);
%             line([x2 x2], [ylims(2)/2, ylims(2)], 'linestyle', ':', 'color', 'g', 'linewidth', 2);
% 
%             title(sprintf('%d, %d', params.clustIds));
%         end

function fNull_t = getFNull_t(refrRange_ms, allFracsRemoved)
    gap_ms = 0.2;
    m_expected = allFracsRemoved(end) / (refrRange_ms(end) - gap_ms);
    c_expected = allFracsRemoved(end) - m_expected * refrRange_ms(end);
    fNull_t = m_expected * refrRange_ms + c_expected;                
end

function n = nUniqueVals(x)
    if isempty(x)
        n = 0; 
        return;
    end
    mx_val = max(x);            
    val_present = false(mx_val, 1); val_present(x) = 1; n = nnz(val_present);

    return;
    nTot = length(x);
    frac = nTot/mx_val;
    if (frac < .01)    
        val_present = false(mx_val, 1); val_present(x) = 1; n = nnz(val_present);            
        n2 = length(unique(x));
    elseif (frac < .02)
        val_present = false(mx_val, 1); val_present(x) = 1; n = nnz(val_present);            
        n2 = length(unique(x));        
    elseif (frac < .03)
        val_present = false(mx_val, 1); val_present(x) = 1; n = nnz(val_present);            
        n2 = length(unique(x));        
    elseif (frac < .04)
        val_present = false(mx_val, 1); val_present(x) = 1; n = nnz(val_present);            
        n2 = length(unique(x));        
    elseif (frac < .05)
        val_present = false(mx_val, 1); val_present(x) = 1; n = nnz(val_present);            
        n2 = length(unique(x));        
    elseif (frac < .06)
        val_present = false(mx_val, 1); val_present(x) = 1; n = nnz(val_present);            
        n2 = length(unique(x));        
    elseif (frac < .07)
        val_present = false(mx_val, 1); val_present(x) = 1; n = nnz(val_present);            
        n2 = length(unique(x));        
    elseif (frac < .08)        
        val_present = false(mx_val, 1); val_present(x) = 1; n = nnz(val_present);            
        n2 = length(unique(x));        
    elseif (frac < .09)
        val_present = false(mx_val, 1); val_present(x) = 1; n = nnz(val_present);            
        n2 = length(unique(x));        
    else
        val_present = false(mx_val, 1); val_present(x) = 1; n = nnz(val_present);            
        n2 = length(unique(x));                        
    end
    assert(isequal(n, n2));
    
end 