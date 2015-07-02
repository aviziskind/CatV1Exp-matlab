function [x_dists, x_dists_lo, x_dists_hi,   y_medians, y_val_lo, y_val_hi,   y_pvalues, y_rand_med, y_rand_lo, y_rand_hi] = getSlidingWindowMediansVsDist(depthDist, X, X_rand, opt)
                    
    nRand = length(X_rand);

    nPairsEach = opt.nPairsEach;
    depthDist_sorted = sort(nonnans(depthDist));
    depthLims = lims(depthDist_sorted);
    nDists = length(depthDist_sorted);
    x_statFunc = @nanmedian;
    y_statFunc = @nanmedian;
    if isfield(opt, 'x_stat_func')
        x_statFunc = opt.x_stat_func;
    end
    if isfield(opt, 'y_stat_func')
        y_statFunc = opt.y_stat_func;
    end

    overlap_f = opt.overlap_f;
    nSteps = ceil( nDists / nPairsEach * overlap_f );

%                     nPairsEach = 2000;
    dLims = linspace(depthLims(1), depthLims(2), nSteps+1);
    [x_dists, x_dists_lo, x_dists_hi, y_medians, y_val_lo, y_val_hi, y_pvalues,  y_rand_med, y_rand_lo, y_rand_hi] = deal( nan(1, nSteps) );
    depth_n = 1: floor(nPairsEach/overlap_f) : nDists;
%                         depth_n = round(linspace(1, length(depthDist_sorted), nSteps+1));
%                     depth_n = round(linspace(1, length(depthDist_sorted), nSteps+1));
%%                  

    prctile_range_vals = 45; %[ 25 - 75];
    prctile_range_rand_vals = 45; %[ 5 - 95];
%                         prctile_range = 45; %[ 5 - 95];

    all_percentiles_val = 50 + [-1, 0, 1]*prctile_range_vals;
    all_percentiles_rand_val = 50 + [-1, 0, 1]*prctile_range_rand_vals;

    vals_thisRange = cell(1, nSteps);
    lastOne = false;
    for i = 1:nSteps
        %%
        idxs = [depth_n(i), depth_n(i)+nPairsEach];
        if idxs(2) > nDists
            idxs(2) = nDists;
            lastOne = true;
        end
        dLims_i = depthDist_sorted(idxs);

        idx_range_i = find( ibetween(depthDist, dLims_i) );

        vals_thisRange{i} = X( idx_range_i );

        randVals_thisRange = cellfun(@(X) X(idx_range_i), X_rand, 'un', 0);

        x_dists(i) = x_statFunc( depthDist( idx_range_i ) );

        [x_dists_lo(i), x_dists_med, x_dists_hi(i)] = dealV( prctile( depthDist( idx_range_i ), all_percentiles_val ) );
%         assert(x_dists_med == x_dists(i));

        y_medians(i) = y_statFunc(  vals_thisRange{i} );
        rand_y_medians = cellfun(y_statFunc, randVals_thisRange );

        [y_val_lo(i), y_val_med, y_val_hi(i)] = dealV( prctile( vals_thisRange{i}, all_percentiles_val ) );
%         assert(y_val_med == y_medians(i));

%                             p_val = signrank(vals_thisRange{i}, opt.bs_median);
        p_val = getRandomizedProb( y_medians(i), rand_y_medians, opt.sigTestTail);
        [y_rand_lo(i), y_rand_med(i), y_rand_hi(i)] = dealV( prctile( rand_y_medians, all_percentiles_rand_val ) );
        y_rand_med(i) = y_statFunc(rand_y_medians);

        y_pvalues(i) = p_val;

%                         depth_idxs = idxs(1) : idxs(2);
%                         xs_depth2(i) = median( depthDist_sorted( depth_idxs ) );
%                             depth_medians(i) = nanmedian(  dX_bcc(samePen_bcc==1 & ibetween(depthDist, dLims_i) ) );

        if lastOne
            break;
        end
    end

    %%
    idx_use = ~isnan(y_medians);

%     x_dists = x_dists(idx_use);
%     y_medians = y_medians(idx_use);
%    
    [x_dists, x_dists_lo, x_dists_hi,  y_medians, y_val_lo, y_val_hi, y_pvalues, y_rand_med, y_rand_lo, y_rand_hi] = deal(...
    x_dists(idx_use), x_dists_lo(idx_use), x_dists_hi(idx_use),  y_medians(idx_use), y_val_lo(idx_use), y_val_hi(idx_use), y_pvalues(idx_use), y_rand_med(idx_use), y_rand_lo(idx_use), y_rand_hi(idx_use));
         
     
     

end