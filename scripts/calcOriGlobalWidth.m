function w_std = calcOriGlobalWidth(r_k, oris_deg, ori_pref_deg, gratingType, ori_diff_max, idx_ori_use)

    r_k = r_k(:);
    oris_deg = oris_deg(:);
    %only include ori diffs within +/-90 of peak.    

%     ori_diffs = mod(oris_deg(:) - ori_pref_deg, 180);
%     ori_diffs = min(ori_diffs, ori_diff_max-ori_diffs);        
    if ~exist('idx_ori_use', 'var')
        idx_ori_use = true(size(r_k));
    end

    if islogical(idx_ori_use)
        idx_ori_unUnused = ~idx_ori_use;
    else
        idx_ori_unUnused = setdiff(1:length(r_k), idx_ori_use);
    end

    
    if ~exist('gratingType', 'var') || strcmp(gratingType, 'flashed')
        assert(nnz(idx_ori_use) == length(r_k))
    elseif strcmp(gratingType, 'drifting') % drifting gratings
        if length(idx_ori_use) < length(oris_deg)
            assert(nnz(idx_ori_use)/length(oris_deg) == 1/2 );
%             assert(sum(r_k(idx_ori_use)) > sum(r_k(idx_ori_unUnused)));
        end
    end
    
    if nargin < 3 || isempty(ori_pref_deg)
        ori_pref_deg = rad2deg( circMean(r_k, deg2rad(oris_deg)) );
    end
        
    ori_diffs = circDist(oris_deg(:), ori_pref_deg, 180);    
            
%     ori_diffs = mod(oris_deg(:) - ori_pref_deg, ori_diff_max);
%     ori_diffs = min(ori_diffs, ori_diff_max-ori_diffs);    
    
    a = sum( r_k(idx_ori_use) .* ori_diffs(idx_ori_use).^2 );
    b = sum( r_k(idx_ori_use) );
    w_std = sqrt(a/b);
    
    dbug = 0;
    if dbug
                
        %%
%         figure(332); clf;
%         ori_diffs_raw = oris_deg(:) - ori_pref_deg;
%         plot(oris_deg, ori_diffs_raw, 'b', oris_deg, ori_diffs, 'g', oris_deg(idx_ori_use), ori_diffs(idx_ori_use), 'r.')
%         set(gca, 'ytick', [-360:45:360]);
%         drawHorizontalLine([-180, 180]);
%         ylim([-360, 360])
        
        figure(333); clf;
        plot(oris_deg, r_k); hold on;
        plot(oris_deg(idx_ori_use), r_k(idx_ori_use), 'r.')
        plot(oris_deg(idx_ori_use(1)), r_k(idx_ori_use(1)), 'ro')
        plot(oris_deg(idx_ori_use(end)), r_k(idx_ori_use(end)), 'rs')
        drawVerticalLine(ori_pref_deg);
        drawVerticalLine(mod(ori_pref_deg+[-90, 90], 360), 'color', 'r');
        drawVerticalLine(mod(ori_pref_deg+[-w_std, w_std], 360), 'color', 'g', 'linestyle', ':' );
        
        3;
    end
    
    
        
end