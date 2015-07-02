function [ori_pref_deg, pval_stats] = getPreferredOriFromOSP(R, oris_deg, spf_pref_idx, gratingType, combine_flag, showWorking) 
    % stimType is 'Flashed Grating' (Default), or 'Drifting Grating'

    showWorking = exist('showWorkingFlag', 'var') && ~isempty(showWorkingFlag);
%     showWorking = 0;
    
    doMonteCarloTest = nargout > 1;
        doAllValsTest = 0;
    
    useAllSpfs = 1;
    combineOppositeDirections_forDriftingGratings = exist('combine_flag', 'var') && isequal(combine_flag, 1);    
    
    fig_offset = combineOppositeDirections_forDriftingGratings*10;
    
    nOri_orig = size(R,1);
    if strcmp(gratingType, 'drifting') && combineOppositeDirections_forDriftingGratings
%         mean(mean(mean(R_full_atSp, 4),3),2)        
        R = R(1:nOri_orig/2, :,:,:) + R(nOri_orig/2+1:nOri_orig, :,:,:);        
        oris_deg = oris_deg(1:nOri_orig/2);
    end        
    oris_deg = oris_deg(:);   
    

    [nOri, nSpf, nPh, nReps] = size(R);
    calc_pval = nReps > 1;
    if nReps > 1        
        R_full = R;
        R = mean(R, 4);
    else
        if nargout > 1
            error('Can''t calculate pval without multiple trials');
        else
            R_full = R;
        end
    end

    R_oriSp = mean(R,3);
    
    
    oris_rad = deg2rad(oris_deg);
    
    % get ori-tuning curve (and full set of responses) at the preferred spatial frequency / all spfs
    if useAllSpfs
        r = mean(R_oriSp, 2);        
        R_full_atSp = mean(R_full, 2);        
    else        
        sp_i = getSelectedSpatialFrequency(R_oriSp, spf_pref_idx)
        r = R_oriSp(:,sp_i);        
        R_full_atSp = R_full(:,sp_i,:,:);
    end
    
    R_trials = reshape(mean(R_full_atSp, 3), [nOri, nReps]);
    
  
    
    % calculate preferred orientation.    
    compressToOE = 0;
    if compressToOE
        %%
        idx_odd = 1:2:nReps;
        idx_even = 2:2:nReps;
        R_trials = [mean(R_trials(:,idx_odd), 2), mean(R_trials(:,idx_even), 2)];
        nReps = 2;
    end    
    
            
    oris_rad_x2 = 2*oris_rad;
    cos_oris_rad_x2 = cos(oris_rad_x2);
    sin_oris_rad_x2 = sin(oris_rad_x2);
    
    [resultant_true, ori_pref_deg, Rij_true] = getOriPrefAndProjections(oris_rad_x2, cos_oris_rad_x2, sin_oris_rad_x2, r, R_trials);    
    if resultant_true == 0
        ori_pref_deg = nan;
    end
    
    [~, pval_stats.proj_ttest] = ttest(Rij_true);

    if doMonteCarloTest
        %%
        LaplaceProb = @(L,N) (L+1)/(N+2);
        rand('state', 0);
        nRand = 10000;
                
        
        [~, rand_idxs] = sort(rand(nOri, nRand), 1); % instead of calling    idx = randperm(nOri); each time
        
        if doAllValsTest  
            [ori_rand, resultant_rand] = deal(  zeros(1, nRand) );
            Rij_rand = zeros(nReps, nRand);
            mean_Rij_rand = zeros(1, nRand);
            
            for i = 1:nRand
                idx = rand_idxs(:,i);
                r_i = r(idx);            
                if doAllValsTest            
                    R_trials_i = R_trials(idx, :);
                    [resultant_rand(i), ori_rand(i), Rij_rand(:,i)] = getOriPrefAndProjections(oris_rad_x2, cos_oris_rad_x2, sin_oris_rad_x2, r_i, R_trials_i);
                else
                    resultant_rand(i) = getOriPrefAndProjections(oris_rad_x2, cos_oris_rad_x2, sin_oris_rad_x2, r_i);
                end
                %{
                if resultant_rand(i) > 300 && showWorking && 0
                    %%
                        oris_deg = rad2deg(oris_rad);
                        nshift = indmin( abs(oris_deg -ori_rand(i)-90));                     
                        figure(104+fig_offset); clf; plot(oris_deg, circshift(R_trials_i, nshift), 'b'); hold on; plot(oris_deg, circshift(r_i, nshift), 'ro-', 'linewidth', 2);  
                    3;
                end
                %}
            end
        
        else
            
            r_permuted = r(rand_idxs);                  
            resultant_rand = getOriPrefAndProjections(oris_rad_x2, cos_oris_rad_x2, sin_oris_rad_x2, r_permuted);
            
        end
        
        
       
        
        L_resultant = nnz(resultant_rand >= resultant_true);
        pval_stats.resultant_prob = LaplaceProb(L_resultant, nRand);
        
        if doAllValsTest
            mean_Rij_true = mean(Rij_true);
            mean_Rij_rand = mean(Rij_rand, 1);
        
            pval_stats.pval_U_allVals = ranksum(Rij_rand(:), Rij_true);
            pval_stats.pval_U_medVals = ranksum(median(Rij_rand, 1), Rij_true);

            [~, pval_stats.pval_T_allVals] = ttest2(Rij_rand(:), Rij_true);

            L_mean_Rij = nnz(mean_Rij_rand >  mean_Rij_true);
            pval_stats.pval_MC_meanVals = LaplaceProb(L_mean_Rij, nRand);                

            L_medVals = nnz(median(Rij_rand, 1) > median(Rij_true));
            pval_stats.pval_MC_medVals = LaplaceProb(L_medVals, nRand);        
        end
        
        %%
        if showWorking
            %%
            if doAllValsTest
                figure(100+fig_offset); clf; hist(mean_Rij_rand, 80); drawVerticalLine(mean_Rij_true, 'color', 'r');
                title(sprintf('Mean vals (resultant): p = %.2g', pval_stats.pval_MC_meanVals));
            end
%             figure(101+fig_offset); clf; hist(Rij_rand(:), 80); drawVerticalLine(Rij_true, 'color', 'r')
%             title(sprintf('All Vals: p_U = %.2g. p_T = %.2g', pval_stats.pval_U_allVals, pval_stats.pval_T_allVals));
    %         figure(102); hist(resultant_rand, 80); drawVerticalLine(resultant_true, 'color', 'r')
    %         profile viewer;
            oris_deg = rad2deg(oris_rad);
            
            %%
    %         pval_stats
            ori_max = roundToNearest(oris_deg(end), 180);
            nshift = indmin( abs(oris_deg -ori_pref_deg-90));                    
            figure(103+fig_offset); clf; 
            subplot(2,1,1);
            plot(oris_deg, circshift(R_trials, nshift), '.-', 'linewidth', 1); hold on; 
            xlim([0 ori_max]); set(gca, 'xtick', [0:45:ori_max])
            subplot(2,1,2);
            r_odd = mean(R_trials(:, 1:2:nReps), 2);
            r_even = mean(R_trials(:, 2:2:nReps), 2);
            plot(oris_deg, circshift(r_odd, nshift), 'ko-', 'linewidth', 1);  hold on;
            plot(oris_deg, circshift(r_even, nshift), 'ro-', 'linewidth', 1);             
            xlim([0 ori_max]); set(gca, 'xtick', [0:45:ori_max])
            3;
        end        
    end
    
    


    if showWorking && 0
        x0 = zeros(size(s));
        figure(101); clf;
        quiver(x0, x0, c, s, 'b');
        hold on;
        hQ = quiver(0, 0, sum(c), sum(s), 'r');
%         set(hQ, 'linewidth', 2)
        axis equal;
        title(sprintf('Preferred ori = %.2f\\circ. pval = %.3f', ori_pref_deg, ori_sel_pval));
        
        figure(102); clf;
        imagesc(n_k_ij);
        3;
    end
        
end


function [resultant, ori_pref_deg, Rij] = getOriPrefAndProjections(oris_rad_x2, cos_oris_rad_x2, sin_oris_rad_x2, r, R_trials)
%     [ori_pref_deg,  = getOriPrefAndProjections(oris_rad, r, 0)
%%
    if size(r, 2) > 1
        s = bsxfun(@times, r, sin_oris_rad_x2);
        c = bsxfun(@times, r, cos_oris_rad_x2);
    else        
        
        s = r .* sin_oris_rad_x2;
        c = r .* cos_oris_rad_x2;    
    end
%     s_rep = R_trials .* sin_2theta(:, ones(1, nReps));
%     c_rep = R_trials .* cos_2theta(:, ones(1, nReps));
    %%
    resultant = hypot(sum(c, 1), sum(s, 1));
    
    if nargout >= 2
        ori_pref_rad = 0.5 * circMean(r, oris_rad_x2);
        ori_pref_deg = rad2deg(ori_pref_rad); %mod(ori_pref_deg, 180);
                
        if ori_pref_deg > 360
            3;
        end
        
    end

    if nargout >= 3
        Rij = pval_projOntoPrefOri(R_trials, oris_rad_x2, ori_pref_rad);
    end
    
    3;
%     if calc_pval && (nargout >= 2)
%     n_k_ij = get_n_k_ij_from_R_trials(R_trials, gratingType);
%     end
end


function [Rij] = pval_projOntoPrefOri(n_k_ij, oris_rad_x2, ori_pref_rad)
    [nOris, nReps] = size(n_k_ij);
%     assert(nOris == length(oris_rad));

%     ori_proj_vector = cos( 2* (oris_rad - ori_pref_rad) );
    ori_proj_vector = cos( oris_rad_x2 - 2*ori_pref_rad);
    ori_proj_vector_ext = ori_proj_vector(:, ones(nReps,1));

    % Rij (a vector) = component of the response of each full-orientation set along the preferred orientation. 
    % each value of R_ij is the "orientation-averaged orientation bias" of one trial (ie. the
    % bias in the direction of the preferred orientation)
    Rij = sum( ori_proj_vector_ext .* n_k_ij, 1 ); % sum over orientations.

%     mean_Rij = mean(Rij);  
end


function n_k_ij = get_n_k_ij_from_R_full(R_full_atSp, gratingType)
    
    [nOri, nSp, nPh, nReps] = size(R_full_atSp);
    assert(nSp == 1);
    switch gratingType
        case 'flashed', phaseAction = 'average';  % phaseAction = 'counterphase'; originally chose this (?), but hard to justify.
        case 'drifting', phaseAction = 'average';
    end
            
            
    % "n_k_ij = the number of spikes elicited by the jth cycle of the ith repetition of the grating of orientation k."  
    % in R_full, cycles & repetitions are just concatenated one after the other, so we just
    % take the entire 2nd dimension (with all trials).

    n_k_ij_nPhXnRep = reshape( R_full_atSp, [nOri, nPh, nReps]);      

    switch phaseAction
        case 'all phases', % take all phases individually
            n_k_ij_allPh = reshape( R_full_atSp, [nOri, nPh*nReps]);     
            n_k_ij = n_k_ij_allPh;

        case 'average',  % average over phases    
            n_k_ij_avPh = reshape( mean(R_full_atSp,3), [nOri, nReps]);  
            n_k_ij = n_k_ij_avPh;

        case 'odd/even',  % average all odd phases, and all even phases
            idx_odd = 1:2:nPh; idx_even = 2:2:nPh;
            n_k_ij_oe = reshape([mean(n_k_ij_nPhXnRep(:, idx_odd, :), 2), mean(n_k_ij_nPhXnRep(:, idx_even, :), 2)], [nOri, nReps*2]);
            n_k_ij = n_k_ij_oe;

        case 'counterphase',  % average each phase with its opposite phase
            idx_cph = arrayfun(@(i) [i, i+nPh/2], 1:nPh/2, 'un', 0);                
            n_k_ij_cPh = cellfun(@(i) mean(n_k_ij_nPhXnRep(:, i, :),2), idx_cph, 'un', 0);
            n_k_ij_cPh = reshape( [n_k_ij_cPh{:}], [nOri, nReps*nPh/2]);
            n_k_ij = n_k_ij_cPh;

        case 'consecutive',  % average each pair of consecutive phases
            idx_consec = arrayfun(@(i) [i, i+1], 1:2:nPh, 'un', 0);
            n_k_ij_consec = cellfun(@(i) mean(n_k_ij_nPhXnRep(:, i, :),2), idx_consec, 'un', 0);
            n_k_ij_consec = reshape( [n_k_ij_consec{:}], [nOri, nReps*nPh/2]);
            n_k_ij = n_k_ij_consec;

    end

    
%         ori_sel_pval(2) = pval_projOntoPrefOri(n_k_ij_oe, oris_rad, ori_pref_rad);    
%         ori_sel_pval(3) = pval_projOntoPrefOri(n_k_ij_cPh, oris_rad, ori_pref_rad);   
%         ori_sel_pval(4) = pval_projOntoPrefOri(n_k_ij_consec, oris_rad, ori_pref_rad); 
%         ori_sel_pval(5) = pval_projOntoPrefOri(n_k_ij_allPh, oris_rad, ori_pref_rad); 
%                 
        
%         ori_proj_vector = cos( 2* (oris_rad - ori_pref_rad) );
%         ori_proj_vector_ext = ori_proj_vector(:, ones(nReps,1));
%         ori_proj_vector_ext_allPh = ori_proj_vector(:, ones(nReps*nPh,1));
% 
%         % Rij (a vector) = component of the response of each full-orientation set along the preferred orientation. 
%         % each value of R_ij is the "orientation-averaged orientation bias" of one trial (ie. the
%         % bias in the direction of the preferred orientation)
%         Rij = sum( ori_proj_vector_ext .* n_k_ij, 1 ); % sum over orientations.
%         Rij_allPh = sum( ori_proj_vector_ext_allPh .* n_k_ij_allPh, 1 ); % sum over orientations.
% 
%         [h,ori_sel_pval] = ttest(Rij);            
%         [h,ori_sel_pval(2)] = ttest(Rij_allPh);           
    
end


function sp_i = getSelectedSpatialFrequency(R_oriSp, spf_pref_idx)
    % code for selecting which spatial frequency (flashed gratings) -- 
%         no longer used: instead, we combine all spatial frequencies.

    if exist('spf_pref_idx', 'var') && ~isempty(spf_pref_idx)
        sp_i  = spf_pref_idx;
    else

        nSpf = size(R_oriSp, 2);
        if (nSpf == 1)
            av_option = 3;
        else
            av_option = 3;
        end

        switch av_option
            case 1,  % option 1: no averaging - pick spf with max response:
                [~, ind_max] = maxElement(R_oriSp);
                [ori_i, sp_i] = elements(ind_max);

            case 2,  % option 2: average over phases - then pick spf with max response:
                [~, ind_max] = maxElement(mean(R_oriSp,1));
                [ori_i, sp_i] = elements(ind_max);

            case 3,  % option 3: average over phases and oris, then spf with max response:
                [~, sp_i] = max(mean(R_oriSp,1));
        end        
    end

end


%   if nargout >= 3
%         magSumRs = norm(double([sum(c), sum(s)]));
%         magRs    = normV([s, c], 2);
%         ori_sel_strength = magSumRs / mean(magRs(magRs >0));
%     end


%{
% persistent ori_rad_1_prev ori_pref_rad_prev ori_proj_vector_ext_prev
% 
%     [nOris, nReps] = size(n_k_ij);
%     
% %     assert(nOris == length(oris_rad));
%     if ~isempty(ori_rad_1_prev) && (oris_rad(1) == ori_rad_1_prev) && (ori_pref_rad == ori_pref_rad_prev) && 
%         ori_proj_vector_ext = ori_proj_vector_ext_prev;
%     else        
%         ori_proj_vector = cos( 2* (oris_rad - ori_pref_rad) );
%         ori_proj_vector_ext = ori_proj_vector(:, ones(nReps,1));
%         
%         ori_rad_1_prev     = oris_rad(1);
%         ori_pref_rad_prev  = ori_pref_rad;
%         ori_proj_vector_ext_prev = ori_proj_vector_ext;        
%     end
%         
%     % Rij (a vector) = component of the response of each full-orientation set along the preferred orientation. 
%     % each value of R_ij is the "orientation-averaged orientation bias" of one trial (ie. the
%     % bias in the direction of the preferred orientation)
%     Rij = sum( ori_proj_vector_ext .* n_k_ij, 1 ); % sum over orientations.
% 
% %     mean_Rij = mean(Rij);  
% end
%}