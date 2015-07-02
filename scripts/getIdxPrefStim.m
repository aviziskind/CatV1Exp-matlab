function idxs = getIdxPrefStim(Gid, R_full)


    % preferred stim with no smoothing.
    R_osp = mean(R_full,4);

    R_avPh_noSm  = mean(R_osp,3);     % take average over phases
    R_maxPh_noSm = max(R_osp, [], 3); % take MAX over phases

    idxs.prefR_avP_noSm  = getOriSpf_idx (  R_avPh_noSm,  1 );   % no smoothing, mean over phases, preferred Stim
    idxs.maxR_avP_noSm   = getOriSpf_idx (  R_avPh_noSm,  0 );   % no smoothing, mean over phases, highest-responding stim
    idxs.prefR_maxP_noSm = getOriSpf_idx (  R_maxPh_noSm, 1 );   % no smoothing, max over phases,  preferred Stim
    idxs.maxR_maxP_noSm  = getOriSpf_idx (  R_maxPh_noSm, 0 );   % no smoothing, max over phases,  highest-responding stim

    
    % preferred stim with smoothing.
    R_osp_smoothed = getSmoothedOSP(Gid, R_osp);

    R_avPh_sm   = mean(R_osp_smoothed,3);
    R_maxPh_sm  = max(R_osp_smoothed, [], 3);
    
    idxs.prefR_avP_sm  = getOriSpf_idx (  R_avPh_sm,  1 );  % smoothing, mean over phases, preferred Stim
    idxs.maxR_avP_sm   = getOriSpf_idx (  R_avPh_sm,  0 );  % smoothing, mean over phases, highest-responding stim
    idxs.prefR_maxP_sm = getOriSpf_idx (  R_maxPh_sm, 1 );  % smoothing, max over phases,  preferred Stim
    idxs.maxR_maxP_sm  = getOriSpf_idx (  R_maxPh_sm, 0 );  % smoothing, max over phases,  highest-responding stim
    
    
end

function idx = getOriSpf_idx (  R_os, prefStim_flag )

    getPreferredStim = isequal(prefStim_flag, 1);
    if getPreferredStim
        ori_pref_idx = indmax( sum(R_os, 2) );              % average over spatial frequencies when estimating preferred orientation
        spf_pref_idx = indmax( R_os(ori_pref_idx(1), :) );  % use preferred orientation when estimating preferred spatial frequency.

    else % simply take max over all stimuli:
        
        [~, ind_max] = max(R_os(:));
        [ori_pref_idx, spf_pref_idx] = ind2sub(size(R_os), ind_max);
%         assert(isequal( [ori_pref_idx, spf_pref_idx], [o_idx, s_idx]));
    end
    idx = [ori_pref_idx, spf_pref_idx];

end


function R_osp_smoothed = getSmoothedOSP(Gid, R_osp)
    % preferred stim after smoothing.
    sd = siteDataFor(Gid);
    stimType = getGratingStimType(Gid);
    switch stimType.gratingType
        case 'drifting',
            isOrientBatch = strcmp(stimType.driftingType, 'Orientation Batch');
            isSpatBatch   = strcmp(stimType.driftingType, 'Spatial Frequency Batch');
        case 'flashed'
            isOrientBatch = 1;
            isSpatBatch = 1; 
    end
  
    R_osp_smoothed = R_osp;
    
    % smooth over phases:
    dw_ph = 1;
    R_osp_smoothed = gaussSmooth(R_osp_smoothed, dw_ph, 1, 3);

    
    % smooth over orientations (if orientation batch):
    if isOrientBatch
        dw_ori = 1;         
        R_osp_smoothed = gaussSmooth(R_osp_smoothed, dw_ori, 1, 1);
    end
    
    % smooth over spatial frequencies (if spatial frequency batch):
    if isSpatBatch
        log_spf = log(sd.spPeriod_pix);
        dw_spf = mean(diff(log_spf));
        R_osp_smoothed = gaussSmooth_nu(log_spf, R_osp_smoothed, dw_spf, 2);    
    end        

    
end