function [PSTHdata, stats] = getPSTHforCell(Gid, cellId)
   
    gt = getGratingStimType(Gid);            

    if ~strcmp(gt.gratingType, 'flashed')
        PSTHdata = [];
        return;
%         error('This function is for flashed grating cells');
    end
    
    frameLength_ms = getFrameLength('Gid', Gid);    
    nBinsPerFrame = round( frameLength_ms / (25/6) );
    

    psthWindow = [-300, 200];    
    [PSTH_bins, PSTH_vals, PSTH_stats] = getCellPSTHvals(Gid, cellId);        

    % 3. Find the window in which cell's response was reproducible. 
%     [timeWindow, windowProfile] = getBestTimeWindowFromPSTH(PSTH_bins,
%     PSTH_vals);
    meanRate = PSTH_stats.meanRate;
    default_bins = [20 120];
    statName = 'cc_p';
    intervalSize = 5;
    nStdTh = 3;
    objTh = 3;
    stimType = getGratingStimType(Gid);     
    allStats = getPSTHwindowData(Gid, cellId, statName, psthWindow);
        
    [L_bin, R_bin, stat_bin_x, stat_val_y, statVals_m, statVals_s] = ...
        getLRbin(Gid, cellId, PSTH_bins, allStats, intervalSize, statName, stimType, nStdTh, objTh, meanRate);
    
    isCellReproducible = (L_bin > 0);
    if ~isCellReproducible
        L_bin = indmin(abs(PSTH_bins-default_bins(1)));
        R_bin = indmin(abs(PSTH_bins-default_bins(2)));
    end    
    if isfield(allStats, 'clustIds')
        allStats = rmfield(allStats, 'clustIds');
    end
    windowStats = structfun(@(X) X(R_bin, L_bin), allStats, 'un', 0);
    if windowStats.cc_p == 0
        3;
    end
      
    %%
    window_start_ok_ms = [10, 85];
    idx_ok = ibetween(PSTH_bins(1:end-nBinsPerFrame+1), window_start_ok_ms);
    diag_stimw = diag(allStats.cc_p, -(nBinsPerFrame-1));
    ind_max = indmax(diag_stimw(idx_ok)) +  find(idx_ok,1)-1;
    [L_bin_stimw, R_bin_stimw] = deal(ind_max, ind_max + nBinsPerFrame-1);
    %%
   
    windowProfile_stimw = PSTH_vals(L_bin_stimw:R_bin_stimw);
    assert(any(windowProfile_stimw))
    
%     [L_bin_stimw, R_bin_stimw] = nBinsPerFrame
    
    stats = struct('statName', statName, 'isRep', isCellReproducible, ...
        'stat_bin_x', stat_bin_x, 'stat_val_y', stat_val_y, 'statVals_m', statVals_m, 'statVals_s', statVals_s, ...
        'allWindowStats', windowStats);

    bw = diff(PSTH_bins(1:2));    
    timeWindow    = PSTH_bins([L_bin, R_bin]) + [-1; 1]*bw/2;
    windowProfile = PSTH_vals(L_bin:R_bin);
    if all(windowProfile==0) 
        windowProfile(:) = 1;
    end
    
    PSTHdata = struct('bins', PSTH_bins, 'vals', PSTH_vals, ...
            'frameLength_ms', frameLength_ms, 'meanRate', PSTH_stats.meanRate, 'bckgRate', PSTH_stats.bckgRate, ...
            'timeWindow_bins', [L_bin, R_bin], 'timeWindow_ms', timeWindow(:)', 'windowProfile', windowProfile(:)', ...
            'timeWindow_stimw',[L_bin_stimw, R_bin_stimw], 'windowProfile_stimw', windowProfile_stimw );    
end

%   Old method of calculating PSTH of short-frame flash-grating cells:
%         [PSTH_bins, PSTH_vals] = calculatePSTHiteratively(Gid, cellId);
