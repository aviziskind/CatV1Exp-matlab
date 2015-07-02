function [R, R_full, uori, usp, uph, utf_Hz, meanRate, bckgSamples] = getOriSpfPhaseProfile_simple(Gid, cellId, PSTH_data, keepAllTrialsFlag)

%     global keepFirstCycle
    % input options: either 
    % (1)  Gid, cellId, spikeWindow and windowProfile (from which can calculate 'relContrOfFrameToSpike'
    %         ie. getOriSpfPhaseProfile(Gid, cellId, spikeWindow, windowProfile, [, meanFiringRate]] )
    %                   OR
    % (2) Gid (so know which stimulus type) and relContrOfFrameToSpike
    %         ie. getOriSpfPhaseProfile(Gid, relContrOfFrameToSpike, [ [bckgRate, meanFiringRate] ])
    correctForCorticalExcitability = 0;
    classicTimeWindow_ms = [35 75];
    keepAllTrials = exist('keepAllTrialsFlag', 'var') && ~isempty(keepAllTrialsFlag) && (keepAllTrialsFlag == true);
%         strcmp(curGroupingType(''), 'cells');
    bckgSkip_ms = 250;
    
    stimType = getGratingStimType(Gid);
    gratingType = stimType.gratingType;% flashedOrDrifting(Gid, 1);                
    
    useSavedData = 1;
%     psthMethod, psthWindow_ms, trialGrouping, splitWindowIfCph, keepFirstDriftingGratingCycle,
    
    [uori, usp, uph, ~, ~, utf_Hz] = dbGetUniqueOriSpPh('Gid', Gid);

%     keepFirstCycle = 0;
%     getHistArgs = {'keepFirstDriftingGratingCycle', keepFirstCycle};
    getHistArgs = {};
%     R_full_t = dbGetStimulusTimes(Gid);    
    
    if keepAllTrials
        R_full_name = 'osp_full';
    else
        R_full_name = 'osp_ph_oe';
    end

    if strcmp(gratingType, 'flashed')  
        useClassicWindow = ~isfield(PSTH_data, 'timeWindow_bins');
        if ~useClassicWindow
            [L_bin, R_bin] = dealV(PSTH_data.timeWindow_bins);
            curPsth = PSTH_data.windowProfile;
        else
            binEdges = binCent2edge( PSTH_data.bins );
            LR_bins = binarySearch(binEdges, classicTimeWindow_ms);
            L_bin = LR_bins(1); R_bin = LR_bins(2);
%             db = diff(bins(1:2))/2;
%             L_bin = indmin( abs(classicTimeWindow_ms(1)+db - bins ));
%             R_bin = indmin( abs(classicTimeWindow_ms(2)-db - bins ));
            curPsth = [];                      
        end 
        meanRate = PSTH_data.meanRate;
        
        
        if useSavedData
            [R, R_full] = getOspDataForPsthWindow(Gid, cellId, [], [], L_bin, R_bin, curPsth, {'osp_ph', R_full_name});
        else
            [R, R_full] = calcOspForPsthWindow(Gid, cellId, L_bin, R_bin, false, curPsth, {'osp_ph', R_full_name}, getHistArgs, meanRate);
        end
            
            % adjust bckgSamples to match windowSize
        nSamplesAv = R_bin-L_bin+1;
        bckgSamples = getBackgroundSpikes(Gid, cellId, bckgSkip_ms, nSamplesAv);       
        if iscell(bckgSamples)
        3;
        end
        
        
    else % drifting gratings
        [L_bin, R_bin] = deal(1);  
        curPsth = [];
        if useSavedData
            [R, R_full, meanRate] = getOspDataForPsthWindow(Gid, cellId, [], [], L_bin, R_bin, curPsth, {'osp_ph', R_full_name, 'meanRate'});
        else
            [R, R_full, meanRate] = calcOspForPsthWindow(Gid, cellId, L_bin, R_bin, [], curPsth, {'osp_ph', R_full_name, 'meanRate'}, getHistArgs);
        end
        
        bckgBinSize_ms = 4+1/6;
        nSamplesAv = round((utf_Hz*1000)/bckgBinSize_ms);        
        bckgSamples = getBackgroundSpikes(Gid, cellId, bckgSkip_ms, nSamplesAv);       
        
    end
%     
%     
%     if keepAllTrials
%         R_full = R_full;
%     else
%         R_full = cat(4, R_odd, R_even);         
%     end    
    
        
    
    assert(~any(isnan(R_full(:))));    
    
    rescaleFactor = meanRate / mean(R(:));
    if strcmp(gratingType, 'drifting')
%         assert(max(rescaleFactor, 1/rescaleFactor) < 1.1)
    else    
        R_full = R_full * rescaleFactor;  % convert from # eff spikes --> #spikes/sec.
        R      = R      * rescaleFactor;  % convert from # eff spikes --> #spikes/sec.
    end    
    
    if 0 && keepAllTrials && correctForCorticalExcitability
        nTrials = size(R_full, 4);
        state_func = getCorticalState(Gid); %
        
        R_full_t = dbGetStimulusTimes(Gid);
        R_full_state = reshape(feval(state_func, R_full_t(:)), size(R_full_t));
        
                
        trialAv_spikes = mean(mean(mean(R_full,1),2),3); trialAv_spikes= trialAv_spikes(:);
%         trialAv_time = mean(mean(mean(R_full_t,1),2),3); trialAv_time = trialAv_time(:);

        trialAv_state = mean(mean(mean(R_full_state,1),2),3); trialAv_state = trialAv_state(:);        
        3;
        state_linfit = polyfit(trialAv_state, trialAv_spikes, 1);
        m = state_linfit(1);
        c = state_linfit(2);
        
        R_full_true = (R_full -c )./(m*R_full_state);
        
        3;  
%         R_full = R_full_true;        
        3;
        idx_odd = 1:2:nTrials;
        idx_even = 2:2:nTrials;
        
        R_odd = R_full(:,:,:,idx_odd);
        R_even = R_full(:,:,:,idx_even);
        R_odd_true = R_full_true(:,:,:,idx_odd);
        R_even_true = R_full_true(:,:,:,idx_even);
        
        R_true = mean(R_full_true, 4);
        R = mean(R_full, 4);
        
        figure(15); 
        [r,p] = corr(R_odd(:), R_even(:));
        [r_true,p_true] = corr(R_odd_true(:), R_even_true(:));
        subplot(1,2,1); plot(R_odd(:), R_even(:), '.');
        subplot(1,2,2); plot(R_odd_true(:), R_even_true(:), '.');
        3;
                
        figure(16);
        subplot(1,2,1); imagesc(mean(R,3));
        subplot(1,2,2); imagesc(mean(R_true,3));
        3;
        
        idx_spf_best = indmax( mean(mean(R, 1),3));
        [nOri, nSpf, nPh, nTrials] = size(R_full);
        R_ori_trials = reshape( mean(R_full(:, idx_spf_best, :, :), 3), [nOri, nTrials]);
        R_ori_trials_true = reshape( mean(R_full_true(:, idx_spf_best, :, :), 3), [nOri, nTrials]);
        
        figure(17); clf
        b = mean( std(R_ori_trials,[], 2)' / mean(R_ori_trials(:)) );
        b_true = mean( std(R_ori_trials_true,[], 2)' / mean(R_ori_trials_true(:)) );
        
        subplot(1,2,1); errorbar(1:nOri, mean(R_ori_trials,2), std(R_ori_trials,[], 2), 'bo-'); title(sprintf('orig : %.2f', b));
        subplot(1,2,2); errorbar(1:nOri, mean(R_ori_trials_true,2), std(R_ori_trials_true,[], 2), 's-'); title(sprintf('adjusted %.2f', b_true));
        3;
        
        
%         cols = jet(102);
%         trialAv_state_norm = (trialAv_state-min(trialAv_state))/diff(lims(trialAv_state));
%         trialAv_state_norm_int = ceil(trialAv_state_norm*100)+1;
%         subplot(1,2,1); hold on;
%         for i = 1:nTrials
%             plot(1:nOri, R_ori_trials(:,i), 'color', cols(trialAv_state_norm_int(i),:) ); 
%         end
%         axis tight;
%         
%         subplot(1,2,2); hold on;
%         for i = 1:nTrials
%             plot(1:nOri, R_ori_trials(:,i), cols(trialAv_state_norm_int,:) ); 
%         end
        3;
        
    end
    
    
%     frameLength_sec = getFrameLength('Gid', Gid, 'sec');    
%     R_full = R_full * / frameLength_sec; % convert from # eff spikes --> #spikes/sec.
%     R = R           / frameLength_sec; % convert from # eff spikes --> #spikes/sec.

%     if compressR_fullToInt
%         f.scale = 255/max(R_full(:));
%         f.idxnan = uint32( find(isnan(R_full(:))) );
%         R_full = uint8(R_full*f.scale);
%     end
        

end



    %{
    switch gratingType, 
        % there are 3 reasons where we want to keep all trials (instead of just odd/even trials)
        % (1) for testing significance of ori tuning (2) testing significance of response to
        % background (3) if we want to reweight by cortical excitability later on. - only need to do
        % these with 'cells' grouping, not with clusters/clustersPruned.
        
        case 'flashed',
            psthWindow = [-300 200];
            retreiveHiddenTrials = 1;  % keep counterphase trials for counterphase flashed gratings exp.
        case 'drifting', 
            psthWindow = [0 8+1/3];
            retreiveHiddenTrials = 0;  % don't keep first cycle (with artifacts)
    end
    
    %     histOpts = struct('psthWindow_ms', psthWindow, 'osp_method', osp_method, 'trialGrouping', 'individual');
%     [bins, stimPSTH_vals, meanRate, bckgSamples] = dbGetCellSpkStimHists(Gid, cellId, histOpts);
    %}
         
%         targetedPSTH = targetedPSTH/sum(targetedPSTH);            
%         if isempty(targetedPSTH) || isnan(targetedPSTH(1)) 
%             error('invalid psth');
%         end
%         stim_tgt_window = stimPSTH_vals_oe(L_bin:R_bin,:, :);
%         stim_tgt_wgt_window = bsxfun(@times, stim_tgt_window, targetedPSTH(:));
%         R_full = reshape( mean ( stim_tgt_wgt_window, 1), [nOri, nSp, nPh, 2]);    
%         R = mean(R_full, 4);    
    
%     else       
%         if nCycles ~= round(nCycles)
%             error('not implemented yet');
%         end
%         nCycles = ceil(nCycles); % for case in which is a fraction
%         nCycTotal = nCycles*nRep;
%         
%         R_full = reshape(stimPSTH_vals_oe, [nOri, nSp, nPh, nCycTotal]);
%                            
%         okTrialsInds = bsxfun(@plus, [0:nRep-1]'*nCycles, 2:nCycles)'; % discard first cycle which may contain artifacts (of sudden contrast change from mean luminance)
%         okTrialsInds = okTrialsInds(:)';
% %         assert(nTrialsMax == max(nTrialsTheseParams(:)));    
%         
%         R = nanmean(R_full(:,:,:,okTrialsInds), 4);
% %         R(nTrialsTheseParams == 0) = 0;        
%         if any(isnan(R(:)))
%             error('Nan detected in R');
%         end        
  