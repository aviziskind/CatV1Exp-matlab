function [PSTH_bins, stimPSTH_vals, meanRate, bckgSamples, spkTsRelToFrame_ms, spike_firstFrameId] = getIndividualPSTHs(Gid, cellId, stimGroupingStyle, psthWindow_ms, bkgrSkip_ms, trialGrouping)
    if nargin < 3
        stimGroupingStyle = 'OSP';
    end
    separateCphOrderTrials_ifCph = true;
    joinBinsForDriftingGrating = true;

    if ~exist('trialGrouping', 'var') || isempty(trialGrouping)
    %     trialGrouping = 'odd/even'; % options: 'odd/even' or 'individual'
        trialGrouping = 'individual'; % options: 'odd/even' or 'individual'
    end
    
    frameLength_ms = getFrameLength('Gid', Gid);

    [nOri, nSpf, nPh, nTrials, isCphFlashed] = getGratingStimType(Gid);
    
    separateCphOrderTrials = separateCphOrderTrials_ifCph && isCphFlashed; % relevant if not separating every trial    
    separateEveryTrial = strncmp(trialGrouping, 'individual', 2);
    separateOddEvenTrials = strncmp(trialGrouping, 'odd/even', 2); 
    
    if ~exist('psthWindow_ms', 'var') || isempty(psthWindow_ms)
        if flashedOrDrifting(Gid) == 1
            [nBinsPerFrame, nFramesPerExtFrame, extFrameLength_ms] = getBinSizeParams(frameLength_ms);          
        else
            [nBinsPerFrame, nFramesPerExtFrame, extFrameLength_ms] = deal(2, 1, frameLength_ms);
        end
        psthWindow_ms = [0 extFrameLength_ms];        
    else
        if length(psthWindow_ms) == 1
            psthWindow_ms = [0 psthWindow_ms];
        end
        [nBinsPerFrame, nFramesPerExtFrame, extFrameLength_ms] = getBinSizeParams(frameLength_ms, diff(psthWindow_ms));          
        psthWindow_ms(2) = psthWindow_ms(1)+extFrameLength_ms;
    end        
    nBinsPerExtFrame = nBinsPerFrame * nFramesPerExtFrame;
    
    if ~exist('bkgrSkip_ms', 'var') || isempty(bkgrSkip_ms)
        bkgrSkip_ms = 200;
    end        

    % 1. Get Position of spike times relative to the frames
    [spkTsRelToFrame_ms, bckgSamples, meanRate, whichPres, spike_firstFrameId] = getParsedSpikes('timing', Gid, cellId, bkgrSkip_ms, psthWindow_ms);     
    
    % 2 Use these spike times to calculate the (extended) PSTHs of individual stimuli 
    frameStimIds = getStimulusFrameSequence(Gid, stimGroupingStyle);
    frameStimIds = frameStimIds(:)';
    [uStimIds, stimIdsIdx] = uniqueList(frameStimIds);
    stimIdsIdx = cellfun(@(x) x(:)', stimIdsIdx, 'un', 0);
    
    nStims = cellfun(@length, stimIdsIdx);
    sameNumAllStims = all(nStims == max(nStims));

    if separateEveryTrial
        N3 = nTrials;
        maxNStims = max(nStims);
        if nTrials ~= maxNStims  
            nWithExtra = nnz( nStims > nTrials);            
            assert(nWithExtra / length(nStims) < .01); % make sure is just a handful of stimuli with extra repetitions N3 = nTrials;
            assert( all(nStims) <= nTrials +1); % make sure are a maximum of 1 over the # reps
            
            %just remove the extra stimuli
            for idx_withExtra = find(nStims > nTrials);
                stimIdsIdx{idx_withExtra} = stimIdsIdx{idx_withExtra}(1:nTrials);
            end
            nStims = cellfun(@length, stimIdsIdx);
            sameNumAllStims = all(nStims == max(nStims));
            
        end
        
    else
        N3 = (separateOddEvenTrials+1)*(separateCphOrderTrials+1);
    end

    % possibilties : cph &  every  xor  o/e
    % Coa
    % cOa
    % COa
    % coa
    % CoA
    % coA
        
    
    nStimIds = length(uStimIds);    
    stimPSTH_vals = nan(nBinsPerExtFrame,  nStimIds, N3, 'single');
    for stim_i = 1:nStimIds
        %%
        frameIndsForStimI = stimIdsIdx{stim_i};        
        nFrames = length(frameIndsForStimI);
        if nFrames == 29;
            3;
        end
        if sameNumAllStims
            assert(N3 <= nFrames);
        end
        N3_here = min(N3, nFrames); % for drifting gratings when not always the same number of presentations of each phase (some end a frame or two early)
        idxs = cell(1,N3);
        idxs(:) = {true(1,nFrames)};
        
        if separateCphOrderTrials  % even if separating every trial, need to do cph in this special way.
            stim_i_cph = mod(stim_i-1 + nStimIds/2, nStimIds)+1;
            frameIndsForStimI_cph = stimIdsIdx{stim_i_cph};
            assert(all(abs(frameIndsForStimI_cph-frameIndsForStimI)==1));            
            idx_1st = frameIndsForStimI < frameIndsForStimI_cph;
            for j = 1:N3/2,     idxs{j} = idxs{j} &  idx_1st;    end
            for j = N3/2+1:N3,  idxs{j} = idxs{j} & ~idx_1st;    end
        end
        
        if separateOddEvenTrials || (separateEveryTrial && separateCphOrderTrials)  %
            if separateCphOrderTrials && xor(idx_1st(1), idx_1st(2))
                idx_odd  = 1:nFrames <= nFrames/2;    % [1 1 0 0]
            else
                idx_odd  = odd(1:nFrames); % [1 0 1 0 1 0 ...]
            end
            for j = 1:2:N3,    idxs{j} = idxs{j} &  idx_odd;    end            
            for j = 2:2:N3,    idxs{j} = idxs{j} & ~idx_odd;   end            
        
        elseif separateEveryTrial

            idxs(:) = {false(1,N3)};
            
            if (flashedOrDrifting(Gid) == 1) || (N3 == nFrames)
                idx_stimPres = 1:N3;                
            else  
                [nOri, nSp, nPh, nCyc, nRep] = dbGetUniqueOriSpPh('Gid', Gid, 'length');  % eg for Gid = 743
                isStim = false(nCyc, nRep);
                [presIds, counts] = uniqueCount(whichPres(frameIndsForStimI));  % not 100% perfect, but adequate.
                for pres_i = 1:length(presIds)
                    isStim(1:counts(pres_i),pres_i) = true;
                end
                idx_stimPres = find(isStim(:));                                
            end
            
            for j = 1:length(idx_stimPres),      
                idxs{idx_stimPres(j)}(j) = true;   
            end      

            
        end
        
        for i3 = 1:N3
            if (nnz(idxs{i3}) == 0), continue; end
            assert(nnz(idxs{i3}) == nFrames/N3_here);
            [PSTH_bins, vals] = calcPSTH( spkTsRelToFrame_ms(frameIndsForStimI( idxs{i3} )), psthWindow_ms, nnz( idxs{i3} ), [], nBinsPerExtFrame);
            if any(isnan(vals))
                beep;
                keyboard;
            end
            stimPSTH_vals(:,stim_i, i3) = vals;
        end
            
        
    end

    bin_idx =  ibetween(PSTH_bins, psthWindow_ms);   
    assert(all(bin_idx));
    if ~all(bin_idx)
        PSTH_bins = PSTH_bins(bin_idx);
        stimPSTH_vals = stimPSTH_vals(bin_idx,:,:);  
    end
    
    if (flashedOrDrifting(Gid) == 2) && (diff(psthWindow_ms) < 10) && joinBinsForDriftingGrating % drifting grating 
        PSTH_bins = mean(PSTH_bins);
        stimPSTH_vals = nanmean(stimPSTH_vals, 1);        
    end
    %%
    
    doReconstruction = 0;
    if doReconstruction
        %%
        nFramesTot = length(spkTsRelToFrame_ms);
        frames_reconstructed = zeros(1, nFramesTot);
        wind = [2, frameLength_ms];
%         wind = [-100, 100];
        windIdxs = binarySearch(PSTH_bins, wind);  
        idx0 = binarySearch(PSTH_bins, 4.16666/2);  
        binIdxs = windIdxs(1) : windIdxs(2);
        nBins = length(binIdxs);
        binIdxs_orig = binIdxs - idx0+1;
        

        for stim_i = 1:length(uStimIds)
            %%
            idxs(:) = {true(1,nFrames)};
            frameIndsForStimI = stimIdsIdx{stim_i};
            stim_i_cph = mod(stim_i-1 + nStimIds/2, nStimIds)+1;
            frameIndsForStimI_cph = stimIdsIdx{stim_i_cph};
            assert(all(abs(frameIndsForStimI_cph-frameIndsForStimI)==1));            
            idx_1st = (frameIndsForStimI < frameIndsForStimI_cph);
%             idx_2nd = (frameIndsForStimI > frameIndsForStimI_cph);
            
%             idx_orig = [idx_1st(1), idx_2nd(1), idx_1st(2),  idx_2nd(2)];
            for j = 1:N3/2,     idxs{j} = idxs{j} &  idx_1st;    end
            for j = N3/2+1:N3,  idxs{j} = idxs{j} & ~idx_1st;    end
                        
            idx_odd  = 1:nFrames <= nFrames/2;    % [1 1 0 0]
            for j = 1:2:N3,    idxs{j} = idxs{j} &  idx_odd;    end            
            for j = 2:2:N3,    idxs{j} = idxs{j} & ~idx_odd;   end            

            idx_orig = cellfun(@find, idxs);

%             idx_1 = find(frameIndsForStimI < frameIndsForStimI_cph);
%             idx_2 = find(frameIndsForStimI > frameIndsForStimI_cph);
            
            
%             [~, idx_orig] = sort(idx_orig_sorted);
            %%
            
            for trial_j = 1:nTrials
                orig_idxs = frameIndsForStimI(idx_orig(trial_j)) + [binIdxs_orig]-1;
                vals = stimPSTH_vals(binIdxs, stim_i, trial_j)';

                idx_use = ibetween( orig_idxs, 1, nFramesTot);
                frames_reconstructed(orig_idxs(idx_use)) = frames_reconstructed(orig_idxs(idx_use)) + vals(idx_use);

            end
        end
        
        frames_reconstructed = frames_reconstructed/mean(frames_reconstructed)*meanRate;

        
        spkCount = cellfun(@length, spkTsRelToFrame_ms);
        spkRate = spkCount / mean(spkCount) * meanRate;
        sig = 50; g = gaussian(-4*sig:4*sig, 0, sig);
        spkRate_sm = conv(spkRate, g, 'same');
        figure(1); clf;
        plot(spkRate_sm);
        hold on;
        frames_reconstructed_sm = conv(frames_reconstructed, g, 'same');
        plot(frames_reconstructed_sm, 'r-');
    end
end


%         if separateCphOrderTrials && separateOddEvenTrials
%             [PSTH_bins, stimPSTH_vals(:,stim_i, 1)] = calcPSTH( spkTsRelToFrame_ms(frameIndsForStimI( idx_1st &  idx_odd)), psthWindow_ms, nnz( idx_1st &  idx_odd), [], nBinsPerExtFrame);
%             [PSTH_bins, stimPSTH_vals(:,stim_i, 2)] = calcPSTH( spkTsRelToFrame_ms(frameIndsForStimI( idx_1st & ~idx_odd)), psthWindow_ms, nnz( idx_1st & ~idx_odd), [], nBinsPerExtFrame);
%             [PSTH_bins, stimPSTH_vals(:,stim_i, 3)] = calcPSTH( spkTsRelToFrame_ms(frameIndsForStimI(~idx_1st &  idx_odd)), psthWindow_ms, nnz(~idx_1st &  idx_odd), [], nBinsPerExtFrame);
%             [PSTH_bins, stimPSTH_vals(:,stim_i, 4)] = calcPSTH( spkTsRelToFrame_ms(frameIndsForStimI(~idx_1st & ~idx_odd)), psthWindow_ms, nnz(~idx_1st & ~idx_odd), [], nBinsPerExtFrame);            
%             
%         elseif separateCphOrderTrials
%             [PSTH_bins, stimPSTH_vals(:,stim_i, 1)] = calcPSTH( spkTsRelToFrame_ms(frameIndsForStimI( idx_1st)), psthWindow_ms, nnz( idx_1st), [], nBinsPerExtFrame);
%             [PSTH_bins, stimPSTH_vals(:,stim_i, 2)] = calcPSTH( spkTsRelToFrame_ms(frameIndsForStimI(~idx_1st)), psthWindow_ms, nnz(~idx_1st), [], nBinsPerExtFrame);
%             
%         elseif separateOddEvenTrials
%             [PSTH_bins, stimPSTH_vals(:,stim_i, 1)] = calcPSTH( spkTsRelToFrame_ms(frameIndsForStimI( idx_odd)), psthWindow_ms, nnz( idx_odd), [], nBinsPerExtFrame);
%             [PSTH_bins, stimPSTH_vals(:,stim_i, 2)] = calcPSTH( spkTsRelToFrame_ms(frameIndsForStimI(~idx_odd)), psthWindow_ms, nnz(~idx_odd), [], nBinsPerExtFrame);
%             
%         else
%             [PSTH_bins, stimPSTH_vals(:,stim_i, 1)] = calcPSTH( spkTsRelToFrame_ms(frameIndsForStimI),           psthWindow_ms, nFrames,       [], nBinsPerExtFrame);
%                                             
%         end


%             elseif strcmp(stimType, '36x10x8(2x8)')
%                 idx_odd  = odd( ceil((1:nFrames)/2) ); % [1 1 0 0 1 1 0 0 ...]
