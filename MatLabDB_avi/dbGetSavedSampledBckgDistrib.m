function bckgSmp_maxed = dbGetSavedSampledBckgDistrib(Gid, cellId, bckgSamples, nStimTot, nPhase, nTrials, nTopStim, nSamples)

    persistent allBckgSamples calc_count

    bckgDistribFileName = [CatV1Path 'MatLabDB_avi' filesep 'allBckgDistribs' curMatchDB('') '.mat'];        

    redo_all = 0;
    redo_current = 0;
    
    fld_name = @(gid, cid) sprintf('bckgSamples_%d_%d', gid, cid);
        
    if isempty(allBckgSamples) 
        
        if exist(bckgDistribFileName, 'file') && ~redo_all        
            S = load(bckgDistribFileName);
            allBckgSamples = S.allBckgSamples;            
        
        else
            allBckgSamples = struct; 
        end    
        calc_count = 0;
    end
    
    if strcmp(Gid, 'save')
        save(bckgDistribFileName, 'allBckgSamples', '-v6');
        calc_count = 0;        
        return;
    end
    
    if ~exist('bckgSamples', 'var') || isempty(bckgSamples)
%         [~, ~, ~, bckgSamples] = dbGetCellSpkStimHists(Gid, cellId);
        error('must supply bckgSamples');
%         if flashedOrDrifting(Gid) == 1
%             % need to use R_bin - L_bin
%         elseif flashedOrDrifting(Gid) == 2
%             [uori, usp, uph, ~, ~, utf_Hz] = dbGetUniqueOriSpPh('Gid', Gid);
%             bckgBinSize_ms = 4+1/6;
%             nSamplesAv = round((utf_Hz*1000)/bckgBinSize_ms);    
%             bckgSamples = getBackgroundSpikes(Gid, cellId, [], nPhase);
%         end

    end
    
%     if ~exist('nPhase', 'var') || isempty(nPhase) || ~exist('nStimTot', 'var') || isempty(nStimTot)
% %         sd = siteDataFor('Gid', Gid); 
%         gt = getGratingStimType(Gid);
%         nPhase = switchh(gt.gratingType, {'flashed', 'drifting'}, [gt.nPh, 1]);
%         nStimTot = gt.nOri * gt.nSpf;        
%     end
%     if ~exist('nTopStim', 'var') || isempty(nTopStim)
%         nTopStim = 5;
%     end
%     if ~exist('nSamples', 'var') || isempty(nSamples)
%         nSamples = 1000;
%     end
        

    curFldName = fld_name(Gid, cellId);
       
    
    if redo_current || ~isfield(allBckgSamples, curFldName) || ...
            any(compSize(allBckgSamples.(curFldName).samples) < [nTrials, max(nTopStim), nSamples]);
      
        bckgSmp_maxed = sampleBckgDistrib(bckgSamples, nStimTot, nPhase, nTrials, nTopStim, nSamples);        

        allBckgSamples.(curFldName) = struct('samples', compress(bckgSmp_maxed), 'nStimTot', nStimTot, 'nPhase', nPhase, 'nTrials', nTrials);
        calc_count = calc_count + 1;
    end
    
    if calc_count >= 500
        save(bckgDistribFileName, 'allBckgSamples', '-v6');
        calc_count = 0;
    end
    
    bckgSmp_S = allBckgSamples.(curFldName);    
    assert(bckgSmp_S.nStimTot == nStimTot);
    assert(bckgSmp_S.nPhase == nPhase);
    assert(bckgSmp_S.nTrials == nTrials);

    bckgSmp_maxed = decompress(bckgSmp_S.samples);
    [nTrialsInSample, nTopStim_available, nSamples_available] = size(bckgSmp_maxed);    
    
    if (nSamples < nSamples_available) || (nTopStim < nTopStim_available)
        idx_samples_use = 1:min(nSamples, nSamples_available);
        idx_nTopStim_use = 1:min(nTopStim, nTopStim_available);
        
        bckgSmp_maxed = bckgSmp_maxed(:, idx_nTopStim_use, idx_samples_use);
    end
    

end


function dims = compSize(s)
    if isfield(s, 'dims')
        dims = s.dims;
    elseif isfield(s, 'vals_idx')        
        dims = size(s.vals_idx);
    elseif isfield(s, 'orig_vals')
        dims = size(s.orig_vals);
    end

end

function bckgSamples_maxed = sampleBckgDistrib(bckgSamples, nStimTot, nPhase, nTrials, nTopStim, nSamples) % = nPhases, nStimTot, nTopStim, nSamples
    rand('state', 0);    
    
    bckgSamples_maxed = zeros(nTrials, nTopStim, nSamples);
    for smp_i = 1:nSamples
        allVals = reshape(randsample(bckgSamples, nStimTot*nPhase*nTrials, true), [nStimTot, nPhase, nTrials]);
        
        % 1. take mean over "phases" 
        bckgSmp_avPh = mean(allVals, 2); 
        
        % 2. take mean over trials, and find top N best stimuli
        bckgSmp_avPh_avTrials = mean(bckgSmp_avPh, 3); 
        bestStim_idxs = indmax_n(bckgSmp_avPh_avTrials, nTopStim);
        
        % 3. collect all the individual trials of the top N best stimuli
        allTrialsOfTopNStim = reshape(bckgSmp_avPh(bestStim_idxs, :, :), [nTopStim, nTrials])';        
        bckgSamples_maxed(:,:,smp_i) = allTrialsOfTopNStim;        
    end
        
        
end

