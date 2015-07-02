function STA_jacks = getJackknifedSTA(r_full_orig, Gid, trial_mode, nJackSegments, jackknifeMethod)
    [nOri, nSpf, nPh, nTrials_orig] = size(r_full_orig);
%     nStim = nOri*nSpf*nPh;

%     STA_scaling = 'none';
    STA_scaling = 'aroundZero';
    
    switch trial_mode
        case 'all',  trial_idxs_use = 1:nTrials_orig;
        case 'odd',  trial_idxs_use = 1:2:nTrials_orig;
        case 'even', trial_idxs_use = 2:2:nTrials_orig;
    end
    
    r_full = r_full_orig(:,:,:, trial_idxs_use);
    nTrials = length(trial_idxs_use);
        
    assert(any(strcmp(jackknifeMethod, {'stimuli', 'trials', 'random_stimuli_and_trials', 'alternating_stimuli_and_trials'})));
    if strcmp(jackknifeMethod, 'trials') && nTrials < nJackSegments % in some cases, have only 4 trials. Then, odd/even only have 2 trials, which is < 4 jackknife segments.
        jackknifeMethod = 'alternating_stimuli_and_trials';
    end
%     
%     if nargin < 5
%         STA_jacks = [];
%         return;
%     end

    
    if strcmp(jackknifeMethod, 'stimuli')   % first average over trials, then jackknife over stimuli (like we do for the MIDs)
                                            % I'm worried this might lead to increased variability, though
        nStims = nOri*nSpf*nPh; % length(r_randOrder);
        r = mean(r_full, 4);
        frameIds = getStimulusFrameSequence(Gid, 'OSP');
        uStimIds = uniqueInOrder(frameIds);

    elseif any(strcmp(jackknifeMethod, {'random_stimuli_and_trials', 'alternating_stimuli_and_trials'} ))  % jackknife over all stimuli & trials

        nStims = nOri*nSpf*nPh*nTrials; % length(r_randOrder);
        r = r_full;

        if strcmp(jackknifeMethod, 'random_stimuli_and_trials')
            if strcmp(trial_mode, 'all')
                frameIds = getStimulusFrameSequence(Gid, 'OSPt');
            elseif any(strcmp(trial_mode, {'odd', 'even'}))
                frameIds = getStimulusFrameSequence(Gid, 'OSP');
                frameIds = reshape(frameIds, [nOri*nSpf*nPh, nTrials_orig]);
                frameIds = frameIds(:,trial_idxs_use);
                frameIds = bsxfun(@plus, frameIds, [0:length(trial_idxs_use)-1]*(nOri*nSpf*nPh) ); % add appropraite offsets
                frameIds = frameIds(:);         
            end
            uStimIds = frameIds;
        elseif strcmp(jackknifeMethod, 'alternating_stimuli_and_trials')
            
            %%
            idx_mark_halfeven_stim = odd( ceil([1:1:nTrials_orig]/2) );
            idx_mark_halfodd_stim = ~idx_mark_halfeven_stim;

            nStim_1trial = nOri*nSpf*nPh;
            stimIds_odd = 1:2:nStim_1trial;
            stimIds_even = 2:2:nStim_1trial;

            idx_alternating = false(nStim_1trial, nTrials);
            idx_alternating(stimIds_even, idx_mark_halfeven_stim) = true;
            idx_alternating(stimIds_odd,  idx_mark_halfodd_stim)  = true;
                %%
            if any(strcmp(trial_mode, {'odd', 'even'}))
                idx_alternating = idx_alternating(:, trial_idxs_use);
                
            end
            
            uStimIds = [find(idx_alternating(:)); find(~idx_alternating(:))];
            
        end
            
            
        
    elseif strcmp(jackknifeMethod, 'trials')  % jackknife only over trials (for some experiments with only 4 trials, don't have enough for even/odd)
        nStims = nTrials;
    end

    jackIdxs = jackknifeIdxs(nStims, nJackSegments);
    
    assert( length( unique( cellfun(@length, jackIdxs))) == 1);


    STA_jacks_C = cell(1, nJackSegments);
    for jack_i = 1:nJackSegments
        if  strcmp(jackknifeMethod, 'trials')
            r_jack_i = r_full(:,:,:,jackIdxs{jack_i});
        else
            r_jack_i = zeros(size(r));
            r_jack_i( uStimIds( jackIdxs{jack_i} ) ) = r( uStimIds( jackIdxs{jack_i} ) );
        end
        STA_jacks_C{jack_i} = getSTAfromOSP(Gid, r_jack_i, STA_scaling);
    end
    STA_jacks = cat(3, STA_jacks_C{:});
    
    
    
    
end