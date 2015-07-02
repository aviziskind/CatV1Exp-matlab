function [frameStimIds, uOri, uSp, uPh, stim_RepeatId, cphOrderId] = getStimulusFrameSequence(Gid, OSPmode, cphOrigOrder_flag)
    
    if nargin < 2
        OSPmode = 'OSP';
    end
    
    [allOri_deg, allSp_pix, allPhase_deg] = getOriSpPhaseForEachStimFrame(Gid);    
    allOri_deg = allOri_deg(:);    
    
%     [uOri, ori_idxs] = uniqueList(allOri_deg(:));
    [uOri, ~, ori_idx] = unique(allOri_deg);
    [uSp, ~, sp_idx]   = unique(allSp_pix);
    [uPh, ~, ph_idx]   = unique(allPhase_deg);
    
    nStimTot = numel(allOri_deg); 
    nUStim = length(uOri)*length(uSp)*length(uPh);
    nTrials = nStimTot/nUStim;
    
    %%
    
    addTrialIdxToStimId = strcmpi(OSPmode(end), 't');
    if addTrialIdxToStimId
        OSPmode(end) = '';
    end    
    
        
%     n = nnz( arrayfun(@(s) any(strcmpi(s, {'o', 's', 'p'})), OSPmode));
    idx = cell(1, length(OSPmode));
    L = zeros(1, length(OSPmode));
    for i = 1:length(OSPmode)
        switch OSPmode(i)
            case 'O', idx{i} = ori_idx; L(i) = length(uOri);
            case 'S', idx{i} = sp_idx;  L(i) = length(uSp);
            case 'P', idx{i} = ph_idx;  L(i) = length(uPh); 
            otherwise, error('Error : must be combination of O,S,P');
        end
    end
    
    switch length(OSPmode)
        case 1, frameStimIds = idx{1};
        case 2, frameStimIds = idx{1} + (idx{2}-1)*L(1);
        case 3, frameStimIds = idx{1} + (idx{2}-1)*L(1) + (idx{3}-1)*L(1)*L(2);
    end

    
    if addTrialIdxToStimId  || (nargout >= 5)
            %%   
        trial_idx = zeros(nStimTot, 1);
            
        cphOrigOrder = exist('cphOrigOrder_flag', 'var') && isequal(cphOrigOrder_flag, 1);
            
        stimType = getGratingStimType(Gid);
        cphGrating = stimType.isCphFlashed;
        if cphGrating && cphOrigOrder
            % Special technique to recover original counterphase flashed
            % grating stimulus order
            frameStimIds = frameStimIds(:)';
            [uStimIds, stimIdsIdx] = uniqueList(frameStimIds); 
            nStimIds = length(uStimIds);
            stimIdsIdx = cellfun(@(x) x(:)', stimIdsIdx, 'un', 0);
            idxs = cell(1,nTrials);
            N3 = nTrials;
            
            for stim_i = 1:nStimIds
                idxs(:) = {true(1,nTrials)};
                frameIndsForStimI = stimIdsIdx{stim_i};
                stim_i_cph = mod(stim_i-1 + nStimIds/2, nStimIds)+1;
                frameIndsForStimI_cph = stimIdsIdx{stim_i_cph};
                assert(all(abs(frameIndsForStimI_cph-frameIndsForStimI)==1));            
                idx_1st = (frameIndsForStimI < frameIndsForStimI_cph);
    %             idx_2nd = (frameIndsForStimI > frameIndsForStimI_cph);

    %             idx_orig = [idx_1st(1), idx_2nd(1), idx_1st(2),  idx_2nd(2)];
                for j = 1:N3/2,     idxs{j} = idxs{j} &  idx_1st;    end
                for j = N3/2+1:N3,  idxs{j} = idxs{j} & ~idx_1st;    end

                idx_odd  = 1:nTrials <= nTrials/2;    % [1 1 0 0]
                for j = 1:2:N3,    idxs{j} = idxs{j} &  idx_odd;    end            
                for j = 2:2:N3,    idxs{j} = idxs{j} & ~idx_odd;   end            

                idx_orig = cellfun(@find, idxs);
                trial_idx(frameIndsForStimI) = idx_orig;
            end
            
        else
            % For all other stimuli (drifting gratings, non-counter-phase flashed gratings)
            % use standard technique to get stimulus repeat index.
           
            stim_count = zeros(nUStim, 1);
           
            for i = 1:nStimTot
                stim_count(frameStimIds(i)) = stim_count(frameStimIds(i))+1;
                trial_idx(i) = stim_count(frameStimIds(i));        
            end
            
        end
                
        if addTrialIdxToStimId
            frameStimIds = frameStimIds + (trial_idx-1)*nUStim;
        end

        if nargout >= 5
            stim_RepeatId = trial_idx;
        end
        
    end
    
        
%%        
            
    
    
            
    
    if nargout >= 6        
        %%% note: this calculation is incorrect!  (?)
        error('! (?)')

        stimType = getGratingStimType(Gid);
        cphGrating = stimType.isCphFlashed;
        if ~cphGrating
            cphOrderId = [];

        else
            cphOrderId = zeros(1, nStimTot);
            cphFlag = 1;

            for frm_i = 1:nStimTot                        
                cphOrderId(frm_i) = cphFlag; 
                cphFlag = cycle(cphFlag, [1, 2]);                        
            end                           
        end
    end

%     if (OSPmode(end) == 't')
%     elseif (OSPmode(end) == 'T')
        
        
    
end


%     if checkForMissingFrame
%         [anyMissingFrames, actualFrameSequence] = getMissedFrames('Gid', Gid);    
%         if anyMissingFrames
%             frameStimIds = frameStimIds(actualFrameSequence);
%         end
%     end


%     nTotalFrames = length(allOri_deg);
%     allStims = zeros(length(allOri_deg), length(OSPmode));
%     for i = 1:length(OSPmode)
%         switch OSPmode(i)
%             case 'O', allStims(:,i) = allOri_deg(:);
%             case 'S', allStims(:,i) = allSp_pix(:);
%             case 'P', allStims(:,i) = allPhase_deg(:);
%             otherwise, error('Error : must be combination of O,S,P');
%         end
%     end
%     % this stimulus list will be sorted by the last column, then by second last...
%     % Thus, the ith column corresponds to the ith dimension of the OSP
%     % array.        
%     uStims = fliplr( unique( fliplr(allStims), 'rows'));   
%     frameStimIds = zeros( nTotalFrames, 1);
%     for stim_i = 1:size(uStims,1);
%          frameStimIds( findRows(uStims(stim_i,:), allStims)  ) = stim_i;
%     end
% 

%{
    tic;    
    t_idx = 1:nTrials;
    trial_idx = zeros(numel(allOri_deg), 1);
    for i = 1:length(stim_idx)
        trial_idx(stim_idx{i}) = t_idx;
    end
    toc;
%}