function r = decorrelatePSTHs(R, frameStimIds, nBinsPerFrame, nFramesPerExtFrame, framesToConsider, nStimsToSkip, eigThreshold)
    % selects the top N PSTHs(the ones with the highest response across the
    % extended frame), and uses the decorrelation procedure only on these
    % frames. (if nStimsToSkip < 1, nStimsToSkip is taken as the fraction of bins to be skipped).
    % 
    % note: size(R) === [nBinsPerFrame, nStims*nFramesPerExtFrame];

    sortBinsAccToRate = true;
    subtractMean = true;
    
    if ~exist('framesToConsider', 'var') || isempty(framesToConsider)
        framesToConsider = 1:nFramesPerExtFrame;
    end
    if ~exist('nStimsToSkip', 'var')
        nStimsToSkip = 0;
    end
    if ~exist('eigThreshold', 'var') || isempty(eigThreshold);
        eigThreshold = 0;
    end
    nStimuli = length(unique(frameStimIds));
    
    
    if (nStimsToSkip == 0)
        A = generateMtxForPSTHdecorrelation(frameStimIds, nFramesPerExtFrame, framesToConsider);
        r = decorrelatePSTHs_helper(R, nBinsPerFrame, framesToConsider, A);
        
%          r = generateMtxAndDecorrelatePSTHs(R, frameStimIds, nBinsPerFrame, nFramesPerExtFrame, framesToConsider);
%     r = generateMtxForPSTHdecorrelation(frameStimIds, nFramesPerExtFrame, nBinsPerFrame, framesToConsider);
%         r = decorrelatePSTHs_helper(R, , framesToConsider, A);
        
        if subtractMean
%             r = r - repmat( mean(r, 2), 1, nStimuli);
        end
        return;
    end

    
    if nStimsToSkip < 1
        nStimsToSkip = round( nStimsToSkip * nStimuli );
    end
    nTopStims = nStimuli - nStimsToSkip;
    
    % Find order 
    if sortBinsAccToRate
        stimIndsInOrder = ord( sum(R,1), 'descend' );
    else
        stimIndsInOrder = 1:nStimuli;
    end

    % Select top N stimuli (with highest responses)
    indsToIgnore = stimIndsInOrder( nTopStims + 1 : end );
    frameTopStimIds = frameStimIds;
    for ind = indsToIgnore(:)'
        frameTopStimIds(frameTopStimIds == ind) = 0;
    end

    % Renumber stimuli from 1 to nTopStims
    frameStimIdsRenumbered = zeros(size(frameStimIds));  
    for i = 1:nTopStims
        frameStimIdsRenumbered( frameStimIds == stimIndsInOrder(i) ) = i;
    end

    % Put PSTHs of top responses into R :
    R_top = R(:,stimIndsInOrder(1:nTopStims));

    % DECORRELATE PSTHS
    A_top = generateMtxForPSTHdecorrelation(frameStimIdsRenumbered, nFramesPerExtFrame, framesToConsider);
    r_top = decorrelatePSTHs_helper(R_top, nBinsPerFrame, framesToConsider, A_top);
    
    if subtractMean
        r_top = r_top - repmat( mean(r_top, 2), 1, nStimuli-nStimsToSkip);
    end

    % Put into r in the correct order.
    r = zeros(size(R));    
    r(:, stimIndsInOrder(1:nTopStims)) = r_top;

    % put PSTHs that were excluded into r
    minBaseline = sum( r(:, stimIndsInOrder(nTopStims) )) / (nBinsPerFrame*nFramesPerExtFrame);
    r(:, stimIndsInOrder(nTopStims+1:end)) = minBaseline;
    
end
