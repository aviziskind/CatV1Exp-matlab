function [PSSH_bins, spikePSSH_vals, meanRate, bckgRate] = getIndividualPSSHs(Gid, cellId, groupingStyle, psthWindow_ms, bkgrSkip_ms)
    
    if nargin < 3
        groupingStyle = 'OSP';
    end
    separateOddEvenTrials = true;
%     trialEncoding = 'odd/even'; %options: 'odd/even', 'trial id'
%     addTrialNumber = true;
    stimType = getGratingStimType(Gid);
    separateCphOrderTrials = stimType.isCphFlashed;
    
    frameLength_ms = getFrameLength('Gid', Gid);

    if ~exist('psthWindow_ms', 'var') || isempty(psthWindow_ms)
        [nBinsPerFrame, nFramesPerExtFrame, extFrameLength_ms] = getBinSizeParams(frameLength_ms);          
        psthWindow_ms = [0 extFrameLength_ms];        
    else
        if length(psthWindow_ms) == 1
            psthWindow_ms = [0 psthWindow_ms];
        end
        [nBinsPerFrame, nFramesPerExtFrame, extFrameLength_ms] = getBinSizeParams(frameLength_ms, diff(psthWindow_ms));          
        psthWindow_ms(2) = psthWindow_ms(1)+extFrameLength_ms;
    end            
    
    if ~exist('bkgrSkip_ms', 'var') || isempty(bkgrSkip_ms)
        bkgrSkip_ms = 200;
    end        

    % 1. Get Position of stimulus times reslative to the spikes
    [spkFrameTimes, frameIds] = getParsedSpikes('spikes', Gid, cellId, psthWindow_ms, [], bkgrSkip_ms);     
    
    frameStimData = cell(size(spkFrameTimes));
    
    % 2. Get ids of the actual stimuli # that were presented (instead of just frame #) 
    [frameStimIds, uOri, uSp, uPh, frameRepId, cphOrderId] = getStimulusFrameSequence(Gid, groupingStyle);
    [uStimIds, stimIdsIdx] = uniqueList(frameStimIds);

    if ~separateOddEvenTrials && separateCphOrderTrials
        error('not implemented');
    end
    
    oddEvenId = zeros(size(frameStimIds));    
    if separateCphOrderTrials        
        nStims = length(uStimIds);
        trial_ids = zeros(1,4);
        for i = 1:nStims   % separate even/odd & cph order simultaneously
            frmIds = stimIdsIdx{i};
            stimCphOrderIds = cphOrderId(frmIds);            
            idx_first = (stimCphOrderIds==1);
            trial_ids(idx_first) = [1 2];
            trial_ids(~idx_first) = [3 4];
            oddEvenId(frmIds) = trial_ids;  
        end
    else  % separate even/odd trials
        oddEvenId = zeros(size(frameStimIds));
        nStims = length(uStimIds);
        for i = 1:nStims
            frmIds = stimIdsIdx{i};            
            idx_odd = odd( frameRepId(frmIds) );
            oddEvenId(frmIds(idx_odd)) = 1;
            oddEvenId(frmIds(~idx_odd)) = 2;
        end        
    end
        
    
    for spk_i = 1:length(spkFrameTimes)       
        frmIds = frameIds{spk_i};
        stimIds = zeros(1,length(frmIds));
        stimOddEvenIds = zeros(1,length(frmIds));        
        
        idx = frmIds > 0;
        stimIds(idx) = frameStimIds(frmIds(idx));
        stimOddEvenIds(idx) = oddEvenId(frmIds(idx));
                
        frameStimData{spk_i} = [stimIds; stimOddEvenIds];
    end
    
    nStimuli = length(uOri)*length(uSp)*length(uPh);
    % 3. "expand" (ie. fill in the gaps: replace frame times with continuum of stimulus Ids.)  
    
    binEdges = linspace(psthWindow_ms(1), psthWindow_ms(2), nBinsPerFrame*nFramesPerExtFrame+1);
    expandedFrameTimes = expandStimulusFrameTimes(spkFrameTimes, frameStimData, binEdges, nStimuli);    
%     expandedFrameTimes = expandPSSH(frameTimes, frameStimData);
    
    N3 = (separateOddEvenTrials+1)*(separateCphOrderTrials+1);
    
    % 4. Histogram over stimuli.
    PSSH_bins = binEdge2cent(binEdges);
    spikePSSH_vals = zeros( length(PSSH_bins) , nStimuli,  N3);
    nBins = length(PSSH_bins);
    
    for i = 1:nBins
        histdata = intHist(expandedFrameTimes(:,i), [1 nStimuli*N3]);            
        spikePSSH_vals(i,:,:) = reshape( histdata, [nStimuli, 1,N3]);
    end    

    [meanRate, bckgRate] = deal([]);
end


  
function n = intHist(data, valRange)        
    if length(valRange) == 1
        valRange = [1 valRange];
    end
    offset = valRange(1)-1;
    numVals = valRange(2)-offset;
    n = zeros(1,numVals);
    
    [uVals, uValCount] = uniqueCount(data);    
    
    idx = find( ibetween(uVals(:)', valRange(1), valRange(2) ) );
    n(uVals(idx)) = uValCount(idx);
%     for i = find()
%         n(idx = uVals(i)-offset;
%         n(id)=n(id)+uValCount(i);
%     end        
end

    
    
