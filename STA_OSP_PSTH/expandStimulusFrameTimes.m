function frameTimes_exp = expandStimulusFrameTimes(frameTimes, stimulusData, binEdges, nStimuli)
    nSpikes = length(frameTimes);
    nBins = length(binEdges)-1;
    
    separateOddEvenTrials = (nargin >= 3);    
%     if isnumeric(
%     
    dbug = false;
    frameTimes_exp = zeros(nSpikes, nBins, 'uint16');
%     
%     stimType = getGratingStimType(Gid);
%     dims = parseStimType(stimType);
%     nStimuli = prod(dims);
    
    for spk_i = 1:nSpikes;
        spkFrameTimes = double(frameTimes{spk_i});
        spkStimIds    = stimulusData{spk_i}(1,:);
        spkStimOddEvenIds = stimulusData{spk_i}(2,:); 
            
%          stimIdx2 = binarySearch(binEdges, spkTimes, -1, -1);
        stimIdx = binarySearch(spkFrameTimes, binEdges, 1, 1);
            
        if dbug
            figure(1012); clf;
            line([1;1]*spkFrameTimes, [0;1] * ones(1, length(spkFrameTimes)), 'color', 'k');
            line([1;1]*binEdges, [1;2] * ones(1, length(binEdges)), 'color', 'r');
            for i = 1:length(spkFrameTimes)-1
                text(spkFrameTimes(i), .5, {num2str(i), num2str(spkStimIds(i))}, 'horiz', 'left')
            end
            for i = 1:nBins
                text(binEdges(i), 1.5, {num2str(i), num2str(stimIdx(i))}, 'horiz', 'left')
            end
            changed = false(1, nBins);
        end
        
        idx_change = diff(stimIdx)'==1;
        for i = find(idx_change)            
            frmT = spkFrameTimes(stimIdx(i+1));            
            if (frmT-binEdges(i)) < (binEdges(i+1)-frmT)
                stimIdx(i) = stimIdx(i) + 1;
%                 if dbug, changed(i) = 1; end
            end                        
        end
        
        if dbug
            for i = find(changed)
                text(binEdges(i), 1.2, num2str(stimIdx(i)), 'horiz', 'left', 'color', 'r');            
            end
        end
    
        stimIds    = spkStimIds(stimIdx(1:nBins));
        stimOddEvenIds = spkStimOddEvenIds(stimIdx(1:nBins));
        
%         stimRepIds_offset = zeros(1, nBins);
        
        if separateOddEvenTrials
            stimRepIds_offset = (stimOddEvenIds-1) .* nStimuli .* (stimIds > 0);            
            stimIds = stimIds + stimRepIds_offset;                    
        end    
        
        frameTimes_exp(spk_i,:) = stimIds;
    end

    
    

end