function [bestOffsets offsetSpikeTimes] = getBestSpikeOffsets(spikeTimes, T, resolution, flag)
    
    offsets = -(T/2):resolution:(T/2);
    n = ceil(spikeTimes(end)/T);  % n = number of presentations
    refSpikes = [];
    weights =   [];
%     weightDec = 0.5;

    if ~exist('flag', 'var') || isempty(flag)
       3; 
    end
    
    
    function s = getOffsetForPres(js)
        
        dists = nan(size(offsets));
        presJsLims = [T*(js(1)-1), T*js(end)];
        
        j_spikes0 = elementsInRange(spikeTimes, presJsLims);
        

        for k = 1:length(offsets)
            j_spikes = sort( mod(j_spikes0 + offsets(k), T) );
%             if exist('flag', 'var')
%                 dists(k) = spikeShiftDistance(refSpikes, j_spikes, T, [], 19, [0 T]);
%             else
            dists(k) = spikeShiftDistance(refSpikes, j_spikes, T, []);
%             end
        end
        idx = indmin(dists); 
        s = offsets(idx);
        
        
            if exist('flag', 'var')
                figure(18);     clf;           
                line([1; 1] * refSpikes(:)', (1 + [-0.5; 0.5]) * ones(1,length(refSpikes)) , 'color', 'b');
                xlim([0 T]);
                figure(19); clf;
                spikeShiftDistance(refSpikes, j_spikes0+s, T, weights, 19, [0 T]);
                figure(20)
                plot(offsets, dists, '.-');
                hold on; plot(offsets(idx), dists(idx), 'ro'); hold off
                3;      
            end
        

%         newRefSpikes = [refSpikes; mod(j_spikes0+s, T)];
%         newWeights   = [weights;   min(weights)*weightDec*ones(size(j_spikes0)) ];
% 
%         [refSpikes, weights] = elements(  sortrows([newRefSpikes, newWeights])  );
        
    end

    presLims = T*[(0:n-1); (1:n)]';
    spikeTimesC = elementsInRange(spikeTimes, presLims);
    nSpikesEachPres = cellfun(@length, spikeTimesC);
    pres_indices = ord(nSpikesEachPres, 'descend');

    % generate the 'reference' presentation to compare all the other ones to
    if median(nSpikesEachPres) == 0
        refPres_idx = 1;                             % very few spikes --> pres with max spike #
    else
        refPres_idx = round(length(pres_indices)/2); % ample spikes --> pres with median spike #
    end
    refPresId = pres_indices(refPres_idx);
    numRefPress = 3;
    refPresLims = [T*(refPresId-1), T*(refPresId-1+numRefPress)];
    refSpikes   = sort(mod( elementsInRange(spikeTimes, refPresLims, 'rel'), T));
    
    weights =   ones(size(refSpikes));
    bestOffsets = zeros(1,n);
    
    nSpikesPerGroup = length(refSpikes);
    numGroups = ceil( sum(nSpikesEachPres)/nSpikesPerGroup);
    groupRanges = [ nSpikesPerGroup*(0:numGroups-1)+1; nSpikesPerGroup*(1:numGroups) ]';
    groupIndices = elementsInRange(cumsum(nSpikesEachPres),  groupRanges, 'index');

    
    for gi = 1:length(groupIndices)
        inds = groupIndices{gi};
        if ~isempty(inds)
            bestOffsets(inds) = getOffsetForPres(inds);
        end
        
    end
    

    function offsetSpikeTimes = getOffsetSpikeTimes(spkTs, offsets)
        allPresLims = [T*(0:n-1); T*(1:n)]';
        spikeTimesC = elementsInRange(spkTs', allPresLims);
        offsetSpikeTimesC = cellfun(@plus, spikeTimesC, num2cell(offsets(:)), 'UniformOutput', false);
        offsetSpikeTimes = sort(  [ offsetSpikeTimesC{:}]  )';
    end
    
    offsetSpikeTimes_temp = getOffsetSpikeTimes(spikeTimes, bestOffsets);

    % realign so that is closer to original alignment:
    originalMedian = median(mod(spikeTimes,T));
    newMedian = median(mod(offsetSpikeTimes_temp,T));
    bestOffsets = bestOffsets - (originalMedian-newMedian);
    
    offsetSpikeTimes = getOffsetSpikeTimes(spikeTimes, bestOffsets);
    

    
    %     whichPres = ceil(spikeTimes / T);

%     spikeTimesWithOffsets = crossOp(spikeTimes(:), @plus, offsets);
%     whichPresWithOffsets =  ceil( spikeTimesWithOffsets / T);
    
        
%         spikeTimesC = cell(1,n);
%     shifts = zeros(1, n);
% 
%     for p_i = 1:n
%         spikeTimesC{p_i} = mod( spikeTimes(whichPres == p_i), T);
%     end        
%     
%     shiftedSpikeTimes = elementsInRange(spikeTimes, 
%     
%     
%     
%     spikeShiftDistance(
% 
%     elementsInRange
    

end