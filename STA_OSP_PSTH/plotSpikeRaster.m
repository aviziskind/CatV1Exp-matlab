function plotSpikeRaster(spikeTimes, T, displayT, col, offsets)

    if iscell(spikeTimes)
        spikeTimesC = spikeTimes;
        n = length(spikeTimesC);
        
    elseif isnumeric(spikeTimes)    % put into cell
        n = ceil(spikeTimes(end)/T);  % n = number of presentations

        presLims = [T*(0:n-1); T*(1:n)]';
        spikeTimesC = elementsInRange(spikeTimes, presLims);
        if exist('offsets', 'var')
%             cumOffsets = offsets(:);
%             presLims = presLims + [cumOffsets, cumOffsets];
            spikeTimesC = cellfun(@plus, spikeTimesC, num2cell(offsets(:)), 'UniformOutput', false);
        end
        spikeTimesC = cellfun(@(x) mod(x,T), spikeTimesC, 'UniformOutput', false);
        
%         spikeTimesC2 = cell(n,1);
%         whichPres = ceil(spikeTimes / T);
%         for p_i = 1:n
%             spikeTimesC2{p_i} = mod( spikeTimes(whichPres == p_i), T);
%         end   
%         chk = cellfun(@minus, spikeTimesC2, spikeTimesC, 'UniformOutput', false);
%         max( vertcat( chk{:} ) );
    end

    if ~exist('displayT', 'var') || isempty(displayT);
        displayT = T;  % sometimes spikeTimes are in normalized units, not absolute units.
    end
    T_ratio = displayT/T;

    if ~exist('col', 'var') || isempty(col);
        col = 'b';  % sometimes spikeTimes are in normalized units, not absolute units.
    end
    
    % 2. Setup axes;
    axis([0 displayT, .5 n+.5])
    axis ij;
    box on;
    hold on;

    % 3. Plot spikes
    for p_i = 1:n
        spikesThisPres = spikeTimesC{p_i};
        if ~isempty(spikesThisPres)
            lineX = (T_ratio * [1; 1]) * spikesThisPres(:)';
            lineY = (p_i + [-0.5; 0.5]) * ones(1,length(spikesThisPres));
            line( lineX, lineY, 'Color', col);
        end
    end

end