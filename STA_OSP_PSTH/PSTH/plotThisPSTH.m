function varargout = plotThisPSTH( arg1, binVals, T, bckgRate, timeWindow_ms, noLabelsArg)
    if (nargin == 1) && isstruct( arg1 ) % allow for input as PSTH struct.
%         [PSTH_bins, PSTH_vals, frameLen_ms, meanRate, bckgRate, timeWindow_ms] = deal(... 
%             S.bins, S.vals, S.frameLength_ms, S.meanRate, S.bckgRate, S.timeWindow_ms, S.windowProfile);
        S = arg1; 
        [binCenters, binVals, bckgRate, timeWindow_ms] = deal(...
            S.bins, S.vals, S.bckgRate, S.timeWindow_ms);
    else
        binCenters = arg1;
    end
     
    noLabels = exist('noLabelsArg', 'var') && ~isempty(noLabelsArg);
    noLabels = true;
    blueOutline = true;
    doBckgRate = false;
    doSmooth = 2;
    doStat     = exist('S', 'var') && isfield(S, 'stats');
    
    sm = @(x) mean(reshape(x(:), 2, length(x)/2),1)';
    if doSmooth > 0
        binCenters = sm(binCenters);        
        binVals    = sm(binVals);
    end
%     h = bar(binCenters, binVals, 1);
    h = stairs(binCenters, binVals);
    if ~exist('T', 'var') || isempty(T)
        T = binCenters(end) + diff(binCenters(1:2))/2;
    end
    xlim([0 T]);    
    ylims = ylim;
    ylim([0 ylims(2)]);
    
%     mean_N = mean(binVals);
%     med_N  = median(binVals);
%     drawHorizontalLine(mean_N, 'LineStyle', ':', 'Color', 'g');
%     drawHorizontalLine(med_N,  'LineStyle', ':', 'Color', 'r');
    if doBckgRate && exist('bckgRate', 'var') && ~isempty(bckgRate)
        drawHorizontalLine(bckgRate(1), 'LineStyle', '-', 'Color', 'c');        
        if length(bckgRate) == 2
            drawHorizontalLine([bckgRate(1) + bckgRate(2)*[-1, 1]], 'LineStyle', '--', 'Color', 'c');
        end            
    end
    if exist('timeWindow_ms', 'var') && ~isempty(timeWindow_ms)
        drawVerticalLine(timeWindow_ms, 'Color', 'r');
    end
    if nargout == 1
        varargout = {h};
    end
    
    if ~noLabels
        xlabel('ms after stimulus onset');
        ylabel('spikes / sec');
        title('PSTH of one cell averaged across stimuli');
    end
    
end