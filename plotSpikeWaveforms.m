function plotSpikeWaveforms(t, wvfm, t_window, plotCol)
    
    if exist('t_window', 'var') && ~isempty(t_window)
        t_window = [t(1), t(end)];
        t_idx = t >= t_window(1) & t <= t_window(2);

        t = t(t_idx);
        wvfm = wvfm(t_idx,:,:);
    end

    
    ax = gca;
    [nT, nChannels, nSpk] = size(wvfm);
    nToPlot = 10000;
%     nToPlot = nSpk;

    nToPlot = min(nToPlot, nSpk);

    wvfm_cat = reshape(wvfm, [nChannels*nT, nSpk]);
    
    t_offset = [t(end)-t(1)]*[0:nChannels-1];
    t_cat = bsxfun(@plus, t(:), t_offset);
    t_cat = t_cat(:);
    
    t_lines = t(1) + t_offset;
    if nargin < 4
        plotCol = '';
    end        
    plot(t_cat(:), wvfm_cat(:, 1:nToPlot)', [plotCol '.-'] );

    if nChannels > 1
        drawVerticalLine( t_lines, 'linestyle', ':')
    
        spc = 0.2;
        t1 = roundToNearest(t(1), spc, 'up');    
        tck1 = t1:spc:t(end);
        m1 = tck1(1)/spc;
        m2 = (t(end)-tck1(end))/spc;
        th1 = .01;
        th2 = .5;
        if (abs(m1 - round(m1)) < th1), tck1(1) = []; end
        if (m2 < th2), tck1(end) = []; end    
    
        showTicksOnAllChannels = true;
        if showTicksOnAllChannels 
            tck_ext = bsxfun(@plus, tck1(:), t_offset);
            tickLabels = cellfun(@num2str, num2cell(tck1), 'un', 0);
            tickLabels = repmat(tickLabels, 1, 4);
            set(gca, 'xtick', tck_ext(:), 'xtickLabel', tickLabels);
        else
            set(gca, 'xtick', tck1(:));
        end        
    end
    
    xlim([t_cat(1), t_cat(end)]);    
    3;
    
%     tck1 = unique(roundT find(t

    
    
end


