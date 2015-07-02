function fig_plotsForPSTHdecorrelation(celldata)

    RandStream.setDefaultStream(RandStream('mt19937ar','seed',6));

    function [hStim, hArr] = drawStimulusWithArrow(mainAxes, stimulus, stimNum, arrowEndsX, arrowColor)
        axes(mainAxes);
        ylims = ylim;

        % draw stimulus
        stimXs = arrowEndsX(1) + diff(arrowEndsX) * [1/4, 3/4];
        stimYs = ylims(1) + diff(ylims)*[.8, .95];

        ll_xy = [stimXs(1), stimYs(1)];
        ur_xy = [stimXs(2), stimYs(2)];

        %         [ll_xy_fig, ur_xy_fig] = dsxy2figxy(mainAxes, ll_xy, ur_xy);
        %         stimBox = corners2box(ll_xy_fig, ur_xy_fig);

        stimBoxDs = corners2box(ll_xy, ur_xy);
        stimBoxPos = dsxy2figxy(mainAxes, stimBoxDs);


        hStim = axes('Position', stimBoxPos);
        imagesc(stimulus);
        colormap('gray')
        axis equal tight;
        set(gca, 'xtick', [], 'ytick', []);
        if (stimNum == 0)
            sub = 'i';
        elseif (stimNum > 0)
            sub = ['i+' num2str(stimNum)];
        elseif (stimNum < 0)
            sub = ['i-' num2str(-stimNum)];
        end

%         xlabel(['$$ s_{' sub '} $$'], 'Rotation', 0, 'interpreter', 'latex');



        arrowY = [1,1] * ylims(1) + diff(ylims)*(.77);
        arrowX = arrowEndsX;
        [arrX, arrY] = dsxy2figxy(mainAxes, arrowX, arrowY);
        hArr = annotation('doublearrow',arrX,arrY, 'Color', arrowColor, 'LineWidth', 2);
    end

    writeBinIdxs = false;

    %     nBinsPerFrame = 3;
    %     nFramesPerExtFrame = 7;
    %     nBinsPerExtFrame = nBinsPerFrame * nFramesPerExtFrame;
    %     [nBinsPerFrame, nFramesPerExtFrame, extFrameLength_ms] = getBinSizeParams(frameLength_ms);
    %       nBinsPerExtFrame = nBinsPerFrame * nFramesPerExtFrame;


    global idealExtFrameLength_ms idealBinSize_ms;    
    idealExtFrameLength_ms = 100;    
    idealBinSize_ms = 1000/(60*3);

    
    useInputPSTH = nargin > 0;
    if useInputPSTH
        PSTH = celldata.PSTH;
        R1 = PSTH.vals;
        frameLength_ms = PSTH.frameLength_ms;
    else        
        frameLength_ms = 30;
    end
%     frameLength_ms = 100;
    frameLength_ms = 1000/60;
    
    [nBinsPerFrame, nFramesPerExtFrame, extFrameLength_ms] = getBinSizeParams(frameLength_ms);
    nBinsPerExtFrame = nBinsPerFrame * nFramesPerExtFrame

    nStimuli = nFramesPerExtFrame;

    Mx = 200;
    % make random PSTHs (for 2nd...nth stimuli)
    r = abs(randn(nBinsPerExtFrame, nStimuli)) * Mx ; %noise
    m = nBinsPerExtFrame/2;
    s = nBinsPerExtFrame/12;
    
    %add some peaks
    for i = 1:nStimuli
        thisPeak = histc( s*randn(1,5000)+m, 1:nBinsPerExtFrame )';
        thisPeak = thisPeak / max(thisPeak) * Mx*rand;
        r(:,i) = r(:,i) + thisPeak;
    end
    
    if useInputPSTH
        r(:,1) = R1(1:nBinsPerExtFrame);
    end
    %         PSTH = celldata.PSTH;
    %         r = PSTH.vals;
    %         r = reshape(r, [nBinsPerExtFrame, nStimuli]);
    

    binStartsForStim  = @(s, nBinsTotal) ([0.5:1:nBinsTotal ]+s*nBinsPerFrame) * (frameLength_ms/nBinsPerFrame);
    binCentersForStim = @(s, nBinsTotal) ([0:  1:nBinsTotal ]+s*nBinsPerFrame) * (frameLength_ms/nBinsPerFrame);

    %     r = r(:,randperm(nStimuli));   % scramble the order

    Gid = 4470;
    gf = getFrameRetrieverFunction(Gid);
    gf('load', Gid);
    frmIds = [3, 5, 7, 10, 11, 15, 18, 21, 84];


    % Figure 1
    figure(1);  clf;
    t_c = binStartsForStim(0, nBinsPerExtFrame);
    t_s = binCentersForStim(0, nBinsPerExtFrame);
    bar( t_c, r(:,1), 1, 'b' );
    aa = -1;
    bb = 1.8;
    if writeBinIdxs
        text(t_c(1)+aa, r(1,1)+bb, 'r_i^{1,1}');
        text(t_c(2)+aa, r(2,1)+bb, 'r_i^{1,2}');
        text(t_c(3)+aa, r(3,1)+bb, 'r_i^{1,3}');
        text(t_c(4)+aa, r(4,1)+bb, 'r_i^{2,1}');
    end

    axes1 = gca;
    set(axes1, 'Color', 'none');
    ylims = ylim;

%     set(axes1, 'xtick', floor([t_s(1):frameLength_ms:t_s(end)]))
%     set(axes1, 'xtick', []);
    set(axes1, 'xtick', roundToDecimalPoint(frameLength_ms*[0:nFramesPerExtFrame], 1));
    xlim(frameLength_ms*[0, nFramesPerExtFrame]);
    ylim([0 200]);
% %     ylim([ylims(1), ylims(2)/.75]);
    xlabel('time (ms)');
    ylabel('spikes/sec');

    drawStimulusWithArrow(axes1, gf(frmIds(1)), 0, [0 frameLength_ms], color_s(1) );
    axes(axes1);
    drawVerticalLine([.2 frameLength_ms], 'Linestyle', ':', 'LineWidth', 3);
    drawVerticalLine([0:nFramesPerExtFrame]*frameLength_ms, 'Linestyle', ':');
    

%     if ~(nStimuli >1)
        return;
%     end

    % Figure 2
    figure(2);  clf;
    t_c = binStartsForStim(-1, nBinsPerExtFrame);
    t_s = binCentersForStim(-1, nBinsPerExtFrame);
    bar( t_c, r(:,2), 1, 'g' );
    if writeBinIdxs
        text(t_c(1)+aa, r(1,2)+bb, 'r_{i-1}^{1,1}');
        text(t_c(2)+aa, r(2,2)+bb, 'r_{i-1}^{1,2}');
        text(t_c(3)+aa, r(3,2)+bb, 'r_{i-1}^{1,3}');
        text(t_c(4)+aa, r(4,2)+bb, 'r_{i-1}^{2,1}');
    end


    axes2 = gca;
    set(axes2, 'Color', 'none');
    xlims = xlim;
    ylims = ylim;

    set(gca, 'xtick', t_s(1):frameLength_ms:t_s(end))
    xlim([t_s(1) t_s(end)]);
    ylim([ylims(1), ylims(2)/.8]);

    xlabel('time (ms)');
    ylabel('spikes/sec');
    matchAxes('Y', axes1, axes2);

    drawStimulusWithArrow(axes2, gf(frmIds(2)), -1, [-frameLength_ms, 0], color_s(2) );
    axes(axes2);
    drawVerticalLine([0 frameLength_ms], 'Linestyle', ':', 'LineWidth', 3);
    drawVerticalLine([-1:nFramesPerExtFrame-1]*frameLength_ms, 'Linestyle', ':');

    % Figure 3  (5 plots)
    figure(3);  clf;

    t_c = binStartsForStim(-nFramesPerExtFrame+1, nBinsPerExtFrame*2 - 1);
    t_s = binCentersForStim(-nFramesPerExtFrame+1, nBinsPerExtFrame*2 - 3);
    R = zeros(nBinsPerExtFrame*2 - 1, nStimuli);
    for st_i = 1:nFramesPerExtFrame
        R([1:nBinsPerExtFrame]+(nStimuli-st_i)*nBinsPerFrame, st_i) =    r(:,st_i);
    end

    hs = bar( t_c, R, 1, 'stacked' );
    for i = 1:nStimuli
        set(hs(i), 'FaceColor', color_s(i))
    end

    mainAxes = gca;
    set(mainAxes, 'Color', 'none');
    xlims = xlim;
    ylims = ylim;

    set(gca, 'xtick', t_s(1):frameLength_ms:t_s(end))
    xlim([t_s(1) t_s(end)]);
    ylim([ylims(1), ylims(2)/.8]);
    drawVerticalLine([0 frameLength_ms], 'Linestyle', ':', 'LineWidth', 2)

    xlabel('time (ms)');
    ylabel('spikes/sec');

    timeWindows = frameLength_ms * [[0:-1:-nFramesPerExtFrame+1]; [1:-1:-nFramesPerExtFrame+2]]';
    for st_i = 1:nFramesPerExtFrame
        drawStimulusWithArrow(mainAxes, gf(frmIds(st_i)), -st_i+1, timeWindows(st_i,:), color_s(st_i) );

    end

end