function [PSTH_bins, PSTH_vals] = calculatePSTHiteratively(Gid, cellId)
    
    extFrameLen_ms = 120;
    timeWindow    = [30 80]; % initial guess
    windowProfile = [1   1]; % initial guess;
    iter = 1;
    maxIterations = 7;

    ccThreshold = 0.98;
    ccBetweenPSTHs = 0;

    currPSTH_bins = [0 extFrameLen_ms/2, extFrameLen_ms];
    currPSTH_vals = [0 0 0];

    if showWorking
        figure(123);
        clf;
        nPlots = 2;  % OSP and PSTH
        nIters2Show = 6;
    end


    while (ccBetweenPSTHs < ccThreshold) && (iter <= maxIterations)

        [OSP, oris, sps, phs] = getOSPFromTimeWindow(Gid, cellId,   timeWindow, windowProfile);
        subplotOSP( OSP, osp_ind );

        [timeWindow, windowProfile,  PSTH_bins, PSTH_vals] = getPSTHwindowFromOSP(Gid, cellId,   OSP, oris, sps, phs,  extFrameLen_ms);
        subplotPSTH( PSTH_bins, PSTH_vals, timeWindow, psth_ind);

        if showWorking
            %             subplot(nPlots, nIters2Show, iter);
            %             imagesc3({oris, sps, phs}, OSP, 3, {'Orientation', 'Sp freq', 'Sp Phase'})
            subplot(nPlots, nIters2Show, nIters2Show+iter);
            plotThisPSTH( PSTH_bins, nextPSTH_vals, extFrameLen_ms, [], timeWindow);
            title(['(' num2str(iter) ')  ' num2str(ccBetweenPSTHs)]);
            ccBetweenPSTHs = crossCorrelation({currPSTH_bins, currPSTH_vals}, {nextPSTH_bins, nextPSTH_vals});
        end

        currPSTH_vals = nextPSTH_vals;
        iter = iter+1;
        if iter > maxIterations
            warning('PSTH:maxIterationsReached', ['Warning: Unable to stabilize PSTH after ' num2str(iter) ' iterations.']);
            break;
        end
        drawnow;

    end
end
