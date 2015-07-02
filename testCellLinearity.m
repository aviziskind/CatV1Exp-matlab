function testCellLinearity

    S = load('cellsGroups_movie_fg');
    movieGroups_fg = S.movieGroups_fg;
    idx = findInStructArray(movieGroups_fg, 'stimType', 'Movie:Flashed_Gratings:36x10x8(16x1)');
    grp = movieGroups_fg(idx(4));
    
    Gid = grp.Gid;
    cellId = grp.cellIds(1);
    
    grpStyle = 'OSP';

    frameLength_ms = getFrameLength('Gid', Gid);
    [nBinsPerFrame, nFramesPerExtFrame, extFrameLength_ms] = getBinSizeParams(frameLength_ms);
      nBinsPerExtFrame = nBinsPerFrame * nFramesPerExtFrame;

    % 1. Get Position of spike times relative to the frames
    [spkTsRelToFrame_ms_ext, bckgRate, meanFiringRate] = getParsedSpikes('timing', Gid, cellId, 200, extFrameLength_ms); 
    [spkTsRelToFrame_ms, bckgRate, meanFiringRate] = getParsedSpikes('timing', Gid, cellId, 200); 
    nSpikesEachExtFrame = cellfun(@length, spkTsRelToFrame_ms_ext);
    nSpikesEachFrame = cellfun(@length, spkTsRelToFrame_ms);
    % 2a. Use these spike times to calculate the (extended) PSTHs of
    % individual stimuli
    [frameStimIds, uori, usp, uph] = getStimulusFrameSequence(Gid, grpStyle);
    
    nStimIds = max(frameStimIds);
    PSTH_vals = zeros(nBinsPerExtFrame,  nStimIds);
    for stim_i = 1:nStimIds
        frameIndsForStimI = find(frameStimIds == stim_i);
        [PSTH_bins, PSTH_vals(:,stim_i)] = calcPSTH( spkTsRelToFrame_ms_ext(frameIndsForStimI), extFrameLength_ms, length(frameIndsForStimI), [], nBinsPerExtFrame);
    end
    stimFiringRate = mean(PSTH_vals, 1);
    stimInds = ord(max(PSTH_vals, [], 1), 'descend');
    meanPSTH = mean( PSTH_vals(:,stimInds(1:10)), 2 );
    [timeWindow] = getBestTimeWindowFromPSTH(PSTH_bins, meanPSTH);
    b1 = floor((timeWindow(1)+.1)/frameLength_ms);
    b2 = ceil((timeWindow(2)-.1)/frameLength_ms);
    
    psth_dots = zeros(1,nStimIds);
    psth_dots_norm = zeros(1,nStimIds); normMeanPSTH = meanPSTH/mag(meanPSTH);
    for si = 1:nStimIds
        psth_i = PSTH_vals(:,si);
        psth_dots(si) = dot(psth_i, meanPSTH);
        psth_dots_norm(si) = dot(psth_i/mag(psth_i), normMeanPSTH);
    end
%     PSTH_vals_norm = 
    
    
    idxStim1 = stimInds(1);
    frmIndsStim1 = find(frameStimIds == idxStim1);
    
    doFigures = 5;
    
    %Sanity check #0: take value of mean of each frames, and compare to
    %predicted based on mean response to that stimulus    
    if any(doFigures == 1)
        figure(1);
        avRates_PSTH = zeros(1,nStimIds);
        avRates_frm = zeros(1,nStimIds);
        for i = 1:nStimIds        
            avRates_PSTH(i) = stimFiringRate(i);
            idx = find(frameStimIds == i);
            avRates_frm(i)    = mean( nSpikesEachExtFrame(idx) ) / (extFrameLength_ms/1000);
        end
        plot(avRates_PSTH, avRates_frm, '.')
        xlabel('average firing rate (PSTH)'); ylabel('average firing rate (frames)')        
    end
    
    %Sanity check #1: take value of mean of each frames, and compare to
    %predicted based on mean response to that stimulus
    if any(doFigures == 2)
        nFrames = length(frameStimIds);    
        avRates = zeros(1,nFrames);
        indivRates_ext = zeros(1,nFrames);
        indivRates_wind = zeros(1,nFrames);
        for i = 1:nFrames
            avRates(i) = stimFiringRate ( frameStimIds(i) );        
            indivRates_ext(i) = nSpikesEachExtFrame(i)/(extFrameLength_ms/1000);        
            idx_frms = i+b1-1: min(i+b2-1, nFrames);
            indivRates_wind(i) = sum(nSpikesEachFrame(idx_frms))/((frameLength_ms*(b2-b1+1))/1000);        
        end
        figure(12); clf;
        plot(avRates+randn(size(avRates)), indivRates_ext+randn(size(avRates))*3, 'bo', 'markersize', 1); 
        hold on;
        figure(13); clf;
        plot(avRates+randn(size(avRates)), indivRates_wind+randn(size(avRates))*4, 'go', 'markersize', 1); 
        plot(avRates, indivRates_wind, 'g.', 'markersize', 2)

        xlabel('average firing rate'); ylabel('individual firing rates')        
        figure(3);
        plot(avRates, indivRates_wind, 'b.')
    end
%     regressionSlopeTtest(avRates, indivRates, .05, '+', 1);
    
    %Take all successive pairs of frames, and compare actual vs predicted
    %based on sum of means for each
    if any(doFigures == 3)
        nPairs = length(frameStimIds)-1;
        actualRates = zeros(1,nPairs);
        actualRatesRnd = zeros(1,nPairs);
        predictedRates_m = zeros(1,nPairs);
        predictedRates_s = zeros(1,nPairs);
        for i = 1:nPairs
    %         idx_frms = i+b1-1: min(i+b2-1, nFrames);
    %         indivRates_wind(i) = sum(nSpikesEachFrame(idx_frms))/((frameLength_ms*(b2-b1+1))/1000);        

            actualRates(i) = nSpikesEachExtFrame(i+1)/(frameLength_ms/1000);        
            actualRatesRnd(i) = nSpikesEachExtFrame(randi(nPairs+1))/(frameLength_ms/1000);
            predictedRates_m(i) = mean( stimFiringRate ( frameStimIds(i:i+1) ) );        
            predictedRates_s(i) = sum( stimFiringRate ( frameStimIds(i:i+1) ) );        
        end
        figure(30);
        plot(predictedRates_m+randn(size(predictedRates_m)), actualRates+randn(size(predictedRates_m)), 'o', 'markersize', 1)
        xlabel('predicted firing rate'); ylabel('actual firing rate')

        devFromMean = zeros(1,nPairs);
        meanOfPrev = zeros(1,nPairs);
        for fi = 1:nFrames-1
    %         idx_frms = i+b1-1: min(i+b2-1, nFrames);
    %         indivRates_wind(i) = sum(nSpikesEachFrame(idx_frms))/((frameLength_ms*(b2-b1+1))/1000);        
            si = frameStimIds(fi+1);
            devFromMean(fi) = nSpikesEachExtFrame(i)/(frameLength_ms/1000)- stimFiringRate ( si );
            si_prev = frameStimIds(fi);
            meanOfPrev(fi) = stimFiringRate ( si_prev );
        end
        figure(40);
        plot(devFromMean, meanOfPrev, 'o', 'markersize', 1)
        xlabel('dev from mean'); ylabel('mean of prev');
    
    end
    
    if any(doFigures == 4)
        
   
        
        
    end
    
%     mean(nSpikesEachExtFrame(frmIndsStim1))
%     max(nSpikesEachExtFrame(frmIndsStim1))
%     mean(PSTH_vals(:,idx1))

%     idx2 = stimInds(2);    
%     inds2 = find(frameStimIds == idx2)
    


end