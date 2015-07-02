function pixelsXY = getMostImportantPixels(Gid, cellId, relContrOfFrameToSpike)

    % fieldName/Value: 'Gid', or 'Did', with corresponding Gid/Did value.
    % frameIds : index of frames

    if (~exist('relContrOfFrameToSpike', 'var') || isempty(relContrOfFrameToSpike))
        relContrOfFrameToSpike = getParsedSpikes('frame',  Gid, cellId, timeWindow );
    end        

    showWorking = true;
    meanIntensity = 0;
    
    getFrame = getFrameRetrieverFunction(Gid);
    getFrame('load', Gid, 'scaling', 'aroundZero');
    nTotalFrames = getFrame('nTotalFrames');
    %     [allFrameIds, spikeFrameIds] = elements(frameIds);

    frameDims = getFrame('size');
    allFramesMean = zeros(frameDims);
    allFramesStd   = zeros(frameDims);
    spikeFramesMean = zeros(frameDims);
    spikeFramesStd   = zeros(frameDims);

    % 	progressBar('init-', nTotalFrames, 50);
    for fi = 1:nTotalFrames
        % progressBar(fi);
        currFrame = getFrame(fi);
        allFramesMean = allFramesMean + currFrame;
        allFramesStd  = allFramesStd  + currFrame .^ 2;
        
        if relContrOfFrameToSpike(fi) > 0
            spikeFramesMean = spikeFramesMean + relContrOfFrameToSpike(fi) * currFrame;
            spikeFramesStd  = spikeFramesStd  + relContrOfFrameToSpike(fi) * (currFrame .^ 2);
        end
    end
    allFramesMean = (allFramesMean / nTotalFrames);
    allFramesStd  = sqrt( (allFramesStd / nTotalFrames) - (allFramesMean).^2 );
    allFramesMean0Sqr = (allFramesMean - meanIntensity) .^2;

    spikeFramesMean = (spikeFramesMean / sum(relContrOfFrameToSpike));
    spikeFramesStd  = sqrt( (spikeFramesStd / sum(relContrOfFrameToSpike)) - (spikeFramesMean).^2 );
    spikeFramesMean0Sqr = abs(spikeFramesMean - meanIntensity);
    varianceDelta = abs(allFramesStd - spikeFramesStd);

    
    desiredFrac = 1/6;
    m0 = 8;
    nPixels = prod(frameDims);

    sigPixels = @(m) ~between( spikeFramesMean0Sqr, allFramesMean0Sqr - allFramesStd/m, allFramesMean0Sqr + allFramesStd/m );

    nExcessSigPixels = @(m) nnz( sigPixels(m) )  - nPixels*desiredFrac;  
    m = fzero(nExcessSigPixels, m0);

    importantPixels = sigPixels(m);
    
    getFrame('close');
    
    if showWorking
        figure(9); clf;
        subplot(2,4,1); graySquareImage(allFramesMean, 'all Frames - mean'); colorbar;
        subplot(2,4,2); graySquareImage(allFramesStd, 'all Frames - var'); colorbar;
        subplot(2,4,3); graySquareImage(allFramesMean0Sqr, 'all Frames - mean^2 around 0'); colorbar;

        subplot(2,4,5); graySquareImage(spikeFramesMean, 'spike Frames - mean'); colorbar;
        subplot(2,4,6); graySquareImage(spikeFramesStd, 'spike Frames - var'); colorbar;
        subplot(2,4,7); graySquareImage(spikeFramesMean0Sqr, 'spike Frames - mean^2 around 0'); colorbar;

        subplot(2,4,4); graySquareImage(importantPixels); colorbar;
        subplot(2,4,8); graySquareImage(varianceDelta, 'var_{all} - var_{spk}'); colorbar;
        
        [pixelsXY(:,1), pixelsXY(:,2)] = find(importantPixels);

    end
    
end