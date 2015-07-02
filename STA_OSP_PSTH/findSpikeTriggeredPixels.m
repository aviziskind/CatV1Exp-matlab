function [movieMean, movieVar, spkTrigMean, spkTrigVar] = findSpikeTriggeredPixels(movieInfo, spikesPerFrame)
    nFrames = length(spikesPerFrame);
    
    xs = movieInfo.xs;
    ys = movieInfo.ys;
    [nrows, ncols] = elements(movieInfo.dims);
    Nrep = movieInfo.Nrep;
    movieFileName = movieInfo.filename;

    
    [xs_grid, ys_grid] = meshgrid(xs, ys);
    hnd = dbOpenExpDb;
    fid = fopen(movieFileName);
    
    movieMean = zeros(nx, ny);
    movieVar = zeros(nx, ny);

    spkTrigMean = zeros(nx, ny);
    spkTrigVar = zeros(nx, ny); 
    
    for fi = 1:nFrames
        thisFrame = (fread(fid, [nrows ncols]) - averIntensity) / (averIntensity);

        movieMean = movieMean + thisFrame;
    	movieVar  = movieVar + thisFrame.^2;

        spkTrigMean = spkTrigMean + thisFrame * spikesPerFrame(fi);
        spkTrigVar  = spkTrigVar  + (thisFrame * spikesPerFrame(fi)).^2; 
        
    end
    fclose(movieFileName);
    
    movieMean = movieMean/nFrames;
    movieVar  = movieVar/nFrames - (movieMean.^2);

    spkTrigMean = spkTrigMean/nFrames;
    spkTrigVar  = spkTrigVar/nFrames - (spkTrigMean.^2);
    
    figure(24); clf(24);
    subplot(2,2,1); imagesc(movieMean'); title('Movie mean');
    subplot(2,2,2); imagesc(movieVar'); title('Movie variance');    
    subplot(2,2,3); imagesc(spkTrigMean'); title('Spike triggered mean');    
    subplot(2,2,4); imagesc(spkTrigVar');  title('Spike triggered variance');    

    
end