function generateNoiseStim(frameDims, nFrames, nGradations, randSeed)
                    
    NoiseRan2('seed', randSeed);
    NoiseRan2;  % generate one random variable after initializing (like the C-version of the program).

    filename = getName('noiseStimFile', randSeed, nGradations);
    fileId = fopen(filename, 'w');
    
    [nrows ncols] = elements(frameDims);
    progressBar('init-', nFrames);
    for fi = 1:nFrames      
        progressBar(fi);
        currFrame = mod( bitshift(NoiseRan2(nrows, ncols), -16), nGradations);
        fwrite(fileId, currFrame, 'uint8');
    end    

    fclose(fileId);
    
end