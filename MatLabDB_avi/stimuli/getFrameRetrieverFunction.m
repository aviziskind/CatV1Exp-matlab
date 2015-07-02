function fhandle = getFrameRetrieverFunction(Gid)
    sd = siteDataFor(Gid);
    stimulusType_full = sd.stimType;
    stimulusType = strtok(stimulusType_full, ':');
    switch stimulusType
        case 'Noise',    fhandle = @getNoiseStimulusFrame;
        case 'Mseq',     fhandle = @getMseqStimulusFrame;  % not yet implemented
        case 'Grating',  fhandle = @getGratingStimulusFrame;
        case 'Movie',    fhandle = @getMovieStimulusFrame;
    end
end