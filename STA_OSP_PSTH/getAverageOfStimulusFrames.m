function framesAverage = getAverageOfStimulusFrames(fieldName, fieldValue, frameIds)

    % fieldName/Value: 'Gid', or 'Did', with corresponding Gid/Did value.
    % frameIds : index of frames 

    getFrame = getFrameRetrieverFunction(Gid);
    getFrame('load', Gid);
	nFramesToRetrieve = length(frameIds);
    framesAverage = zeros(getFrame('size'));

    % progressBar('init', nFramesToRetrieve, 50);
    for fi = 2:nFramesToRetrieve
        % progressBar(fi);
        framesAverage = framesAverage + getFrame(frameIds(fi));
    end
    framesAverage = framesAverage / nFramesToRetrieve;
    
    getFrame('close');
end