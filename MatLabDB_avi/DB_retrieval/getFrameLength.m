function frameLength = getFrameLength(idType, idVal, outputType)
    if nargin < 3
        outputType = 'ms'; % default;
    end
    
    fps = dbGetFramesPerSecond(idType, idVal);
    frameLength_ms = 1000/fps;
    if ~strcmp(outputType, 'ms')
        Did = dbLookup('Did', idType, idVal);
        frameLength = dbConvertTimeMeasures(Did, frameLength_ms, 'ms', outputType);
    else
        frameLength = frameLength_ms;
    end
    
%     Did = dbLookup('Did',  idType, idVal);
%     frameLength = dbConvertTimeMeasures(Did, 1, 'frame', outputType);
end
