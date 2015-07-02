function t2 = dbConvertTimeMeasures(Did, t1, fromType, toType)

%     hnd = dbOpenExpDb;
    
    if strcmp(fromType, toType) 
        t2 = t1;
        return;
    end    
    if any(strcmp('tick', {fromType, toType}))
        sd = siteDataFor('Did', Did);
        ticksPerSecond = sd.dataFileInfo.samplingRateHz;        
%         ticksPerSecond = getFieldsFromDatabaseTable(hnd, 'DBL_SAMPLING_RATE_HZ', 'TBL_DATA_FILES', {'DATAFILE_ID', Did}, [], 1, 1);; 
    end
    if any(strcmp('frame', {fromType, toType}))
        framesPerSecond = dbGetFramesPerSecond('Did', Did); 
    end
    msPerSecond = 1000;
        
    

    %%%%%%% DO THE TIME-CONVERSION %%%%%%%    
        
    % convert from 'fromType' to 'seconds'
    switch fromType
        case 'sec'
            timeSec = t1;
            
        case 'ms'
            timeSec = t1 / msPerSecond;
            
        case 'frame'              
            timeSec = t1 / framesPerSecond;

        case 'tick'
            timeSec = t1 / ticksPerSecond;
            
        otherwise
            error('unknown input type');
            
    end

    % convert from 'seconds' to 'toType'
    switch toType
        case 'sec'
            t2 = timeSec;

        case 'ms'
            t2 = timeSec * msPerSecond;

        case 'frame'
            t2 = timeSec * framesPerSecond;

        case 'tick'            
            t2 = round(timeSec * ticksPerSecond);
    
        otherwise
            error('unknown output type')
    end

end


% when converting from 'ms'/'sec'/'tick' to 'frame', we want to know exactly
% which frame (and the relative position in the frame) that that ms was in --> use binarySearch with (0.5).
%        syncTimes_sec = dbGetSyncs('Did', Did, 'sec');
%        relPosInFrame = binarySearch(syncTimes_sec, timeSec, 0, 0.5);
%        t2 = relPosInFrame;    
