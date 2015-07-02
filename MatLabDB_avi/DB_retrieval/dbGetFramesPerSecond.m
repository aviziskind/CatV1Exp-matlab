function fps = dbGetFramesPerSecond(idType, idVal)
%     fps = dbGetFramesPerSecond('Gid', Gid)
%     fps = dbGetFramesPerSecond('Did', Did)

%     hnd = dbOpenExpDb;
%     tableName = getDatabaseTableForDid(Did);
%     monitorFramesPerSecond = getFieldsFromDatabaseTable(hnd, 'DBL_FRAME_RATE_HZ', tableName, {'DATAFILE_ID', Did}, [], 1, 1);
%     
%     % (2) get  # frames / updates
%     if dbDoesFieldExist(hnd, 'LNG_FRAMES_PER_UPDATE', tableName)
%         framesPerUpdate = getFieldsFromDatabaseTable(hnd, 'LNG_FRAMES_PER_UPDATE', tableName, {'DATAFILE_ID', Did}, [], 1, 1);
%     else
%         framesPerUpdate = 1;  % for gratings
%     end
%     fps = (monitorFramesPerSecond / framesPerUpdate);
% 
%     return;
%     
    persistent allDidFrameRates allGidFrameRates;    
%     mlock;
    
    if isempty(allDidFrameRates)    
        frameRateFileName = [CatV1Path 'MatLabDB_avi' filesep 'dbFrameRates.mat'];    
    
        if exist(frameRateFileName, 'file');
            load(frameRateFileName);
        end
        
        if isempty(allDidFrameRates)
            tic; fprintf('Frame rate file doesn''t exist. Creating ... \n');         
            hnd = dbOpenExpDb;
            allDids = unique(getFieldsFromDatabaseTable(hnd, 'DATAFILE_ID', 'TBL_DATA_FILES'));                
            
            allDidFrameRates = zeros(max(allDids), 1);
            progressBar('init-', length(allDids), 30)
            for i = 1:length(allDids)
                progressBar;
                Did = allDids(i);
                Gid = dbLookup('Gid',  'Did', Did);
                fps = getFramesPerSecond(hnd, Did);            
                allDidFrameRates(Did) = fps;
                allGidFrameRates(Gid) = fps;
            end
            progressBar('done');
            save(frameRateFileName, 'allDidFrameRates', 'allGidFrameRates');
        end
    end
    
    switch idType
        case 'Gid',  
            fps = allGidFrameRates(idVal);
        case 'Did',  
            fps = allDidFrameRates(idVal);
        otherwise, 
            Did = dbLookup('Did',  idType, idVal);
            fps = allDidFrameRates(Did);
    end
    

end
   

function ups = getFramesPerSecond(hnd, Did)
    % DBL_FRAME_RATE_HZ = fixed monitor rate. (usually 120 Hz)
    % we really want the 'update rate', the rate that new images are
    % shown to the screen, which is DBL_FRAME_RATE_HZ / LNG_FRAMES_PER_UPDATE
    
    % (1) get  # frames / second
    tableName = getDatabaseTableForDid(Did);
    monitorFramesPerSecond = getFieldsFromDatabaseTable(hnd, 'DBL_FRAME_RATE_HZ', tableName, {'DATAFILE_ID', Did}, [], 1, 1);
    
    % (2) get  # frames / updates
    if dbDoesFieldExist(hnd, 'LNG_FRAMES_PER_UPDATE', tableName)
        framesPerUpdate = getFieldsFromDatabaseTable(hnd, 'LNG_FRAMES_PER_UPDATE', tableName, {'DATAFILE_ID', Did}, [], 1, 1);
    else
        framesPerUpdate = 1;  % for gratings
    end
    ups = (monitorFramesPerSecond / framesPerUpdate);
end