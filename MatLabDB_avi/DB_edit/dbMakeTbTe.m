function dbMakeTbTe(hnd, Did)
    % This is my approach to generating the begin/end ticks for
    % presentations. I think it is a bit more robust
    % (using my dbParseSyncs function), and it is more
    % flexibile & adaptable when syncs are not perfect.
    global DidsNeedTbTeFixing;
    dbug = true;

    [stimulusType subType] = getStimulusTypeForDid(Did);
    stimTypes = {'movie', 'noise', 'grating', 'mseq'};
    stimulusTypeId = find(strcmpi(stimulusType, stimTypes));
    
    fprintf(['[' stimulusType ':' subType ']']);
    tableName = getDatabaseTableForDid(Did);

	syncs = dbGetSyncs('Did', Did, 'tick', 1);
    if isempty(syncs)
        fprintf('No syncs');
        return;
    end
    if ~isempty(syncs) && ~(any(syncs(1:5) < 10))   % add beginning sync.
        syncs = [0; syncs]; 
    end
    
    [nFramesPerPres_Runs, endStatus, runs] = dbParseSyncs(Did, syncs);
    if (endStatus < 0)
        fprintf('Couldn''t parse syncs. ');
        return;
    end
    emptyPresInds = find(nFramesPerPres_Runs == 0);
    
    nFramesPerPres_Runs(emptyPresInds) = [];   
    presRuns = [ones(size(nFramesPerPres_Runs)), nFramesPerPres_Runs]';    
    
    C1 = cumsum(runs(:,2))+1;
    C = cumsum(presRuns(:))+1;
    indStarts = C(1:2:end-1);
    indEnds   = C(2:2:end);
    tb = syncs(indStarts);
    te = syncs(indEnds);

    [tbInDB, teInDB, presNum] = getFieldsFromDatabaseTable(hnd, {'LNG_START_TICK', 'LNG_END_TICK', 'LNG_PRESENT_NO'}, tableName, {'DATAFILE_ID', Did}, 'LNG_PRESENT_NO');
    tbInDB(emptyPresInds) = [];
    teInDB(emptyPresInds) = [];
    presNum(emptyPresInds) = [];
    
    output = [{'BeginTick: Now', 'BeginTick: DB', 'diff', 'EndTick: Now', 'EndTick: DB', 'diff'}; ...
        num2cell([tb, tbInDB, tb-tbInDB, te, teInDB, te-teInDB])];
%     disp(output);

    discrep = sum(abs([tb-tbInDB;te-teInDB]));

    if dbug
        if (discrep == 0)
            fprintf(' MATCHES DATABASE: good:');
%             disp(runs);
        else
            disp(' MISMATCH: editing database:');
%             disp(runs);
            disp(output);
            if ~all(isnan(tbInDB))
%                 beep; keyboard;
            end
            
            for i = 1:length(presNum)
%                 DidsNeedTbTeFixing{stimulusTypeId} = unique([DidsNeedTbTeFixing{stimulusTypeId}, Did]);
                updateValueInDatabaseTable(hnd, tb(i), 'LNG_START_TICK', tableName, {'DATAFILE_ID', Did; 'LNG_PRESENT_NO', presNum(i)});
                updateValueInDatabaseTable(hnd, te(i), 'LNG_END_TICK',   tableName, {'DATAFILE_ID', Did; 'LNG_PRESENT_NO', presNum(i)});
            end

        end
    end


end


