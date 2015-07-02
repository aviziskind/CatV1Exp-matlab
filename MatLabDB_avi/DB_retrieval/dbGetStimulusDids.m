function allStimDids = dbGetStimulusDids(stimType)
    % returns all the DatafileIds for stimulus type 'stimType'
    % stimType can be 'movie', 'grating', 'noise', or 'mseq'. (or a cell
    % array of 1 or more of these)

%     removeNonCatV1Experiments = false;
    
    if nargin == 0
        stimType = {'movie', 'grating', 'noise', 'mseq'};
    end

    if iscellstr(stimType)
        allDids = cellfun(@dbGetStimulusDids, stimType, 'un', 0);
        allStimDids = unique(cat(1,allDids{:}));
        return;
    end    
    
    didsFileName = [CatV1Path 'MatLabDB_avi' filesep 'stimulusDids.mat'];
    fieldname = [lower(stimType) 'Dids'];

    
    didFileExists = exist(didsFileName, 'file');        
    if didFileExists
        vars = who('-file', didsFileName);
    end
    if didFileExists && any(strcmp(fieldname, vars))
        S = load(didsFileName, fieldname);            
        allStimDids = S.(fieldname);
        
    else
        hnd = dbOpenExpDb;
        tic; fprintf(['List of Datafile Ids for ' stimType ' stimuli doesn''t exist. Creating ... ']);         
        allStimDids = unique(getFieldsFromDatabaseTable(hnd, 'DATAFILE_ID', ['TBL_' stimType '_PRES']));                        
        appendOption = iff(didFileExists, {'-append'}, {});
        eval([fieldname ' = allStimDids;']);
        save(didsFileName, fieldname, appendOption{:});
        fprintf(' done'); toc;
    end
end