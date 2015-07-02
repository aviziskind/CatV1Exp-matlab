function allStimGids = dbGetStimulusGids(stimType, keepAllFlag)
    % returns all the GroupIds for stimulus type 'stimType'
    % stimType can be 'movie', 'grating', 'noise', or 'mseq'. (or a cell
    % array of 1 or more of these)
    if (nargin == 0) || isempty(stimType)
        stimType = {'movie', 'grating', 'noise', 'mseq'};
    end

    removeNonCatV1Experiments = ~exist('keepAllFlag', 'var') || isempty(keepAllFlag);

    if iscellstr(stimType)
        if removeNonCatV1Experiments
            allGids = cellfun(@dbGetStimulusGids, stimType, 'un', 0);
        else
            allGids = cellfun(@(s) dbGetStimulusGids(s, 1), stimType, 'un', 0);
        end
            
        allStimGids = unique(cat(1,allGids{:}));
        return;
    end
        
    gidsFileName = [CatV1Path 'MatLabDB_avi' filesep 'stimulusGids.mat'];
    
    fieldname = [lower(stimType) 'Gids'];
    fieldname_all = [fieldname '_all'];
    fieldname_return = iff(removeNonCatV1Experiments, fieldname, fieldname_all);            
            
    ghostGids = [134 470 622 801 802];

    redo = 0;
    gidFileExists = exist(gidsFileName, 'file');        
    if gidFileExists
        S = load(gidsFileName);
        flds = fieldnames(S);         
    else
        S = struct;
    end
    
    if ~gidFileExists || ~any(strcmp(fieldname_return, flds)) || redo
        hnd = dbOpenExpDb;
        tic; fprintf(['List of Group Ids for ' stimType ' stimuli doesn''t exist. Creating ... ']);         
        allStimDids_all = dbGetStimulusDids(stimType);        
        allStimGids_all = sort(getFieldsFromDatabaseTable(hnd, 'GROUP_ID', 'TBL_GROUPS', {'DATAFILE_ID', allStimDids_all(:)}));        
        allStimGids = allStimGids_all;
%         [uDids, idx_ofGid, idx_ofDid] = unique(correspondingDids);
        
        % remove non-CatV1 experiments for 'allStimGids' variable
        if true
        
            % there are some that are labeled as from LGN, or are from Rat or a Diode (test)
            locData = dbGetLocationData('Gid', allStimGids);
            idx_nonCatV1 = arrayfun(@(s) ~strcmp(s.AnimalType, 'Cat') || ~strcmp(s.brainStruct, 'V1'), locData);

            % there are some that appear to be photodiode tests 
            % (<4 channels, and spike times are regular), but are not
            % labeled as such. we remove them here.        
            haveAtLeast4TotChannels = 1;
            haveAtLeast4UsedChannels = 1;  % some groups with just 3 active channels seem like real experiments... (?)

            idxTooFewChannels = false(length(allStimGids),1);
            if haveAtLeast4TotChannels
                correspondingDids = arrayfun(@(Gid) dbLookup('Did', 'Gid', Gid), allStimGids); % some datafiles have more than 1 group associated with them.
                nChannelsTot = arrayfun(@(Did) getFieldsFromDatabaseTable(hnd, 'LNG_N_CHANNELS', 'TBL_DATA_FILES', {'DATAFILE_ID', Did}), correspondingDids);        
                idxTooFewChannels = idxTooFewChannels | (nChannelsTot < 4);
            end
            if haveAtLeast4UsedChannels
                nChannelFields = 16;
                channelsUsed = cell(1,nChannelFields);
                channelsUsedFieldnames = arrayfun(@(ch_id) ['BLN_CH_' num2str(ch_id, '%02d')], 1:nChannelFields, 'un', 0);
                nChannelsUsed = zeros(length(allStimGids),1);
                for gi = 1:length(allStimGids)                
                    [channelsUsed{1:nChannelFields}] = getFieldsFromDatabaseTable(hnd, channelsUsedFieldnames, 'TBL_GROUPS', {'GROUP_ID', allStimGids(gi)});
                    nChannelsUsed(gi) = nnz([channelsUsed{:}]);
                end                                        
                idxTooFewChannels = idxTooFewChannels | (nChannelsUsed < 4);
            end

            idx_remove = idx_nonCatV1 | idxTooFewChannels;
            fprintf('%d nonCatV1, %d test, %d total to remove, out of %d total\n', nnz(idx_nonCatV1), nnz(idxTooFewChannels), nnz(idx_remove), nnz(allStimGids) );
            allStimGids(idx_remove) = [];
            allStimGids = setdiff(allStimGids(:), ghostGids(:));        
        end
        
%         appendOption = iff(gidFileExists, {'-append'}, {});
        S.(fieldname) = allStimGids;
        S.(fieldname_all) = allStimGids_all;        
        save(gidsFileName, '-struct', 'S');
        fprintf(' done'); toc;
    end
        
    allStimGids = S.(fieldname_return);


    
end


