function checkFrameLengths
    hnd = dbOpenExpDb;
    fid1 = 1;
    fid = fopen('checkFrameLengths.log', 'wt');
    fprintf2 = @(s) [fprintf(s), fprintf(fid, s)];
    xifne = @(x) iff(isempty(x), x);
    pad = @(s, N) [repmat(' ', 1, N-length(s)), s];
            
    stimTypes = {'movie', 'noise', 'grating', 'mseq'};
    typesToDo = [1];

for ti = typesToDo
    stimType = stimTypes{ti};
    fprintf2(['DOING ' upper(stimType) ' GROUPS\n-----------------------\n']);
    
    S = load(['cellsGroups_' stimType '_DB']);
    groupData = S.([stimType 'Groups']);
    eval([stimType 'Groups = groupData;']);
    allStimGids = [groupData.Gid];
    
    
    
    
    allGids = getFieldsFromDatabaseTable(hnd, 'GROUP_ID', 'TBL_GROUPS', [], 'GROUP_ID');
    allStimIds = zeros(size(allGids));
    for Gid_i = 1:length(allGids)
        allStimIds(Gid_i) = getStimulusIdForDid( dbLookup('Did', 'Gid', allGids(Gid_i)) );
    end
    
    allFrameLen_ms_Sync = zeros(size(allGids));
    allFrameLen_ms_Defn = zeros(size(allGids));
    allDiscreps = zeros(size(allGids));

    load('frameLengthCalculations.mat');
    

    for Gid_i = 1:length(allGids);
        Gid = allGids(Gid_i);
        
        fprintf2(['Gid = ' pad(num2str(Gid), 4),  ': ']);

        Did = dbLookup('Did',  'Gid', Gid);

        % Movie/Noise/Grating IDs,  begin/end ticks of presentations
        frameLen_ms_Sync = convertTimeMeasuresTemp(Did, 1, 'frame', 'ms', 1);
        frameLen_ms_Defn = convertTimeMeasuresTemp(Did, 1, 'frame', 'ms', 2);

        allFrameLen_ms_Sync(Gid_i) = xifne(frameLen_ms_Sync);
        allFrameLen_ms_Defn(Gid_i) = xifne(frameLen_ms_Defn);

        if ~isempty(frameLen_ms_Sync) && ~isempty(frameLen_ms_Defn)
            if (frameLen_ms_Sync == frameLen_ms_Defn)
                fprintf2(['Good: Both are ' num2str(frameLen_ms_Defn) '.\n' ]);
            else
                discrep = abs(frameLen_ms_Defn - frameLen_ms_Sync) / (frameLen_ms_Sync+eps) * 100;
                fprintf2([' -- ERR: ( ' pad( num2str(discrep, '%3.2f'),5) ' %% )     From syncs: ' pad( num2str(frameLen_ms_Sync, '%4.2f'), 7) '.  From Defn: ' pad(num2str(frameLen_ms_Defn, '%10.3f'),7) '\n']);
            end

            allDiscreps(Gid_i) = discrep;
        else
            fprintf2([' -- ERR: ( ' pad( num2str(0, '%3.2f'),5) ' %% )     From syncs: ' pad( num2str(xifne(frameLen_ms_Sync), '%4.2f'), 7) '.  From Defn: ' pad(num2str(xifne(frameLen_ms_Defn), '%10.3f'),7) '\n']);
            allDiscreps(Gid_i) = NaN;

        end

    end

    figure;
    hist(allDiscreps);
    save('frameLengthCalculations.mat', 'all*');
    keyboard;
%     fclose(fid);
end