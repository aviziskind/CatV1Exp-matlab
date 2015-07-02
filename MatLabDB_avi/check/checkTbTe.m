% function checkTbTe
hnd = dbOpenExpDb;
fileOutput = 'checkTbTe.log';

stimTypes = {'movie', 'noise', 'grating', 'mseq'};
typesToDo = [1];
fid = fopen(fileOutput, 'wt');
fprintf2 = @(s) [fprintf(s), fprintf(fid, s)];

for ti = typesToDo
    stimType = stimTypes{ti};
    fprintf2(['DOING ' upper(stimType) ' GROUPS\n-----------------------\n']);
    
    S = load(['cellsGroups_' stimType '_DB']);
    groupData = S.([stimType 'Groups']);
    eval([stimType 'Groups = groupData;']);
    allStimGids = [groupData.Gid];
    
    gidsWithTbTeNaNs = [];
    gidsWithNoTbTe = [];
    gidsWithNoSyncs = [];
    gidsWithOffsetsFromSync1 = [];
    gidsWithDiscrepFromDB = [];
    
    %     Gid_start = find(allGids == 1030, 1);
    tableName = getDatabaseTableForDid(dbLookup('Did',  'Gid', allStimGids(1)));
    for Gid_i = 1:length(allStimGids);
        cmpSyncs = true;
        Gid = allStimGids(Gid_i);
        ok = true;
        
        Did = dbLookup('Did',  'Gid', Gid);
        fprintf2(['(' outOf(Gid_i, length(allStimGids)) ')  Gid = ' num2str(Gid),  ', Did = ' num2str(Did) ' : ']);
        
        [beginTicks, endTicks] = dbGetTbTe(hnd, Did);
        if isempty(beginTicks)
            gidsWithNoTbTe = [gidsWithNoTbTe Gid]; %#ok<AGROW>
            fprintf2('   No ticks');
            ok = false;            
        end
        if isnan(beginTicks)
            gidsWithTbTeNaNs = [gidsWithTbTeNaNs Gid]; %#ok<AGROW>
            fprintf2('   Ticks are NaNs');
            ok = false;
        end
        
        syncs = dbGetSyncs('Did', Did, 'tick', 1);
        if isempty(syncs)
            gidsWithNoSyncs = [gidsWithNoSyncs Gid]; %#ok<AGROW>
            fprintf2('   No Syncs');
            ok = false;
        end
        
        if ok
            % compare with syncs(1) and syncs(end)
            nPres = length(beginTicks);
            initialOffset = beginTicks(1) - syncs(1);
            finalOffset =   endTicks(end) - syncs(end);
            % fprintf2([' initial offset = ' num2str(initialOffset), '.  finalOffset = ' num2str(finalOffset)]);
            if (initialOffset ~= 0) || (finalOffset ~= 0)
                gidsWithOffsetsFromSync1 = [gidsWithOffsetsFromSync1 Gid]; %#ok<AGROW>
                fprintf2(['  offsets of [' num2str(initialOffset), ', ' num2str(finalOffset) ']']);
                ok = false;
            else
                fprintf2('  (no offset)');
            end
        end
        
        if ok
            % compare with db fields
            [nSustainedFrames] = getFieldsFromDatabaseTable(hnd, {'LNG_N_SUSTAINED_FRM'}, tableName, {'DATAFILE_ID', Did});
            if ~isempty(nSustainedFrames)
                nFrames_fromSyncs = length(syncs)-nPres;  % each pres causes an extra sync to be used on an inter-presentation spike. And I added 2 above (at beginning & end)
                nFrames_fromDB = sum(nSustainedFrames);
                if (nFrames_fromDB ~= nFrames_fromSyncs)
                    fprintf2(['  nFramesMismatch (syncs: ' num2str(nFrames_fromSyncs) ', db: ' num2str(nFrames_fromDB) ')']);
                    gidsWithDiscrepFromDB = [gidsWithDiscrepFromDB Gid]; %#ok<AGROW>
                    ok = false;
                else
                    fprintf2(' (db match) ');
                end
            end

        end
        
        if ok
            fprintf2('ok');
        end
        fprintf2('\n');
        
    end
    
    
    
end
fclose(fid);
    % end
    
    %
    %
    
    
    
    
    
    
    
    
    
    
    
    %       gidsWithTbTeNaNs = [];
    %         gidsWithNoTbTe = [500 890 891 896 897 909 910 930];
    %         gidsWithNoSyncs = [801 802 1024];
    %         gidsWithOffsetsFromSync1 = [133 134 437 470 606 838 839 1019 1020];
    %         gidsWithDiscrepFromDB = [];
    %