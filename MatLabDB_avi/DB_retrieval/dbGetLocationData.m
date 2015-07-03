function loc_data = dbGetLocationData(idtype, idval)
    
    persistent allLocationData;

%     nStepsEdits = getLocationChanges();
 
    saveMultipleInstances = 0;  % a couple places with double tetrodes -- just pick 1
    ii = 1;
%     allLocationData = [];
    redo = 0;
    
    if isempty(allLocationData) || redo
        locationsFileName = [CatV1Path 'MatLabDB_avi\locationNames.mat'];
        if exist(locationsFileName, 'file') && ~redo
            S = load(locationsFileName);
            allLocationData = S.allLocationData;
        else    
            maxGid = 5336;     %2450 Gids
%             maxDid = 3822;     %2384 Dids
            
            allLocationData = repmat(struct, maxGid,1);

%             allLocationData.Gid = cell(maxGid,1);
%             allLocationData.Did = cell(maxDid,1);
            hnd = dbOpenExpDb;

%             allDids = getFieldsFromDatabaseTable(hnd, 'DATAFILE_ID', 'TBL_DATA_FILES');
            allGids = dbGetStimulusGids([], 1); % get location data for all sites.
            
            fieldnames = {'TBL_DATA_FILES.DATAFILE_ID', 'TBL_GROUPS.GROUP_ID', 'TBL_ANIMALS.TXT_LAB_NAME', ...
                'TBL_ANIMALS.ANIMAL_ID','TBL_PENETRATIONS.PENETRATION_ID', 'TBL_LOCATIONS.LOCATION_ID',  'TBL_LOC_DEPTHS.LOC_DEPTH_ID', ...
                'TBL_LOCATIONS.TXT_LOCATION_NAME', 'TBL_PENETRATIONS.TXT_PENETRATION_NAME', ...
                'TBL_ELECTRODE_TYPES.ELECTRODE_TYPE_ID', 'TBL_PENETRATIONS.BRAIN_STRUCT_TYPE_ID', ...
                'TBL_ELECTRODES.TXT_ELECTRODE_NAME', 'TBL_ELECTRODES.ELECTRODE_ID', 'TBL_ELECTRODES.MEM_ELECTRODE_NOTES', ...
                ...
                'TBL_PENETRATIONS.BLN_HEMISPHERE_L', 'TBL_PENETRATIONS.BLN_HEMISPHERE_R',  ...
                'TBL_PENETRATIONS.DBL_AP', 'TBL_PENETRATIONS.DBL_ML', ...
                'TBL_AP_ML_ZERO.DBL_AP_ZERO', 'TBL_AP_ML_ZERO.DBL_ML_ZERO', 'TBL_LOC_DEPTHS.DBL_DEPTH_MANIP_STEPS', ...
                ...
                'TBL_GROUPS.DTM_CREATED', 'TBL_DATA_FILES.DTM_CREATED', 'TBL_PENETRATIONS.DTM_CREATED', ...
                'TBL_LOCATIONS.DTM_CREATED', 'TBL_AP_ML_ZERO.DTM_CREATED', ...
                };
% 'TBL_LOCATION.DBL_LEFT_OPTIC_DISK_X', 'DBL_LEFT_OPTIC_DISK_Y', 'DBL_RIGHT_OPTIC_DISK_X
            
            T1 = {'TBL_AP_ML_ZERO', 'ELECTRODE_ID', 'TBL_ELECTRODES'};
            T2 = {T1, 'TBL_AP_ML_ZERO', 'AP_ML_ZERO_ID', 'TBL_PENETRATIONS'};               
            T3 = {'TBL_ANIMALS', 'ANIMAL_ID', T2, 'TBL_ELECTRODES'};
            T4 = {T3, 'TBL_PENETRATIONS', 'PENETRATION_ID', 'TBL_LOCATIONS'};
            T5 = {'TBL_DATA_FILES', 'DATAFILE_ID', 'TBL_LOCS_FILES_LINKS'};
            T6 = {'TBL_LOC_DEPTHS', 'LOC_DEPTH_ID',  T5, 'TBL_LOCS_FILES_LINKS'};
            T7 = {T4, 'TBL_LOCATIONS', 'LOCATION_ID', T6, 'TBL_LOC_DEPTHS'};
            T8 = {T7, 'TBL_DATA_FILES', 'DATAFILE_ID', 'TBL_GROUPS'};
            T9 = {'TBL_ELECTRODE_TYPES', 'ELECTRODE_TYPE_ID', T8, 'TBL_ELECTRODES'};
            joinedTables = T9;
            
                        
%             criterea = {'TBL_ELECTRODES.ELECTRODE_TYPE_ID', 2};
            criterea = {};
            tic;
            [Did_tbl, Gid_tbl, CatName_tbl, AnimalID_tbl, PenetrID_tbl, LocID_tbl, LocDepthID_tbl, LocName_tbl, PenName_tbl, ...
                electTypeId_tbl, brainStruct_tbl, electrode_name_tbl, electrode_id_tbl, electrode_notes_tbl, ...
                ...
                l_hemi_tbl, r_hemi_tbl, ap_tbl, ml_tbl, ap_zero_tbl, ml_zero_tbl, depth_steps_tbl, ...
                date_grp, date_df, date_pen, date_loc, date_apml] = ...
                getFieldsFromDatabaseTable(hnd, fieldnames, joinedTables, criterea);

            toc;
          
            isThomas = cellfun(@(str) ~isempty(strfind(str, 'Thomas')), electrode_notes_tbl);
            isNiCr = cellfun(@(str) ~isempty(strfind(lower(str), 'nic')), electrode_notes_tbl);
            isTungsten12_7 = cellfun(@(str) ~isempty(strfind(str, 'Tungsten 12.7')), electrode_notes_tbl)';
            isTungsten7_6 = cellfun(@(str) ~isempty(strfind(str, 'Tungsten wire, 7.6um')), electrode_notes_tbl)';
            isTungsten25 = cellfun(@(str) ~isempty(strfind(str, 'Tungsten wire, center; 25.4 um')), electrode_notes_tbl)';
%%
            if 0
                %%
                Did_tbl(idx)
                Gid_tbl(idx),
                CatName_tbl(idx),
                AnimalID_tbl(idx), 
                PenetrID_tbl(idx), 
                LocID_tbl(idx), 
                LocName_tbl(idx), 
                PenName_tbl(idx), ...
                elect_tbl(idx), 
                
                brainStruct_tbl(idx), 
                electrode_name_tbl(idx), 
                electrode_notes_tbl(idx), ...
                ...
                %%
                l_hemi_tbl(idx), r_hemi_tbl(idx), 
                %%
                ap_tbl(idx), ml_tbl(idx), ap_zero_tbl(idx), ml_zero_tbl(idx), depth_steps_tbl(idx) 
                
            end

            %%
            electrode_types_tbl = cell(size(electrode_notes_tbl));
            electrode_types_tbl(isThomas) = {'Thomas'};
            electrode_types_tbl(isNiCr) = {'NiCr 12.7'};
            electrode_types_tbl(isTungsten12_7) = {'Tungsten 12.7'};
            electrode_types_tbl(isTungsten7_6) = {'Tungsten 7.6'};
            electrode_types_tbl(isTungsten25) = {'Tungsten 25'};
            
            [nStepsEdits_tblCols, nStepsEdits_tbl] = getLocationChanges();
        %%
        
%             printGidRecords(2500);
                        
            
%             showIdx = @(i) [{'DatafileId', 'GroupID', 'AnimalId', 'PenetrId', 'LocId', 'ElectrodeId'};
%                             num2cell([Did_tbl(i(:)), Gid_tbl(i(:)), AnimalID_tbl(i(:)), PenetrID_tbl(i(:)), LocID_tbl(i(:)), elect_tbl(i(:))])] ;
            
                        
            3;
%             l_idx = zeros(size(allGids));
%             progressBar('init-', length(allGids), 30);                
            for Gid_i = 1:length(allGids)
%                 progressBar(Gid_i);
                Gid = allGids(Gid_i);
                
                if Gid == 5
                    3;
                end
                
                idx = find(Gid_tbl == Gid);                
%                  n_idx(Gid_i) = length(idx);
%                 l_idx(Gid_i) = length(idx);
                idx_orig = idx;
                
                if length(idx) > 1
                    electrodeTypeIds = electTypeId_tbl(idx);
                    nTetrodes = nnz(electrodeTypeIds == 2);
                    idx2 = find(electrodeTypeIds == 2); % 2 == Tetrode.
                    
                    if nTetrodes == 0
                        idx = idx(1);
                        
                    elseif nTetrodes == 1
                        idx = idx(electrodeTypeIds == 2);  % pick the tetrode one
                                            
                    elseif nTetrodes == 2
                        
% %                         if 
%                         try
%                             sd = siteDataFor('Did', Did_tbl(idx(1))); nCh = sd.dataFileInfo.nChannels;
%                         catch
%                             nCh = 0;
%                         end
%                         fprintf('Gid_i = %d. Gid = %d. Did = %d. (nCh = %d) types = %d/%d, ids = %d/%d. names= %s/%s\n', ...
%                             Gid_i, Gid_tbl(idx(1)), Did_tbl(idx(1)), nCh, electTypeId_tbl(idx),  electrode_id_tbl(idx), electrode_name_tbl{idx(1)}, electrode_name_tbl{idx(2)});

                        idx = idx(1); % not sure why there are 2 in many cases.... for the purposes of i just picked 1
                    end
                end
                
%                     else
%                         Did = dbLookup('Did',  'Gid', Gid);
%                         [st, sbt] = getStimulusTypeForDid(Did);
%                         cts = CatName_tbl(idx);
%                         pns = PenetrID_tbl(idx);
%                         lcs = LocID_tbl(idx);
%                         els = elect_tbl(idx);
%                         fprintf('Duplicates: Gid = %d [%s:%s]. Cats: %s. Pens = %s. Locs = %s, Electrodes = %s.\n', Gid, st, sbt, cellstr2csslist(cts), vec2csslist(pns), vec2csslist(lcs), vec2csslist(els));
%                     end
%                 end                        
                
%                 if isempty(idx)
%                     Did = dbLookup('Did',  'Gid', Gid);
%                     [st, sbt] = getStimulusTypeForDid(Did);
%                     fprintf('No entry: Gid = %d [%s : %s].\n', Gid, st, sbt);
%                     continue;
%                 end
                
%                 if length(idx) > 1
%                     3;
%                     Did = dbLookup('Did',  'Gid', Gid);
%                     [st, sbt] = getStimulusTypeForDid(Did);
%                     cts = CatName_tbl(idx);
%                     pns = PenetrID_tbl(idx);
%                     lcs = LocID_tbl(idx);
%                     els = elect_tbl(idx);
%                     fprintf('Duplicates: Gid = %d [%s:%s]. Cats: %s. Pens = %s. Locs = %s, Electrodes = %s.\n', Gid, st, sbt, cellstr2csslist(cts), vec2csslist(pns), vec2csslist(lcs), vec2csslist(els));
%                 end

                if length(idx) == 1
                    3;                    
                end
                
                if Gid == 2196
                    3;
                end
%                 Gid = dbLookup('Gid',  'Did', Did);
%                 loc_depth_id = getFieldsFromDatabaseTable(hnd, 'LOC_DEPTH_ID', 'TBL_LOCS_FILES_LINKS', {'DATAFILE_ID', Did});
%                 loc_id = getFieldsFromDatabaseTable(hnd, 'LOCATION_ID', 'TBL_LOC_DEPTHS', {'LOC_DEPTH_ID', loc_depth_id});
%                 loc_name = getFieldsFromDatabaseTable(hnd, 'TXT_LOCATION_NAME', 'TBL_LOCATIONS', {'LOCATION_ID', loc_id});                            
%                 
%                 allLocationData.Did(Did) = {loc_name};
%                 allLocationData.Gid(Gid).CatName = CatName_tbl(idx);
%                 allLocationData.Gid(Gid).PenetrId = PenetrID_tbl(idx);
%                 allLocationData.Gid(Gid).LocId = LocID_tbl(idx);



                if length(idx) > 1
                    3;
                end

                CatName = CatName_tbl{idx};
                switch CatName(1),
                    case 'K',        AnimalType = 'Cat'; 
                    case {'M', 'R'}, AnimalType = 'Rat'; 
                    case 'T',        AnimalType = 'Test';
                end                                
                brainStruct = switchh(unique(brainStruct_tbl(idx)), [1,2,100], {'LGN', 'V1', 'Diode'});
                
                electrode_types = electrode_types_tbl(idx);
%                 CatId = str2double( strtok(CatName, 'K') );
%                 if isnan(CatId)
%                     3;
%                 end   

                isUnique = @(x) length(unique(x)) == 1;

                if length(idx_orig) > 1
                    Animal_orig = AnimalID_tbl(idx_orig);
                    Pen_orig = PenetrID_tbl(idx_orig);
                    Loc_orig = LocID_tbl(idx_orig);
                    Gid_orig = Gid_tbl(idx_orig);
                    Did_orig = Did_tbl(idx_orig);
                    
%                     if ~isUnique(Animal_orig);
%                         fprintf('(%d) Animal Id : %d \n', ii, length(unique(Animal_orig))); 
%                     end
%                     if ~isUnique(Pen_orig);
%                         fprintf('(%d) Pen Id : %d \n', ii, length(unique(Pen_orig))); 
%                     end
%                     if ~isUnique(Loc_orig);
%                         fprintf('(%d) Loc Id : %d \n', ii, length(unique(Loc_orig))); 
%                     end
%                     if ~isUnique(Gid_orig);
%                         fprintf('(%d) Gid : %d \n', ii, length(unique(Gid_orig))); 
%                     end
%                     if ~isUnique(Did_orig);
%                         fprintf('(%d) Did : %d \n', ii, length(unique(Did_orig))); 
%                     end                   
%                     
%                     if ~isUnique(Animal_orig) || ~isUnique(Pen_orig) || ~isUnique(Loc_orig) || ~isUnique(Gid_orig) || ~isUnique(Did_orig);
%                         ii=ii+1;
%                     end                    
                end
                if saveMultipleInstances
                    idx_save = idx_orig;
                else
                    idx_save = idx;
                end

                haveElecVals = ~cellfun(@isempty, electrode_types_tbl(idx_orig));
                if nnz(haveElecVals) == 0
                    electrode_type = '';
                elseif (length(idx_orig) == 1) || (nnz(haveElecVals) == 1) || isUnique( electrode_types_tbl(idx_orig) )
                    electrode_type = electrode_types_tbl{idx_orig(haveElecVals)};
                else
                    electrode_type = ['? ' cellstr2csslist(electrode_types_tbl(idx_orig), ' / ')];
                end
                
                allLocationData(Gid).AnimalType  = AnimalType;
                allLocationData(Gid).brainStruct = brainStruct;
                allLocationData(Gid).ElectrodeType = electrode_type; %unique( electrode_types_tbl(idx_orig) );                
                allLocationData(Gid).AnimalId    = unique( AnimalID_tbl(idx_save) );
                allLocationData(Gid).PenetrId    = unique( PenetrID_tbl(idx_save) );
                allLocationData(Gid).LocId       = unique( LocID_tbl(idx_save) );
                allLocationData(Gid).LocDepthId  = unique( LocDepthID_tbl(idx_save) );
                allLocationData(Gid).CatName     = unique( CatName_tbl(idx_save) );
                allLocationData(Gid).PenName     = unique( PenName_tbl(idx_save) );
                allLocationData(Gid).LocName     = unique( LocName_tbl(idx_save) );
                allLocationData(Gid).Gid         = unique( Gid_tbl(idx_save) );
                allLocationData(Gid).Did         = unique( Did_tbl(idx_save) );
                     
                isLeft = unique(l_hemi_tbl(idx_save)) ~= 0;
                isRight = unique(r_hemi_tbl(idx_save)) ~= 0;
                if (isLeft ~= isRight)
                    LR_str = iff(isLeft, 'L', 'R');
                else
                    LR_str = '?';
                    isLeft = nan;
                    isRight = nan;
                end
                
               
                
                [nStepsEdits_tblCols, nStepsEdits_tbl] = getLocationChanges();
                 idx_loc_depth_id = find(strcmp(nStepsEdits_tblCols, 'LocDepthID'),1);
                idx_location_id = find(strcmp(nStepsEdits_tblCols, 'LocationID'),1);
                idx_curSteps = find(strcmp(nStepsEdits_tblCols, 'CurrentSTEPS'),1);
                idx_newSteps = find(strcmp(nStepsEdits_tblCols, 'newSTEPS'),1);
                
               
                depth = depth_steps_tbl(idx_save);
                 
                if length(idx_save) == 1
                    idx_replace = find(  nStepsEdits_tbl(:,idx_loc_depth_id) == LocDepthID_tbl(idx_save) );
                    if ~isempty(idx_replace)
                        assert(length(idx_replace) == 1);
                        assert(nStepsEdits_tbl(idx_replace, idx_location_id) == LocID_tbl(idx_save) );
                        assert(nStepsEdits_tbl(idx_replace, idx_curSteps) == depth_steps_tbl(idx_save) );
                        newNSteps = nStepsEdits_tbl(idx_replace, idx_newSteps);
                        
                        depth = newNSteps;
                        
                    end
                    
                end
                
               
                if any(depth == [0, 999, 9999, 2222222])
                    depth = nan;
                end
%                 if 
                
            
                allLocationData(Gid).hemi       = LR_str;
                allLocationData(Gid).hemiId     = double(isRight);
                allLocationData(Gid).AP         = unique( ap_tbl(idx_save) );
                allLocationData(Gid).ML         = unique( ml_tbl(idx_save) );
                allLocationData(Gid).AP_ZERO    = unique( ap_zero_tbl(idx_save) );
                allLocationData(Gid).ML_ZERO    = unique( ml_zero_tbl(idx_save) );
                allLocationData(Gid).depth      = depth;

                allLocationData(Gid).Pen_date  = unique( date_pen(idx_save) );
                allLocationData(Gid).Loc_date  = unique( date_loc(idx_save) );
                allLocationData(Gid).APML_date = unique( date_apml(idx_save) );
                allLocationData(Gid).Df_date   = unique( date_df(idx_save) );
                allLocationData(Gid).Grp_date  = unique( date_grp(idx_save) );
                
                
                3;
                
            end
%             progressBar('done');
            3;
            save(locationsFileName, 'allLocationData');
        end
        
    end
    
    if nargin == 0
        loc_data = allLocationData;
        return;
    end
    
    if ~strcmp(idtype, 'Gid')
        Gid = dbLookup('Gid', idtype, idval);
    else
        Gid = idval;
    end
    loc_data = allLocationData(Gid);
%     switch lower(idtype)
%         case 'did'
%             loc_name = allLocationData.Did{idval};
%         case 'gid'
%             loc_name = allLocationData.Gid{idval};
%         otherwise
%             error('Idtype must be "Gid" or "Did"');
%     end
    
end

function s = trunc(s, n)
    if length(s) > n
        s = s(1:n);
    end        
    s(s==13) = '/';
    s(s==10) = '';
end

function printGidRecords(Gid)

    T1 = {'TBL_AP_ML_ZERO', 'ELECTRODE_ID', 'TBL_ELECTRODES'};
    T2 = {T1, 'TBL_AP_ML_ZERO', 'AP_ML_ZERO_ID', 'TBL_PENETRATIONS'};
    T3 = {'TBL_ANIMALS', 'ANIMAL_ID', T2, 'TBL_ELECTRODES'};
    T4 = {T3, 'TBL_PENETRATIONS', 'PENETRATION_ID', 'TBL_LOCATIONS'};
    T5 = {'TBL_DATA_FILES', 'DATAFILE_ID', 'TBL_LOCS_FILES_LINKS'};
    T6 = {'TBL_LOC_DEPTHS', 'LOC_DEPTH_ID',  T5, 'TBL_LOCS_FILES_LINKS'};
    T7 = {T4, 'TBL_LOCATIONS', 'LOCATION_ID', T6, 'TBL_LOC_DEPTHS'};
    T8 = {T7, 'TBL_DATA_FILES', 'DATAFILE_ID', 'TBL_GROUPS'};
    T9 = {'TBL_ELECTRODE_TYPES', 'ELECTRODE_TYPE_ID', T8, 'TBL_ELECTRODES'};
    joinedTables = T9;

    

    hnd = dbOpenExpDb;
    criterea = {'TBL_GROUPS.GROUP_ID', Gid};
    [recs, fld_names] = getFieldsFromDatabaseTable(hnd, '*', joinedTables, criterea);
    rec_diff = zeros(size(fld_names));
    for i = 1:length(recs)
        rec_diff(i) = ~isequalwithequalnans(recs(1,i), recs(2,i));
    end
    [fld_names, recs', num2cell(rec_diff)]
    nTrunc = 40;
    for i = 1:length(recs)
        if strcmp(fld_names{i}, 'MEM_SPIKETIMES_MTX')
            3;
        end
        fprintf('%20s : %40s, %40s, [%d]\n', trunc(fld_names{i}, nTrunc), trunc(var2str(recs{1,i}), nTrunc), trunc(var2str(recs{2,i}), nTrunc), rec_diff(i))                 
    end
3;

end


function [nStepsEdits_tblCols, nStepsEdits_tbl] = getLocationChanges()
    nStepsEdits_tblCols = {'LocDepthID', 'LocationID', 'CurrentSTEPS', 'newSTEPS'};
    nStepsEdits_tbl = [
        43  35, -1, 8379
        492, 453, 0, 13754
        507, 469, 0, 19375
        509, 471, 0, 19640
        511, 473, 19640, 20345
        513, 475, 0, 20926
        515, 477, 0, 21602
        517, 479, 0, 22262
        519, 481, 0, 23336
        521, 483, 0, 24000
        523, 485, 0, 8924
        525, 487, 0, 9440
        1004, 944, 0, 1
        1005, 945, 0, 1];
    
end
