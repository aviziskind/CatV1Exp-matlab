function dbCorrectIpsiContraEye
    hnd = dbOpenExpDb;

   
    [Did, bln_ipsi, bln_contra] = getFieldsFromDatabaseTable(hnd, {'DATAFILE_ID', 'BLN_CONTRA_EYE', 'BLN_IPSI_EYE'}, 'TBL_DATA_FILES');
    
    idx_not_ok = find( (bln_ipsi ~= 0) & (bln_contra ~= 0)); % are 52 of these 
    Did_notOK = Did(idx_not_ok);
    
    
    allSD = siteDataFor(getAllGids);
    allDid_use = [allSD.Did];
    
    [~, ia, ib] = intersect(Did_notOK, allDid_use);
    idx_not_ok_use = idx_not_ok(ia);
    
%     allLocData = [allSD.locationData];
%     allLocIds = [allLocData.LocId];
    
    loc_data = dbGetLocationData;
    loc_data = loc_data(~cellfun(@isempty, {loc_data.Gid}));
    
    loc_Dids = [loc_data.Did];
    loc_LocIds = arrayfun(@(s) s.LocId(end), loc_data ); % [loc_data.LocId];    
%     loc_PenIds = arrayfun(@(s) s.PenetrId(end), loc_data ); % [loc_data.LocId];
        
    
    
    
    locId_prev = nan;
    for i = 1:length(idx_not_ok_use)
        %%
        % find all datafiles from this location
        Did_i = Did(idx_not_ok(i));
        
        locId = loc_LocIds( loc_Dids == Did_i );
        if locId_prev == locId
            continue;
        end
        locId_prev = locId;
%         PenId = loc_PenIds( find(loc_Dids == Did_i) );        

%         allDids_atThisPen = loc_Dids ( find(loc_PenIds == PenId) );        
        allDids_atThisLoc = loc_Dids ( loc_LocIds == locId );
           
        Df_tbl_idx = binarySearch(Did, allDids_atThisLoc);
        
        if length(Df_tbl_idx) == 1
            continue
        end
        
        loc_ipsi   = bln_ipsi(Df_tbl_idx);
        loc_contra = bln_contra(Df_tbl_idx);
        
        fprintf('locId = %d \n', locId);
        [allDids_atThisLoc', loc_ipsi, loc_contra]
        3;
    end
    
    
    %%
    
    for i = 1:0
        pen_right_i = pen_right(cat_idx{i});
        pen_left_i = pen_left(cat_idx{i});
        if ~ all( xor(pen_right_i, pen_left_i))
            continue;
        end
       
        if any(pen_right_i) && any(pen_left_i)
            figure(i); clf; hold on;            
            
            mls = ml( cat_idx{i} );
            aps = ap( cat_idx{i} );
            
            idx_R = (pen_right_i ~= 0);
            idx_L = (pen_left_i ~= 0);

            ids = binarySearch(ap_ml0_id_key, ap_ml0_id( cat_idx{i} ));
            ml0s = ml0(ids); 
            ap0s = ap0(ids);
            
            mls = mls-ml0s;
            aps = aps-ap0s;
            
%             plot(ml0s, ap0s, 'k*');
            3;
            plot(mls(idx_L), aps(idx_L), 'bo', mls(idx_R), aps(idx_R), 'ro'); hold on;
            
            
            
            
            xlabel('ML'); ylabel('AP');
            drawHorizontalLine(0); drawVerticalLine(0);
            3;
        end
        
    end
    
    3;
    %%
    return;
    if ~isempty(idx_not_ok)
       
        for i = 1:length(idx_not_ok)
            switch pen_id(i)
                case 32, % (neither selected) first 9/10 are in L hemisphere. assume this one is left as well
                    L = -1;
                    R =  0;
                case 163 % (L and R selected) all 11 are on the L hemisphere. but for #5, both are selected)
                    L = -1;
                    R =  0;                    
                case 178  % (L and R selected) all 3 are on the L hemisphere. but for 35, both are selected)
                    L = -1;
                    R =  0;                                        
                case 293 % (L and R selected) all 8 are on the L hemisphere. but for #6, both are selected)
                    L = -1;
                    R =  0;                                        
            end
        
%             updateValueInDatabaseTable(hnd, L, 'BLN_HEMISPHERE_L', 'TBL_PENETRATIONS', {'PENETRATION_ID', pen_id(i)} );
%             updateValueInDatabaseTable(hnd, R, 'BLN_HEMISPHERE_R', 'TBL_PENETRATIONS', {'PENETRATION_ID', pen_id(i)} );            
        end
    
    end
end







