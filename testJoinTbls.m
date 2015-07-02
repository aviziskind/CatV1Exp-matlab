function testJoinTbls



% joinedT
% {
% {'TBL_DATA_FILES', 'DATAFILE_ID', 'TBL_LOCS_FILES_LINKS'}

strFrm = FROM TBL_ELECTRODE_TYPES ...
           INNER JOIN ((((TBL_ANIMALS ...             %#ok<*NOPRT,*NOPTS>
             INNER JOIN ((TBL_AP_ML_ZERO  
               INNER JOIN TBL_ELECTRODES ON TBL_AP_ML_ZERO.ELECTRODE_ID = TBL_ELECTRODES.ELECTRODE_ID) ...
                 INNER JOIN TBL_PENETRATIONS ON TBL_AP_ML_ZERO.AP_ML_ZERO_ID = TBL_PENETRATIONS.AP_ML_ZERO_ID) ON TBL_ANIMALS.ANIMAL_ID = TBL_ELECTRODES.ANIMAL_ID) ...
                   INNER JOIN TBL_LOCATIONS ON TBL_PENETRATIONS.PENETRATION_ID = TBL_LOCATIONS.PENETRATION_ID) ...
                      INNER JOIN (TBL_LOC_DEPTHS ...
                        INNER JOIN (TBL_DATA_FILES 
                          INNER JOIN TBL_LOCS_FILES_LINKS ON TBL_DATA_FILES.DATAFILE_ID = TBL_LOCS_FILES_LINKS.DATAFILE_ID) ON TBL_LOC_DEPTHS.LOC_DEPTH_ID = TBL_LOCS_FILES_LINKS.LOC_DEPTH_ID) ON TBL_LOCATIONS.LOCATION_ID = TBL_LOC_DEPTHS.LOCATION_ID) ...
                            INNER JOIN TBL_GROUPS ON TBL_DATA_FILES.DATAFILE_ID = TBL_GROUPS.DATAFILE_ID) ON TBL_ELECTRODE_TYPES.ELECTRODE_TYPE_ID = TBL_ELECTRODES.ELECTRODE_TYPE_ID ;

                        %#ok<*NOPRT,*NOPTS>
strFrm = FROM T9 ;
                                                
   T1 = {'TBL_AP_ML_ZERO', 'ELECTRODE_ID', 'TBL_ELECTRODES'};               
   T2 = {T1, 'TBL_AP_ML_ZERO', 'AP_ML_ZERO_ID', 'TBL_PENETRATIONS'};               
   T3 = {'TBL_ANIMALS', 'ANIMAL_ID', T2, 'TBL_ELECTRODES'}               
   T4 = {T3, 'TBL_PENETRATIONS', 'PENETRATION_ID', 'TBL_LOCATIONS'}               
   T5 = {'TBL_DATA_FILES', 'DATAFILE_ID', 'TBL_LOCS_FILES_LINKS'};
   T6 = {'TBL_LOC_DEPTHS', 'LOC_DEPTH_ID',  T5, 'TBL_LOCS_FILES_LINKS'};
   T7 = {T4, 'TBL_LOCATIONS', 'LOCATION_ID', T6, 'TBL_LOC_DEPTHS'};
   T8 = {T7, 'TBL_DATA_FILES', 'DATAFILE_ID', 'TBL_GROUPS'};
   T9 = {'TBL_ELECTRODE_TYPES', 'ELECTRODE_TYPE_ID', T8, 'TBL_ELECTRODES'};
               
                
                  
                        

function from_txt = joinTables(joinTxt)
    [tbl1,fld,tbl2] = deal(joinTxt{:});
    if iscell(tbl1)
        tbl1 = joinTables(tbl1);
    end
    if iscell(tbl2)
        tbl2 = joinTables(tbl2);
    end
    from_txt = [tbl1 ' INNER JOIN ' tbl2 ' ON ' tbl1 '.' fld ' = ' tbl2 '.' fld ];
end


hnd = dbOpenExpDb;
stimId = getFieldsFromDatabaseTable(hnd, 'STIMULUS_TYPE_ID', 'TBL_DATA_FILES', {'DATAFILE_ID', {@lt, 76}})


tbl_str = {'TBL_DATA_FILES', 'DATAFILE_ID', 'TBL_LOCS_FILES_LINKS'}
% tbl_str2 = 'TBL_DATA_FILES, TBL_LOCS_FILES_LINKS  WHERE TBL_DATA_FILES.DATAFILE_ID=TBL_LOCS_FILES_LINKS.DATAFILE_ID';
tbl_str_expl = 'TBL_DATA_FILES INNER JOIN TBL_LOCS_FILES_LINKS ON TBL_DATA_FILES.DATAFILE_ID=TBL_LOCS_FILES_LINKS.DATAFILE_ID'; 
[stimId, loc_depth_id] = getFieldsFromDatabaseTable(hnd, {'STIMULUS_TYPE_ID', 'LOC_DEPTH_ID'}, tbl_str_expl, {'TBL_DATA_FILES.DATAFILE_ID', {@lt, 76}})


FROM TBL_DATA_FILES, TBL_LOCS_FILES_LINKS  

ON  TBL_DATA_FILES.DATAFILE_ID=TBL_LOCS_FILES_LINKS.DATAFILE_ID 

'LOC_DEPTH_ID', 'TBL_LOC_DEPTHS'};



'LOC_DEPTH_ID', 'TBL_LOC_DEPTHS'};

tbl_str_expl = '(TBL_DATA_FILES INNER JOIN TBL_LOCS_FILES_LINKS ON TBL_DATA_FILES.DATAFILE_ID=TBL_LOCS_FILES_LINKS.DATAFILE_ID) INNER JOIN TBL_LOC_DEPTHS ON TBL_LOCS_FILES_LINKS.LOC_DEPTH_ID = TBL_LOC_DEPTHS.LOC_DEPTH_ID'; 

data_links = {'TBL_DATA_FILES', 'DATAFILE_ID', 'TBL_LOCS_FILES_LINKS'};
tbl_str = {data_links, 'TBL_LOCS_FILES_LINKS', 'LOC_DEPTH_ID',   'TBL_LOC_DEPTHS'};
[stimId, loc_depth_id, loc_id] = getFieldsFromDatabaseTable(hnd, {'STIMULUS_TYPE_ID', 'TBL_LOCS_FILES_LINKS.LOC_DEPTH_ID', 'LOCATION_ID'}, tbl_str, {'TBL_DATA_FILES.DATAFILE_ID', {@lt, 76}})



end