function tableData = getTableData(tableName, allFieldNames)
    tableFileName = [CatV1Path 'MatLabDB_avi' filesep 'DB_matlabcopy' filesep tableName '.mat'];
%     tableName = [lower(stimType) 'PresTable'];

    tableFileExists = exist(tableFileName, 'file');
    if tableFileExists
        allFieldsInTable = who('-file', tableFileName);        
                
        whichFieldsPresent = cellfun(@(fld) any(strcmp(fld, allFieldsInTable)), allFieldNames);        
        
%         tableData = load(tableFileName, tableName);
        
        if all(whichFieldsPresent);
            tableData = load(tableFileName, allFieldNames{:});
            return;
        else
            disp( [ 'Fields : ' allFieldNames{~whichFieldsPresent} ' not present. Reloading to table']);
            tableData = load(tableFileName, allFieldNames{whichFieldsPresent});
        end
        fieldsToLoad = setdiff(allFieldNames, allFieldsInTable);
    else
        fieldsToLoad = allFieldNames;
    end
    
    tic;
    fprintf(['Loading fields from table ' tableName ' ... ']);
    hnd = dbOpenExpDb;
    [newTableData{1:length(fieldsToLoad)}] = getFieldsFromDatabaseTable(hnd, fieldsToLoad, tableName);
    for fld_i = 1:length(fieldsToLoad)
        eval([fieldsToLoad{fld_i} ' = newTableData{fld_i};']);    % put into current working space
        tableData.(fieldsToLoad{fld_i}) = newTableData{fld_i};   % put into struct for return argument
    end    
    fprintf('done\n');
    toc;
    
    appendArgs = iff(tableFileExists, {'-append'}, {});
    save(tableFileName, fieldsToLoad{:}, appendArgs{:}, '-v6');

end
