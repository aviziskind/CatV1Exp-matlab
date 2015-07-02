function tf = dbDoesFieldExist(hnd, fieldName, tableName)

    tableFields = dbGetTableFields(hnd, tableName);
    
    if iscell(fieldName) % ie. multiple field names
        tf = cellfun(@(s) any(strcmp(s, tableFields)), fieldName );
    else
        tf = any(strcmp(fieldName, tableFields));
    end
    
end
