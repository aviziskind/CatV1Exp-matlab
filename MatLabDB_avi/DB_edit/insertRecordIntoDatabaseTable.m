function insertRecordIntoDatabaseTable(hnd, tableName, values)
    
    insTxt = ['INSERT INTO ' tableName ];

    fieldNames = dbGetTableFields(hnd, tableName);
%     if length(fieldNames) ~= length(values)
%         error(['Number of values must be the same as the number of fields ( ' num2str(length(fieldNames)) ')']);
%     end
    if length(fieldNames) ~= length(values)
        fieldNames = fieldNames(1:length(values));
    end
    
    fieldListStr = [' (' cellstr2csslist(fieldNames) ')'];
    
    if iscell(values)
        valuesC = values;
    else
        [m,n] = size(values);
        valuesC = mat2cell(values, m, n);
    end
    valueStrs = var2str(valuesC);
    valueListStr = [' VALUES ( ' cellstr2csslist(valueStrs) ')'];
    
    insStr = [insTxt fieldListStr valueListStr];

    hRst = actxserver('ADODB.Recordset');
    hRst.Open( insStr, hnd, 3, 1, 1);

end

