function result = setFieldsInDatabaseTable(hnd, data, fieldNames, tableName, criteria, orderField)
    
    nInputsSpecified = iff(iscell(fieldNames), length(fieldNames), 1);
    nInputsProvided = iff(ischar(data), 1, size(data,2));
    if (nInputsSpecified ~= nInputsProvided) 
%         error(['Number of field names (' num2str(nInputsSpecified) ') does not match data size ( ' num2str(nInputsProvided) ')']);
    end

    fieldNamesList = cellstr2csslist(fieldNames);
    sel = ['SELECT ' fieldNamesList];
    frm = [' FROM ' tableName];
    
    whr = [];
    if (exist('criteria', 'var') && ~isempty(criteria))
        whr = [' WHERE ( ' makeCriteriaList(criteria) ' )']; 
    end

    ord = [];
    if (exist('orderField', 'var') && ~isempty(orderField))
        ord = [' ORDER BY ' orderField ';'];
    end
    
    
    editStr = [sel frm whr ord];
    
    if ~iscell(data) 
        if ischar(data)
            dataC = cellstr(data); 
        else
            dataC = num2cell(data); 
        end
    else
        dataC = data;
    end    
    
    result = edbRunEditQry(hnd, editStr, dataC);

end

