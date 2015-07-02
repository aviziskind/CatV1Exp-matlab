function updateValueInDatabaseTable(hnd, data, fieldName, tableName, criteria, dateTimeFlag)
    
    if (~ischar(data) && (numel(data) > 1)) || (iscell(fieldName))
        error('Currently, only one data element at a time supported');
    end

    selTxt = ['UPDATE ' tableName ];
    if exist('dateTimeFlag', 'var') && ~isempty(dateTimeFlag)
        setTxt = ['   SET ' fieldName '=CAST(' var2str(data) ' AS TIMESTAMP)'];
    else
        setTxt = ['   SET ' fieldName '=' var2str(data)];
    end

    whrTxt = [];
    if (exist('criteria', 'var') && ~isempty(criteria))
        whrTxt = [' WHERE ( ' makeCriteriaList(criteria) ' )']; 
    end
    
    editStr = [selTxt setTxt whrTxt];

    hRst = actxserver('ADODB.Recordset');
    hRst.Open( editStr, hnd, 3, 1, 1);
        
end

