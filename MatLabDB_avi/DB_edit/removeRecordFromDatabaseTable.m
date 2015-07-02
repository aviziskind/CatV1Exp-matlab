function removeRecordFromDatabaseTable(hnd, tableName, criteria)
    
    delTxt = ['DELETE FROM ' tableName ];

    whrTxt = [];
    if (exist('criteria', 'var') && ~isempty(criteria))
        whrTxt = [' WHERE ( ' makeCriteriaList(criteria) ' )']; 
    end
    
    removeStr = [delTxt whrTxt];

    hRst = actxserver('ADODB.Recordset');
    hRst.Open( removeStr, hnd, 3, 1, 1);

end

