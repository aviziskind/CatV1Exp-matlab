function tableFields = dbGetTableFields(hnd, tableName)

    % I'm not sure how to do this quickly using SQL. 
    % So I make sure I only need to do this once for each table. A list of
    % all the fields in each table is kept in a data file in the path
    % below. This makes accessing the field names extremely quick.
    
    path = getName('MatlabDB_path');
    filename = 'dbTableFields.mat';
    if exist([path filename], 'file')
        load([path filename]);
    end
    varname = [tableName '_fields'];
    
    if ~exist(varname, 'var')
        sqlstr = ['SELECT * FROM ' tableName];    
        tableFields = edbGetColumnNames(hnd, sqlstr); 

        eval([varname ' = tableFields;']);
        save([path filename], '*_fields');
    else
        tableFields = eval(varname);
    end
        
end
