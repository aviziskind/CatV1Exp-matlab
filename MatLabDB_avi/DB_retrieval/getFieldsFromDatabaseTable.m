function varargout = getFieldsFromDatabaseTable(hnd, fieldNames, tableName, criteria, orderField, nmin, nmax)
    
    nFieldsRequested = iff(iscell(fieldNames), length(fieldNames), 1);
    nOutputs = iff(nargout > 0, nargout, 1);    
    getWholeRecords = (ischar(fieldNames) && strcmp(fieldNames, '*')); 
    if (nFieldsRequested ~= nOutputs)  && ~getWholeRecords
        error(['Number of field names (' num2str(nFieldsRequested) ') does not match number of outputs (' num2str(nOutputs) ')']);
    end

    if isempty(hnd)
        hnd = dbOpenExpDb;
    end
    
    fieldNamesList = cellstr2csslist(fieldNames);
    sel = ['SELECT ' fieldNamesList];
    
    if iscell( tableName )
        tableName = joinTables(tableName);
    end
    frm = [' FROM ' tableName];

    whr = [];
    if (exist('criteria', 'var') && ~isempty(criteria))
        whr = [' WHERE ( ' makeCriteriaList(criteria) ' )']; 
    end
    
    ord = [];
    if (exist('orderField', 'var') && ~isempty(orderField))
        ord = [' ORDER BY ' orderField ';'];
    end
    
    
    queryStr = [sel frm whr ord];
    data = edbRunSelQry(hnd, queryStr);
    
    if exist('nmin', 'var') && ~isempty(nmin) && (size(data, 1) < nmin)
        error(['There are fewer than ' num2str(nmin) ' results for the query']);
    end
    

    convertToStringIfOneString = (exist('nmax', 'var') && ~isempty(nmax)) && (nmax == 1);
    
    if isempty(data)
        varargout = cell(1, nargout);

    elseif getWholeRecords
        nFields = size(data,2);
        if nargout == 1
            varargout{1} = data;
        elseif nargout == 2
            varargout{1} = data;
            
            wholeRecordFieldnames = edbGetColumnNames(hnd, queryStr); 
            varargout{2} = wholeRecordFieldnames;
            
        elseif nargout == nFields
            for i = 1:nFields
                varargout{i} = [data{:,i}]';
            end
        else 
            error(['Number of fields (' num2str(nFields) ') does not match number of outputs (' num2str(nargout) ')']);
        end
            

    else  % assign field columns to individual output arguments
        nRecords = size(data,1);
        if (~exist('nmax', 'var') || isempty(nmax)) || (nmax > size(data,1))
            nmax = nRecords;
        end
        
        for field_i = 1:size(data,2)
            sampleFromField = data{1,field_i};
            
            if isnumeric(sampleFromField)  % if numeric type, convert to double format
                varargout{field_i} = double( vertcat( data{1:nmax,field_i}) );

            elseif ischar(sampleFromField)
                if convertToStringIfOneString && (nRecords == 1)
                    varargout{field_i} = data{:,field_i};  % if just 1 string or nmax = 1, just output the string.
                else
                    varargout{field_i} = data(1:nmax,field_i);  % if cell array of strings - output the cell array without converting to matrix.
                end
            end
        end    
    end

end


%--------------------------------------------------------------------
function from_txt = joinTables(joinTxt)
    
    % syntax1: joinTxt = {tbl1, fld, tbl2}
    % syntax2: joinTxt = {Tbls1, tbls1_fld, fld, tbl2} or {tbl1, fld, Tbls2, tbls2_fld}
    % syntax3: joinTxt = {Tbls1, tbls1_fld, fld, Tbls2, tbls2_fld} 
    
    possibleSyntaxes = {{'char', 'char', 'char'};
                       {'cell', 'char', 'char', 'char'}; 
                       {'char', 'char', 'cell', 'char'};
                       {'cell', 'char', 'char', 'cell', 'char'}};
    synt = cellfun(@class, joinTxt, 'un', 0);
    synt_matches = cellfun(@(s) isequal(synt, s), possibleSyntaxes);    
    if ~any(synt_matches)
        error('unknown syntax');
    end
    switch find(synt_matches)
        case 1,
            [tbl1_fld,fld,tbl2_fld] = deal(joinTxt{:});   
            tbl1_txt=tbl1_fld;
            tbl2_txt=tbl2_fld;
        case 2,
            [Tbls1Grp, tbl1_fld, fld, tbl2_fld] = deal(joinTxt{:});   
            tbl1_txt = ['(' joinTables(Tbls1Grp) ')'];
            tbl2_txt = tbl2_fld;
        case 3,
            [tbl1_fld, fld, Tbls2Grp, tbl2_fld] = deal(joinTxt{:});   
            tbl1_txt = tbl1_fld;
            tbl2_txt = ['(' joinTables(Tbls2Grp) ')'];            
        case 4
            [Tbls1Grp, tbl1_fld, fld, Tbls2Grp, tbl2_fld] = deal(joinTxt{:});   
            tbl1_txt = ['(' joinTables(Tbls1Grp) ')'];
            tbl2_txt = ['(' joinTables(Tbls2Grp) ')'];
        otherwise
            error('Invalid syntax');
    end
            
    from_txt = [tbl1_txt ' INNER JOIN ' tbl2_txt ' ON ' tbl1_fld '.' fld ' = ' tbl2_fld '.' fld ];            
    
end






% function ranges = convertDiscreteIntegersToRanges(vals)
%     grps = continuousGroupings(vals);
%     starts = cellfun(grps, @min);
%     ends   = cellfun(grps, @max);
%     ranges = [starts(:), ends(:)]
%     
% end
% 
% 