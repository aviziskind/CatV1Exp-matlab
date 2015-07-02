function allCriteria = makeCriteriaList(criteriaPairs)
    if ischar(criteriaPairs)
        allCriteria = criteriaPairs;
        return;
    end
    
    if size(criteriaPairs, 2) > 2
        warning('getFields:criteria', 'criteriaPairs variable not constructed correctly. Some criteria may be omitted. \n Be sure to put each set of criteria on a separate row');
    end
    for ci = 1:size(criteriaPairs,1)
        name = criteriaPairs{ci,1};
        opStr = '=';
        values = criteriaPairs{ci,2};
        % 'values' can either be: a single number; a vector of numbers,   [[ if given a row vector, it is automatically flipped to a column vector ]];
        % a single string or a cell array of strings.
        if (isnumeric(values) || iscell(values)) && isrow(values)
            values = values';  % must be a column vector when convert to string.
        end
        if ~ischar(values) && size(values, 2) > 2
            warning('getFields:criteria','criteriaPairs variable not constructed correctly (Alternate values must be in columns). Some criteria may be omitted');
        end
        if isempty(values)
            error('Criteria value field cannot be empty if you specify a criteria name');
        end
        if iscell(values)   % allow for {'>', N} or {'<=', N} or {@ge, N}, or {@lt, N}
            opStr = values{1};
            if isa(opStr, 'function_handle')
                opStr = op2str(opStr, 'SQL');
            end
            values = values{2};
        end
        values = var2str(values);  % need to feed a string to the query.

        for vi = 1:size(values, 1)
            if iscell(values)
                thisPart = ['((' name ')' opStr values{vi} ')'];
            else
                thisPart = ['((' name ')' opStr values(vi,:) ')']; % add all
            end
            if vi == 1
                thisCriterion = ['(' thisPart];
            else
                thisCriterion = [thisCriterion ' OR ' thisPart];   %#ok<AGROW>
            end
        end
        thisCriterion = [thisCriterion ')'];  %#ok<AGROW>

        %  add this criterion to the list of all criteria
        if (ci == 1)
            allCriteria = thisCriterion;
        else
            allCriteria = [allCriteria ' AND ' thisCriterion]; %#ok<AGROW>
        end
    end
end
