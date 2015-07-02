function tf_out = curMatchDB(tf)

    persistent curMatchDB_val

    filename = [CatV1Path 'MatLabDB_avi' filesep 'curMatchDB.mat'];    

    setType = (nargin == 1) && ~isempty(tf);

    if setType % set current grouping type 
        
        curMatchDB_val = logical(tf); 
        save(filename, 'curMatchDB_val');
        
    else
        if isempty(curMatchDB_val)
            if exist(filename, 'file')
                S = load(filename);
                curMatchDB_val = S.curMatchDB_val;
            else
                error('Grouping Type file does not exist');
            end
        end        
        tf_out = curMatchDB_val; 
        
        if exist('tf', 'var') && isempty(tf)
            tf_out = iff(tf_out, '_DB', '');
        end
        
    end

end
