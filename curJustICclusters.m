function tf_out = curJustICclusters(tf)

    persistent curJustICclusters_val

    filename = [CatV1Path 'MatLabDB_avi' filesep 'curJustICclusters.mat'];    

    setType = (nargin == 1) && ~isempty(tf);

    if setType %
        
        curJustICclusters_val = logical(tf); 
        save(filename, 'curJustICclusters_val');
        
    else
        if isempty(curJustICclusters_val)
            if exist(filename, 'file')
                S = load(filename);
                curJustICclusters_val = S.curJustICclusters_val;
            else
                error('Grouping Type file does not exist');
            end
        end        
        tf_out = curJustICclusters_val; 
        
        if exist('tf', 'var') && isempty(tf)
            tf_out = iff(tf_out, '_DB', '');
        end
        
    end

end
