function tf_out = curIdealPruning(tf)

    persistent curIdealPruning_val

    filename = [CatV1Path 'MatLabDB_avi' filesep 'curIdealPruning.mat'];    

    setType = (nargin == 1) && ~isempty(tf);

    if setType % set current grouping type 
        
        curIdealPruning_val = logical(tf); 
        save(filename, 'curIdealPruning_val');
        
    else
        if isempty(curIdealPruning_val)
            if exist(filename, 'file')
                S = load(filename);
                curIdealPruning_val = S.curIdealPruning_val;
            else
                error('Ideal Pruning file does not exist');
            end
        end        
        tf_out = curIdealPruning_val; 
        
        if exist('tf', 'var') && isempty(tf)
            tf_out = iff(tf_out, '_ideal', '');
        end
        
    end

end
