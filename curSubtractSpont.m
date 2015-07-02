function tf_out = curSubtractSpont(tf)

    persistent curSubtractSpont_val

    filename = [CatV1Path 'MatLabDB_avi' filesep 'curSubtractSpont.mat'];    

    setType = (nargin == 1) && ~isempty(tf);

    if setType % set current grouping type 
        
        curSubtractSpont_val = logical(tf); 
        save(filename, 'curSubtractSpont_val');
        
    else
        if isempty(curSubtractSpont_val)
            if exist(filename, 'file')
                S = load(filename);
                curSubtractSpont_val = S.curSubtractSpont_val;
            else
                error('Grouping Type file does not exist');
            end
        end        
        tf_out = curSubtractSpont_val; 
        
        if exist('tf', 'var') && isempty(tf)
            tf_out = iff(tf_out, '_SS', '');
        end
        
    end

end
