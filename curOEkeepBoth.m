function tf_out = curOEkeepBoth(tf)

    persistent curOEkeepBoth_val

    filename = [CatV1Path 'MatLabDB_avi' filesep 'curOEkeepBoth.mat'];    

    setType = (nargin == 1) && ~isempty(tf);

    if setType % set current grouping type 
        
        curOEkeepBoth_val = logical(tf); 
        save(filename, 'curOEkeepBoth_val');
        
    else
        if isempty(curOEkeepBoth_val)
            if exist(filename, 'file')
                S = load(filename);
                curOEkeepBoth_val = S.curOEkeepBoth_val;
            else
                error('Grouping Type file does not exist');
            end
        end        
        tf_out = curOEkeepBoth_val; 
        
        if exist('tf', 'var') && isempty(tf)
            tf_out = iff(tf_out, '', '_av');
        end
        
    end

end
