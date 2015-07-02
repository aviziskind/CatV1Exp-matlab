function tf_out = curSameAnimal(tf)
    % controls the option of only using sites at the same penetration in the between
    % site-distribution

    persistent curSameAnimal_val

    filename = [CatV1Path 'MatLabDB_avi' filesep 'curSameAnimal.mat'];    

    setType = (nargin == 1) && ~isempty(tf);

    if setType % set whether want to only consider pref spf / spf width at the same penetration
        
        curSameAnimal_val = logical(tf); 
        save(filename, 'curSameAnimal_val');
        
    else
        if isempty(curSameAnimal_val)
            if exist(filename, 'file')
                S = load(filename);
                curSameAnimal_val = S.curSameAnimal_val;
            else
                error('SameAnimal file does not exist');
            end
        end        
        tf_out = curSameAnimal_val; 
        
        if exist('tf', 'var') && isempty(tf)
            tf_out = iff(tf_out, '_sameAnimal', '');
        end
        
    end

end
