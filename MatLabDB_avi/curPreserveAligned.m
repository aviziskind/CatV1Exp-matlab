function tf_out = curPreserveAligned(tf)
    % controls the option of preserving # of simple/complex at a site when
    % computing the within-site randomized control distribution

    persistent curPreserveAligned_val

    filename = [CatV1Path 'MatLabDB_avi' filesep 'curPreserveAlignedFile.mat'];    

    setType = (nargin == 1) && ~isempty(tf);

    if setType % set whether want to only consider pref spf / spf width at the same penetration
        
        curPreserveAligned_val = logical(tf); 
        save(filename, 'curPreserveAligned_val');
        
    else
        if isempty(curPreserveAligned_val)
            if exist(filename, 'file')
                S = load(filename);
                curPreserveAligned_val = S.curPreserveAligned_val;
            else
                error('PreserveAligned file does not exist');
            end
        end        
        tf_out = curPreserveAligned_val; 
        
        if exist('tf', 'var') && isempty(tf)
            tf_out = iff(tf_out, '_pAB', '');
        end
        
    end

end
