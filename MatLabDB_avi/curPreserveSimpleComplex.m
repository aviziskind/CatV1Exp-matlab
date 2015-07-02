function tf_out = curPreserveSimpleComplex(tf)
    % controls the option of preserving # of simple/complex at a site when
    % computing the within-site randomized control distribution

    persistent curPreserveSimpleComplex_val

    filename = [CatV1Path 'MatLabDB_avi' filesep 'curPreserveSimpleComplexFile.mat'];    

    setType = (nargin == 1) && ~isempty(tf);

    if setType % set whether want to only consider pref spf / spf width at the same penetration
        
        curPreserveSimpleComplex_val = logical(tf); 
        save(filename, 'curPreserveSimpleComplex_val');
        
    else
        if isempty(curPreserveSimpleComplex_val)
            if exist(filename, 'file')
                S = load(filename);
                curPreserveSimpleComplex_val = S.curPreserveSimpleComplex_val;
            else
                error('PreserveSC file does not exist');
            end
        end        
        tf_out = curPreserveSimpleComplex_val; 
        
        if exist('tf', 'var') && isempty(tf)
            tf_out = iff(tf_out, '_pSC', '');
        end
        
    end

end
