function tf_out = curSamePen(tf)
    % controls the option of only using sites at the same penetration in the between
    % site-distribution

    persistent curSamePen_val

    filename = [CatV1Path 'MatLabDB_avi' filesep 'curSamePen.mat'];    

    setType = (nargin == 1) && ~isempty(tf);

    if setType % set whether want to only consider pref spf / spf width at the same penetration
        
        curSamePen_val = logical(tf); 
        save(filename, 'curSamePen_val');
        
    else
        if isempty(curSamePen_val)
            if exist(filename, 'file')
                S = load(filename);
                curSamePen_val = S.curSamePen_val;
            else
                error('SamePen file does not exist');
            end
        end        
        tf_out = curSamePen_val; 
        
        if exist('tf', 'var') && isempty(tf)
            tf_out = iff(tf_out, '_samePen', '');
        end
        
    end

end
