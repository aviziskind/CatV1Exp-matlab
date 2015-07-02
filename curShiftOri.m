function tf_out = curShiftOri(tf)

    persistent curShiftOri_val

    filename = [CatV1Path 'MatLabDB_avi' filesep 'curShiftOri.mat'];    

    setType = (nargin == 1) && ~isempty(tf);

    if setType % set current grouping type 
        
        curShiftOri_val = logical(tf); 
        save(filename, 'curShiftOri_val');
        
    else
        if isempty(curShiftOri_val)
            if exist(filename, 'file')
                S = load(filename);
                curShiftOri_val = S.curShiftOri_val;
            else
                error('Shift ori file does not exist');
            end
        end        
        tf_out = curShiftOri_val; 
        
        if exist('tf', 'var') && isempty(tf)
            tf_out = iff(tf_out, '_shiftOri', '');
        end
        
    end

end
