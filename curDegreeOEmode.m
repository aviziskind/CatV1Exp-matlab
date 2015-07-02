function varargout = curDegreeOEmode(s)

    persistent DegreeOEmode

    filename = [CatV1Path 'MatLabDB_avi' filesep 'curDegreeOEmodeFile.mat'];
    
    allDegreeOEmodes = {'aa', 'oe_same', 'oe_diff', 'oe_same_av', 'oe_diff_av'};
    
    setType = (nargin == 1) && ~isempty(s);
    
    if setType    % set current grating type 
    
        idx = find(strcmp(s, allDegreeOEmodes));
        if isempty(idx)
            error('Invalid DegreeOE mode : %s', s);
        end
        DegreeOEmode = allDegreeOEmodes{idx};

%         DegreeOEmode_out = DegreeOEmode;
        save(filename, 'DegreeOEmode');
        
    else
        if isempty(DegreeOEmode)            
            if exist(filename, 'file')
                S = load(filename);
                DegreeOEmode = S.DegreeOEmode;
            else
                error('DegreeOE mode file does not exist');
            end
        end
%         DegreeOEmode_out = DegreeOEmode;
        
%         DegreeOEmode_id_out = find(strcmp(DegreeOEmode, allDegreeOEmodes)); 
        
%         if (nargin==1) && isempty(s)
%             DegreeOEmode_id_out = DegreeOEmode_out;
%         end                
    end
    
    if nargout <= 1
        varargout = {DegreeOEmode};
    else
        if strcmp(DegreeOEmode, 'aa')
            oe_vs_aa = 'aa';
            diff_vs_same = '';
            keep_vs_av = '';
        else
            oe_vs_aa = 'oe';
            diff_vs_same = switchh(DegreeOEmode(1:7), {'oe_same', 'oe_diff'}, {'same', 'diff'});
            keep_vs_av = iff(strcmp(DegreeOEmode(end-2:end), '_av'), 'average', 'keepBoth');
%             oe_action = switchh(PhaseOEmode(9:end), {'av', 'k1', 'k2', 'k12'}, {'average', 'keep1', 'keep2', 'keepBoth'});
        end
        varargout = {oe_vs_aa, diff_vs_same, keep_vs_av};            
    end
    

end
