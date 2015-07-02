function varargout = curPhaseOEmode(s)

    persistent PhaseOEmode

    filename = [CatV1Path 'MatLabDB_avi' filesep 'curPhaseOEmodeFile.mat'];
    
    allPhaseOEmodes = {'aa', ...
                      'oe_same_av', 'oe_same_k1', 'oe_same_k2', 'oe_same_k12', ...
                      'oe_diff_av', 'oe_diff_k1', 'oe_diff_k2', 'oe_diff_k12'};
    
    
    setType = (nargin == 1) && ~isempty(s);
    
    if setType    % set current grating type 
    
        idx = find(strcmp(s, allPhaseOEmodes));
        if isempty(idx)
            error('Invalid PhaseOE mode : %s', s);
        end
        PhaseOEmode = allPhaseOEmodes{idx};

%         PhaseOEmode_out = PhaseOEmode;
        save(filename, 'PhaseOEmode');
        
    else
        if isempty(PhaseOEmode)            
            if exist(filename, 'file')
                S = load(filename);
                PhaseOEmode = S.PhaseOEmode;
            else
                error('PhaseOE mode file does not exist');
            end
        end
%         PhaseOEmode_out = PhaseOEmode;
        
%         PhaseOEmode_id_out = find(strcmp(PhaseOEmode, allPhaseOEmodes)); 
        
%         if (nargin==1) && isempty(s)
%             PhaseOEmode_id_out = PhaseOEmode_out;
%         end                
    end
    
    if nargout <= 1
        varargout = {PhaseOEmode};
    else
        nice_str = '';
        if strcmp(PhaseOEmode, 'aa')
            oe_vs_aa = 'aa';
            diff_vs_same = '';
            oe_action = '';
            nice_str = '(all trials)';
        else
            oe_vs_aa = 'oe';
            diff_vs_same = switchh(PhaseOEmode(1:7), {'oe_same', 'oe_diff'}, {'same', 'diff'});
            oe_action = switchh(PhaseOEmode(9:end), {'av', 'k1', 'k2', 'k12'}, {'average', 'keep1', 'keep2', 'keepBoth'});
            
            nice_str = ['odd/even trials' iff(strcmp(diff_vs_same, 'diff'), '', '(same)'), iff(strcmp(oe_action, 'keepBoth'), '', [';' oe_action])];
        end
        varargout = {oe_vs_aa, diff_vs_same, oe_action, nice_str};            
    end

end
