function [bccType_out, bccType_out_str] = curBccType(s)

    persistent bccType

    filename = [CatV1Path 'MatLabDB_avi' filesep 'curBccTypeFile.mat'];
    
    setType = (nargin == 1) && ~isempty(s);

    allBccTypes = {'samePen', 'diffPen', 'sameAnimal', 'diffAnimal', 'diffPen-sameAnimal', 'full'};
    allBccTypes_abbrev = {'WP',  'BP',   'WA',         'BA',         'BP-WA',              ''};
%     bccType_out_str = switchh(bccType, , ...
                                           
    
    if setType    % set current bcc type 
        idx = [find(strcmpi(s, allBccTypes), 1),   find(strcmpi(s, allBccTypes_abbrev), 1)];
        if isempty(idx)            
            error('invalid bcc type: options: %s\n', cellstr2csslist(allBccTypes))  % enter "p[enetration]" or "a[nimal]" or "f[ull]"/("")');
        end
        bccType = allBccTypes{idx};
        bccType_out = bccType;
        save(filename, 'bccType');
        
    else
        if isempty(bccType)            
            if exist(filename, 'file')
                S = load(filename);
                bccType = S.bccType;
            else
                error('BccType Type file does not exist');
            end
        end
        bccType_out = bccType;
        idx = find(strcmpi(bccType, allBccTypes), 1);
        bccType_out_str = allBccTypes_abbrev{idx};
        if ~isempty(bccType_out_str) % ie. not full bcc ("")
            bccType_out_str = ['_' bccType_out_str];
        end
        
        if (nargin==1) && isempty(s)
            bccType_out = bccType_out_str;
        end                
    end

end
