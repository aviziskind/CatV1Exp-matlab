function [cmpType_id_out, cmpType_out] = curCmpType(s)
    
    persistent cmpType

    filename = [CatV1Path 'MatLabDB_avi' filesep 'curComparisonTypeFile.mat'];
    cmpTypes = {'phase', 'degree', 'clusters', 'psth'};    

    setType = (nargin == 1) && ~isempty(s);
    
    if setType   % set current Comparison type 
    
        if isequal(s, 1) || strncmpi(s, 'phase', 1)
            cmpType = 'phase';
        elseif isequal(s, 2) || strncmpi(s, 'degree', 1)
            cmpType = 'degree';
        elseif isequal(s, 3) || strncmpi(s, 'clusters', 1)
            cmpType = 'clusters';
        elseif isequal(s, 4) || strncmpi(s, 'psth', 1)
            cmpType = 'psth';
        else
            error('invalid Comparison type: enter "p[hase uning]", "d[egree of tuning]", "c[lusters]", or "p[sth]"');
        end
        save(filename, 'cmpType');
        
    else   % retrieve current Comparison type
        if isempty(cmpType)
            if exist(filename, 'file')
                load(filename);
            else
                error('Comparison Type file does not exist');
            end
        end
        cmpType_out = cmpType;
        cmpType_id_out = find(strcmp(cmpType, cmpTypes)); 
        
        if (nargin==1) && isempty(s)
            cmpType_id_out = cmpType;
        end        

        
    end

end
