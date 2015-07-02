function [gratingType_id_out, gratingType_out] = curGratingType(s)

    persistent gratingType

    filename = [CatV1Path 'MatLabDB_avi' filesep 'curGratingTypeFile.mat'];
    
    setType = (nargin == 1) && ~isempty(s);
    
    if setType    % set current grating type 
    
        if isequal(s, 1) || strncmp(s, 'flashed', 1)
            gratingType = 'flashed';
        elseif isequal(s, 2) || strncmp(s, 'drifting', 1)
            gratingType = 'drifting';
        else
            error('invalid grating type: enter "f[lashed]" or "d[rifting]"');
        end
        gratingType_out = gratingType;
        save(filename, 'gratingType');
        
    else
        if isempty(gratingType)            
            if exist(filename, 'file')
                S = load(filename);
                gratingType = S.gratingType;
            else
                error('Grating Type file does not exist');
            end
        end
        gratingType_out = gratingType;

        allGratingTypes = {'flashed', 'drifting'};
        gratingType_id_out = find(strcmp(gratingType, allGratingTypes)); 
        
        if (nargin==1) && isempty(s)
            gratingType_id_out = gratingType_out;
        end                
    end

end
