function [respType_id_out, respType_out] = curResponseType(s)
    
    persistent respType

    filename = [CatV1Path 'MatLabDB_avi' filesep 'curResponseTypeFile.mat'];
    respTypes = {'raw', 'gainCorrected'};    

    setType = (nargin == 1) && ~isempty(s);
    
    if setType   % set current Comparison type 
    
        if isequal(s, 1) || strncmpi(s, 'raw', 1)
            respType = 'raw';
        elseif isequal(s, 2) || strncmpi(s, 'gainCorrected', 1)
            respType = 'gainCorrected';
        else
            error('invalid Response type: enter "r[aw]", "g[ainCorrected]", "m[ID predicted]"');
        end
        save(filename, 'respType');
        
    else   % retrieve current Comparison type
        if isempty(respType)
            if exist(filename, 'file')
                load(filename);
            else
                error('Response type file does not exist');
            end
        end
        respType_out = respType;
        respType_id_out = find(strcmpi(respType, respTypes)); 
        
        if (nargin==1) && isempty(s)
            respType_id_out = respType;
        end        
        
    end

end
