function str = getResponseTypeStr(responseTypeId, short_flag)

    respTypes_str = {'raw', 'gainCorrected'};
    respTypes_short_str = {'', '_GC'};    
    
    if ischar(responseTypeId) && ~isempty(responseTypeId);
       responseTypeId = find(strcmp(responseTypeId,  respTypes_str));
       assert(~isempty(responseTypeId))
    end
    
    useShort = nargin > 1 && isequal(short_flag, 1);
    if useShort
        str = respTypes_short_str{responseTypeId}; 
    else
        str = respTypes_str{responseTypeId}; 
    end
    
    
end