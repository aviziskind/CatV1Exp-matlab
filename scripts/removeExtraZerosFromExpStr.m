function str = removeExtraZeros(str)
    idx_e = strfind(lower(str), 'e');
    if isempty(idx_e)
        return;
    end
    
%     while true        
%         N = length(str);
    while strcmp(str(idx_e+2), '0')
        str(idx_e+2) = '';
    end
            

end
