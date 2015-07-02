function s2 = rmds(s)
    % remove double spaces in the string s
    s2 = strrep(s, '  ', ' ');
    while length(s) ~= length(s2)
        s = s2;
        s2 = strrep(s2, '  ', ' ');        
    end
end

