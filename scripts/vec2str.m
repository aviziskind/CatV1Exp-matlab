function s = vec2str(v, sep)
    if nargin < 2
        sep = '-';
    end            
    c = [cellfun(@num2str, num2cell(v(:)), 'un', false), repmat({sep}, length(v),1)]';
    c{end} = [];
    s = [c{:}];
end
