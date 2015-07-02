function y = anyNonNanIn(x)
    x = x(~isnan(x));
    if isempty(x)
        y = nan;
        return;
    elseif length(x) == 1
        y = x;
    else
        y = x(randi(length(x)));
    end

end