function timeWindow_str = getTimeWindowStr(timeWindow)
%     timeWindow_str = iff(isnumeric(timeWindow), sprintf('__%d_%d', timeWindow), iff(strcmp(timeWindow, 'best'), '', ['__' timeWindow] ));
    if isnumeric(timeWindow)
        timeWindow_str = sprintf('__%d_%d', timeWindow);
    elseif strcmp(timeWindow, 'best') || isempty(timeWindow)
        timeWindow_str = '';
    end
%         || isempty(timeWindow)
%         
%         , iff(strcmp(timeWindow, 'best'), '', ['__' timeWindow] ));
end