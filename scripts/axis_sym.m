function axis_sym(ax)
    ch = get(ax, 'children');
    hasXY = arrayfun(@(h) isprop(h, 'xdata'), ch);
    ch = ch(hasXY);
    allX = arrayfun(@(h) get(h, 'xdata'), ch, 'un', 0);  allX = [allX{:}];
    allY = arrayfun(@(h) get(h, 'ydata'), ch, 'un', 0);  allY = [allY{:}];
    
    xy_max = max([abs(allX) abs(allY)])*1.05;

    set(ax, 'xlim', [-xy_max, xy_max]);
    set(ax, 'ylim', [-xy_max, xy_max]);    
end