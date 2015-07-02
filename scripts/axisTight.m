function axisTight(ax)

    ch = get(ax, 'children');
    ch_visible = strcmp(get(ch, 'visible'), 'on');
        
    xdata = get( ch( ch_visible ), 'xdata');
    ydata = get( ch( ch_visible ), 'ydata');

    if ~isempty(xdata)
        if iscell(xdata)
            xdata = [xdata{:}];
            ydata = [ydata{:}];
        end            
        xlims = lims(xdata, .02);
        ylims = lims(ydata, .02);

        set(ax, 'xlim', xlims, 'ylim', ylims);
    end

end