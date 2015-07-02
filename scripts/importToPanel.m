function importToPanel(h_fig, h_panel)
    h_ch = get(h_fig, 'children');
    h_ax_all = h_ch( strcmp( get(h_ch, 'type'), 'axes') );
    h_ax_cbar = h_ax_all ( strcmp( get(h_ax_all, 'tag'), 'Colorbar') );
    h_ax_plot = setdiff(h_ax_all, h_ax_cbar);

    for j = 1:length(h_ax_plot)        
        set(h_ax_plot(j), 'XLimMode', 'manual', 'ylimMode', 'manual')
        pause(.1);
        set( h_ax_plot(j), 'parent', h_panel );
    end
    
    %%
    if ~isempty(h_ax_cbar)
        axis(h_ax_plot(1))
        p = get(h_ax_plot(1), 'position');
        colorbar off;
        h_new_cbar = colorbar('EastOutside');
        set(h_ax_plot(1), 'position', p);
        set(h_new_cbar, 'yTick', [0:45:135])
    end
    
end