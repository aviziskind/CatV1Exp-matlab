% f
rawFigTitles = {'sampleOriMapB', 'diffOriHist', 'median_vs_sigma', 'pct_outliers_vs_sigma', 'distSigma', 'distPctOutliers', 'medianDiff_vs_dist', 'propOutliers_vs_dist'}; % 'oriMapSample', 
nSubplots = 8;

    xlabel_fsize = 12;
    ylabel_fsize = 11;
    title_fsize = 11;
    legend_fsize = 8;


fig_id = 14;
figure(fig_id); clf;
% set(fig_id, 'color', 'w', 'position', [1000, 200, 780, 900]);
figureFolder = [CatV1Path 'Figures' filesep 'DegreePaper' filesep];

%%
fig_filename = sprintf('%sFigure%d_simulatedOriMaps.pdf', figureFolder, fig_id);
set(fig_id, 'color', 'w', 'windowstyle', 'normal',  'colormap', hsv(250));
refresh(fig_id)
row_height_pix = 320;
set(fig_id, 'position', [1000, 150, 700, row_height_pix*4]);
    
%%
subM = 4;
subN = 2;
h = [];

subSpcM = [0 0 0];
subSpcN = [0 0, 0];


for i = 1:length(rawFigTitles)
    if isempty(rawFigTitles{i}) 
        continue;
    end
    charLabel = 'A' + i -1;
    subFigName = sprintf('%sFig15%s_simOri_%s.fig', figureFolder, charLabel, rawFigTitles{i});
    h_saved(i) = hgload(subFigName);
    pause(.1);

    %%
    subplotInfo_C = {subM, subN, i, [], subSpcM, subSpcN};
    p_subplot = getNormPosition(subplotInfo_C{:});

    h_panel(i) = uipanel( fig_id, 'units','normal','pos', p_subplot, 'BackgroundColor', 'w', 'bordertype', 'none' );
%     set(h_panel(i)
%     ohs=hgload( 'untitled.fig');
%%
    h_ch = get(h_saved(i), 'children');
    h_ax_all = h_ch( strcmp( get(h_ch, 'type'), 'axes') );
    h_ax_cbar = h_ax_all ( strcmp( get(h_ax_all, 'tag'), 'Colorbar') );
    h_ax_plot = setdiff(h_ax_all, h_ax_cbar);

    for j = 1:length(h_ax_plot)
        
        set(h_ax_plot(j), 'XLimMode', 'manual', 'ylimMode', 'manual')
        pause(.1);
        set( h_ax_plot(j), 'parent', h_panel(i) );
        
        h_xlab = get(h_ax_plot(j), 'xLabel');  set(h_xlab, 'fontsize', xlabel_fsize);
        h_ylab = get(h_ax_plot(j), 'yLabel');  set(h_ylab, 'fontsize', ylabel_fsize);
        
            %%
    end
    
    %%
    if ~isempty(h_ax_cbar)
        axis(h_ax_plot(1))
        p = get(h_ax_plot(1), 'position');
        colorbar off;
        h_new_cbar = colorbar('EastOutside');
        set(h_ax_plot(1), 'position', p);
        set(h_new_cbar, 'yTick', [0:45:135])
        
        extendCbarCdata = 1;
        if extendCbarCdata
            h_cbar_im = get(h_new_cbar, 'children');
            cbar_cdata = get(h_cbar_im, 'cdata');
            cbar_cdata_ext = repmat(cbar_cdata, [1, 6]);
            assert(numel(cbar_cdata_ext) > 1005);
            set(h_cbar_im, 'cdata', cbar_cdata_ext);
            set(h_new_cbar, 'xlim', [0 1]);
        end
        
        
    end

    figure(fig_id)
    addSubplotLetter(subM, subN, i, [],   subSpcM, subSpcN, char('A'+(i)-1), []);
    3;
    %%
    close(h_saved(i))

    
%     figure(fig_id);
%     copyFigureToSubplot(h_saved(i), fig_id, subplotInfo_C);
    3;

    
end
%%
export_fig(fig_id, '-pdf', sprintf('%s%s', figureFolder, 'Figure15_simOriMaps_rgb.pdf'))
export_fig(fig_id, '-pdf', '-cmyk', sprintf('%s%s', figureFolder, 'Figure15_simOriMaps_cmyk.pdf'))

3;