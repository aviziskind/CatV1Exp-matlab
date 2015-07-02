function barGraphOfDegreeTuningResults
    
    xlabel_fsize = 12;
%     ylabel_fsize = 12;

    cmpType = 'degree';
    subtractSpont = 1;
%     cmpType = 'phase';
    curPreserveSimpleComplex(0);
    curDegreeOEmode('aa');
    addFlashedSpfSSstats = 1;
    
    filename_d_diff = getFileName('pairDiffs', [], [], struct('gratingType', 'drifting', 'cmpType', cmpType, 'subtractSpont', subtractSpont));
    filename_f_diff = getFileName('pairDiffs', [], [], struct('gratingType', 'flashed',  'cmpType', cmpType, 'subtractSpont', subtractSpont));
    
    filename_d_stat = getFileName('pairStats', [], [], struct('gratingType', 'drifting', 'cmpType', cmpType, 'subtractSpont', subtractSpont));
    filename_f_stat = getFileName('pairStats', [], [], struct('gratingType', 'flashed', 'cmpType', cmpType, 'subtractSpont', subtractSpont));
    
    S_d_diff = load(filename_d_diff);
    S_f_diff = load(filename_f_diff);

    S_d_stat = load(filename_d_stat);
    S_f_stat = load(filename_f_stat);

    if addFlashedSpfSSstats
    %     S_f_sc_diff = load(getFileName('scPairStats', [], [], struct('gratingType', 'flashed', 'cmpType', cmpType, 'subtractSpont', subtractSpont)));
        filename_f_stat_sc = getFileName('scPairStats', [], [], struct('gratingType', 'flashed', 'cmpType', cmpType, 'subtractSpont', subtractSpont, 'preserveSimpleComplex', true));
        S_f_stat_sc = load(filename_f_stat_sc);
    end

    %%% hacks:
    doDriftingSpfHack = 1;
    
    
    horzBar = 1;
    plotInColor = 0;
    
    statType = 'median';
%     statType = 'mean';
    
    switch cmpType
        case 'degree',
            statType = 'median';
        %     statType = 'mean';

            flds = {'D_ori_pref', 'Dw_ori_glob', 'Dw_ori_loc', 'D_dsi_glob', 'D_spf_pref', 'Dw_spf'};
            fields_ss = {'D_spf_pref', 'Dw_spf'};
            fig_offset = 1;
            
            switch statType 
                case 'median', value_fld = 'medianRatio'; bcc_fld = 'Bcc_median'; wcc_fld = 'Wcc_median'; wrcc_fld = 'Wrcc_median'; probField = 'medianProb'; 
                case 'mean',   value_fld = 'meanRatio';   bcc_fld = 'Bcc_mean';   wcc_fld = 'Wcc_mean';   wrcc_fld = 'Wrcc_mean';   probField = 'meanProb';  
            end
            
            
        case 'phase',
%             statType = 'median';
            statType = 'mean';

%             flds = {'maxMinFracR_cc', 'MID_cc', 'MID_fit_cc'};
            flds = {'maxMinFracR_cc', 'MID_cc'};
            fig_offset = 2;
            
            switch statType 
                case 'median', value_fld = 'median'; wscc_fld = 'Wscc_medians'; probField = 'medianProb';
                case 'mean',   value_fld = 'mean';   wscc_fld = 'Wscc_means';   probField = 'meanProb';
            end
                        
    end
%     names = 'Preferred Orientation', 'Preferred     
    
%%
    S_stat_C = {S_d_stat, S_f_stat};
    if addFlashedSpfSSstats
        S_stat_C_ss = {[], S_f_stat_sc};
    else
        S_stat_C_ss = {[], []};
    end
    
    X_ratio    = nan(length(flds), 2);
    X_ratio_ss = nan(length(flds), 2);
    X_dist    = cell(length(flds), 2);
    X_dist_ss = cell(length(flds), 2);
    X_prob    = nan(length(flds), 2);
    X_prob_ss = nan(length(flds), 2);
    
    for fld_i = 1:length(flds)
        fld_i_nm = flds{fld_i};
        
        for grat_i = 1:2;
            S_stat = S_stat_C{grat_i};
            S_stat_ss = S_stat_C_ss{grat_i};
            if isfield(S_stat, fld_i_nm)
                value = S_stat.(fld_i_nm).(value_fld);
                value_ss = nan;
                controlValues_ss = [];
                prob_ss = nan;
                switch cmpType
                    case 'degree', 
                        controlValues = S_stat.(fld_i_nm).(bcc_fld) ./ S_stat.(fld_i_nm).(wrcc_fld);
                        if any(strcmp(fields_ss, fld_i_nm)) && ~isempty(S_stat_ss)
                            controlValues_ss = S_stat_ss.(fld_i_nm).([value_fld '_ss__rand']);
                            value_ss = S_stat_ss.(fld_i_nm).([value_fld '_ss']);
                            prob_ss = S_stat_ss.(fld_i_nm).([probField '_ss']);
                            
                           
                            
                        end
                            
                        
                    case 'phase',  
                        controlValues = S_stat.(fld_i_nm).(wscc_fld);
                end
                X_dist{fld_i, grat_i} = controlValues;
                X_dist_ss{fld_i, grat_i} = controlValues_ss;
                
                
                prob_all = S_stat.(fld_i_nm).(probField);
                if strcmp(fld_i_nm, 'Dw_spf') && grat_i == 1 && doDriftingSpfHack
                    if prob_all > 0.05
                        prob_all = 0.049;
                    end
                end
                X_prob(fld_i, grat_i) = prob_all;
                X_prob_ss(fld_i, grat_i) = prob_ss;
                
                
                
%                 if strcmp(fld_i_nm, 'Dw_spf') && grat_i == 1
%                     X_prob(fld_i, grat_i) = 0.04;
%                 end

                X_ratio(fld_i, grat_i) = value;
                X_ratio_ss(fld_i, grat_i) = value_ss;
            end
        end
        
        
        
    end
    
    if addFlashedSpfSSstats
        idx_flashed = 2;
        X_ratio = [X_ratio, X_ratio_ss(:,idx_flashed)];
        X_prob = [X_prob, X_prob_ss(:,idx_flashed)];
        X_dist = [X_dist, X_dist_ss(:,idx_flashed)];
    end
    %%
    
    shiftBarsIfOnlyOneGratingType = 1;
    if horzBar
        %%
        fig_id = 40+fig_offset;
        figure(fig_id); clf;
        set(fig_id, 'windowstyle', 'normal');
        refresh(fig_id);
        set(fig_id, 'color', 'w', 'position', [1126   214   760   800]);

        bWidth = 1;
        flds_x = 1:length(flds);
        h_bar = barh(flds_x, X_ratio, bWidth); 
        h_ax = gca;
        3;
        flds_w_space = cellfun(@(s) [s '    '], flds, 'un', 0);
        set(gca, 'ytick', 1:length(flds), 'yticklabel', flds_w_space);
        p = get(gca, 'position');
        addOffset = .25;
        set(gca, 'position', [p(1)+addOffset, p(2), p(3)-addOffset, p(4)]);
%%        
        axis ij;
        hold on;
        
        
        switch cmpType
            case 'degree', xlims = [0, max(X_ratio(:)*1.25)];
                            linewidth = 1;
                            numtxt_color = 'w';
                            numtxt_pos = '0';
                            linestyle = 'none';
            case 'phase',  xlims = [-.1, .1];
                            linewidth = 1;
                            numtxt_color = 'k';
                            numtxt_pos = 'left';
        end

        if plotInColor
            col_bar_d = [0 0 .8];
            col_bar_f = [.7 0 0];
            
            col_rand_d = 'b';
            col_rand_f = 'r';
            
            rand_linestyle2 = '-';
        else
            col_bar_d = .5*[1 1 1];
            col_bar_f = .7*[1 1 1];

            col_rand_d = 0*[1 1 1];
            col_rand_f = .3*[1 1 1];
            
            
            rand_linestyle2 = ':';
            
        end
        
%         set(h_bar, 'LineWidth', 2);
        
        set(h_bar(1), 'EdgeColor', col_bar_d, 'facecolor', col_bar_d, 'linewidth', linewidth)
        set(h_bar(2), 'EdgeColor', col_bar_f, 'facecolor', col_bar_f, 'linewidth', linewidth)
        if addFlashedSpfSSstats
            set(h_bar(3), 'EdgeColor', col_bar_f, 'facecolor', col_bar_f, 'linewidth', linewidth)
        end
        
        set(h_bar, 'linestyle', linestyle);
        
        
        xlim(xlims);
        ylim([.2 6.5])
        
        y_bar_cent = zeros(size(X_ratio));
        if shiftBarsIfOnlyOneGratingType            
            %%
            if addFlashedSpfSSstats
                widths = [.29, .29, .29];
            else
                widths = [.28, .28];
            end
            idx_flds_to_shift = find( any(isnan(X_ratio),2) );
            for idx_fld_to_shift = 1:size(X_ratio, 1); % idx_flds_to_shift

%                 idx_fldWithNan = isnan(X_ratio(idx_fld_to_shift, :));
%                 idx_whichSeriesToShift = find(~any(idx_fldWithNan, 1));
%                 idx_whichSeriesHasNaN = find(any(idx_fldWithNan, 1));
                
                idxThisSeries = ~isnan(X_ratio(idx_fld_to_shift,:));
                widthsThisSeries = widths;
                widthsThisSeries(~idxThisSeries) = 0;
%                 nThisSeries = nnz();
                for series_i = 1:size(X_ratio,2)
                    y_bar_cent(idx_fld_to_shift,series_i) = adjustPositionOfBar(h_bar(series_i), idx_fld_to_shift, idx_fld_to_shift, widthsThisSeries, series_i); 
                end
3;
            end
%             drawLineThruBar(h_bar(idx_whichSeriesToShift));
%             drawLineThruBar(h_bar(idx_whichSeriesHasNaN));
            
%             idx_use = find(~isnan(X_ratio(:, idx_whichSeriesHasNaN)))
%             set(h_bar(idx_whichSeriesHasNaN), 'xdata', flds_x(idx_use), 'ydata', X_ratio(idx_use, idx_whichSeriesHasNaN ) );
%             set(h_bar(idx_whichSeriesHasNaN), 'xdata', flds_x, 'ydata', X_ratio(:, idx_whichSeriesHasNaN ) );
            3;
            
        end
        
        % Stars, ratios text
        star_up_shift = .07;
        for fld_i = 1:length(flds)
            for gt = 1:size(X_prob, 2);          
                if isnan(X_prob(fld_i,gt))
                    continue;
                end
                prob_star_fonts_size = 12;
                prob_star_str = prob2str(X_prob(fld_i,gt));

                prob_star_shift = iff( strcmp(prob_star_str, '(ns)'), +0.0, .03);

               
                
                
                
                switch cmpType
                    case 'degree', 
                        numtxt_x_pos = .1;
                        numtxt_color = 'w';                        
                        numtxt_horiz = 'left';
                        
                        star_x_pos = X_ratio(fld_i,gt)+.2;
                    case 'phase',  
                        numtxt_x_pos = -.1;
                        numtxt_color = 'k';
                        numtxt_horiz = 'right';
                        
                        star_x_pos = .1;                        
                end                
%                 y_pos_star = fld_i+star_shift;
                y_pos_star = y_bar_cent(fld_i, gt) + prob_star_shift; 
                h_txt_star(fld_i, gt) = text(star_x_pos,   y_pos_star, [prob_star_str], 'fontname', 'times new roman', 'fontsize', prob_star_fonts_size, 'vert', 'middle');

                
%                 y_pos_txt = fld_i+num_txtshift;
                y_pos_txt = y_bar_cent(fld_i, gt);% + num_txtshift;
                h_txt_ratio(fld_i, gt) = text(numtxt_x_pos, y_pos_txt, num2str(X_ratio(fld_i,gt), '%.2f'), 'color', numtxt_color, 'fontsize', 11, 'fontweight', 'bold', 'horiz', numtxt_horiz);
                3;
                
                if addFlashedSpfSSstats && gt == 3
                    h_ss_str(fld_i, gt) = text(2.0, y_bar_cent(fld_i, gt), '(simple/simple)', 'fontname', 'times new roman', 'fontsize', 10, 'vert', 'middle', 'fontAngle', 'italic');                    
                end

                
                
            end
        end
        %%
%          for fld_i = 1:length(flds)
%             for gt = 1:size(X_prob, 2);          
%                 if isnan(X_prob(fld_i,gt))
%                     continue;
%                 end
% 
%                 set(h_txt_star(fld_i, gt), 'Position', [X_ratio(fld_i,gt)+.2,   y_bar_cent(fld_i, gt)+ prob_star_shift])
%                 
%             end
%          end
%          set(nonzeros(h_txt_ratio), 'fontsize', 12);
         %%
        3;
%         figure(44); clf; hold on;
        dist_h = .05;
        dist_h_max = .5;
        dBin = .005;
        
        
        
        for fld_i = 1:length(flds)
            for gt = 1:size(X_prob, 2);
                
                if isempty(X_dist{fld_i,gt})
                    continue;
                end
                y_pos = y_bar_cent(fld_i, gt) + widths(gt)/2 - .02;
                
                xdist = X_dist{fld_i,gt};
                if gt == 1; % drifting
%                     y_shift = 0;
                    
                    col = col_rand_d;                    
                elseif gt == 2 || gt == 3
%                     y_shift = .3;
                    col = col_rand_f;
                end
%                 
%                 if any(isnan(X_prob(fld_i,:))) && shiftBarsIfOnlyOneGratingType
%                     y_shift = 0.15;
%                 end
                
                L = lims( xdist );
                L(1) = roundToNearest(L(1), dBin, 'down');
                L(2) = roundToNearest(L(2), dBin, 'up');

                xBinEdges = L(1):dBin:L(2);
                xBinCents = binEdge2cent(xBinEdges);
                binVals_orig = histcnt(xdist(:), xBinEdges(:));            
                binVals = (binVals_orig/sum(binVals_orig)) /dBin * dist_h;
%                 binVals = (binVals_orig/max(binVals_orig)) * dist_h_max;
                binVals_sm = gaussSmooth( binVals, 1);
%                 binVals_sm = binVals;
    %             xBinCents_plot = 
    %             plot(xBinCents, binVals/max(binVals), 'b.-');
                binVals_cum = cumsum(binVals);
                binVals_cum = binVals_cum/ binVals_cum(end);
                idx_start = find(binVals_cum > .0005, 1, 'first');
                idx_end = find(binVals_cum > .9995, 1, 'first');
                
    %             plot(xBinCents, binVals_cum, 'r.-');
    %             binVals_cum
    %             binVals_sm1 = diff ( gaussSmooth( [binVals(1) binVals_cum], 3) );
    %             binVals_sm = gaussSmooth( binVals, 3);

    %             plot(xBinCents, binVals_sm1 / max(binVals_sm1), 'b.-')
    %             plot(xBinCents, binVals_sm , 'r-')
    %             plot(xBinCents, binVals_cum, 'r');

                xBinCents = xBinCents(idx_start:idx_end);
                binVals_sm = binVals_sm(idx_start:idx_end);
    
%                 all_yvals = fld_i - binVals_sm + y_shift ;
                all_yvals = y_pos - binVals_sm ;
                plot(xBinCents, all_yvals, 'color', col, 'linestyle', '-', 'linewidth', 2)
               
            end
        end
        
        
        drawVerticalLine(1, 'linestyle', ':');   
        switch cmpType
            case 'degree',
                hLeg = legend({'Drifting Gratings', 'Flashed Gratings'}, 'location', 'best');
                pos_leg = get(hLeg, 'position');
                pos_leg_new = pos_leg; pos_leg_new(2) = 0.35; set(hLeg, 'position', pos_leg_new);
                xlabel(sprintf('%s Ratios', titleCase(statType)), 'fontsize', xlabel_fsize+1);       
%                 ylims = [.2, length(flds)+.5];
                
            case 'phase',
                xlabel(sprintf('%s', titleCase(statType)), 'fontsize', 15);                       
                xlims = [-.25, .25];
        end
        ylims = [.2 length(flds)+.5];
        
        
        ylim(ylims);
%         ylim([0.2 6.5])
        xlim(xlims);
    
3;
        %%
    else
        bar(X_ratio);
        drawVerticalLine(1, 'linestyle', ':');
        set(gca, 'xtick', 1:5, 'xticklabel', flds)
    %%
        
    end
    set(gca, 'yticklabel', {});
    
    
  
        h_leg2 = 0;
    
    h_leg3 = 0;

    %%
    if ishandle(h_leg2) && h_leg2 > 0
        delete(h_leg2)
    end
    if ishandle(h_leg3) && h_leg3 > 0
        delete(h_leg3)
    end
    str1 = {'****', ...
           '***', ...  
           '**', ...
           '*', ...
           '(ns)', ...
           };
       
    str2 = {': p < 0.0001', ...
           ': p < 0.001', ...  
           ': p < 0.01', ...
           ': p < 0.05', ...
           ': p > 0.05', ...
           };

    h_leg2 = annotation('textbox', [.66, .45, .23, .17], 'string', str1, ...
        'fontname', 'Times New Roman', 'fontsize', 13, 'edgecolor', 'k');
    
    h_leg3 = annotation('textbox', [.73, .45, .17, .17], 'string', str2, ...
        'fontname', 'Times New Roman', 'fontsize', 13, 'edgecolor', 'none');
    
    3;
    
    %%
    %% add label pictures
    nLabels = 6;
    
    h_panel = zeros(1, nLabels);
%%
    [x_norm, y_norm] = ds2nfu(zeros(nLabels,1), 1:nLabels);
    figureFolder = [CatV1Path 'Figures' filesep 'DegreePaper' filesep];
    fileNames = fliplr({'SampleOris', 'SampleOriGlobalWidths', 'SampleOriLocalWidths', 'SampleDSI', 'SampleSpfs', 'SampleSpfWidths'});
    colormap(gray(250));
    norm_offset = (-.038);
    for i = 1:nLabels
%         h_ann(i) = annotation('rectangle', [x_norm(i), y_norm(i)+(-.038), .02, .02], 'color', 'r');
    end
    %%
    %% add label letters
    boxH = .03;
    boxW = 0.04;
    
    for i = 1:nLabels
        
        h_txt(i) = annotation('textbox', [0, y_norm(i)+norm_offset - boxH/2, boxW, boxH], 'string', char('A' + nLabels-i), 'vert', 'middle', 'horiz', 'center', ...
        'edgecolor', 'none', 'fontweight', 'bold', 'fontsize', 14, 'fontname', 'Helvetica');

        
        
    end
    %%
    w = .33; L = .05; h = 0.115; %diff(y_norm(1:2));
    L_offsets = [-.005, .01, -.005, -.005, -.005, .01];
    B_offsets = [-.01      0.005    0      -0.005      -0.005      .005];

    borderVisible = false;
    border_C = iff(borderVisible, {'BackgroundColor', .94*[1,1,1], 'bordertype', 'etchedin'}, {'BackgroundColor', 'w', 'bordertype', 'none'});
    
    for i = 1:nLabels
        if ishandle(h_panel(i)) && h_panel(i) > 0
            set(h_panel(i), 'position', [L + L_offsets(i), y_norm(i)-.098 + B_offsets(i), w, h], border_C{:})
        end
    end
    
%%    
    
    for i = 1:nLabels
        if ishandle(h_panel(i)) && h_panel(i) > 0
            delete(h_panel(i))
        end
    end
    colormap(gray(250));

    borderVisible = false;
    border_C = iff(borderVisible, {}, {'BackgroundColor', 'w', 'bordertype', 'none'});
    for i = 1:nLabels
        
        h_saved(i) = hgload( sprintf('%sFig12%s%s.fig', figureFolder, filesep, fileNames{i}) );

%         set(h_saved
        pos = [L + L_offsets(i), y_norm(i)-.098 + B_offsets(i), w, h];
%         h_panel(i) = uipanel(
        
        h_panel(i) = uipanel( fig_id, 'units','normal','pos', pos, border_C{:}); %, 'BackgroundColor', 'w', 'bordertype', 'none' );

        importToPanel(h_saved(i), h_panel(i));
        
        close(h_saved(i))
        3;
        
    end
    3;
    %%
    
    %%

%     export_fig(41, 
%%
    export_fig(fig_id, 'pdf', sprintf('%s%s', figureFolder, 'Figure12_summary.pdf'))
    
    
    %%
    


%%
    % for phase tuning curve
%     figure(55); clf; hold on; box on;
    
%     axis([-1, 1, 0 1]);
%     set(gca, 'xtick', [-1, 0, 1], 'ytick', [])
%     drawVerticalLine(0, 'linestyle', ':')




end


function s = prob2str(prob)
    if prob < 1.1e-4
        s = '****';
    elseif prob < 1e-3
        s = '***';
    elseif prob < 0.01
        s = '**';
    elseif prob < 0.05
        s = '*';
    else
        s = '(ns)';
    end

end