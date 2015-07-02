function plotMeanAllCCs

    gratingType = curGratingType('');

    showGratingTypeInTitle = 0;
    xlim_max_flashed = 100; %#ok<NASGU>
    xlim_max_drifting = 36; %#ok<NASGU>
    
%     plotMode = 'justAA';
    plotMode = 'allLines';

    allCCs_filename = [CatV1Path 'allCCs_vs_rank' curMatchDB('') '.mat'];
    
    S = load(allCCs_filename);
    fn_oe_same = [gratingType '_oe_same'];
    fn_oe_diff = [gratingType '_oe_diff'];
    fn_aa = [gratingType '_aa'];
    
    if ~isfield(S, fn_oe_same),  error('%s field is missing', fn_oe_same); end
    if ~isfield(S, fn_oe_diff),  error('%s field is missing', fn_oe_diff); end
    if ~isfield(S, fn_aa),  error('%s field is missing', fn_aa); end
    
    aa_name = 'All-All';
    oe_same_name = 'Even-Even, Odd-Odd';
    oe_diff_name = 'Even-Odd, Odd-Even';        
    
    
    figure(600+curGratingType); clf;
    hold on; box on;
    legLocation = 'best';
    
    switch plotMode
        case 'justAA',
            
            allWCCs_aa = S.(fn_aa);
            col = 'b';
            face_fill = {'markerfacecolor', col};
    
            xlim_max = eval(['xlim_max_' gratingType]);
            switch gratingType
                case 'flashed', ylims = [-0.05, 0.15];
                case 'drifting', ylims = [-0.02, 0.06];
            end
                
            meanCCs_aa = nanmean(allWCCs_aa, 1);                
            h = plot(meanCCs_aa, [col 'o-'], 'markersize', 3, face_fill{:});
            x_lim2 = min(xlim_max, size(allWCCs_aa, 2));
            xlim([0 x_lim2]);
            ylim(ylims);
                
            drawHorizontalLine(0, 'linestyle', ':')
            xlabel('Stimulus Rank #'); ylabel('Mean CC');
            
            n = size(allWCCs_aa, 1);
            if showGratingTypeInTitle                    
                title(sprintf('%s gratings (N = %d)', titleCase(gratingType), n))
            else
                title(sprintf('(N = %d)', n))
            end
            
            
        case 'allLines',
            plot_types = {'aa', 'oe_diff', 'oe_same'};
            plot_Names = {aa_name, oe_diff_name, oe_same_name};
            allWCCs_C = cellfun(@(s) S.([gratingType '_', s]), plot_types, 'un', 0);  % {S.([gratingType '_', plot_order{2}]), S.([gratingType '_', plot_order{3}]) };
            
            nLines = length(plot_types);
            
            
            for i = 1:nLines
                        
                allWCCs = allWCCs_C{i};
                switch plot_types{i}
                    case 'aa', 
                        col = 'm';                         
                    case 'oe_same', 
                        col = 'b';                        
                    case 'oe_diff',                             
                        col = 'r';                    
                end
                face_fill = {'markerfacecolor', col};
    
                xlim_max = eval(['xlim_max_' gratingType]);
                switch gratingType
                    case 'flashed', ylims = [-0.1, 0.15];
                    case 'drifting', ylims = [-0.04, 0.06];
                end
                
                meanCCs = nanmean(allWCCs, 1);                
                h(i) = plot(meanCCs, [col 'o-'], 'markersize', 3, face_fill{:});
                x_lim2 = min(xlim_max, size(allWCCs, 2));
                xlim([0 x_lim2]);
                ylim(ylims);
                
                if i == 1;
                    drawHorizontalLine(0, 'linestyle', ':')
                    xlabel('Stimulus Rank #'); ylabel('Mean CC');
                    n = size(allWCCs_C{1}, 1);
                    if showGratingTypeInTitle                    
                        title(sprintf('%s gratings (N = %d)', titleCase(gratingType), n))
                    else
                        title(sprintf('(N = %d)', n))
                    end
                    h_tmp = plot(0,0, 0,0, 0,0);
                    legend(plot_Names, 'location', legLocation);
                    pos_final = get(gca, 'position');                    
                    delete(h_tmp);
                end

                legend(h, plot_Names(1:i), 'fontsize', 9, 'location', legLocation);
                set(gca, 'position', pos_final);
                3;
                
            end
               

            
            3;
    end

    3;
                
end