redo_cmps = 0;
redo_stats = 1;

% OE_keep_both = 1;
gratingTypes = {'drifting'};
% gratingTypes = {'flashed', 'drifting'}; 

test_mode = 'minR';
% test_mode = 'maxP';
all_min_r_oe_test = sort([0:.15:.75]);
% all_min_r_oe_test = [0.25, 0.5, 0.75, 0.8];
% all_max_p_oe_test = [0.05, .01, .001];


plotOverlay = 0;
plotDotted = 0;
test_oe_mode = 'hoe';

% plotOverlay = 1;
% plotDotted = 0;
% test_oe_mode = 'hoe';


criteria = struct('n_phases', struct('op', @eq, 'value', 60));

switch test_mode
    case 'minR', vals =  all_min_r_oe_test; set_func = @curMinR_oe;  curMaxP_oe(nan);
    case 'maxP', vals =  all_max_p_oe_test; set_func = @curMaxP_oe;  curMinR_oe(nan);
end
N = length(vals);

% curOEkeepBoth(OE_keep_both);


    
%     for i = 1:length(all_min_r_oe_test)
    for i = 1:N

%         min_r_oe = all_min_r_oe_test(i);
%         curMinR_oe(min_r_oe);

        val_i = vals(i);
        out_str = set_func(val_i, test_oe_mode);


        fprintf('*************(%d/%d) : min_r = %.2f, mode = %s ("%s") ... \n', i, N, val_i, test_oe_mode, out_str)

        if any(strcmp(gratingTypes, 'flashed'))
            curGratingType(1);
            if ~exist(getFileName('comparisons'), 'file') || redo_cmps
                generateGratingComparisonsDatafile;
            end
            if ~exist(getFileName('pairDiffs'), 'file') || ~exist(getFileName('pairStats'), 'file') || redo_stats
                printPhaseTuningComparisons([1 2], criteria);
            end
        end
        
        if any(strcmp(gratingTypes, 'drifting'))
            curGratingType(2);
            if ~exist(getFileName('comparisons'), 'file') || redo_cmps
                generateGratingComparisonsDatafile;
            end
            if ~exist(getFileName('pairDiffs'), 'file') || ~exist(getFileName('pairStats'), 'file') || redo_stats
                printPhaseTuningComparisons([1 2], criteria);
            end
        end
        
    end



%%
Nroe = length(vals);

% measures = {'cc', 'dphi'};
measures = {'cc'};
nGratTypes = length(gratingTypes);
 nMeasures = length(measures);
[meds, meds_ctrl, med_probs, means, stderrs, means_ctrl, mean_probs, KS_stats, KS_probs, nCells, nPairs] = deal( zeros(Nroe, nMeasures, nGratTypes) );
location = 'maxMinFracR';
% location = 'maxR1xR2';
% statsToPlot = {'mean', 'median', 'ks-stat', 'nPairs'};
statsToPlot = {'mean', 'nPairs'};
nStatsPlot = length(statsToPlot);

%%
for r_oe_i = 1:Nroe
%     min_r_oe = all_min_r_oe_test(r_oe_i);
%     curMinR_oe(min_r_oe);
    val_i = vals(r_oe_i);
    set_func(val_i, test_oe_mode);


    for gi = 1:nGratTypes
        pairDiffsFile = getFileName('pairDiffs', [], [], struct('gratingType', gratingTypes{gi}));
        pairStatsFile = getFileName('pairStats', [], [], struct('gratingType', gratingTypes{gi}));
        pairDiffs_S = load(pairDiffsFile);
        pairStats_S = load(pairStatsFile);
    
        for mi = 1:nMeasures
            %%
            pairDiffs_Wcc  = pairDiffs_S.([location '_' measures{mi} '_Wcc']);
            pairDiffs_Wscc = pairDiffs_S.([location '_' measures{mi} '_Wscc']);
            pairStats      = pairStats_S.([location '_' measures{mi}]);

            
%             pairStr = pairDiffs_Wcc.N;
%             nPairs(r_oe_i, mi, gi) = sscanf(pairStr, '%d Pairs');
            nPairs(r_oe_i, mi, gi)    = pairDiffs_Wcc.N_pairs;

            meds(r_oe_i, mi, gi)      = pairDiffs_Wcc.vals_median;
                                        assert(isequal(meds(r_oe_i, mi, gi),  pairStats.median ) );            
            meds_ctrl(r_oe_i, mi, gi) = pairDiffs_Wscc.vals_median;
            med_probs(r_oe_i, mi, gi) = pairStats.medianProb;
            
            means(r_oe_i, mi, gi)     = pairDiffs_Wcc.mean;
                                        assert(isequal(means(r_oe_i, mi, gi), pairStats.mean) );                        
            stderrs(r_oe_i, mi, gi)   = pairDiffs_Wcc.std / sqrt(pairDiffs_Wcc.N_pairs);
            means_ctrl(r_oe_i, mi, gi)= pairDiffs_Wscc.mean;
            mean_probs(r_oe_i, mi, gi)= pairStats.meanProb;
            KS_stats(r_oe_i, mi, gi)  = pairStats.ksStat;            
            KS_probs(r_oe_i, mi, gi)  = pairStats.ksProb;
    
        end
        
    end
    
end

KS_probs(KS_probs == 0) = 1e-4;
med_probs(med_probs == 0) = 1e-4;
mean_probs(mean_probs == 0) = 1e-4;

    %%
    

    switch test_mode
        case 'minR', x_vals = vals; 
        case 'maxP', x_vals = -log10(vals);
    end
    
    nullMeasures = [0, 90];
    xlims = lims(x_vals, .02);
    M = nStatsPlot;
    
    if ~plotOverlay
        h_ax = zeros(M,2);
    end
    [h1, h2] = deal(zeros(M,1));
    %%
    for gi = 1:nGratTypes
        
        for mi = 1:nMeasures
            
            figure(100*gi + mi*10); 
            if ~plotOverlay
                clf;            
            end
            
            sub_i = 1;

            for si = 1:nStatsPlot
            
                if (si > 1) && plotOverlay
                    subplotGap(M,1,si);
                end
                if strcmp(statsToPlot{si}, 'mean')
                    if ~plotOverlay
                        [h_ax(si,:), h1(si), h2(si)] = plotyy(x_vals, means(:,mi,gi), x_vals, -log10(mean_probs(:,mi, gi) )); 
                        set(h_ax(si,:), 'nextPlot', 'add');
                    else
                        h1(si) = plot(h_ax(si,1), x_vals, means(:,mi,gi));
                        h2(si) = plot(h_ax(si,2), x_vals, -log10(mean_probs(:,mi, gi) ));                         
                    end
%                     h3(si) = plot(h_ax(si,1), x_vals, means_ctrl(:,mi,gi), 'o-', 'color', [0 .7 0]);
                    h4(si) = errorbar(h_ax(si,1), x_vals, means(:,mi,gi), stderrs(:,mi,gi));
                    if plotDotted
                        set([h1(si), h2(si)], 'linestyle', ':', 'linewidth', 1);
                    else
                        set([h1(si), h2(si)], 'linestyle', '-', 'linewidth', 2);
                    end
                    
                    title(sprintf('%s gratings; %s', gratingTypes{gi}, measures{mi}));
                    ylabel('mean');
                    set(h_ax(si,:), 'xlim', xlims);
                    %%
                    mean_ylims1 = [-0.01, 0.13];  rel0_1 = abs(mean_ylims1(1))/diff(mean_ylims1);
                    ylims2(2) = 4.2;
                    ylims2(1) = -ylims2(2)*rel0_1/(1-rel0_1);
                    mean_ylims2 = [ylims2];
%                     set(h_ax(si,1), 'ylim', lims([means(:,mi,gi) + stderrs(:,mi,gi); means_ctrl(:,mi,gi)], .05) );
                    set(h_ax(si,1), 'ylim', mean_ylims1, 'ytick', [0:.02:.1]);
                    set(h_ax(si,2), 'ylim', mean_ylims2, 'ytick', [0:1:4]);
                    
                    set(h_ax(si,1), 'box', 'off')
                    set(h_ax(si,2), 'box', 'off')
                    if ~plotOverlay
                        ylabel(h_ax(si,2), '-log_{10} p-value');
                        drawHorizontalLine(nullMeasures(mi), 'linestyle', ':', 'color', 'k');
                    end
                    3;
                
                elseif strcmp(statsToPlot{si}, 'median')                        
                    [h_ax(si,:), h1(si), h2(si)] = plotyy(x_vals, meds(:,mi,gi), x_vals, -log10(med_probs(:,mi, gi) )); 
                    set(h_ax(si,1), 'nextPlot', 'add');
                    h3(si) = plot(h_ax(si,1), x_vals, meds_ctrl(:,mi,gi), 'o-', 'color', [0 .7 0]);
                    drawHorizontalLine(nullMeasures(mi), 'linestyle', ':', 'color', 'k');
                    ylabel('median');
                    set(h_ax(si,:), 'xlim', xlims);
                    set(h_ax(si,1), 'ylim', lims([meds(:,mi,gi); meds_ctrl(:,mi,gi)], .05) );
                    set(h_ax(si,:), 'outerposition', getNormPosition(M, 1, 2, 1));

                elseif strcmp(statsToPlot{si}, 'ks-stat')
            
                    [h_ax(si,:), h1(si), h2(si)] = plotyy(x_vals, KS_stats(:,mi,gi), x_vals, -log10(KS_probs(:,mi, gi) ));             
        %             drawHorizontalLine(nullMeasures(mi), 'linestyle', ':', 'color', 'k');
                    ylabel('KS stat');
                    set(h_ax(si,:), 'xlim', xlims);            
                    set(h_ax(si,:), 'outerposition', getNormPosition(M, 1, 3, 1));
                elseif strcmp(statsToPlot{si}, 'nPairs')

                    h_ax(si,1) = gca;
                    h1(si) = plot(x_vals, nPairs(:,mi,gi) );             
        %             drawHorizontalLine(nullMeasures(mi), 'linestyle', ':', 'color', 'k');
                    ylabel('N Pairs');
                    set(h_ax(si,1), 'xlim', xlims);
            
                end
                if si == nStatsPlot
%                     xlabel('min r ')
                    %%
                    xlabel('Minimum required CC between even/odd trials')
                end
                if ~plotOverlay
%                     set(h_ax(si,:), 'outerposition', getNormPosition(M, 1, si, 1));
                end
                
            end
            3;
            
            
            
            col1 = [0 0 1]; col2 = [1 0 0];
            set(h1, 'color', col1, 'marker', '.')
            set(h_ax(:,1), 'ycolor', col1)
            set(nonzeros(h2), 'color', col2, 'marker', 's')            
            set(nonzeros(h_ax(:,2)), 'ycolor', col2)
            
            doLegend = 1;
            if doLegend
                %%
                legend([h1(si), h2(si)], {'mean CC', '–log(p-value)'}, 'location', 'NW')
            end
            3;
        end
    end
    3;
    
%%    






