% function slides
%     load indivCells_movie_fg;
%

% makeSlides = 1012;
% makeSlides = 1;
makeSlides = 'SetsOfGratings';

global showWorking;

saveToFile = true;
figureFolder = [CatV1Path 'Figures' filesep 'DegreePaper' filesep];

if any(strcmp(makeSlides, 'SetsOfGratings'))   % gratings for OSP axes 

    nPixSide = 1024;
    nPixSide_disp = 64;

    dims = [nPixSide nPixSide];
    dims_show = [nPixSide_disp, nPixSide_disp];
    oris = [10:65:175];            ori_forSpf = 115; ori_forPhs = 45;
%     spfs = nPixSide*[1/2, 1/3, 1/4, 1/8]; spf_forOri = nPixSide/2; spf_forPhs = nPixSide;
    spfs = nPixSide*[1/2, 1/4, 1/8]; spf_forOri = nPixSide/2; spf_forPhs = nPixSide;
    phis = [0:90:270];             phs_forOri = 0;   phs_forSpf = 0; 
%     params = struct('O', oris, 'S', spfs, 'P', phis);
    params = struct('O', {{oris, spf_forOri*ones(size(oris)), phs_forOri*ones(size(oris))}}, ...
                    'S', {{ori_forSpf*ones(size(spfs)), spfs, phs_forSpf*ones(size(spfs))}}, ...
                    'P', {{ori_forPhs*ones(size(phis)), spf_forPhs*ones(size(phis)), phis}});
    subplotArrangement = {'horiz', 'horiz', 'stacked'};
    doDirectionArrows = false;
    spacingFactor = 1.2;
    LB_0 = [10, 10];
    doAxes = [1, 2];

            getSubplotArrangement = @(arr, i) switchh(arr, {'vert', 'horiz', 'stacked'}, {[20 -140] + [0, i*nPixSide_disp*spacingFactor], ...
                                                                                          [-140 10] + [i*nPixSide_disp*spacingFactor, 0 ], ...
                                                                                          [-80 -20] + [i*nPixSide_disp*rx, i*nPixSide_disp*ry] });
            getSubplotArrangement = @(arr, i) switchh(arr, {'vert', 'horiz'}, {LB_0 + [0, (i-1)*nPixSide_disp*spacingFactor], ...
                                                                               LB_0 + [(i-1)*nPixSide_disp*spacingFactor, 0 ], ...
                                                                           });
                                                                       
                                                                       
    axOrder = 'SOP';

    if any(doAxes == 1)
        figure(1); clf; hold on;
        [oris, spfs, phis] = elements(params.(axOrder(1)));
        n = length(oris);
        for i = 1:n;
            ax1(i) = subplot(n, 1, i); %#ok<AGROW>
%             ax1(i) = subplot(1, n, i); %#ok<AGROW>
        end         
        for i = 1:n
            axes(ax1(i));            
            imagesc(generateGratingFrame(dims, oris(i), spfs(i), phis(i)));
            
            colormap('gray'); axis equal tight; set(gca, 'xtick', [], 'ytick', []);  
            LB = getSubplotArrangement('horiz', i);
            
            set(ax1(i), 'Units', 'Pixels', 'Position', [LB, dims_show]);        
        end
        %%
        w = LB_0(1)*2 + n*nPixSide_disp + (n-1)*(spacingFactor-1)*nPixSide_disp;
        h = LB_0(2)*2 + nPixSide_disp;
        set(1, 'position', [1000, 250, w, h], 'color', 'w')
        
        %%
        if saveToFile
            %%
            set(1, 'windowstyle', 'normal');
            refresh;
%             set(1, 'color', 'w', 'position', )
%             hgsave(figureFolder
        end
        
    end
    
    if any(doAxes == 2)
        figure(2); clf;
        [oris, spfs, phis] = elements(params.(axOrder(2)));
        n = length(oris);
        for i = 1:n
            ax2(i) = subplot(n, 1, i); %#ok<AGROW>
        end       
        for i = 1:n
            axes(ax2(i));       
            imagesc(generateGratingFrame(dims, oris(i), spfs(i), phis(i) ));
            
            LB = getSubplotArrangement('horiz', i);

            colormap('gray'); axis equal tight; set(gca, 'xtick', [], 'ytick', []);      
            set(ax2(i), 'Units', 'Pixels', 'Position', [LB, dims_show]);        
            
            if doDirectionArrows
                 X0 = [nPixSide_disp*(.5); nPixSide_disp*(.5)];
                 r = nPixSide_disp*(.3);
                 theta = deg2rad( oris(i) )+pi/2;

                 Xa = X0 - r* [cos(theta); sin(theta)];
                 Xb = X0 + r* [cos(theta); sin(theta)];
                 U = Xb-Xa;
                 hold on;
                 h = quiver(Xa(1), Xa(2), U(1), U(2), 'color', 'r', 'linewidth', 10, 'maxheadsize', 58);
                 hold off;
            end
            
        end
        
        
        w = LB_0(1)*2 + n*nPixSide_disp + (n-1)*(spacingFactor-1)*nPixSide_disp;
        h = LB_0(2)*2 + nPixSide_disp;
        set(2, 'position', [1000, 350, w, h], 'color', 'w', 'windowstyle', 'normal')
        3;
    end
    
    % axis 3 (z): orientation
    if any(doAxes == 3)
        [oris, spfs, phis] = elements(params.(axOrder(3)));
        n = length(oris);
        figure(3); clf; hold on;
        for i = 1:n;
            ax3(i) = subplot(n, 1, i); %#ok<AGROW>
        end    
        for i = n:-1:1  % so that axes in lower left are in the front
            axes(ax3(i));
            imagesc( generateGratingFrame(dims, oris(i), spfs(i), phis(i)) );
            LB = getSubplotArrangement('horiz', i);

            colormap('gray'); axis equal tight; set(gca, 'xtick', [], 'ytick', []);      
            rx = .8;
            ry = (.8)/3;
            set(ax3(i), 'Units', 'Pixels', 'Position', [LB, dims]);
            if doDirectionArrows
                 X0 = [nPixSide_disp*(.5); nPixSide_disp*(.5)];
                 r = nPixSide_disp*(.3);
                 theta = deg2rad( oris(i) )+pi/2;

                 Xa = X0 - r* [cos(theta); sin(theta)];
                 Xb = X0 + r* [cos(theta); sin(theta)];
                 U = Xb-Xa;
                 hold on;
                 h = quiver(Xa(1), Xa(2), U(1), U(2), 'color', 'g', 'linewidth', 7, 'maxheadsize', 12);
                 hold off;
            end
            
        end
        
            
            
    end
    
    
end


sampleFigWH = [250, 85];
if any(makeSlides == 101) % illustration of orientation / spatial frequency tuning widths
    
    
    %% show two gaussians with different widths    
    doGrayScale = 1;
    if doGrayScale
        col1 = [0 0 0]; col2 = [.6 .6 .6];
        linestyle1 = '-'; linestyle2 = '--';
    else
        col1 = 'b'; col2 = 'r';
        linestyle1 = '-'; linestyle2 = '-';
    end
    
    oris_deg = -90:90;
    sig1 = 15;
    sig2 = 35;
    off1 = .09;
    off2 = .18;
    expfun = @(x, sig, offset) exp(-x.^2/(2*sig^2)) + offset;
    g1 = expfun(oris_deg, sig1, off1);
    g2 = expfun(oris_deg, sig2, off2);
    %%
    xticks = -60:30:60;
%     xticks = -90:15:90;
    xlims = [-90, 90];
    
    nSub = 7;
    figure(101); clf;    
    set(101, 'windowStyle', 'normal', 'position', [1200, 550 sampleFigWH], 'color', 'w')
    h_arrow = [];

    h_ax1 = subplotGap(nSub,1,1:nSub-1);
    h_l = plot(oris_deg, g1, oris_deg, g2);    
    set(h_l(1), 'color', col1, 'linestyle', linestyle1);
    set(h_l(2), 'color', col2, 'linestyle', linestyle2);
    set(h_l, 'linewidth', 3);
    set(h_ax1, 'ytick', []);
    xlim(xlims);      
    ylim([0, max([g1, g2])*1.1])
    hold on;
    
    
    %%
    arrow_type = 'local';
%     arrow_type = 'global';
    switch arrow_type
        case 'global',
            xLR = [-45, 45];
            idx = binarySearch(oris_deg, xLR);
            yLR = g1(idx)+.25;
            
        case 'local',
            %%
            xLR = [-sig1, sig1];
            idx = binarySearch(oris_deg, xLR);
            yLR = g1(idx) - .15;
            
    end
    
% %      hh = plot(xLR, yLR, 'ro-');
%     [x_norm, y_norm] = ds2nfu(xLR, yLR);
%     h_arrow = annotation('doublearrow', x_norm, y_norm);
%     arrowSize = 15;
%     arrowWidth = 5;
%     %             set(h_arrow, 'linewidth', 3, 'head1Style', 'cback1')
%     set(h_arrow, 'linewidth', arrowWidth, 'head1Style', 'plain', 'head2Style', 'plain', ...
%         'head1Length', arrowSize, 'head1Width', arrowSize, 'head2Length', arrowSize, 'head2Width', arrowSize);
%     
%     h_arrow = arrowh(xLR, yLR, 'r', [500, 300], 100);
    
    for i = 1:length(h_arrow)
        if ishandle(h_arrow(i))
            delete(h_arrow(i));
        end
    end

    arrowSize = 300;
    arrowHeadWidth = 600;
    arrowLineWidths = 300;
    h_arrow(1) = plot(xLR*.9, yLR, 'k-', 'linewidth', 3);
    h_arrow(2) = arrowh(xLR, yLR, 'k', [arrowSize, arrowHeadWidth], 100);
    h_arrow(3) = arrowh(fliplr(xLR), yLR, 'k', [arrowSize, arrowHeadWidth], 100);

%%
    
    
    %%
    set(gca, 'xtick', xticks, 'xticklabel', {''}, 'ytick', [])
    

    %%
    % orientation bars for x-axis.
%     figure(102);    
    p1_norm = get(h_ax1, 'position');
  
    h_ax2_offset = -.25;
    h_ax2_height = .25;
    p2 = [p1_norm(1), p1_norm(2) + h_ax2_offset, p1_norm(3), h_ax2_height];

    h_ax2 = axes('position', p2);
    set(h_ax2, 'units', 'pixels', 'xlim', xlims);
    p2_pix = get(h_ax2, 'position');
    %%
    ar = p2_pix(3)/p2_pix(4);
    
    h_width = 5;
    h_height = h_width * ar;
    spacing_factor = 2;
    
    yrange = diff(xlims)/ar;
    ylims = yrange*[-.5, .5];
    set(h_ax2, 'ylim', ylims, 'visible', 'off');
    %%
    
%     set(h_ax2, 'xlim', xlims, 'ylim', ylims2); 
    
    %%
%     subplotGap(nSub-1, 1, nSub-1);
%%
    L = (yrange/2) * .8;
%     axis equal;
%     ylims = [-L*1.5, L*1.5];
%     axis([xlims, ylims]);    
    h = zeros(size(xticks));
    for i = 1:length(xticks);
        xc = xticks(i);
        theta = deg2rad(xticks(i));
        X = [xc - cos(theta)*L, xc + cos(theta)*L];
        Y = [-sin(theta)*L, sin(theta)*L];
        h(i) = line(X, Y);        
    end
    %%
    set(h, 'color', 'k', 'lineWidth', 3);
%     set(gca, 'xtick', [], 'ytick', []);
    box on;
    
    %%
    
    export_fig(101, 'pdf', [figureFolder 'Fig12\OriGlobalWidthSample.pdf'])
    %%
    
    3;
end

if any(makeSlides == 1011)
    %% show two skewed lognormals with different widths    

    spfs = 0.1:.01:4;
    w1 = 0.6;
    w2 = 1;
    fopt = 1;
    sln = @(x, w) skewLogNormal(x, fopt, 1, 0, w, 0);
        
    sln1 = sln(spfs, w1);
    sln2 = sln(spfs, w2);
    
    nTicks = 4;
    xticks = linspace(0.5, 3.5, nTicks);
    xticks_cpd = logspace(log10(0.5), log10(2.5), nTicks);
    
    nSub = 7;

    figure(103); clf;   
    
    set(103, 'windowStyle', 'normal', 'position', [1200, 250 sampleFigWH], 'color', 'w')
    h_arrow = [];

    h_ax1 = subplotGap(nSub,1,1:nSub-1);
    
    
    h_l = plot(spfs, sln1, spfs, sln2);    
    set(h_l(1), 'color', col1, 'linestyle', linestyle1);
    set(h_l(2), 'color', col2, 'linestyle', linestyle2);
    
    set(h_l, 'linewidth', 3)
    xlim(spfs([1, end]));        
%     set(gca, 'xtick', xticks, 'ytick', [])
    set(h_ax1, 'xtick', xticks, 'xticklabel', {''}, 'ytick', []);
    ylim([0, max([sln2, sln1])*1.1]);
    hold on;
    
    
    %%
    maxVal = max(sln1);
    idxL = find(sln1 > maxVal/2, 1, 'first');
    idxR = find(sln1 > maxVal/2, 1, 'last');
    
    for i = 1:length(h_arrow)
        if ishandle(h_arrow(i))
            delete(h_arrow(i));
        end
    end

    xLR = spfs([idxL, idxR]) + [.1, -.15];
    yLR = sln1([idxL, idxR]);
            
    arrowSize = 300;
    arrowHeadWidth = 600;
    arrowLineWidths = 300;
    h_arrow(1) = plot(xLR+[.1, -.15], yLR, 'k-', 'linewidth', 3);
    h_arrow(2) = arrowh(xLR, yLR, 'k', [arrowSize, arrowHeadWidth], 100);
    h_arrow(3) = arrowh(fliplr(xLR), yLR, 'k', [arrowSize, arrowHeadWidth], 100);

    
    %%
        
    % sine waves for x-axis.
%     figure(104); clf; hold on;

    p1_norm = get(h_ax1, 'position');
  
    h_ax2_offset = -.2;
    h_ax2_height = .2;
    p2 = [p1_norm(1), p1_norm(2) + h_ax2_offset, p1_norm(3), h_ax2_height];

    h_ax2 = axes('position', p2, 'color', 'r', 'nextplot', 'add');
    set(h_ax2, 'units', 'pixels', 'xlim', xlims);
    p2_pix = get(h_ax2, 'position');
    %%
    ar = p2_pix(3)/p2_pix(4);
    
    h_width = 5;
    h_height = h_width * ar;
    spacing_factor = 2;
    
    yrange = diff(xlims)/ar;
    ylims = yrange*[-.5, .5];
    set(h_ax2, 'ylim', ylims); %, 'visible', 'off');

%%
    L = 0.2;
    
    xlims = spfs([1, end]);
%     axis equal;
%     axis([xlims, -L, L]);    
    
    
    phis = linspace(-2*pi, 2*pi, 100);
    
    w = diff(xlims)/(nTicks*2);
    scl = w/(phis(end)-phis(1));
%     w = 1/(2*pi*10);    
    
    for i = 1:length(xticks);
        xc = xticks(i);
        cpd = xticks_cpd(i);
        X = [xc + phis*scl];
        Y = (L*.6)*sin((phis-phis(1)) *cpd);
        h(i) = plot(X, Y, 'k', 'linewidth', 1);
    end
%     set(h, 'color', 'k', 'lineWidth', 4);
%     set(gca, 'xtick', [], 'ytick', [], 'color', 'none');
%     box on;
    %%
    set(h_ax2, 'visible', 'off')
    3;
    %%
     export_fig(103, 'pdf', [figureFolder 'Fig12\SpfWidthSample.pdf'])
    
    
    
    
end


if any(makeSlides == 1012)
    %%
    figure(105); clf;   
    
    set(105, 'windowStyle', 'normal', 'position', [1200, 350 sampleFigWH], 'color', 'w')
    h_arrow = zeros(2,2);
    h_line = zeros(2,2);
    h = subplot(1,1,1);
    %%
    box on; hold on;
    set(gca, 'xtick', [], 'ytick', []);
    axis([0 1, 0 1]);
%     set(h, 'tic
    arr_lengths = [.32, .32; 
                   .22, .32];
    gap = .1;
    ys = [0.75, 0.25,];
    cols = {[1, 1, 1]*.5, [0, 0, 0]};
    
    arrowSizes      = [1100, 1100; 
                       500, 1100];
    
    arrowHeadWidths = [280, 280; 
                       260, 280];
    
    arrowLineWidths = [7, 7; 
                       3, 7];
    
    for cell_i = 1:2
      
        y = ys(cell_i);
        for lr = 1:2;
            
            m = 0.64;
            if lr == 1
                x = (0.5 - gap/2) + [0, -arr_lengths(cell_i, lr)];
                x_line = (0.5 - gap/2) + [0, -arr_lengths(cell_i, lr)*m];
            elseif lr == 2;
                x = (0.5 + gap/2) + [0, +arr_lengths(cell_i, lr)];
                x_line = (0.5 + gap/2) + [0, +arr_lengths(cell_i, lr)*m];
            end

            if ishandle(h_arrow(cell_i, lr)) && h_arrow(cell_i, lr) > 0;
                delete(h_arrow(cell_i, lr))
            end
            if ishandle(h_line(cell_i, lr)) && h_line(cell_i, lr) > 0;
                delete(h_line(cell_i, lr))
            end

            
            h_line(cell_i, lr) = plot(x_line, y*[1, 1], 'color', cols{cell_i}, 'linewidth', arrowLineWidths(cell_i, lr));
            
            h_arrow(cell_i, lr) = arrowh(x, y*[1, 1], 'k', [arrowSizes(cell_i, lr), arrowHeadWidths(cell_i, lr)], 100);
            set(h_arrow(cell_i, lr), 'facecolor', cols{cell_i}, 'edgecolor', cols{cell_i})
            
        end

    end
    3;
%          x = [arr_l LR = spfs([idxL, idxR]) + [1, -1]*.05;
%     yLR = sln1([idxL, idxR]);
%             
%     arrowSize = 300;
%     arrowHeadWidth = 300;
%     arrowLineWidth = 300;
%     h_arrow(1) = plot(xLR+[.1, -.1], yLR, 'k-', 'linewidth', 3);
%     h_arrow(2) = arrowh(xLR, yLR, 'k', [arrowSize, arrowHeadWidth], 100);
%     h_arrow(3) = arrowh(fliplr(xLR), yLR, 'k', [arrowSize, arrowHeadWidth], 100);
%     end
    
    
end


if any(makeSlides == 102) % 4 sparse firing rate vectors;
    %%
    
    wrp = @(x) [x, x(1)];
    
    nPh = 72;
    phi = linspace(0,360, nPh+1);
    v0 = zeros(1, nPh+1);
%     s = 1;  n = .5;
%     s = .2;  n = .6;    
    s = .3;  n = .9;    
    nCells = 2;
    
    nTrials = 4;
    showAverage = 0;
    fig_id = 103;
    
    ph_cell1 = 90; ph_cell2 = 160;
    ph_cells = {ph_cell1, ph_cell2};
    rand('state', 4);
    odd_idx = 1;
    even_idx = 2;
    
    func = 'sparse_noise';
    sparseness = 5;
    gauss_sig = 30;

    showNoise = 1;
    
%     justShowTrials = [1 2];
%     trials_toShow = [1];
%     trials_toShow = {[1], [1];
%                      [2], [2];
%                      [3], [3];
%                      [4], [4]};
    trials_toShow = {[3,4], [2,4]};
    
%     set(gcf, 'position', [826   426   244    91])
    
    cells_toShow = [2];
    showTrialsIn = 'columns';
    
    
    v = cell(2,3);
    if nTrials == 2
        trial_noiseSpks = {[45, 270], [225, 325]};
        trial_noiseMag = [1, 2/3];
    elseif nTrials == 4
%         trial_noiseSpks = {[45, 270], [225, 325]};
        trial_noiseSpks = {[45, 270], [45, 270], [225, 325], [225, 325]};
        trial_noiseMag = [1, 1, 2/3, 2/3];        
    end    
    
    trial_noiseIdxs = cell(size(trial_noiseSpks));
    for ti = 1:nTrials
        trial_noiseIdxs{ti} = binarySearch(phi, trial_noiseSpks{ti});
    end
%     idx_225 = indmin( abs(phi - 225) );
%     idx_270 = indmin( abs(phi - 270) );
%     idx_45  = indmin( abs(phi - 45 ) );
%     idx_300 = indmin( abs(phi - 325 ) );
%     func = 'cosine';
%     func = 'gaussian';
    gs = @(x,mu,sig) exp( -(x-mu).^2./(2*sig.^2));    
            
    for cell_i = 1:nCells
        for ti = 1:nTrials
            switch func
                case 'gaussian'
                    v_i = gs( phi, ph_cells{cell_i}, gauss_sig )*s ; 
                case 'cosine'
                    v_i = rectified( cos( deg2rad( phi - ph_cells{cell_i} ) )*s );
                case 'sparse_noise',
                    v_i = v0 + double(rand(1, nPh+1) < 1/sparseness)*s;
            end
            v{cell_i, ti} = v_i;
        end
    end
                
        
    v_wn = v;
%     [v1_odd_wn, v2_odd_wn, v1_even_wn, v2_even_wn] = deal(...
%         v1_odd, v2_odd, v1_even, v2_even);


    for ti = 1:nTrials
        for spk_i = 1:length(trial_noiseIdxs{ti})
            j = trial_noiseIdxs{ti}(spk_i);
            for cell_i = 1:nCells
                for tj = 1:nTrials
                    v{cell_i,tj}([j-2:j+2]) = 0;                
                    v_wn{cell_i,tj}([j-2:j+2]) = 0;                
                end
            end
        end
    end

    for ti = 1:nTrials
        for spk_i = 1:length(trial_noiseIdxs{ti})
            j = trial_noiseIdxs{ti}(spk_i);
            for cell_i = 1:nCells
                v_wn{cell_i,ti}(j) = v_wn{cell_i,ti}(j) + n*trial_noiseMag(ti);                
%                 v_wn{cell_i,ti}(j) = n*trial_noiseMag(ti);                                
            end
        end
    end    
    
    for cell_i = 1:nCells
        v{   cell_i,nTrials+1} = mean(vertcat(v{   cell_i, 1:nTrials }));
        v_wn{cell_i,nTrials+1} = mean(vertcat(v_wn{cell_i, 1:nTrials }));
    end
    
    if nTrials == 4
        trialLabels = arrayfun(@(i) sprintf('Trial %d', i), 1:nTrials, 'un', 0);
    elseif nTrials == 2;
        trialLabels = {'Odd # Trials', 'Even # Trials'};
    end
    
    nTrialsToShow = size(trials_toShow, 1) + showAverage;
    nCellsToShow = length(cells_toShow);
    M_spc = [0.0, 0, 0.0];
%     N_spc = [.05, .05, .05];
    N_spc = [0.0, 0, 0];
            
    cols = {'b', [0 .7 0]};
    h_ax = zeros(nCellsToShow, nTrialsToShow);
    
    switch showTrialsIn
        case 'columns', [M, N] = deal(nTrialsToShow, nCellsToShow);
        case 'rows',    [M, N] = deal(nCellsToShow, nTrialsToShow);
    end
    
    figure(fig_id); clf;
        
    for cell_i = 1:nCellsToShow            
        for ti = 1:nTrialsToShow
            switch showTrialsIn
                case 'columns', 
                    h_ax(cell_i, ti) = subplotGap(nTrialsToShow, nCellsToShow, ti, cell_i, M_spc, N_spc);
                case 'rows',    
                    h_ax(cell_i, ti) = subplotGap(nCellsToShow, nTrialsToShow, cell_i, ti, M_spc, N_spc);
            end
            
            trial_idxs = trials_toShow{ti, cells_toShow(cell_i)};            
            
            v_wn_i = mean( vertcat( v_wn {cell_i, trial_idxs } ), 1);
            v_i = mean( vertcat( v { cell_i, trial_idxs } ), 1);
            
            if showNoise
                h_p(cell_i, ti) = plot(phi, v_wn_i, 'r-', 'linewidth', 2);   
%                 h_p(cell_i, ti) = plot(phi, v_wn{cell_i,ti}, 'r-', 'linewidth', 2);   
                hold on;
            end
%             h_p1(cell_i, ti) = plot(phi, v{cell_i,ti}, '-', 'linewidth', 2, 'color', cols{cell_i});    
            h_p1(cell_i, ti) = plot(phi, v_i, '-', 'linewidth', 2, 'color', cols{cells_toShow(cell_i)});    

            set(h_ax(cell_i, ti), 'xtick', [0:90:360], 'xlim', [0, 360], 'ylim', [0 1.3], 'ytick', []); 
            if (ti <= nTrials) 
                
%                 h_ylab = ylabel('Trial 2+22');
%                 set(h_ylab, 'rotation', 0, 'horiz', 'right', 'vert', 'middle')
%                 p = get(gca, 'position');
                if length(trial_idxs) == 1
                    trialLabel_i = sprintf('Trial %d', trial_idxs);
                else
                    trialLabel_i = sprintf('Trials %s', cellstr2csslist( arrayfun(@num2str, trial_idxs, 'un', 0), '+'));
                end                
%                 set(h_ylab, 'string', trialLabel_i);
%                 set(gca, 'position', p); 
                h_txt = text(90, 1.2, trialLabel_i, 'horiz', 'left', 'vert', 'top', 'fontsize', 12);
                
                
                
                
                
%                 ylabel(trialLabels{ti}); 
            else
                ylabel('Average'); 
            end
%             xlabel(' '); 
            if ti == 1
%                 title(sprintf('Cell %d', cell_i), 'fontsize', 13);  
            end
            if ti == 3;
%                 xlabel('Phase');
            end
        end

    end
    3;
        
%     else
%         figure(103); clf;
%         
%         
%         for cell_i = 1:2
%             
%             h_ax(ti, cell_i) = subplotGap(,1,cell_i, 1);
%             if showNoise
%                 h_p(ti, cell_i) = plot(phi, v_wn{cell_i,trials_toShow(cell_i)}, 'r-', 'linewidth', 2);   
%                 hold on;
%             end
%             h_p1(ti, cell_i) = plot(phi, v{cell_i,trials_toShow(cell_i)}, '-', 'linewidth', 2, 'color', cols{cell_i});    
% 
%             if cell_i == 1
% %                 set(gca, 'xtick', []);
%             end
%             if (ti < 3) 
% %                 ylabel(sprintf('%s # Trials', trialSetNames{ti})); 
%             else
% %                 ylabel('Average'); 
%             end
% %             xlabel(' '); 
%             if ti == 1
% %                 title(sprintf('Cell %d', cell_i), 'fontsize', 13);  
%             end
%             if ti == 3;
% %                 xlabel('Phase');
%             end
%         end
%         
%     end
    set(h_ax, 'xtick', [0:90:360], 'xlim', [0, 360], 'ylim', [0 1.3], 'ytick', []); 
%     corr_noNoise = corr(v{1,3}(:), v{2,3}(:));
%     corr_withNoise = corr(v_wn{1,3}(:), v_wn{2,3}(:));
%     fprintf('\n. cc no noise: %.2f. with noise: %.2f\n', corr_noNoise, corr_withNoise)
    3;
%     h_ax(1,1) = subplotGap(3,2,1,1, M_spc, N_spc); 
%     h_ax(1,2) = subplotGap(3,2,1,2, M_spc, N_spc); h_p(1,2) = plot(phi, v2_odd, 'gs-');   ylabel('Trial 1'); xlabel(' '); title('Cell 2'); 
%     h_ax(2,1) = subplotGap(3,2,2,1, M_spc, N_spc); h_p(1,2) = plot(phi, v1_even, 'bs-');   ylabel('Trial 2'); xlabel(' ');
%     h_ax(2,2) = subplotGap(3,2,2,2, M_spc, N_spc); h_p(1,2) = plot(phi, v2_even, 'gs-');   ylabel('Trial 2'); xlabel(' ');
%         
%     h_ax(3,1) = subplotGap(3,2,3,1, M_spc, N_spc); h_p(1,2) = plot(phi, (v1_even+v1_odd)/2, 'bs-', 'linewidth', 3);   xlabel('Phase'); ylabel('Average'); 
%     h_ax(3,2) = subplotGap(3,2,3,2, M_spc, N_spc); h_p(1,2) = plot(phi, (v2_even+v2_odd)/2, 'gs-', 'linewidth', 3);   xlabel('Phase'); ylabel('Average');
    
    
    3;
    
    
    
end

if any(makeSlides == 2)  % absolute vs relative spatial phase.
    
    figure(1); clf;% show clustering of relative spatial phase: 2 gabors with same cosine parameter, but offset
    
    s = 2*512;
    x = 1:s;
    y = 1:s;
    theta_deg = 35;
    [x_grid, y_grid] = meshgrid(x,y);
    theta = deg2rad(theta_deg);
    mu1 = [.5*s;.7*s];
    sig_x = .17*s; sig_y = .1*s;
    k = 20/s;
    d = .4*s;
    mu2 = mu1      - d*[cos(pi/2-theta); sin(pi/2-theta)];
    
    mu2 = mu2+(pi/k)*[cos(theta); -sin(theta)];
%         (A, mu_x, sig_x, mu_y, sig_y,     k, phi, theta, const)
    g1p = [1, mu1(1), sig_x, mu1(2),  sig_y, k, pi/2, -theta, 0];
    g2p = [1, mu2(1), sig_x, mu2(2),  sig_y, k, -pi/2, -theta, 0];
    
    f1 = gabor(g1p, {x_grid, y_grid});
    f2 = gabor(g2p, {x_grid, y_grid});
    imagesc(x, y, f1+f2); axis equal tight; axis xy; colormap('jet')
%     colorbar; 
    set(gca, 'xtick', [], 'ytick', []);
    caxis([-1 1]);
    hold on
    plot(mu1(1), mu1(2), 'ko', 'markersize', 7);
    plot(mu2(1), mu2(2), 'k+', 'markersize', 7);
%     plot(mu3(1), mu3(2), 'ko');

    hold off
    figure(2); clf;

    frm = generateGratingFrame([s 1.5*s], 90-theta_deg, 2*pi/k, pi/2);
    imagesc(frm); axis equal tight; axis xy; colormap('gray'); %colorbar; 
%     caxis([-2 2]);
%     set(h2, 'Position', [.1 .1, .5 .5], 'color', 'none');
    set(gca, 'xtick', [], 'ytick', []);
    
        
%         set(ax1(oi), 'Units', 'Pixels', 'Position', [[-100 20] + [oi*s*1.2, 0 ], dims]);        
    
    
    
    
end

if any(makeSlides == 3)  % find a couple of nice OSP profiles / phase tuning curves
        
    % nice drifting grating OSP: 1145:0
    % nice flashed grating OSP: 4476,0    
    
%     N = 8000;
%     figure(1);
%    inds = findInStructArray(movieGroups_fg, 'nSpikes', [], @(x) any(x(2:end) > N));
%    inds2 =findInStructArray(movieGroups_fg, 'spPh_deg', [], @(x) length(x) > 5);
%    inds = intersect(inds, inds2);
%    for i = 1:length(inds)
% %        grp = movieGroups_fg(inds(i));
%        Gid = grp.Gid;
%        cells = grp.cellIds(  grp.nSpikes > N );

    GC1 = [4476, 0];
    GC2 = [1145, 1];
    
    [Gid, cellId] = elements(GC2);    
    
%     celldata = eval(getName('celldata', Gid, cellId));
%     OSP = celldata.OSP;
    if flashedOrDrifting(Gid) == 1
        S = load('flashedGratingCells_DB_all.mat');
    else
        S = load('driftingGratingCells_DB_all.mat');
    end        
    idx = find(arrayfun(@(s) s.Gid == Gid & s.cellId == cellId, S.allCells), 1);
    cell_s = S.allCells(idx);    
%     OSP = cell_s. OSP;
    
    figure(1); clf;
    hs1 = imageOSP(cell_s, 'subplots:diagonal', 'SPO', 'noLabels'); colormap('jet');
            
    figure(2); clf;
    hs2 = imageOSP(cell_s, 'mean:ph');
    matchAxes('C', [hs1 hs2]);
    colormap('jet');
    
%     disp( [num2str(i) ')' num2str(Gid) ' / ' num2str(ci)] );
%     input('');
            
            
    
% 2911/2
% 2919/2*
% 
% 4466/2
% 4470/2
% 4476/0*   
end

if any(makeSlides == 4)  % show how you get a phase tuning curve 
%     Gid = 2478, [1 2];
% 2312 - 1,2,4

%%% 2761, 1, 3
%%% 2771, 1, 3
    x1 = [2761, 1, 3];
    xt = [1875, 0, 1];
    xt2 = [1877, 0, 1];
    
    y1 = [4538, 0]; % flashed grating
    y2 = [2771, 9];
    y3 = [4482, 0];
    y4 = [4494, 0];

    z1 = [1145, 0];
    z2 = [1140, 0];
    z3 = [1171, 0];
    z4 = [4047, 0];
    z5 = [1131, 4];% (very sparse)
    [Gid, cellId] = elements(y1);
    
    stimtype = 'flash';
    
    if strcmp(stimtype, 'flash')
        resultsFilename = 'indivCells_movie_fg';
    else
        resultsFilename = 'indivCells_grating';
    end
        
%     load cellsGroups_movie_fg;    
%     grp = movieGroups_fg( findInStructArray(movieGroups_fg, 'Gid', Gid));
    nm = getName('celldata', Gid, cellId);
    S = load(resultsFilename, nm);
    C(1) = S.(nm);
    
%     C(2) = eval(getName('celldata', Gid, cellId2));
        
%     ph = [];
    for i = 1:1
        
%         C(i).OSP.sp = linspace(0, max(C(i).OSP.sp), length(C(i).OSP.sp));
        OSP = C(i).OSP;
        [R, oris, spfreqs, phases] = deal(OSP.R, OSP.ori, OSP.sp, OSP.ph);
        OS = mean(R, 3);
        [tmp, inds] = maxElement(OS);
        [ori_ind, sp_ind] = elements(inds);        
        
        figure(10*i+1); clf;
        h = imageOSP(C(i).OSP, 'subplots:diagonal', 'SPO', 'nolabels');
        return;
        [x1, x2, y1, y2] = drawSquareAroundBlock(gca, spfreqs, oris, sp_ind, ori_ind, 'Color', 'w', 'LineWidth', 2);
        lb = [x1 y1]; ur = [x2 y2];        
        dx = diff(spfreqs(1:2));
        [fx, fy] = dsxy2figxy(get(h, 'Parent'),  x1+dx *.5, 175-mean([y1 y2]));
        annotation('arrow', [fx-.2 fx], [fy fy], 'Color', 'k', 'LineWidth', 3);

        % re-order axes;
        figure(10*i+2); clf;
        hs = imageOSP(C(i).OSP);
        for j = 1:length(hs)
            ax(j) = get(hs(j), 'parent'); 
            axes(ax(j));
%             drawSquare(lb, ur, 'Color', 'k', 'LineWidth', 2);
        end
        
        % draw arrow at pref __
        dx = (x2-x1)/10;
        for j = 1:length(hs)-2
            [fx1, fy1] = dsxy2figxy(ax(j),   x2+dx, 175-mean([y1 y2]));
            [fx2, fy2] = dsxy2figxy(ax(j+1), x1-dx, 175-mean([y1 y2]));
%             annotation('line', [fx1 fx2], [fy1 fy2], 'Color', 'k', 'LineWidth', 3);
        end

        figure(10*i+3); if i == 1, clf; hold on; end
        clrs = 'bg';
        ph(i,:) = squeeze(R(ori_ind, sp_ind, :)); 
        plot(phases, ph(i,:), ['s-' clrs(i)]);
        set(gca, 'xtick', phases);
        xlim([0 315]);
        xlabel('Spatial Phase'); ylabel('firing rate');
    end
    
    

    
    
end

if any(makeSlides == 5)  % show the different measures of phase-tuning curve differences.

    Gid = 2500;
    cellId1 = 4;
    cellId2 = 5;
%     load cellGroups_movie_fg;
%     grp = movieGroups_fg( findInStructArray(movieGroups_fg, 'Gid', Gid));
    C(1) = eval(getName('celldata', Gid, cellId1));
    C(2) = eval(getName('celldata', Gid, cellId2));
    ph_cell1 = getPhaseTuningCurveFromOSP(C(1).OSP);
    ph_cell2 = getPhaseTuningCurveFromOSP(C(2).OSP);
    phases = C(1).OSP.ph;

    measures = 2;    
    
    if any(measures == 1)         % Peak-to-peak difference in phase


        figure(2); clf;
        xlims = [0 360];
        plot(phases, ph_cell1, ['s-'], phases, ph_cell2, ['s-']);    set(gca, 'xtick', phases); xlim(xlims);
        xlabel('Spatial Phase'); ylabel('firing rate'); 
        i1 = indmax(ph_cell1);
        i2 = indmax(ph_cell2);
        x = phases([i1, i2]);
        drawVerticalLine(x, 'LineStyle', ':');
        y = mean([ph_cell1(i1), ph_cell2(i2)]);
        [fx, fy] = dsxy2figxy(gca, x, [y y]);
        annotation('doublearrow', fx, fy);

            % example with wrapping around 360
        figure(3); clf;
        shf1 = 2;
        shf2 = 2;
        xlims = [0 360];
        ph3 = ph_cell1([shf1+1:end, 1:shf1]);
        ph4 = ph_cell2([end-shf2+1:end, 1:end-shf2]);
        plot(phases, ph3, ['s-'], phases, ph4, ['s-']);    set(gca, 'xtick', phases); xlim(xlims);
        xlabel('Spatial Phase'); ylabel('firing rate'); 
        i1 = indmax(ph3);
        i2 = indmax(ph4);
        [x1, x2] = elements(phases([i1, i2]));
        drawVerticalLine([x1 x2], 'LineStyle', ':');
        y = mean([ph3(i1), ph4(i2)]);
        [fx1, fy1] = dsxy2figxy(gca, [x1 x2], [y y]-.5);
        [fx2a, fy2a] = dsxy2figxy(gca, [xlims(1) x1], [y y]+.5);
        [fx2b, fy2b] = dsxy2figxy(gca, [xlims(2) x2], [y y]+.5);
        annotation('doublearrow', fx1, fy1, 'color', [.4 .4 .4]);
        annotation('arrow', fx2a, fy2a);
        annotation('arrow', fx2b, fy2b);
    
    end
    
    if any(measures == 2)  % correlation coefficient
    
        figure(1);
        v1 = [2 -3 5];
        v1_odd = [1 5  -4];
        quiver3([0 0], [0 0], [0 0], [v1(1) v1_odd(1)], [v1(2) v1_odd(2)], [v1(3) v1_odd(3)]);

        % Two phase curves - non mean subtracted
        figure(16);
        plot(phases, ph_cell1, ['s-'], phases, ph_cell2, ['s-']);    set(gca, 'xtick', phases); xlim([0 315]);
        xlabel('Spatial Phase'); legend('cell 1', 'cell 2'); title('Actual firing rates');
        h1 = gca;
        zeroaxes;

        % Two phase curves, mean subtracted
        figure(17);
        plot(phases, ph_cell1-mean(ph_cell1), 's-', phases, ph_cell2-mean(ph_cell2), 's-');
        set(gca, 'xtick', phases); xlim([0 315]);
        xlabel('Spatial Phase'); legend('cell 1', 'cell 2'); title('Mean-subtracted firing rates');
        h2 = gca;
        zeroaxes;
        matchAxes('Y', [h1, h2])
        
        
        % show curves with many zeros, show how is bad measure.
        ph3 = [0 1 .4 2 10 21 1 1.2];
        ph4 = [3 15 2 1.1 0 .7 .2 2.5];

        % Two phase curves - non mean subtracted
        figure(18);
        plot(phases, ph3, ['s-'], phases, ph4, ['s-']);    
        set(gca, 'xtick', phases); xlim([0 315]);
        xlabel('Spatial Phase'); title('Actual firing rates');
        h3 = gca;
        zeroaxes;

        % Two phase curves, mean subtracted
        figure(19);
        plot(phases, ph3-mean(ph3), 's-', phases, ph4-mean(ph4), 's-');
        set(gca, 'xtick', phases); xlim([0 315]);
        xlabel('Spatial Phase'); title('Mean-subtracted firing rates');
        h4 = gca;
        zeroaxes;
        matchAxes('Y', [h3, h4])
        
        
        
    end

    if any(measures == 3)  % rho
    
        % Two phase curves - non ranked
        figure(18);
        plot(phases, ph, ['s-']);    set(gca, 'xtick', phases); xlim([0 315]);
        xlabel('Spatial Phase'); ylabel('firing rate'); legend('cell 1', 'cell 2');
        title('  ');

        % Two phase curves - ranked
        figure(19);
        rank1 = tiedrank(ph_cell1);
        rank2 = tiedrank(ph_cell2);
        plot(phases, rank1, 's-', phases, rank2, 's-');
        set(gca, 'xtick', phases); xlim([0 315]);
        xlabel('Spatial Phase'); ylabel('rank'); legend('cell 1', 'cell 2');

    end
    
    if any(measures == 4)  % %%% Kendall's tau
    
        figure(21); clf;
        n = length(rank1);
        plot(phases, ph_cell1, 'k.-', phases, ph_cell2, 'b.-'); hold on
        xlabel('Spatial Phase'); ylabel('firing rate'); legend('cell 1', 'cell 2');
        hc_i = plot(0,0, 'o', 'color', 'k', 'markersize', 6, 'lineWidth', 2);
        hc_j = plot(0,0, 's', 'color', 'b', 'markersize', 6, 'lineWidth', 2);
        hold off;
        ht = text(200, 10, ' ', 'fontsize', 15);
        xlim([-5, 315+5]); 

        ph_cell2(2) = ph_cell2(1);
        for i = 1:2
            for j = i+1:n;
                x_i = phases(i); x_j = phases(j);
                y_i1 = ph_cell1(i); y_i2 = ph_cell2(i);
                y_j1 = ph_cell1(j); y_j2 = ph_cell2(j);
                si = sign(y_i1-y_j1); sj = sign(y_i2-y_j2);
                switch si * sj
                    case 0, col = 'k'; txt = 'neither';
                        if (si == 0) && (sj ~= 0), txt = '"extra x" pair'; end
                        if (si ~= 0) && (sj == 0), txt = '"extra y" pair'; end

                    case 1, col = 'g'; txt = 'concordant pair';
                    case -1, col = 'r'; txt = 'discordant pair';
                end
                set(hc_i, 'xdata', [x_i x_j], 'ydata', [y_i1 y_j1], 'markerfacecolor', col);
                set(hc_j, 'xdata', [x_i x_j], 'ydata', [y_i2 y_j2], 'markerfacecolor', col);

                set(ht, 'string', txt, 'color', col);
    %             input('');
            end
        end
        
    end
    

% 2911/2
% 2919/2*
% 
% 4466/2
% 4470/2
% 4476/0*   
end

if any(makeSlides == 6)  % show the different ways to pick ori/sp for phase tuning curve.

    x1 = [2478, 2 3]; %#ok<*NASGU> % bad
    x2 = [2480, 2 3]; % *bad
    x3 = [2312, 2 4]; 
    
    x4 = [2478, 0 3];    
    x5 = [2787, 0,2];
    x6 = [2911, 0,2];
    
    [Gid, cellId1, cellId2] = elements(x5);

%     load cellGroups_movie_fg;
%     grp = movieGroups_fg( findInStructArray(movieGroups_fg, 'Gid', Gid));
    C(1) = eval(getName('celldata', Gid, cellId1));
    C(2) = eval(getName('celldata', Gid, cellId2));
    
    C(1).OSP.sp = linspace(0, max(C(1).OSP.sp), length(C(1).OSP.sp));
    C(2).OSP.sp = linspace(0, max(C(2).OSP.sp), length(C(2).OSP.sp));
    [R, oris, spfreqs, phases] = elements(C(1).OSP);

    mR1 = mean(C(1).OSP.R, 3); fR1 = mR1/max(mR1(:));
    mR2 = mean(C(2).OSP.R, 3); fR2 = mR2/max(mR2(:));   
    
    oriSpInds = zeros(4,2);
    [tmp, oriSpInds(1,:)] = maxElement(mR1); % will ignore 
    [tmp, oriSpInds(2,:)] = maxElement(mR2);
          oriSpInds(3,:)  = round(  mean([oriSpInds(1,:); oriSpInds(2,:)], 1 ));
    [tmp, oriSpInds(4,:)] = maxElement(fR1 .* fR2);
    [tmp, oriSpInds(5,:)] = maxElement(min(fR1, fR2));
%         figure(10); imageOSP(mOSP1 .* mOSP2, 'mean:Phase');
    mxOSP = C(1).OSP;
    mxOSP.R = fR1 .* fR2;
    
%         for j = 1:4
    % 1. max of 1
    fsize = 14;
    w = 2.5;
    figure(1); clf; 
    g = [.0, 1, .0];
    subplot(1,3,1); imageOSP(C(1).OSP, 'mean:Phase', 'OSP'); title('Cell 1 profile', 'fontsize', fsize);
    hold on; colormap('cool');
    [x1, x2, y1, y2] = drawSquareAroundBlock(gca, spfreqs, oris, oriSpInds(1,2), oriSpInds(1,1), 'Color', 'k', 'LineWidth', w);                    
    [x1, x2, y1, y2] = drawSquareAroundBlock(gca, spfreqs, oris, oriSpInds(2,2), oriSpInds(2,1), 'Color', 'k', 'LineWidth', w);
    [x1, x2, y1, y2] = drawSquareAroundBlock(gca, spfreqs, oris, oriSpInds(4,2), oriSpInds(4,1), 'Color', 'r', 'LineWidth', w);

    % subplot 2:
    subplot(1,3,2); imageOSP(C(2).OSP, 'mean:Phase', 'OSP'); title('Cell 2 profile', 'fontsize', fsize); hold on;
    [x1, x2, y1, y2] = drawSquareAroundBlock(gca, spfreqs, oris, oriSpInds(1,2), oriSpInds(1,1), 'Color', 'k', 'LineWidth', w);                    
    [x1, x2, y1, y2] = drawSquareAroundBlock(gca, spfreqs, oris, oriSpInds(2,2), oriSpInds(2,1), 'Color', 'k', 'LineWidth', w);
    [x1, x2, y1, y2] = drawSquareAroundBlock(gca, spfreqs, oris, oriSpInds(4,2), oriSpInds(4,1), 'Color', 'r', 'LineWidth', w);

    subplot(1,3,3); imageOSP(mxOSP, 'mean:Phase'); title('norm(1) x norm(2)', 'fontsize', fsize); hold on;
    [x1, x2, y1, y2] = drawSquareAroundBlock(gca, spfreqs, oris, oriSpInds(4,2), oriSpInds(4,1), 'Color', 'r', 'LineWidth', w);
    

%     % 2. max of 2
% %     figure(2); clf;
%     subplot(1,3,1); imageOSP(C(1).OSP.R, 'mean:Phase'); title('Cell 1 profile');
%     [x1, x2, y1, y2] = drawSquareAroundBlock(gca, spfreqs, oris, oriSpInds(2,2), oriSpInds(2,1), 'Color', 'w', 'LineWidth', w);                    
%     subplot(1,3,2); imageOSP(C(2).OSP.R, 'mean:Phase'); title('Cell 2 profile');
%     [x1, x2, y1, y2] = drawSquareAroundBlock(gca, spfreqs, oris, oriSpInds(2,2), oriSpInds(2,1), 'Color', 'w', 'LineWidth', w);
% 
%     % 4. max 1x2
%     figure(4); clf;
%     subplot(1,3,1); imageOSP(C(1).OSP.R, 'mean:Phase');  title('Cell 1 profile');
%     [x1, x2, y1, y2] = drawSquareAroundBlock(gca, spfreqs, oris, oriSpInds(4,2), oriSpInds(4,1), 'Color', 'w', 'LineWidth', w);
%     subplot(1,3,2); imageOSP(C(2).OSP.R, 'mean:Phase');  title('Cell 2 profile');
%     [x1, x2, y1, y2] = drawSquareAroundBlock(gca, spfreqs, oris, oriSpInds(4,2), oriSpInds(4,1), 'Color', 'w', 'LineWidth', w);
%     
%     subplot(1,3,3); imageOSP(mOSP1 .* mOSP2, 'mean:Phase'); title('Product: 1 x 2');
%     [x1, x2, y1, y2] = drawSquareAroundBlock(gca, spfreqs, oris, oriSpInds(4,2), oriSpInds(4,1), 'Color', 'w', 'LineWidth', w);
%         
    
%     [fx, fy] = dsxy2figxy(gca,  x2, 175-mean([y1 y2]));
%     annotation('arrow', [fx-.2 fx], [fy fy], 'Color', 'k', 'LineWidth', 2);
%     
%     % 3. midpoint
%     figure(3); clf;
%     subplot(1,3,1); imageOSP(C(1).OSP.R, 'mean:Phase');  title('Cell 1 profile');
%     [x1, x2, y1, y2] = drawSquareAroundBlock(gca, spfreqs, oris, oriSpInds(3,2), oriSpInds(3,1), 'Color', 'w', 'LineWidth', w);
%     subplot(1,3,2); imageOSP(C(2).OSP.R, 'mean:Phase');  title('Cell 2 profile');
%     [x1, x2, y1, y2] = drawSquareAroundBlock(gca, spfreqs, oris, oriSpInds(3,2), oriSpInds(3,1), 'Color', 'w', 'LineWidth', w);
%     
    return;
    % 4. max 1x2
    figure(4); clf;
    subplot(1,3,1); imageOSP(C(1).OSP.R, 'mean:Phase');  title('Cell 1 profile');
    [x1, x2, y1, y2] = drawSquareAroundBlock(gca, spfreqs, oris, oriSpInds(4,2), oriSpInds(4,1), 'Color', 'w', 'LineWidth', w);
    subplot(1,3,2); imageOSP(C(2).OSP.R, 'mean:Phase');  title('Cell 2 profile');
    [x1, x2, y1, y2] = drawSquareAroundBlock(gca, spfreqs, oris, oriSpInds(4,2), oriSpInds(4,1), 'Color', 'w', 'LineWidth', w);
    
    subplot(1,3,3); imageOSP(mOSP1 .* mOSP2, 'mean:Phase'); title('Product: 1 x 2');
    [x1, x2, y1, y2] = drawSquareAroundBlock(gca, spfreqs, oris, oriSpInds(4,2), oriSpInds(4,1), 'Color', 'w', 'LineWidth', w);
    
    colormap('cool');
%     subplot(1,3,1); [x1, x2, y1, y2] = drawSquareAroundBlock(gca, spfreqs, oris, oriSpInds(3,2), oriSpInds(3,1), 'Color', 'w', 'LineWidth', w);

    
    
end

if any(makeSlides == 7)    % show difference between absolute spatial & temporal response phase
    
    doFig = [1];
    s = 450;
    
    if any(doFig == 1)
        
        figure(1); clf; % two flashed grating stimuli
        
        sp = 20;        
        k = 3;
        mu1 = 2*(2*pi/k);
        mu2 = mu1+1.5*pi/k;
        s_x = 1.2;
        expsqr = @(x, mu, sigma)  exp( -((x-mu).^2)./(2*sigma^2));    
        G = @(x, A, mu_x, sig_x, k, phi) A * expsqr(x, mu_x, sig_x) .* cos( k * x + phi ) ;
        f1y = 2;
        f2y = 0;
        x_max = 10;
        x = [0:.01:x_max];
        lb = [20 20];

        f1 = G(x, 1, mu1, s_x, k, 0);
        f2 = G(x, 1, mu2, s_x, k, 0);
        rf_ax = axes; hold on;
        plot(x, f1+f1y, 'b', 'linewidth', 4); hold on;
        plot(x, f2+f2y, 'g', 'linewidth', 4);
        set(gca, 'xtick', [], 'ytick', []);
        ylim([-1 5]);
        drawHorizontalLine([f1y, f2y], 'linestyle', ':', 'color', 'k', 'linewidth', 2);
        axis square; box on;
        set(gca, 'units', 'pixels');
        p = get(gca, 'position');
        set(gca, 'position', [lb, s s]);

        r = .8;
        spPer = (2*pi/k)*(s/x_max);
        fg_ax = axes;        
        set(fg_ax, 'units', 'pixels', 'position', [lb(1), lb(2)+s*r, s, s*(1-r)]);

        imagesc(generateGratingFrame([s*(1-r), s], 90, spPer, +90 )); 
        set(gca, 'xtick', [], 'ytick', [])
        colormap('gray'); axis xy;
    end
    
    
    if any(doFig == 2)
        % show clustering of relative spatial phase: 2 gabors with same cosine parameter, but offset

        s = 450;
        x = 1:s;
        y = 1:s;
        nlevels = 8;
        line_w = 5;
        theta_deg = -35;
        [x_grid, y_grid] = meshgrid(x,y);
        theta = deg2rad(theta_deg);
        mu1 = [.45*s;.48*s];
        sig_x = .1*s; sig_y = .07*s;
        k = 20/s;
        d = .35*s;
        mu2 = mu1      - d*[cos(pi/2-theta); sin(pi/2-theta)];

        mu2 = mu2+(pi/k)*[cos(theta); -sin(theta)];
    %         (A, mu_x, sig_x, mu_y, sig_y,     k, phi, theta, const)
        g1p = [1, mu1(1), sig_x, mu1(2),  sig_y, k, -pi/2, -theta, 0];
        g2p = [1, mu2(1), sig_x, mu2(2),  sig_y, k, +pi/2, -theta, 0];
        
        n = 100;
        lin = linspace(0,1,n)';
        cub = .3*lin.^2;
        zer = zeros(n,1);
        blue_map =  [cub, cub, lin];
        green_map = [cub, lin, cub];
                
        f1 = gabor(g1p, {x_grid, y_grid});
        f2 = gabor(g2p, {x_grid, y_grid});
        th = .001;
        f1_small = abs(f1) < th;
        f1( f1_small ) = 0; 
        f2_small = abs(f2) < th;
        f2( f2_small ) = 0; 
        
        lb = [20 20];
        figure(2); clf; 
        h1 = axes; contour(f1, nlevels, 'linewidth', line_w); 
        axis([0 s 0 s]);
        axis xy;
        set(h1, 'units', 'pixels'); h1p = get(h1, 'position');
        set(h1, 'color', 'none'); 
        set(h1, 'Position', [lb, s, s]);
        colormap(blue_map);
        set(gca, 'xtick', [], 'ytick', []);

        figure(3); clf;
        h2 = axes; contour(f2, nlevels, 'linewidth', line_w);
        axis([0 s 0 s]);
        axis xy;
        set(h2, 'color', 'none'); hold on;        
        set(h2, 'units', 'pixels', 'Position', [lb, s, s]);
        colormap(green_map);
        set(gca, 'xtick', [], 'ytick', []);
        
        figure(4); clf;
        r = .7;
        h3 = axes; 
        spPer = (2*pi/k);       
        imagesc(generateGratingFrame([s*(1-r), s], 90-theta_deg, spPer, +90 )); colormap('gray');
        axis xy;
        set(h3, 'units', 'pixels', 'position', [lb(1), lb(2) + s*r, s, s*(1-r)]);
    
%         frm = generateGratingFrame([s 1.5*s], 90-theta_deg, 2*pi/k, pi/2 );
%         imagesc(frm); axis equal tight; axis xy; colormap('gray'); %colorbar; 
    %     caxis([-2 2]);
    %     set(h2, 'Position', [.1 .1, .5 .5], 'color', 'none');
        set(gca, 'xtick', [], 'ytick', []);
    
    end        
    
    
    
end
    
if any(makeSlides == 8)    % show histogram of different frame lengths
    load cellsGroups_movie_fg
    makeRed = 1;
    frmLng = [movieGroups_fg.frameLength_ms];
    [uLngths, cnt] = uniqueCount(frmLng);
    bar(uLngths, cnt); hold on;
    cnt2 = cnt; cnt2(setdiff([1:length(cnt)], makeRed)) = 0;
    bar(uLngths, cnt2, 'r'); 
    set(gca, 'xtick', floor(uLngths), 'ytick', []);   
    set(gca, 'xticklabel', arrayfun(@(l) [num2str(l) ' ms'], floor(uLngths), 'un', 0))
    for i = 1:length(uLngths)
        text(uLngths(i), cnt(i), num2str(cnt(i)), 'horiz', 'center', 'vert', 'bottom');
    end
    title('Length of each flashed grating frame');    
    ylabel('# of experiments');
end
    
if any(makeSlides == 9)    % show swivelling PSTH
%     c4470_2 = celldata__Group_4470__Cell_2;
    PSTH_bins = c2502_1.PSTH.bins';
    PSTH_vals = c2502_1.PSTH.vals';
    bEdge = binCent2edge(PSTH_bins);
    xlim2 = 120;%bEdge(end);
    tw = getBestTimeWindowFromPSTH(PSTH_bins, PSTH_vals);
%     figure(1);      % just axes, without any graph contents
%     h1 = bar(PSTH_bins, PSTH_vals, 1);
%     ylim([0 600]);
%     axis(axis);
%     set(h1, 'visible', 'off')
    
    figure(1);     % original PSTH
    h1 = bar(PSTH_bins, PSTH_vals, 1);
    set(h1, 'edgecolor', 'b', 'facecolor', 'b');
    ylim1 = 125;
    axis([0 xlim2 0 ylim1]);
    xlabel('time (ms)'); ylabel('spikes/sec');    
    drawVerticalLine(tw, 'color', 'r');    
%     xlabel(' '); ylabel(' ');
%     set(gca, 'xtick', 0, 'xticklabel', ' ', 'ytick', 0, 'yticklabel', ' ', 'color', 'none'); 
        
    figure(2);     % PSTH with stuff outside window set to 0
    PSTH_vals_flat = PSTH_vals; PSTH_vals_flat([1:5, 11:end]) = 0;
    h2 = bar(PSTH_bins, PSTH_vals_flat, 1);    
    set(h2, 'edgecolor', 'b', 'facecolor', 'b');
    axis([0 xlim2 0 ylim1]);
    drawVerticalLine(tw, 'color', 'r');
    xlabel('time (ms)'); ylabel('spikes/sec');    
%     xlabel(' '); ylabel(' ');
%     set(gca, 'xtick', 0, 'xticklabel', ' ', 'ytick', 0, 'yticklabel', ' ', 'color', 'none'); 
    
    figure(3);     % normalized PSTH: forwards
    PSTH_vals_flat_norm = PSTH_vals_flat/sum(PSTH_vals_flat);
    h3 = bar(PSTH_bins, PSTH_vals_flat_norm, 1);
    set(h3, 'edgecolor', 'b', 'facecolor', 'b');
    ylim2 = max(PSTH_vals_flat_norm)/(max(PSTH_vals)/ylim1);
    axis([0 xlim2 0 ylim2]);
    xlabel('time (ms)');  ylabel('P(spike | stim at 0)');    
    drawVerticalLine(tw, 'color', 'r');
                

    figure(4);     % normalized & flipped
%     PSTH_vals_flat_norm_rev = fliplr(PSTH_vals_flat_norm);
%     PSTH_bins_rev = fliplr(PSTH_bins);
    h4 = bar(xlims(2)-PSTH_bins, PSTH_vals_flat_norm, 1);
    set(h4, 'edgecolor', 'b', 'facecolor', 'b');    
    axis([0 xlim2 0 ylim2]);
    xlabel('time (ms)');  ylabel('P(stim | spike at 0)');
    xticks = get(gca, 'xtick'); newxticks = fliplr(-xticks); newxticks(end) = 0;
    set(gca, 'xticklabel', newxticks);
    drawVerticalLine(xlims(2)-tw, 'color', 'r');
    
%     figure(7);     % normalized & flipped PSTH: actual graph
%     h7 = bar(PSTH_bins, PSTH_vals_flat_norm_rev, 1);
%     xlims = xlim;
%     ylim([0 ylim2]);
%     xlabel(' '); ylabel(' ');
%     xlabel(' '); ylabel(' ');
%     set(gca, 'xtick', 0, 'xticklabel', ' ', 'ytick', 0, 'yticklabel', ' ', 'color', 'none'); 
%     box off;
%     rectangle('position', [0, 0 .001, ylim2], 'edgecolor', 'w')



end    
% end

if any(makeSlides == 10)    % show examples of cell profiles & their reproducibility plots
    
   
% good examples:
allcells = [...
...2346,1;
...4946,0; % (not tuned, not reproducible)
...4878,4; % (semi- tuned, semi- reprodubible)
...2484,2; % (tuned & reproducoble)
...2292,7; % xx (not that well tuned, but reproducible)
4718,2; % xx (not that well tuned, but reproducible)
2861,2]; % tuned, but not reproducible, 


    for i = 1 :size(allcells,1);

        showWorking = 'oriSpf_rep';

        gid = allcells(i,1); cellid = allcells(i,2);
        nm = getName('celldata', gid, cellid);
        vr = eval(nm);
        redoStats_id = [gid cellid];
        figure(1);
        imageOSP(vr.OSP, 'mean:ph');
        set(gca, 'xtick', []);
        title(sprintf('%d : %d', gid, cellid));
        redoStats;
        3;

    end    
    
end

if any(makeSlides == 11)   % draw a couple of cosines for diagrammatic representation of cells & sites
    
    t = 0:.01:2*pi;
   for i = 1:10
       figure(i); clf;
       f = cos(t+rand*2*pi);
       plot(t,f, 'linewidth', 10); axis([0 2*pi -1.1 1.1])
       set(gca, 'position', get(gca, 'outerposition'))
       set(gca, 'visible', 'off');
   end
    
    
end

if any(makeSlides == 12)   % show directions that were picked for drifting grating example
    
    figure(1); clf;
    id = 13;
%     dirs_deg = [260   265   270   275   280   285   290   295   300   305   310   315   325   330   335   340   345];
    dirs_deg = gratingGroups_dFree(id).ori_deg;

    dirs = deg2rad(dirs_deg);
    n = length(dirs);
    r_q = .7;
    r_c = .65;
    quiver(zeros(1,n), zeros(1,n), r_q*cos(dirs), r_q*sin(dirs));
    hold on;
    plot(r_c*cos(dirs), r_c*sin(dirs), 'o')    
    axis equal;
    axis([-1 1, -1 1]);    
    set(gca, 'xtick', [], 'ytick', []);
    drawCircle(.7, [0;0]);
    r_txt = .8;
    for i = 1:n
        text(r_txt*cos(dirs(i)), r_txt*sin(dirs(i)), [num2str(dirs_deg(i)) '\circ'], 'hor', 'cent', 'vert', 'middle', 'fontsize', 12);
    end 
    set(gca, 'visible', 'off')    
    
end


%{
 version with multiple locations (unneeded
if any(makeSlides == 102) % 4 sparse firing rate vectors;
    %%
    
    wrp = @(x) [x, x(1)];
    
    nPh = 72;
    phi = linspace(0,360, nPh+1);
    v0 = zeros(1, nPh+1);
    s = 1;  n = .5;
    s = .2;  n = .6;
    nLocations = 1;
    nCells = 2;
    
    nTrials = 4;
%     showAverage = 
        
    
    ph_cell1 = [90 270]; ph_cell2 = [160, 200];
    ph_cells = {ph_cell1, ph_cell2};
%     rand('state', 0);
    odd_idx = 1;
    even_idx = 2;
    
    func = 'sparse_noise';
    sparseness = 5;
    gauss_sig = 30;

    showNoise = 1;
    
%     justShowTrials = [1 2];
    justShowTrials = [];
    showTrialsInRows = 0;
    
    
    v = cell(2,3);
    if nTrials == 2
        trial_noiseSpks = {[45, 270], [225, 325]};
        trial_noiseMag = [1, 2/3];
    elseif nTrials == 4
%         trial_noiseSpks = {[45, 270], [225, 325]};
        trial_noiseSpks = {[45, 270], [45, 270] , [225, 325], [225, 325]};
        trial_noiseMag = [1, 1, 2/3, 2/3];        
    end    
    
    trial_noiseIdxs = cell(size(trial_noiseSpks));
    for ti = 1:nTrials
        trial_noiseIdxs{ti} = binarySearch(phi, trial_noiseSpks{ti});
    end
%     idx_225 = indmin( abs(phi - 225) );
%     idx_270 = indmin( abs(phi - 270) );
%     idx_45  = indmin( abs(phi - 45 ) );
%     idx_300 = indmin( abs(phi - 325 ) );
%     func = 'cosine';
%     func = 'gaussian';
    gs = @(x,mu,sig) exp( -(x-mu).^2./(2*sig.^2));    
            
    for cell_i = 1:nCells
        for loc_i = 1:nLocations
            for ti = 1:nTrials
                switch func
                    case 'gaussian'
                        v_i = gs( phi, ph_cells{cell_i}(loc_i), gauss_sig )*s ; 
                    case 'cosine'
                        v_i = rectified( cos( deg2rad( phi - ph_cells{cell_i}(loc_i) ) )*s );
                    case 'sparse_noise',
                        v_i = v0 + double(rand(1, nPh+1) < 1/sparseness)*s;
                end
                v{cell_i,loc_i, ti} = v_i;
            end
        end
    end
                
        
    v_wn = v;
%     [v1_odd_wn, v2_odd_wn, v1_even_wn, v2_even_wn] = deal(...
%         v1_odd, v2_odd, v1_even, v2_even);


    for ti = 1:nTrials
        for loc_i = 1:nLocations
            for spk_i = 1:length(trial_noiseIdxs{ti})
                j = trial_noiseIdxs{ti}(spk_i);
                for cell_i = 1:nCells
                    v_wn{cell_i,loc_i,ti}(j) = v_wn{cell_i,loc_i,ti}(j) + n*trial_noiseMag(ti);                
                end
            end
        end
    end

    for cell_i = 1:nCells
        for loc_i = 1:nLocations
            v{   cell_i,loc_i,nTrials+1} = mean(vertcat(v{   cell_i,loc_i,[1:nTrials]}));
            v_wn{cell_i,loc_i,nTrials+1} = mean(vertcat(v_wn{cell_i,loc_i,[1:nTrials]}));
        end
    end
    
    if nTrials == 4
        trialLabels = arrayfun(@(i) sprintf('Trial %d', i), 1:nTrials, 'un', 0);
    elseif nTrials == 2;
        trialLabels = {'Odd # Trials', 'Even # Trials'};
    end
    
    M_spc = [0.0, 0, 0.0];
    N_spc = [.05, .05, .05];
    M = nTrials+1;
    cols = {'b', [0 .7 0]};
    h_ax = zeros(nCells, nLocations, nTrials);
    if isempty(justShowTrials)
        figure(102); clf;
        for cell_i = 1:nCells
            for loc_i = 1:nLocations;
            
                for ti = 1:nTrials+1
                    h_ax(cell_i, loc_i, ti) = subplotGap(M,2*nLocations, ti, (cell_i-1)*nLocations + loc_i, M_spc, N_spc);
                    if showNoise
                        h_p(cell_i, loc_i, ti) = plot(phi, v_wn{cell_i,loc_i,ti}, 'r-', 'linewidth', 2);   
                        hold on;
                    end
                    h_p1(cell_i, loc_i, ti) = plot(phi, v{cell_i,loc_i,ti}, '-', 'linewidth', 2, 'color', cols{cell_i});    

                    if (ti <= nTrials) 
                        ylabel(trialLabels{ti}); 
                    else
                        ylabel('Average'); 
                    end
                    xlabel(' '); 
                    if ti == 1
        %                 title(sprintf('Cell %d', cell_i), 'fontsize', 13);  
                    end
                    if ti == 3;
        %                 xlabel('Phase');
                    end
                end
            end
        end
        
    else
        figure(103); clf;
        for cell_i = 1:2
            
            h_ax(ti, cell_i) = subplotGap(2,1,cell_i, 1);
            if showNoise
                h_p(ti, cell_i) = plot(phi, v_wn{cell_i,justShowTrials(cell_i)}, 'r-', 'linewidth', 2);   
                hold on;
            end
            h_p1(ti, cell_i) = plot(phi, v{cell_i,justShowTrials(cell_i)}, '-', 'linewidth', 2, 'color', cols{cell_i});    

            if cell_i == 1
%                 set(gca, 'xtick', []);
            end
            if (ti < 3) 
%                 ylabel(sprintf('%s # Trials', trialSetNames{ti})); 
            else
%                 ylabel('Average'); 
            end
%             xlabel(' '); 
            if ti == 1
%                 title(sprintf('Cell %d', cell_i), 'fontsize', 13);  
            end
            if ti == 3;
%                 xlabel('Phase');
            end
        end
        
    end
    set(h_ax, 'xtick', [0:90:360], 'xlim', [0, 360], 'ylim', [0 1.3], 'ytick', []); 
%     corr_noNoise = corr(v{1,3}(:), v{2,3}(:));
%     corr_withNoise = corr(v_wn{1,3}(:), v_wn{2,3}(:));
%     fprintf('\n. cc no noise: %.2f. with noise: %.2f\n', corr_noNoise, corr_withNoise)
    3;
%     h_ax(1,1) = subplotGap(3,2,1,1, M_spc, N_spc); 
%     h_ax(1,2) = subplotGap(3,2,1,2, M_spc, N_spc); h_p(1,2) = plot(phi, v2_odd, 'gs-');   ylabel('Trial 1'); xlabel(' '); title('Cell 2'); 
%     h_ax(2,1) = subplotGap(3,2,2,1, M_spc, N_spc); h_p(1,2) = plot(phi, v1_even, 'bs-');   ylabel('Trial 2'); xlabel(' ');
%     h_ax(2,2) = subplotGap(3,2,2,2, M_spc, N_spc); h_p(1,2) = plot(phi, v2_even, 'gs-');   ylabel('Trial 2'); xlabel(' ');
%         
%     h_ax(3,1) = subplotGap(3,2,3,1, M_spc, N_spc); h_p(1,2) = plot(phi, (v1_even+v1_odd)/2, 'bs-', 'linewidth', 3);   xlabel('Phase'); ylabel('Average'); 
%     h_ax(3,2) = subplotGap(3,2,3,2, M_spc, N_spc); h_p(1,2) = plot(phi, (v2_even+v2_odd)/2, 'gs-', 'linewidth', 3);   xlabel('Phase'); ylabel('Average');
    
    
    3;
    
    
    
end
%}