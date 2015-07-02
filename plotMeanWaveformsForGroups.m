
%%
Gids = getAllGids('f');
% Gids = Gids(1:100);
nGids = length(Gids);

t_ms = [-0.9:.05:1.2];
t_ms_rep = repmat(t_ms, [1, 4]);
idx_0_cat = find(t_ms_rep == 0);

% col_order = [0,0,1;  0,.5,0;  1,0,0;   0, 0.75, 0.75;  0.75, 0, 0.75;  0.75, 0.75, 0;  0.25, 0.25, 0.25];
col_order = [1,0,0;  0,.5,0;   0,0,1;  0, 0.75, 0.75;  0.75, 0, 0.75;  0.75, 0.75, 0;  0.25, 0.25, 0.25];

% allCC_neg = cell(1,nGids);
% allCC_wvfm = cell(1,nGids);
dists_pca = cell(1,nGids);
dists_neg = cell(1,nGids);

% mean_ratio_dists = zeros(1,nGids);
mean_ratio_dists = cellfun(@(x,y) median(x./y), dists_pca, dists_neg);


curVoltageFilter('default');
%%
% progressBar('init-', 120);
% for gid_i = 1:120
%     progressBar(gid_i);
    % gid_i = gid_i+1;
    gid_i = 116;
    Gid = Gids(gid_i);

    [uCells, cellSpkIdxs] = getCellSorting(Gid, 'cells');
    S_prop = load(getFileName('properties', Gid));
    spkIds = getCellSorting(Gid, 'cells');
    idx_use = uCells > 0;
%     if nnz(idx_use) <= 2
%         continue;
%     end
    cellSpkIdxs = cellSpkIdxs(idx_use);
    
    
    spks = getSpikeWaveforms(Gid, [], 0, 0, 0)/4;
    spks_ccw = getSpikeWaveforms(Gid, [], 0, 1, 0)/4;
    nC = length(cellSpkIdxs);

    [nT, nCh, nSpk] = size(spks);
    spks_cat = reshape(spks, [nT*nCh, nSpk]);
    spks_cat_ccw = reshape(spks_ccw, [nT*nCh, nSpk]);

    mWvfm_C = cellfun(@(ii) mean(spks_cat(:,ii), 2), cellSpkIdxs, 'un', 0);
    mWvfms = [mWvfm_C{:}];
    
    mNegAmps = mWvfms(idx_0_cat,:);

    % all_negAmps = cellfun(@(ii) spks_cat(idx_0_cat,ii), cellSpkIdxs, 'un', 0);
    all_negAmps = spks_cat(idx_0_cat,:)';
    all_negAmps_ccw = spks_cat_ccw(idx_0_cat,:)';

    negAmps_pca = doPCA(all_negAmps_ccw, 2)';

    cell_negAmps_pca = cellfun(@(ii) negAmps_pca(ii,:), cellSpkIdxs, 'un', 0);
    cell_negAmps_pca_M = cellfun(@mean, cell_negAmps_pca, 'un', 0);
    cell_negAmps_pca_C = cellfun(@cov, cell_negAmps_pca, 'un', 0);    
    MD_neg = cellfun(@(X) mahal(X, X), cell_negAmps_pca, 'un', 0);

    wvfm_pca = getGroupWaveformCoefficients('PCA', Gid, 2, 'concat', 0, 1, 0);
    cell_wvfm_pca = cellfun(@(ii) wvfm_pca(ii,:), cellSpkIdxs, 'un', 0);
    cell_wvfm_pca_M = cellfun(@mean, cell_wvfm_pca, 'un', 0);
    cell_wvfm_pca_C = cellfun(@cov, cell_wvfm_pca, 'un', 0);
    MD_wvfm = cellfun(@(X) mahal(X, X), cell_wvfm_pca, 'un', 0);
    
    pr_idxs = nchoosek(1:nC, 2);

    dists_pca_i = zeros(1, length(pr_idxs));
    dists_neg_i = zeros(1, length(pr_idxs));
    for pr_i = 1:size(pr_idxs, 1)   
       i1 = pr_idxs(pr_i,1);
       i2 = pr_idxs(pr_i,2);
       dists_neg_i(pr_i) = -quadProdGaussians(cell_negAmps_pca_M{i1}, cell_negAmps_pca_C{i1},...
                                             cell_negAmps_pca_M{i2}, cell_negAmps_pca_C{i2}, 'log', [], [], 1);
       dists_pca_i(pr_i) = -quadProdGaussians(cell_wvfm_pca_M{i1}, cell_wvfm_pca_C{i1},...
                                             cell_wvfm_pca_M{i2}, cell_wvfm_pca_C{i2}, 'log', [], [], 1);

    end
    dists_pca{gid_i} = dists_pca_i;
    dists_neg{gid_i} = dists_neg_i;

    mean_ratio_dists(gid_i) = median(dists_pca_i ./ dists_neg_i);

    negAmps = mWvfms(idx_0_cat,:);
    maxAbs = max(abs(negAmps), [], 1);

% end
3;
% cc_neg = corr(negAmps);
% cc_wvfm = corr(mWvfms);

% idx_prs = find( tril(ones(nC), -1) );

% if any(  (cc_neg(:) < .2) & (cc_wvfm(:) < .5) )
%    beep; 
% end


figure(54); clf;
L = lims(mWvfms(:), .05);
for i = 1:nC
    subplotGap(nC,1,i,1);
    plot(mWvfms(:,i));
    set(gca, 'xtick', [], 'ylim', L);    
    drawHorizontalLine(-200, 'linestyle', ':', 'color', 'k')
end
%%
figure(55);
cell_idxs = [2,3];
nSpk = 4;
for j = 1:2
    subplot(2,1,j);
    plot(1:length(t_ms_rep), spks_cat(:, cellSpkIdxs{cell_idxs(j)}(1:nSpk)  ));    
end

%%
cell_idxs_use = 1:nC; %[1, 2,3,4];
% idxs = [1, 3];
figure(56); clf;
subplot(2,1,1);
hold on; box on;
nStd = 2;
nPts = 50;
% cols = get(gca, 'colorOrder');
for ci = cell_idxs_use
    h_neg_scat(ci) = plot(cell_negAmps_pca{ci}(:,1), cell_negAmps_pca{ci}(:,2), '.', 'color', col_order(ci-cell_idxs_use(1)+1,:), 'markersize', 1); 
end
h_cur_neg(1) = plot(0,0, 'visible', 'off');
h_cur_neg(2) = plot(0,0, 'visible', 'off');
for ci = cell_idxs_use
    [x_el, y_el] = ellipsoidFromCov(cell_negAmps_pca_M{ci}, cell_negAmps_pca_C{ci}, nStd, nPts);
    h_neg_elps(ci) = plot(x_el, y_el, '-', 'color', 'k', 'linewidth', 2);
end

axis tight square;
xlabel('PCA 1 (Negative Amplitudes)'); ylabel('PCA 2 (Negative Amplitudes)');
 
subplot(2,1,2);
hold on; box on;
for ci = cell_idxs_use
    h_wvfm_scat(ci) = plot(cell_wvfm_pca{ci}(:,1), cell_wvfm_pca{ci}(:,2), '.', 'color', col_order(ci-cell_idxs_use(1)+1,:), 'markersize', 1);
end
h_cur_pca(1) = plot(0,0, 'visible', 'off');
h_cur_pca(2) = plot(0,0, 'visible', 'off');
for ci = cell_idxs_use
    [x_el, y_el] = ellipsoidFromCov(cell_wvfm_pca_M{ci}, cell_wvfm_pca_C{ci}, nStd, nPts);
    h_wvfm_elps(ci) = plot(x_el, y_el, '-', 'color', 'k', 'linewidth', 2);
end
axis tight square;
xlabel('PCA 1 (Full Waveform)'); ylabel('PCA 2 (Full Waveform)');

%%
% figure(60);
% clear idxs negAmps wvfms spk_idxs;
cell_idxs = [3, 2, 1, 4];

[~, idx_srt1] = sort( MD_wvfm{cell_idxs(1)}, 'descend' );
[~, idx_srt2] = sort( MD_wvfm{cell_idxs(2)}, 'descend' );

nSpkEachCell = cellfun(@length, cellSpkIdxs);
rng(0);
i_rand1 = randperm(nSpkEachCell(cell_idxs(1)));
i_rand2 = randperm(nSpkEachCell(cell_idxs(2)));

% spk_idxs{1} = [i_rand1(1:30)]; 
% spk_idxs{2} = [i_rand2(1:30)];
spk_idxs{1} = [1:30]; 
spk_idxs{2} = [1:30];


% spk_idxs{1} = [idx_srt1(1:30)]; 
% spk_idxs{2} = [idx_srt2(1:30)];
% nSpks = length(spk_idxs{1});


% negAmps_raw = 

idxs{1} = cellSpkIdxs{cell_idxs(1)}(spk_idxs{1});
idxs{2} = cellSpkIdxs{cell_idxs(2)}(spk_idxs{2});

negAmps{1} = S_prop.negAmps(idxs{1},:)/4;
negAmps{2} = S_prop.negAmps(idxs{2},:)/4;

wvfms{1} = reshape( spks_cat(:,idxs{1}), [nT, nChannels, nSpks]);
wvfms{2} = reshape( spks_cat(:,idxs{2}), [nT, nChannels, nSpks]);
wvfms_cat{1} = spks_cat(:,idxs{1});
wvfms_cat{2} = spks_cat(:,idxs{2});
ylims = lims([wvfms{1}(:); wvfms{2}(:)], .05);

%%



%%

neg_amp_fig = 91;
wvfm_fig = 92;
figure(neg_amp_fig); clf; hold on;  box on;
h_neg = plot(0,0);
h_neg_amp_ax = gca;

figure(wvfm_fig); clf; hold on;  box on;
for ch_i = 1:nChannels
    h_wvfm(ch_i) = plot(0,0);
end
h_wvfm_ax = gca;
set([h_neg_amp_ax, h_wvfm_ax], 'xtick', []);

% set(h_neg, 'parent', h_wvfm_ax);


nChannels = 4;
% t_ms_ext = [t_ms, [1 : length(t_ms) * (nChannels-1)]*diff(t_ms(1:2)) + t_ms(end)];
t_ms_ext_M = bsxfun(@plus, t_ms(:), [0:nChannels-1]*diff(t_ms([1, end])) );
t_ms_ext = t_ms_ext_M(:);
t_neg_ext = t_ms_ext(idx_0_cat);
%%
% h_ax(1) = axes(getNormPosition( ));

clusters_show = [1 2 3 4];
show_ellipses = 1;
cell_i = 4;
spk_i = 1;
col = arrayfun(@(i) col_order(i,:), cell_idxs, 'un', 0);
% col{1} =  col_order(cell_idxs(1),:);
% col{2} =  col_order(cell_idxs(2),:);

col_dark{1} = col{1}*.5;
col_dark{2} = col{2}*.5;


showMean = 1;

if showMean
    ylims = lims(mWvfms(:), .02);
end

for spk_i = 4;% 1:1; %1:5; %1:10; %1:1; %1:1
    if showMean
        
        negAmps_show = mNegAmps(:,cell_idxs(cell_i));
        wvfms_show = reshape(mWvfms(:, cell_idxs(cell_i)), 43, 4);
    else
        wvfms_show = wvfms{cell_i}(:,:,spk_i);
        negAmps_show = negAmps{cell_i}(spk_i, :);
    end
    
    set(h_neg, 'xdata', 1:4, 'ydata', negAmps_show, 'marker', 'o', 'linestyle', 'none',...
         'color', col{cell_i}, 'markerfacecolor', col{cell_i});
...        'color', 'r', 'markerfacecolor', 'none');
        
    set(h_neg_amp_ax, 'xlim', lims([1, 4], .1), 'xtick', [], 'ylim', ylims);
    
    for ch_i = 1:nChannels
        set(h_wvfm(ch_i), 'xdata', t_ms_ext_M(:,ch_i), 'ydata', wvfms_show(:,ch_i), ...
            'marker', '.', 'color', col{cell_i})    
    end
    set(h_wvfm_ax, 'xlim', lims(t_ms_ext), 'ylim', ylims);

    
    if ~showMean
        negAmps_xy = cell_negAmps_pca{cell_idxs(cell_i)}(spk_idxs{cell_i}(spk_i),:);
        pca_xy = cell_wvfm_pca{cell_idxs(cell_i)}(spk_idxs{cell_i}(spk_i),:);
        %%
        set(h_cur_neg, 'xdata', negAmps_xy(1), 'ydata', negAmps_xy(2), 'visible', 'on');
        set(h_cur_pca, 'xdata', pca_xy(1),     'ydata', pca_xy(2), 'visible', 'on');


        set([h_cur_neg(1), h_cur_pca(1)], 'marker', '*', 'markersize', 12, 'markerfacecolor', col_dark{ cell_i}, 'linewidth', 2, 'color', 'k');

    %     set(h_cur_neg(2), 'marker', 's', 'markersize', 8, 'color', 'k', 'linewidth', 1, 'visible', 'off');

    %     set(, 'marker', '*', 'markersize', 12, 'color', 'k');

        set([h_cur_neg(2), h_cur_pca(2)], 'xdata', pca_xy(1),     'ydata', pca_xy(2), 'visible', 'off');

    %     set(h_cur_pca(2), 'marker', '*', 'markersize', 12, 'color', 'k', 'visible', 'off');

        set([h_neg_scat(clusters_show), h_neg_elps(clusters_show), h_wvfm_scat(clusters_show), h_wvfm_elps(clusters_show)], 'visible', 'on');
        clusters_hide = setdiff(1:nC, clusters_show);
        set([h_neg_scat(clusters_hide), h_neg_elps(clusters_hide), h_wvfm_scat(clusters_hide), h_wvfm_elps(clusters_hide)], 'visible', 'off');
    end
    if ~show_ellipses
        set([h_neg_elps, h_wvfm_elps], 'visible', 'off');
    end
    drawnow;
    pause(.2);
    
    set([h_cur_neg, h_cur_pca], 'visible', 'off');
    
    
end
        %%
                
        
        set([h_ax, h_ax2], 'xlim', x([1, end]), 'ylim', ylims);










%%
%{
figure(101); clf;
spk_i = 5;
plot(t_ms_ext, spks_cat(idx_0_cat,spk_i), '.-');
hold on;
plot(t_ms_ext(idx_0_cat), S_prop.negAmps(spk_i, :)/4, 'ro');

% clf;
% plot(spks_cat(idx_0_cat,1:70000)', S_prop.negAmps(1:70000, :)/4, '.');
% hold on;
% fplot(@(x) x, xlim, 'k:')
%}