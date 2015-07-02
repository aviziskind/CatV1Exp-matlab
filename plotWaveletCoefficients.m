function plotWaveletCoefficients

    redo = true;
    coefMode = 'pca';
    % coefMode = 'haar';

    channelActions = {'average', 'concat', 'separate'};
    channelAction = channelActions{3};

    showAveragedWaveforms = strcmp(channelAction, 'average');
    showSeparateWaveforms = any(strcmp(channelAction, {'concat', 'separate'}));
        showSepOnSamePlot = 0;
        showSepOnDifferentPlots = 1;

    nStdShow = 4;
    nStdCutoff = 2.5;

    matchSpiker = 0;    
    
    if ~exist('wvfms', 'var') || redo
        Gid = 4470;
%         Gid = 1822;
    %     Gid = 1699;
        grp = siteDataFor('Gid', Gid, 1);
        cellIds = grp.cellIds;
        % tWindow = [-.25, 1.25];


        [wvfms, t_ms] = getSpikeWaveforms(Gid, [], 'unnorm', 'raw', matchSpiker);
    %     if alignSpikes
    %         [wvfms, t_ms_algn] = alignSpikeWaveforms(wvfms, t_ms);
    %     end


    end
    [nT, nChannels, nSpk] = size(wvfms);
    nCoefs = nT;

    % idx = ord(ks_stat, 'descend');


    drawWaveFormOnHist = 1;
    %     waveFormLocation = 'corner';
        waveFormLocation = 'fill';
    %     drawWaveFormInCorner = true;
    %     drawWaveFormInCorner = true;
    warning('off', 'stats:lillietest:OutOfRangeP')




    [coeffs, V] = getWaveformCoefficients(coefMode, wvfms, [], channelAction);        
    meanWaveform = mean(wvfms, 3);

    M = floor(sqrt(nCoefs));
    N = ceil(nCoefs/M);

    M = floor(3);
    N = ceil(4);


    if showAveragedWaveforms

        idx = 1:nCoefs;    
        V = haar(nCoefs);    


        for i = 1:nCoefs
            x = wvfm_m_haar(i,:);
            all_m(i) = mean(x);
            all_s(i) = std(x);
            idx_include = x > all_m(i)-all_s(i)*nStdCutoff & x < all_m(i)+all_s(i)*nStdCutoff;    
            [h,p(i),ks_stat(i)] = lillietest(x(idx_include));
        end    


        figure(1); clf;
        for i = 1:nCoefs
            hnd_ax(1,i) = subplot(M,N,i);
            j = idx(i);
            x = wvfm_m_haar(j,:);
            m = all_m(j);
            s = all_s(j);

            nBin = 30;
            x_edges_plot = linspace(m-s*nStdShow, m+s*nStdShow, nBin);
            x_cent_plot = binEdge2cent(x_edges_plot);            
            n = histcnt(x, x_edges_plot);
            bar(x_cent_plot, n, 1);    

            title(sprintf('(%d)', idx(i)))
            xlim([x_edges_plot(1), x_edges_plot(end)]);    
            set(gca, 'ytick', [], 'yticklabel', []);
%             set(gca, 'xtick', [])
            ylim([0 max(n)*1.1]);
        %     drawVerticalLine([m - s*nStdCutoff, m + s*nStdCutoff], 'linestyle', ':');

            if drawWaveFormOnHist
                p = get(hnd_ax(1,i), 'position');
                hnd_ax(2,i) = axes('position', p);
                wvfm = V(j,[1:end, end]);
                wvfm = wvfm/max(abs(wvfm));
                hnd_plot(i) = stairs(hnd_ax(2,i), 1:nCoefs+1, wvfm, 'r');
                set(hnd_ax(2,i), 'color', 'none', 'xtick', [], 'ytick', [], 'xlim', [1, nCoefs+1], 'ylim', [-1.2, 1.2], 'box', 'off');
            end
            3;

        end

        if drawWaveFormOnHist
            drawnow;
            for i = 1:nCoefs
                set(hnd_ax(2,i), 'position', get(hnd_ax(1,i), 'position'));
            end
        end
    end


    if showSeparateWaveforms    

        [all_m, all_s] = deal(zeros(nCoefs, nChannels));
        for co_i = 1:nCoefs
            for ch_i = 1:nChannels
    %             if strcmp(channelAction, 'concat')
    %             elseif strcmp(channelAction, 'separate')
    %         end
                x = coeffs(co_i,ch_i,:); 
                x = x(:);
                m = mean(x);
                s = std(x);
                idx_include = x > m-s*nStdCutoff & x < m+s*nStdCutoff;    
                [h,p,ks_stat(co_i,ch_i)] = lillietest(x(idx_include));             %#ok<AGROW>
                all_m(co_i,ch_i) = m;
                all_s(co_i,ch_i) = s;
            end
        end


        idx = ord(ks_stat(:), 'descend');
    %     idx = 1:nCoefs*nChannels;

        cols = 'bcgr';
        ext = @(x) [x, x(end)+diff(x(1:2))];
        rep = @(x) x([1:end, end]);

        if showSepOnSamePlot        
            nChannelsToPlot = 1;

        elseif showSepOnDifferentPlots        

            glob_i = 1;
            nChannelsToPlot = nChannels;
            nChannelsToPlot = 1;
        end


        nBin = 40;    

        for fig_i = 1:nChannelsToPlot
            figure(fig_i); clf;

            for co_i = 1:nCoefs
                if co_i > M*N
                    break;
                end
                hnd_ax(fig_i,co_i,1) = subplot(M,N,co_i); hold on; %#ok<AGROW>

                if showSepOnSamePlot                
                    rngs = [1;1]*all_m(co_i,:) + [-1;1] * all_s(co_i,:)*nStdShow;        
                    x1 = min(rngs(:));
                    x2 = max(rngs(:));

                    x_edges_plot = linspace(x1, x2, nBin);
                    x_cent_plot = binEdge2cent(x_edges_plot);            
                    maxval = 0;
                    for j = 1:nChannels
                        x = coeffs(co_i,j,:); x = x(:);           
                        n = histcnt(x, x_edges_plot);
            %             bar(x_cent_plot, n, 1);    
                        stairs( x_edges_plot, rep(n), 'color', cols(j), 'linewidth', 2);            
                        maxval = max([maxval; n(:)]);
                    end
                    title(sprintf('(%d)', idx(co_i)))


                elseif showSepOnDifferentPlots                                
                    j = idx(glob_i);
                    set(hnd_ax(fig_i,co_i,1), 'UserData', j)
                    [coef_i, ch_i] = ind2sub([nCoefs, nChannels], j);

                    x = coeffs(coef_i, ch_i,:); x = x(:);
                    m = all_m(j);  s = all_s(j);

                    nBin = 40;
                    x_edges_plot = linspace(m-s*nStdShow, m+s*nStdShow, nBin);
                    x_cent_plot = binEdge2cent(x_edges_plot);            
                    n = histcnt(x, x_edges_plot);
                    bar(x_cent_plot, n, 1, 'facecolor', cols(ch_i), 'edgecolor', cols(ch_i));    

        %             title(sprintf('(%d) %d : %d', glob_i, ch_i, coef_i )); %, ks_stat(coef_i,ch_i)
                    title(sprintf('(%d) %d : %d [%.2g]', glob_i, ch_i, coef_i, ks_stat(coef_i,ch_i) ));


                    maxval = max(n);                
    %                 drawVerticalLine([m - s*nStdCutoff, m + s*nStdCutoff], 'linestyle', ':');

                    glob_i = glob_i + 1;
                end

                xlim([x_edges_plot(1), x_edges_plot(end)]);    
                set(gca, 'ytick', [], 'yticklabel', []);
%                 set(gca, 'xtick', []);
                ylim([0 maxval*1.1]);

            end

            if drawWaveFormOnHist
                drawnow;
                for co_i = 1:nCoefs
                    if co_i > M*N,
                        break;
                    end
                    p = get(hnd_ax(fig_i,co_i,1), 'position');
                    switch waveFormLocation
                        case 'corner',
                            p(2) = p(2) + p(4)*(2/3);
                            p(3:4) = p(3:4)/3;
                        case 'fill', % do nothing                            
                    end

                    hnd_ax(fig_i,co_i,2) = axes('position', p);
                    gi = get(hnd_ax(fig_i,co_i,1), 'UserData');
                    [coef_i, ch_i] = ind2sub([nCoefs, nChannels], gi);
                    wvfm_m = meanWaveform(:,ch_i);
                    wvfm_m = wvfm_m/norm(wvfm_m);

                    switch coefMode
                        case 'haar',                                      
                            ylims = [min(wvfm_m), max(wvfm_m)]; ylims = ylims + diff(ylims)*[-1, 1]/20;
                            haar_wvlt = rep(V(coef_i,:));
                            haar_wvlt = haar_wvlt/max(abs(haar_wvlt)) * .95;
                            haar_wvlt = (diff(ylims)/2)*haar_wvlt + mean(ylims);

                            stairs(hnd_ax(fig_i,co_i,2), 1:nCoefs+1, rep(wvfm_m), 'color', [.7, .7, .7]); hold on;
                            stairs(hnd_ax(fig_i,co_i,2), 1:nCoefs+1, haar_wvlt, 'k');

                            set(hnd_ax(fig_i,co_i,2), 'color', 'none', 'xtick', [], 'ytick', [], 'xlim', [1, nCoefs+1], 'ylim', ylims, 'box', 'off');
                        case 'pca',

                            pca_wvfm = V{ch_i}(coef_i,:);
                            pca_wvfm = pca_wvfm / norm(pca_wvfm); %* norm(wvfm_m) *.7;

                            plot(hnd_ax(fig_i,co_i,2), 1:nCoefs, wvfm_m, '.-', 'color', [.7, .7, .7]); hold on;
                            plot(hnd_ax(fig_i,co_i,2), 1:nCoefs, pca_wvfm(:), 'k.-');
                            set(hnd_ax(fig_i,co_i,2), 'color', 'none', 'xtick', [], 'ytick', [], 'xlim', [1, nCoefs], 'box', 'off');                        
                            3;
                    end

                end
            end


        end




    end

end

%{
if showSeparateWaveforms && showSepOnDifferentPlots && 0;
    
    cols = 'brgk';
    
    glob_i = 1;
    for i1 = 1:nChannels
        figure(i1); clf;
        for ci = 1:nCoefs
            hnd_ax(i1,ci) = subplot(M,N,ci);
            j = idx(glob_i);
            [coef_i, ch_i] = ind2sub([nCoefs, nChannels], j);
            
            x = coeffs(coef_i, ch_i,:); x = x(:);
            m = all_m(j);
            s = all_s(j);

            nBin = 40;
            x_edges_plot = linspace(m-s*nStdShow, m+s*nStdShow, nBin);
            x_cent_plot = binEdge2cent(x_edges_plot);            
            n = histcnt(x, x_edges_plot);
            bar(x_cent_plot, n, 1, 'facecolor', cols(ch_i), 'edgecolor', cols(ch_i));    

%             title(sprintf('(%d) %d : %d', glob_i, ch_i, coef_i )); %, ks_stat(coef_i,ch_i)
            title(sprintf('(%d) %d : %d [%.2g]', glob_i, ch_i, coef_i, ks_stat(coef_i,ch_i) ));
             
            
            if drawWaveFormOnHist
                p = get(hnd_ax(1,i), 'position');
                hnd_ax(2,i) = axes('position', p);
                wvfm = V(j,[1:end, end]);
                wvfm = wvfm/max(abs(wvfm));
                hnd_plot(i) = stairs(hnd_ax(2,i), 1:nCoefs+1, wvfm, 'r');
                set(hnd_ax(2,i), 'color', 'none', 'xtick', [], 'ytick', [], 'xlim', [1, nCoefs+1], 'ylim', [-1.2, 1.2], 'box', 'off');
            end

            3;

            glob_i = glob_i + 1;
        end
    end
%     drawnow;
%     for i = 1:nCoefs
%         set(hnd_ax(2,i), 'position', get(hnd_ax(1,i), 'position'));
%     end
end
%}

    %{
    switch channelActions, 

        case 'average';

            switch coefMode
                case 'haar',
                    wvfms_m = squeeze(mean(wvfms, 1));
                    wvfm_m_haar = haarTransform(wvfms_m);
                    nCoefs = size(wvfm_m_haar,1);
                case 'pca'

            end


        case 'concat',

            switch coefMode
                case 'haar'
                    wvfms_cat = reshape(wvfms, [nT, nChannels*nSpk]);
                    coeffs = haarTransform(wvfms_cat, 1);            
                    nCoefs = size(coeffs, 1);
                case 'pca'
                    [nT, nChannels, nSpk] = size(wvfms);            
                    spk_wvfm_cat = reshape(wvfms, [nT*nChannels, nSpk])';
                    meanWaveform = mean(spk_wvfm_cat,1);
                    spk_wvfm_cov = cov(spk_wvfm_cat); % covariance matrix
                    [wvfm_eigvc, wvfm_ev] = eig(spk_wvfm_cov);
                    wvfm_ev = diag(wvfm_ev);
                    [tmp, idx_top_eigs] = sort(wvfm_ev, 'descend');
                    pca_comps = wvfm_eigvc(:, idx_top_eigs );
                    spk_wvfm_cat_msub = bsxfun(@minus, spk_wvfm_cat, meanWaveform);
                    coeffs = [pca_comps' * spk_wvfm_cat_msub'];            
                    nCoefs = nT;
            end                
            coeffs = reshape(coeffs, [nCoefs, nChannels, nSpk]);   

        case 'separate',

            switch coefMode


                case 'pca'
                    [nT, nChannels, nSpk] = size(wvfms);            
                    coeffs = zeros(nT, nChannels, nSpk, 'single');
                    [V, meanWaveform] = deal( cell(1,nChannels) );
                    for ch_i = 1:nChannels
                        wvfms_i = reshape(wvfms(:,ch_i,:), [nT, nSpk])';
                        meanWaveform{ch_i} = mean(wvfms_i,1);
                        spk_wvfm_cov = cov( wvfms_i ); % covariance matrix
                        [wvfm_eigvc, wvfm_ev] = eig(spk_wvfm_cov);
                        wvfm_ev = diag(wvfm_ev);
                        [tmp, idx_top_eigs] = sort(wvfm_ev, 'descend');
                        pca_comps = wvfm_eigvc(:, idx_top_eigs );
                        spk_wvfm_msub = bsxfun(@minus, wvfms_i, meanWaveform{ch_i});
                        coeffs_i = [pca_comps' * spk_wvfm_msub'];            
                        coeffs(:,ch_i,:) = coeffs_i; 
                        V{ch_i} = pca_comps;
                    end
                    nCoefs = nT;
            end    



    end
    %}
