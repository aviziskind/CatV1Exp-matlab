function plotCellData(cellData, fig_id)
    global addendum;
    if nargin < 2
        figure;
    else
        figure(fig_id); clf;
    end
    
    showMID_versions = 0;
    
    plotPSTH = 0 ;%~isempty(cellData.PSTH);
    if ~showMID_versions
        plotOSP  = ~isempty(cellData.OSP);    
        plotSTA  = ~isempty(cellData.STAs);
    else
        plotOSP  = ~isempty(cellData.OSP); %&& isfield(cellData.OSP, 'R_pred');    
        plotSTA  = ~isempty(cellData.STAs);% && isfield(cellData.STAs, 'MID');
    end
%     frameDims = size(STAs(:,:,1));
    
    if plotSTA
        numSTAs = size(cellData.STAs.STA, 3);
    else
        numSTAs = 0;
    end
    numPlots = plotPSTH + (plotOSP) + numSTAs;
    
    ad = iff(isempty(addendum), '', [', [' addendum ']']);
    suptitle_2( sprintf(' Group %d ; Cell %d  (Frame length = %3.0f ms) %s',  cellData.Gid, cellData.cellId, getFrameLength('Gid', cellData.Gid), ad) );
    addendum = [];
    
    plot_i = 1;
    if plotPSTH
        subplot(1,numPlots,plot_i);
        S = cellData.PSTH;
        [PSTH_bins, PSTH_vals, frameLen_ms, meanRate, bckgRate, timeWindow] = deal(... 
            S.bins, S.vals, S.frameLength_ms, S.meanRate, S.bckgRate, S.timeWindow_ms);
                    
        tw = iff(exist('timeWindow', 'var'), {timeWindow}, {});
        plotThisPSTH(PSTH_bins, PSTH_vals, [], bckgRate, tw{:});
        
        axis square;
        xlabel('ms after frame start'); ylabel('spikes / sec');
        plot_i = plot_i + 1;
    end

    if plotOSP
        subplot(1,numPlots,plot_i);
        OSP = cellData.OSP;
        if showMID_versions && isfield(OSP, 'R_pred');
            OSP.R = OSP.R_pred;
        end
        imageOSP(OSP, 'mean:Ph');
%         [OSP, oris, sps, phs] = elements(cellData.OSP);
%         imagesc3({oris, sps, phs}, OSP, 3, {'Orientation', 'Spatial Period', 'Phase'});
        
%         colormap('jet')
%         axis square xy;
%         colorbar        
        plot_i = plot_i + 1;

%         subplot(1,numPlots,plot_i);
%         [mx, inds] = maxElement(osp);
%         [ori_max_ind, sp_max_ind] = elements(inds);
%             k = 2*pi / sps(sp_max_ind);
%             theta = deg2rad( oris(ori_max_ind) );
%         frame = generateGratingFrame(frameDims, k, inf, theta );
%         imagesc(frame);
%         axis square xy;
%         plot_i = plot_i + 1;
        
    end
    
    
    if plotSTA
        STA_S = cellData.STAs;
        timeWindows =  STA_S.timeWindow_ms;
        if showMID_versions
            STAs = STA_S.MID;
        else
            STAs = STA_S.STA;
        end
            
        [nx,ny] = size(STAs(:,:,1));
        if isfield(STA_S, 'degPerPix');            
            xs = [1:nx]*STA_S.degPerPix;
            ys = [1:ny]*STA_S.degPerPix;
        else
            xs = 1:nx;
            ys = 1:ny;
        end
            
        for j = 1:numSTAs
            subplot(1,numPlots,plot_i);
            STA_j = STAs(:,:,j);
            imagesc(xs, ys, STA_j);
            axis square ij;
            c = caxis; cL = max(abs(c));
            caxis([-cL, cL]);
%             colormap('gray');
%             colorbar;
            plot_i = plot_i + 1;
            
            if exist('timeWindows', 'var')
                xlabel( sprintf(' %.2f - %.2f ms ', timeWindows(j,1), timeWindows(j,2)) ); 
            end
        end
    end
    
    drawnow;
    
end




%         ospSTA = zeros(frameDims);
%         for ori_i = 1:length(oris)
%             for sp_i = 1:length(sps)
%                 k = 2*pi / sps(sp_i);
%                 theta = oris(ori_i));
%                 frame = generateGratingFrame(frameDims, theta, 1/k, 0);
%                 ospSTA = ospSTA + osp(ori_i, sp_i) * frame;
%             end
%         end
%         toc;
%         imagesc(ospSTA);
%         axis square xy;
%         plot_i = plot_i + 1;
