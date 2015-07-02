function fig_locations(x_in)


    function [OSP1, OSP2, OSP_mu] = getOSPs(GCC)
        [Gid, cellId1, cellId2] = elements(GCC);
        S = load('flashedGratingCells_all');
        allCells = S.allCells;
        gids = [allCells.Gid];
        cellids = [allCells.cellId];
        id1 = find(gids == Gid & cellids == cellId1, 1);
        id2 = find(gids == Gid & cellids == cellId2, 1);
        idmu = find(gids == Gid & cellids == 0, 1);
        OSP1 = allCells(id1);
        OSP2 = allCells(id2);        
        OSP_mu = allCells(idmu);        
    end

    function R = shiftOSP(R, shft)
        R = R([shft+1:end, 1:shft],:,:);        
    end

    wrp = @(x) x([1:end, 1]);
    
    l_maxR1xR2 = 1;
    l_maxMinFrac = 2; %#ok<*NASGU>
    l_maxMU = 3;

    x1 = [2478, 2 3]; % ok - peaks are quite separate
    x2 = [2480, 2 3]; % *peaks separate
    x3 = [2312, 2 4];     
    x4 = [2478, 0 3];    
    x5 = [2787, 0,2]; % (Both same)
    x6 = [2911, 0,2];  %ok.
    x7 = [1792, 1,3];  %ok. but need to shift
    x8 = [1792, 5,3];  %ok. but need to shift

    shft = [10];
    
    doFigs = l_maxMinFrac;
    
    if nargin == 1
        xx = x_in;
    else
        xx = x7;
    end
    doLocPlots(xx);
    
    function doLocPlots(GCC)
%     if any(doFigs == l_maxR1xR2)
    
%        celldata1 = allCells(
        [OSP1, OSP2, OSP_mu] = getOSPs(GCC);                
        R1 = OSP1.R;           R2 = OSP2.R;            Rmu = OSP_mu.R;
        mR1 = mean(R1, 3);     mR2 = mean(R2, 3);      mRmu = mean(Rmu, 3);  
        fR1 = mR1/max(mR1(:)); fR2 = mR2/max(mR2(:));  fRmu = mRmu/max(mRmu(:));      
        
%         ori1_idx = max(mR1);
%         
%         ori1 = getPreferredOriFromOSP(R1);
%         ori2 = getPreferredOriFromOSP(R2);
%         mm = round mean([ori1, ori2]);
%         if  <
        
        if ~isempty(shft)
            R1 = shiftOSP(R1, shft);
            R2 = shiftOSP(R2, shft);
            Rmu = shiftOSP(Rmu, shft);
        end
        
        mR1 = mean(R1, 3);     mR2 = mean(R2, 3);      mRmu = mean(Rmu, 3);  
        fR1 = mR1/max(mR1(:)); fR2 = mR2/max(mR2(:));  fRmu = mRmu/max(mRmu(:));      
        OSP1.R = fR1;          OSP2.R = fR2;           OSPmu.R = fRmu;
        
        [OSP_1x2, OSP_min] = deal(OSP1);
        R1x2 = fR1 .* fR2;
        Rmin = min(fR1, fR2); 
        [tmp, idx_maxmin] = maxElement(Rmin);
             
        showAll = false;
%         if showAll

       col_min = 'r'; col_min_nm = 'red';
       col_1x2 = 'm'; col_1x2_nm = 'magenta';
       col_mu = 'k';  col_mu_nm = 'black';

       figure(1); clf;
       fsize1 = 15;
       fsize2 = 13;
       h_ax1 = subplot(2,3,1); imageOSP(OSP1, 'mean:ph', 'SOP', 'nolabels'); title('\color{blue}\bfR1_{ }', 'fontsize', fsize1);
       h_ax2 = subplot(2,3,2); imageOSP(OSP2, 'mean:ph', 'SOP', 'nolabels'); title('\color{darkgreen}\bfR2_{ }', 'fontsize', fsize1);
       h_ax3 = subplot(2,3,3); imageOSP(OSP1, 'mean:ph', 'SOP', 'nolabels'); 
       h_tit3 = title(['\color{' col_min_nm '}\bfmin(R1, R2)_{ }'], 'fontsize', fsize1);       
       h_im3 = get(h_ax3, 'children'); set(h_im3, 'cdata', R1x2');
%        if showAll  % make an extra plot for minR1,R2
%            h_ax4 = subplot(2,4,4); imageOSP(OSP1, 'mean:ph', 'SOP', 'nolabels'); 
%            title('min(R1,R2)');
%            caxis(h_ax4, [0 1]);
%            h_im4 = get(h_ax4, 'children'); 
%            set(h_im4, 'cdata', Rmin');
%        else
           p3 = get(h_ax3, 'position');
           h_cb3 = colorbar;
           set(h_ax3, 'position', p3);           
%        end
       caxis(h_ax1, [0 1]);
       caxis(h_ax2, [0 1]);
       caxis(h_ax3, [0 1]);       

       ylabel(h_ax1, 'Spat. Freq', 'fontsize', fsize2);
       xlabel(h_ax1, 'Orientation', 'fontsize', fsize2);       
       
%        [tmp, oriSp_ind] = maxElement(mean(OSP_1x2.R, 3));
%        [ori_ind, sp_ind] = elements(oriSp_ind);
       
       ax = 'ori';
       switch ax
           case 'ori'
               xs = 0:5:180;
               x_tick = 0:45:180;
               w = @(x) wrp(x);
               dim = 2;
           case 'spf'
               xs = 1:10;
               x_tick = 1:10;
               w = @(x) x;
               dim = 1;
       end
       
%        y1 = fR1(:,sp_ind); y2 = fR2(:,sp_ind); 
        idx = idx_maxmin(2);
%        y1 = mean(fR1, dim); y2 = mean(fR2, dim); 
       y1 = fR1(:,idx); y2 = fR2(:,idx); 
       
       y1x2 = y1.*y2; ymin = min(y1, y2)-.03;
       yMU  = max(fRmu,[], dim);
       [y1x2_max, ind_y1x2max] = max(y1x2);
       [ymin_max, ind_yminmax] = max(ymin);
       [yMU_max, ind_yMUmax] = max(yMU);
       

       
       h_axB = subplot(2,6,8:11); 
       h_l = plot(xs, w(y1), 'b.-', ...
                  xs, w(y2), 'g.-', ...
                  xs, w(ymin), [col_min ':'], ...
                  xs(ind_yminmax), ymin_max, [col_min 's'], ...
                  xs, w(y1x2), [col_1x2 ':'], ...
                  xs(ind_y1x2max), y1x2_max, [col_1x2 'o'], ...
                  xs, w(yMU),  [col_mu ':'], ...
                  xs(ind_yMUmax), yMU_max, [col_mu 'd']);
        idx_min = [3,4];
        idx_1x2 = [5,6];
        idx_mu  = [7,8];
%        set(h_l(6), 'markerfacecolor', 'darkgreen');
       set(h_l(idx_min(2)), 'markerfacecolor', col_min);
       set(h_l(idx_1x2(2)), 'markerfacecolor', col_1x2);
       set(h_l(idx_mu(2)), 'markerfacecolor', col_mu);
       set(h_axB, 'activePositionProperty', 'outerposition');
       xlabel('Orientation (at preferred spatial frequency)',  'fontsize', fsize2);       
       
        set(gca, 'xtick', x_tick);
        xlim([xs(1) xs(end)]);
        p = get(h_axB, 'position');
        hLeg = legend({'R1', 'R2'}, 'location', 'bestoutside');
        set(h_axB, 'position', p);            
        ylim([-.1, 1]);
        % fix in place
%         set([h_ax1, h_ax2, h_ax3, h_axB], 'units', 'pixels');
        

        if showAll
            
            
        else
            % first figure (no third plot)
            set([h_l(3:end); h_ax3; h_im3], 'visible', 'off')

            % figure for min(R1,R2)
            set(h_im3, 'cdata', Rmin');
            caxis(h_ax3, [0 1]);

            set([h_l(idx_min); h_ax3; h_im3], 'visible', 'on');
            set(h_l(idx_min(1)), 'linewidth', 3);  set(h_l(idx_min(2)), 'markersize', 8);            
            legend({'R1', 'R2', 'min(R1,R2)', 'max min(R1, R2)'}, 'location', 'bestoutside');
            3;
            
            return;
            
%             set(h_l(3:end), 'linewidth', 1);  %set(h_l(6), 'markersize', 1);            

            % figure for R1xR2
            set(h_l(idx_min), 'visible', 'off');
            
            set([h_l(idx_1x2); h_ax3; h_im3], 'visible', 'on');
            set(h_l(3:end), 'linewidth', 1);  set(h_l(6:end), 'markersize', 3);
            set(h_l(idx_1x2(1)), 'linewidth', 3);  set(h_l(idx_1x2(2)), 'markersize', 8);            
            set(h_tit3, 'string', ['\color{' col_min_nm '}\bfR1xR2_{ }']);                        
            legend({'R1', 'R2', 'R1 x R2', 'min(R1,R2)'}, 'location', 'bestoutside');
            3;

        
            % figure for R_mu
            set(h_im3, 'cdata', fRmu');
            caxis(h_ax3, [0 1]);
            set(h_l([5,8]), 'visible', 'on');
            set(h_l(3:end), 'linewidth', 1);  %set(h_l(7), 'markersize', 1);            
            set(h_l(5), 'linewidth', 3);  set(h_l(8), 'markersize', 8);
            set(h_tit3, 'string', ['\color{' col_mu_nm '}\bfR_{MU}']);           
            legend({'R1', 'R2', 'R1 x R2', 'min(R1,R2)', 'R_{MU}'}, 'location', 'bestoutside');
            3;
        end
%         set(h_l(5), 'linewidth', 3);  set(h_l(6), 'markersize', 8);
%         caxis(h_ax3, [0 1]);                

        
        
        
    end








end



%{
AA = [...
        2671        2671           1           3
        2697        2697           2           3
        2771        2771           1           3
        2771        2771           2           3
        2771        2771           2           4
        2771        2771           3           5
        2771        2771           4           5
        2801        2801           2           4
        2817        2817           2           3
        2829        2829           1           7
        2829        2829           2           7
        2829        2829           3           7
        2861        2861           2           8
        2861        2861           3           8
        2883        2883           1           2
        2883        2883           1           3
        2883        2883           2           4
        2919        2919           1           2
        2919        2919           1           3
        2919        2919           1           4
        2919        2919           2           4
        2919        2919           3           5
        3057        3057           2           6
        3057        3057           2           8
        3057        3057           4           8
        3057        3057           6           8
        3085        3085           3           5
        3085        3085           3           6
        3085        3085           3           9
        3085        3085           5           6
        3085        3085           6           7
        3085        3085           6           9
        3085        3085           8           9
        4462        4462           1           4
        4462        4462           2           3
        4462        4462           2           4
        4470        4470           1           4
        4470        4470           1           5
        4470        4470           2           5
        4470        4470           3           4
        4476        4476           2           4
        4476        4476           2           5
        4482        4482           1           2
        4482        4482           1           5
        4482        4482           1           6
        4482        4482           2           6
        4494        4494           2           3
        4498        4498           1           5
        4498        4498           2           5
        4498        4498           4           5
        4538        4538           1           4
        4538        4538           2           4
        4712        4712           1           2
        4716        4716           4           8
        4722        4722           2           3
        4722        4722           3           4
        4784        4784           3           7
        4784        4784           4           7
        4874        4874           2           3
        4990        4990           2           4
        4990        4990           3           4
        4990        4990           4           5
        4998        4998           1           3
        5048        5048           3           4
        5114        5114           1           4
        5114        5114           1           6
        5114        5114           3           6
        5158        5158           1           3
        5158        5158           1           6
        5158        5158           3           5
        5158        5158           3           7
        5158        5158           4           6
        5170        5170           2           3
        5176        5176           2           4
        5214        5214           3           4
        5226        5226           1           3
        5226        5226           2           3
        5236        5236           2           5
        5256        5256           2           5
        5256        5256           3           5
        5262        5262           1           4
        5262        5262           3           4
        5262        5262           3           9
        5262        5262           4           9
        5270        5270           2           4
        5276        5276           1           3
        5294        5294           2           5
        5300        5300           1           3
        2022        2022           1           2
        2246        2246           3           4
        2346        2346           1           5
        2346        2346           2           4
        2346        2346           2           5
        2346        2346           2           6
        2346        2346           4           6
        2486        2486           1           2
        2486        2486           1           6
        2486        2486           1           8
        2486        2486           2           4
        2486        2486           2           6
        2486        2486           2           8
        2486        2486           6           8
        2548        2548           1           4
        2548        2548           1           5
        2548        2548           1           7
        2548        2548           2           6
        2548        2548           2           7
        2548        2548           4           7
        1648        1648           1           2
        1699        1699           1           3
        1699        1699           2           5
        1699        1699           3           4
        1707        1707           2           3
        1707        1707           2           4
        1725        1725           3           4
        1735        1735           3           4
        1735        1735           3           5
        1792        1792           1           4
        1792        1792           2           3
        1792        1792           2           4
        1792        1792           2           5
        1792        1792           3           4
        1792        1792           4           5
        1822        1822           1           2
        2256        2256           1           4
        2256        2256           2           4
        2270        2270           1           2
        2270        2270           1           3
        2270        2270           1           4
        2270        2270           1           5
        2270        2270           1           7
        2270        2270           2           3
        2270        2270           2           5
        2270        2270           2           7
        2270        2270           2           8
        2270        2270           3           4
        2270        2270           3           5
        2270        2270           3           7
        2270        2270           4           5
        2270        2270           4           7
        2270        2270           5           8
        2292        2292           1           4
        2292        2292           3           4
        2292        2292           3           5
        2292        2292           4           5
        2292        2292           4           6
        2292        2292           4           7
        2304        2304           1           3
        2304        2304           1           4
        2304        2304           2           4
        2334        2334           1           3
        2334        2334           1           4
        2334        2334           1           5
        2334        2334           1           6
        2334        2334           1           7
        2334        2334           2           7
        2334        2334           3           5
        2334        2334           3           7
        2462        2462           1           3
        2514        2514           2           3
        2514        2514           2           4
        2514        2514           2           5
        2514        2514           2           6
        2514        2514           2           7
        2514        2514           2           9
        2514        2514           3           6
        2514        2514           5           9
        2432        2432           1           2
        2432        2432           3           4
        2500        2500           1           4
        2500        2500           1           7
        2500        2500           2           7
        2500        2500           3           4
        2500        2500           3           6
        2500        2500           3           7
        2500        2500           4           5
        2500        2500           5           6
        2500        2500           6           7
        ];

for i = 1:length(AA)
    B = AA(i,2:4);
    fig_locations(B);
end

%}