
%%

doPlot = 1;

allG = [S.allCells.Gid] ;
allC = [S.allCells.cellId] ;

if doPlot == 1
    for i = 18%1:length(GC)
        Gid = GC(i,1);
        cellId1 = GC(i,2);
        cellId2 = GC(i,3);

        ci1 = find(allG == Gid & allC == cellId1);
        ci2 = find(allG == Gid & allC == cellId2);

        R1 = S.allCells(ci1).R;
        R2 = S.allCells(ci2).R;

        R1p = mean(R1,3);  R1p = R1p/max(R1p(:));
        R2p = mean(R2,3);  R2p = R2p/max(R2p(:));
        [m1, i1] = maxElements(min(R1p, R2p), ntop);

        figure(1); clf;
        imageOSP(R1, 'mean:ph'); colorbar off;
        title(sprintf('Site # %d, cell # %d', Gid, cellId1))

        figure(2); clf;
        imageOSP(R2, 'mean:ph'); colorbar off;
        title(sprintf('Site # %d, cell # %d', Gid, cellId2))

        %%
        wrp = @(x) x([1:end, 1]);
        for j = 2;
            y1 = squeeze(R1(i1(j,1), i1(j,2), :)); 
            y2 = squeeze(R2(i1(j,1), i1(j,2), :)); 

            ph = 0:45:360;
            figure(3); clf; h3_orig = plot(ph, wrp(y1), 'bo-'); xlim([0, 360]); set(gca, 'xtick', 0:90:360); title('Cell #1'); hold on; 
            figure(4); clf; h4_orig = plot(ph, wrp(y2), 'bo-'); xlim([0, 360]); set(gca, 'xtick', 0:90:360); title('Cell #2'); hold on;
            3;
        end
%%
        set([h3_orig, h4_orig], 'linewidth', 3, 'marker', 'o');
        
        3;
        cnt = 0;
        while true
            shft1s = randperm(8-1); %[1  -1 2  4];
            shft2s = randperm(8-1); %[-1, 3 1 -2];

            
            d1 = [diff(shft1s), shft1s(1)-shft1s(7)];
            d2 = [diff(shft2s), shft2s(1)-shft2s(7)];            
            d12 = shft1s - shft2s;
            
            n12 = nnz(d12==0);
            nd11 = nnz(abs(d1) == 1);
            nd21 = nnz(abs(d1) == 1);
            ndd1 = nnz(diff(shft1s, 2) == 0);
            ndd2 = nnz(diff(shft2s, 2) == 0);
            
            n_same = nnz(diff(n12) == 0);
            
            ok = (n12 < 1) && (nd11 < 1) && (nd21 < 1) && (ndd1 < 1) && (ndd2 < 1) && (n_same < 1);
            cnt= cnt+1;
            if ok
                break;
            end            
            
        end        
        
        nshft = length(shft1s);
        set([h3_orig, h4_orig], 'linewidth', 1, 'linestyle', ':');
        figure(3); h3_shift = plot(ph, wrp(y1), 'ro-', 'linewidth', 3);
        figure(4); h4_shift = plot(ph, wrp(y2), 'ro-', 'linewidth', 3);
        for j = 1:nshft
            y1_shft = circshift(y1, shft1s(j));
            y2_shft = circshift(y2, shft2s(j));
            set(h3_shift, 'ydata', wrp(y1_shft)); xlim([0, 360]); set(gca, 'xtick', 0:90:360)
            set(h4_shift, 'ydata', wrp(y2_shft)); xlim([0, 360]); set(gca, 'xtick', 0:90:360)
            3;
        end
        
            
        3;
    end
end

if doPlot == 2

    Gid_c = 5158;
    cellId_c = 6;    
    ci_c = find(allG == Gid_c & allC == cellId_c);
    R_c = S.allCells(ci_c).R;
    Rp_c = mean(R_c, 3);
    Rp_c = Rp_c/max(Rp_c(:));

    allM = zeros(1, length(allG));
    for ci = 1:length(allG)
        
        Gid = S.allCells(ci).Gid;
        cellId = S.allCells(ci).cellId;
                
        R = S.allCells(ci).R;
        
        if size(R,3) ~= 8
            continue;
        end

        Rp = mean(R, 3);
        Rp = Rp / max(Rp(:));
        
        [m1, i1] = maxElements(min(Rp_c, Rp));
        allM(ci) = m1;        
        if m1 > .07
            continue;
        end

        figure(10); clf;
        imageOSP(R, 'mean:ph'); colorbar off;
        title(sprintf('Site # %d, cell # %d', Gid, cellId))
        
        
        y1 = squeeze(R(i1(1), i1(2), :)); y1 = y1([1:end, 1]);
        y2 = squeeze(R_c(i1(1), i1(2), :)); y2 = y2([1:end, 1]);
        

        ph = 0:45:360;
        figure(3); plot(ph, y1, 'bs-'); xlim([0, 360]); set(gca, 'xtick', 0:90:360)
        figure(4); plot(ph, y2, 'bs-'); xlim([0, 360]); set(gca, 'xtick', 0:90:360)
        3;
    end

    3;
    
end


    