function testCellMUdistribution

    % generate nD random vectors, all with a certain angle from a central
    % vector. compute correlation between the random vectors.
    fmt = '%.3f';
    str_mn_std = @(x) sprintf( [fmt ' \\pm ' fmt ' (%d)'], mean(x), stderr(x), length(x));
    nDim = 3;
    
    doIntermediatePlots = true;
    evenDistrib = false;
    nBins = 30;
    subtractMean = false;
    
    cc_fun = @(x,y) sum(x(:).*y(:));   % assume all vectors are of length 1, so don't need to normalize
    
    ccMUs = [0:.05:1];
%     ccMUs = .5;
    nccMUs = length(ccMUs);


    nCellsMax = 200;
    B = round(nthroot(nCellsMax, nDim-2));
    nCells = B^(nDim-2);
    
    pairs = nchoosek(1:nCells, 2); nPairs = size(pairs,1);
    fprintf('nCells = %d, nPairs = %d\n', nCells, nPairs);

    
    if evenDistrib
        siz = ones(1,nDim-2)*B;
        cell_phis = zeros(nDim-2,nCells);
        phis = linspace(0, 2*pi, B+1);
        thetas = linspace(0, pi, B+1);
        
        for ci = 1:nCells
            v = ind2subV(siz, ci);
            for ti = 1:nDim-3
                cell_phis(ti,ci) = thetas(v(ti));
            end
            cell_phis(nDim-2,ci) = phis(v(nDim-2));
        end
        
    end
    

    
    progressBar('init-', nccMUs);
    ccCells_m = zeros(1, nccMUs);
    ccCells_std = zeros(1, nccMUs);
    ccCells_sem = zeros(1, nccMUs);
    
    for ci = 1:nccMUs
        progressBar(ci);
        theta = acos(ccMUs(ci));
    
        MU_c = [1, zeros(1,nDim-1)];
        if subtractMean
            MU_c = MU_c - mean(MU_c);
            MU_c = MU_c/mag(MU_c);
        end
        
        MU_p = NDcart2pol(MU_c);
    
        cell_coords = zeros(nDim,nCells);

        for cell_i = 1:nCells
            c1_p = MU_p;
            c1_p(2) = c1_p(2)+theta;
            if ~evenDistrib
                c1_p(3:nDim-1) = c1_p(3:nDim-1)+rand(1,length(3:nDim-1))*pi;
                c1_p(nDim) = c1_p(nDim)+rand*2*pi;
            else
                c1_p(3:nDim) = c1_p(3:nDim)+cell_phis(:,cell_i);                
            end
            c1_c = NDpol2cart(c1_p);
            if subtractMean
                c1_c = c1_c - mean(c1_c);
                c1_c = c1_c/mag(c1_c);
            end            
            cell_coords(:,cell_i) = c1_c(:);                            

            theta2 = angleBetweenVectors(c1_c', MU_c');
            assert(abs(mag(c1_c)-1) < 1e-10);

            assert(abs(theta2-theta) < 1e-10);
        end
                    
        
        cell_ccs = zeros(1,nPairs);
%         cell_ccs2 = zeros(1,nPairs);
        for i = 1:nPairs
            cell1 = cell_coords(:,pairs(i,1)); cell2 = cell_coords(:,pairs(i,2));
            cell_ccs(i) = cc_fun(cell1, cell2);
%             cell_ccs2(i) = cc_fun2(cell1, cell2);
%             assert(abs(cell_ccs(i) - cell_ccs2(i))<1e-10);
        end

        if doIntermediatePlots
            figure(1); clf;
            plot3([0, MU_c(1)], [0, MU_c(2)], [0, MU_c(3)], 'b-'); hold on;
            for i = 1:nCells
                plot3([0, cell_coords(1,i)], [0, cell_coords(2,i)], [0, cell_coords(3,i)], 'r-');
            end
            axis equal;
            zeroaxes;
            title(sprintf('cc = %.3f', ccMUs(ci)));

            figure(2); clf;
            hist(cell_ccs, nBins);
            title( str_mn_std (cell_ccs));
            xlim([-1 1]);
            drawnow;
        end
        
        ccCells_m(ci) = mean(cell_ccs);
        ccCells_std(ci) = std(cell_ccs);        
        ccCells_sem(ci) = stderr(cell_ccs);        
    end
    progressBar('done');
    if nccMUs > 1

        figure(3); clf;
        
        errorbar(ccMUs, ccCells_m, ccCells_sem, 'b');  hold on
        errorbar(ccMUs, ccCells_m, ccCells_std, 'g'); 
%         h1 = plot(ccMUs, ccMUs, 'k:');
        title(sprintf('Ndim = %d', nDim));
        xs = [-1:.01:1];
        h2 = plot(xs, xs.^2, 'r:');
        axis([-1 1, -1 1]);
%         xlabel('cc between cell and MU');
%         ylabel('<cc> between cells');
    end
end    
    
%     x1 = r*sin(t1)*sin(t2)*sin(t3)*sin(t4);
%     x2 = r*sin(t1)*sin(t2)*sin(t3)*cos(t4);    
%     x3 = r*sin(t1)*sin(t2)*cos(t3);
%     x4 = r*sin(t1)*cos(t2);
%     x5 = r*cos(t1);
%     X = [x1 x2 x3 x4 x5];

% end

% function R = NDpol2cart(x1, x2, x3, x4)  % 4d 
% 
%     r = norm([x1, x2, x3, x4])
%     
% 
%     x1 = r*sin(t1)*sin(t2). sin(tn-1)sin(tn)
%     x2 = r*sin(t1)*sin(t2). sin(tn-1)cos(tn)
%     
%     x3 = r*sin(t1)*sin(t2). cos(tn-1)
%     x2 = r*sin(t1)*cos(t2)
%     x1 = r*cos(t1)
% 
% end


%     x1 = r*sin(t1)*sin(t2). sin(tn-1)sin(tn)
%     x2 = r*sin(t1)*sin(t2). sin(tn-1)cos(tn)
%     
%     x3 = r*sin(t1)*sin(t2). cos(tn-1)
%     x2 = r*sin(t1)*cos(t2)
%     x1 = r*cos(t1)
