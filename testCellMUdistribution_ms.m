function testCellMUdistribution_ms

    % generate nD random vectors, all with a certain angle from a central
    % vector. compute correlation between the random vectors.
    fmt = '%.3f';
    str_mn_std = @(x) sprintf( [fmt ' \\pm ' fmt ' (%d)'], mean(x), stderr(x), length(x));
    nDim = 5;
    
    doIntermediatePlots = false;
    evenDistrib = false;
    nBins = 30;
    subtractMean = false;
    
    cc_fun = @(x,y) sum(x(:).*y(:));   % assume all vectors are of length 1, so don't need to normalize
    
    ccMUs_init = [0:.01:1];
%     ccMUs_init = .5;
    nccMUs = length(ccMUs_init);

    nCells = 50;
    cc_pairs = nchoosek(1:nCells, 2); nPairs = size(cc_pairs,1);
    fprintf('nCells = %d, nPairs = %d\n', nCells, nPairs);
    
    nBiases = 10;
    MUbiases = linspace(0, pi/4, nBiases);
    as = zeros(1,nBiases);
    cis = zeros(2,nBiases);
    
    progressBar('init-', nccMUs*nBiases);
    
    for bi = 1:nBiases
    
        MUbias = MUbiases(bi);

        corrCellMU_m     = zeros(1, nccMUs);
        corrCellMU_sem   = zeros(1, nccMUs);
        corrCellCell_m   = zeros(1, nccMUs);
        corrCellCell_std = zeros(1, nccMUs);
        corrCellCell_sem = zeros(1, nccMUs);

        for ci = 1:nccMUs
            progressBar;
            theta = acos(ccMUs_init(ci));

            MU_c = [1, zeros(1,nDim-1)];
            if subtractMean
                MU_c = MU_c - mean(MU_c);
                MU_c = MU_c/mag(MU_c);
            end

            MU_p = NDcart2pol(MU_c);
            MU_p_biased = MU_p;
            MU_p_biased(2) = MU_p_biased(2) + MUbias;
            MU_c_biased = NDpol2cart(MU_p_biased);

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

                if ~subtractMean
                    assert(abs(theta2-theta) < 1e-10);
                end
            end


            corrCellCells = zeros(1,nPairs);
            corrCellMUs = zeros(1,nCells);
    %         corrCellCells2 = zeros(1,nPairs);
            for i = 1:nPairs
                cell1 = cell_coords(:,cc_pairs(i,1)); cell2 = cell_coords(:,cc_pairs(i,2));
                corrCellCells(i) = cc_fun(cell1, cell2);
            end        
            for i = 1:nCells
                cell1 = cell_coords(:,i);
                corrCellMUs(i) = cc_fun(cell1, MU_c_biased);            
            end

            if doIntermediatePlots
                figure(1); clf;
                plot3([0, MU_c(1)], [0, MU_c(2)], [0, MU_c(3)], 'b-'); hold on;
                for i = 1:nCells
                    plot3([0, cell_coords(1,i)], [0, cell_coords(2,i)], [0, cell_coords(3,i)], 'r-');
                end
                axis equal;
                zeroaxes;
                title(sprintf('cc = %.3f', ccMUs_init(ci)));

                figure(2); clf;
                hist(corrCellCells, nBins);
                title( str_mn_std (corrCellCells));
                xlim([-1 1]);
                drawnow;
            end

            corrCellMU_m(ci) = mean(corrCellMUs);
            corrCellMU_sem(ci) = stderr(corrCellMUs);        
            corrCellCell_m(ci) = mean(corrCellCells);        
            corrCellCell_std(ci) = std(corrCellCells);        
            corrCellCell_sem(ci) = stderr(corrCellCells);        
        end
        

        if exist('MUbias', 'var')
            p = polyfit(corrCellMU_m, corrCellCell_m, 2);
            [a,resid,J,sigma]= nlinfit(corrCellMU_m, corrCellCell_m, @(beta, x) beta * x.^2, 1);
            ci = nlparci(a,resid,'covar',sigma, 'alpha', .33);
            if nBiases > 1
                as(bi) = a;
                cis(:,bi) = ci;
            end        
            
        end        
        if (nccMUs > 1) % && (nBiases == 1)

            figure(3); clf;
            xyerrorbar(corrCellMU_m, corrCellCell_m, corrCellMU_sem, corrCellCell_sem, 'bo');  hold on;
            
            fprintf('y = %.2f*x.^2\n', a);
            fprintf('y = %.2f*x.^2 + %.2f *x + %.2f\n', p(1), p(2), p(3));

            %         errorbar(ccMUs_init, corrCellCell_m, corrCellCell_sem, 'b');  hold on
    %         errorbar(ccMUs_init, corrCellCell_m, corrCellCell_std, 'g'); 
    %         h1 = plot(ccMUs_init, ccMUs_init, 'k:');
            title(sprintf('Ndim = %d', nDim));
            xs = [-1:.01:1];
            h2 = plot(xs, xs.^2, 'r:');
            h3 = plot(xs, a*xs.^2, 'g:');
            axis([-1 1, -1 1]);
            xlabel('cc between cell and MU');
            ylabel('<cc> between cells');
            drawnow;
        end
        
    end
    progressBar('done'); 
    
    figure(10);
    errorbar(MUbiases, as, cis(1,:), cis(2,:));
    
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
