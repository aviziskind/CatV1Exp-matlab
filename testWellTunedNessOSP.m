function testWellTunedNessOSP

    S = load('flashGratingData.mat');
    
    allOSPs = S.allOSPs;
    
    N = length(allOSPs);
%     N = min(N, 100);
    Us1 = zeros(1, N);
    Us2 = zeros(1, N);
    wrapDims = 1;
    
%     progressBar('init-', N);
%     for i = 1:N
%         progressBar(i);
%         sOSP = sum(allOSPs(i).OSP, 3);
%         [U1, U1max] = gravitationalPotential(sOSP, wrapDims, 1);
%         Us1(i) = U1 / U1max;
%         [U2, U2max] = gravitationalPotential(sOSP, wrapDims, 2);
%         Us2(i) = U2 / U2max;
%     end
%     save U.mat Us1 Us2;
    load U.mat;

    U = Us1;
    
    N = round(length(Us)/2);
    U = U(1:N);

    U_scaled = zeros(size(U));
    for i = 1:N
        U_scaled(i) = U(i) * sqrt( max(allOSPs(i).OSP(:)) );
    end
    
    U = U_scaled;
    
    inds = ord(U);
    gridSubPlot(3,4, 100);
    for j = inds
        S = allOSPs(j);
        sOSP = sum(S.OSP, 3);
        gridSubPlot;
        imagesc(sOSP); colorbar;
        set(gca, 'xtick', [], 'ytick', []);
        title(sprintf('(%d) [%d, %d], %d', j, S.GroupId, S. cellId, U(j)) );        
    end
    
        
end
    
    
    