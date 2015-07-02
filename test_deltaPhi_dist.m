function test_deltaPhi_dist(vals)

    nBins = 11;
    if nargin == 0 
        nPh = 60; %[40, 44, 60];
        ph = linspace(0, 360, nPh+1); ph = ph(1:nPh);

        T = 60000;
%         dphs = zeros(1,T);
%         ph1 = 0;
%         ph2_ind = randi(nPh,1,T);
%         for t = 1:T;
%             dphs(t) = circDist(ph1, ph(ph2_ind(t)), 360);
%         end 

        [vals, count] = deltaPhiNull(nPh, T);
        dphs = arrayfun(@(v,c) ones(1,c)*v, vals, count, 'un', 0);
        dphs = [dphs{:}];
        
        
    else
        nPh = 60;%length(unique(vals));
        ph = linspace(0, 360, nPh+1); ph = ph(1:nPh);

        dphs = vals;
    end
        
%     binCenters = linspace(0, 180, nPh/6+1); diff(binCenters(1:2))    
%     binEdges = binCent2edge(binCenters);
    
    binCent_indv = linspace(0, 180, nPh+1); 
    binEdges_indv = binCent2edge(binCent_indv);
    bin_n_indv = histcnt(dphs, binEdges_indv);
%     set(gca, 'xtick', binEdges)
    
    dphiBinE = @(nBin) binCent2edge( linspace(0, 180, nBin) );
    
    binEdges = dphiBinE(nBins);
    binCenters = binEdge2cent(binEdges);
%     binEdges(end) = binEdges(end)+1;
    binW = diff(binCenters(1:2))    

    figure(555);
    subplot(2,1,1);
    bin_n = histcnt(dphs, binEdges);
    bar(binCenters, bin_n, 1)
    xlim([binEdges(1), binEdges(end)])
    set(gca, 'xtick', binCenters)

    subplot(2,1,2);
    bar(binCent_indv, bin_n_indv, 1)
    xlim([binEdges_indv(1), binEdges_indv(end)])
    xlim([binEdges(1), binEdges(end)])
    set(gca, 'xtick', binEdges)
    
    
    % set(gca, 'xtick', binCenters);
    mean(dphs)

end