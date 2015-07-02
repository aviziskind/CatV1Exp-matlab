function [xbins nbins] = getBestBinsForTwoDistributions(x1, x2)
    % we choose the number of bins to maximize the mutual information
    % between the two distributions.
    nbinsRange = [5:100];
    nbins = 20;
   
    infos = zeros(1,    length(nbinsRange));
    dists = zeros(1,    length(nbinsRange));
    for T = 1:2
        for bi = 1:length(nbinsRange)
            nbins = nbinsRange(bi);
            [tmp xbins] = hist([x1; x2], nbins);
            P1 = hist(x1, xbins);
            P2 = hist(x2, xbins);

            P1 = P1 / sum(P1);
            P2 = P2 / sum(P2);

            if T == 2
                figure(8);
                bar(xbins, P1, 'EdgeColor', 'b', 'FaceColor', 'none'); hold on
                bar(xbins, P2, 'EdgeColor', 'g', 'FaceColor', 'none'); hold off
                legend('P1', 'P2');
                set(h, 'xdata', nbinsRange(bi), 'ydata', infos(bi));
                pause(.5);
            end

            infos(bi) = mutualInfo(P1, P2);
            dists(bi) = 
        end
        if T == 1
            figure(9); clf;
            plot(nbinsRange, infos, '.-'); hold on
            h = plot(nbinsRange(1), infos(1), 'ro');
        end
    end
        
end