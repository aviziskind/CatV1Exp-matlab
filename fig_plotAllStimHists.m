% run this while in debug mode of  "calcPSTHforFlashGratingCells"
Gid = 2502;
cellId = 1;
bkgrSkip_ms=200;
firstLast = 1; %1 = first, 2 = last;

getfrm = @getMovieStimulusFrame;
getfrm('load', 'Gid', Gid);
[frameStimIds] = getStimulusFrameSequence(Gid, 'OSP');


[PSTH_bins, allPSTH_vals, spkTsRelToFrame_ms, bckgRate, meanFiringRate] = getIndividualPSTHs(Gid, cellId, 'OSP', bkgrSkip_ms);

% 2c. Take top N most preferred stimuli, and use their PSTHs to get
% a good estimate of the cell response profile.
N = 10; 

%     stimIndsInOrder = ord(max(allPSTH_vals, [], 1), 'descend');        
stimIndsInOrder1 = ord(mean(allPSTH_vals, 1).*max(allPSTH_vals, [], 1), 'descend');        
PSTH1 = allPSTH_vals(:,stimIndsInOrder1(1));

ylim2 = roundToNearest(max(PSTH1), 25, 'up');




if firstLast == 1    
    figure(2); clf;
    dots = sum(bsxfun(@times, PSTH1, allPSTH_vals),1);
    stimIndsInOrder = ord( dots, 'descend');        

    idx = [stimIndsInOrder(1:20 )];
    m = 4; n = 5;
elseif firstLast == 2
    figure(3); clf;
    stimIndsInOrder = ord(mean(allPSTH_vals, 1), 'descend');        
    idx = [stimIndsInOrder(2501:2505)];
    m = 1; n = 5;    
end
hs = zeros(m,n);
ind = 1;
bEdge = binCent2edge(PSTH_bins);
for i = 1:m
    for j = 1:n
%         subplot(m,n,(i-1)*n+j)
        subplot(m,n,ind)
        hs(i,j) = bar(PSTH_bins, allPSTH_vals(:,idx(ind)),1);
        set(hs(i,j), 'edgecolor', 'b', 'facecolor', 'b');        
        xlim([bEdge(1), bEdge(end)])
        set(gca, 'xtick', [])
        ylim([0 ylim2]);
        set(gca, 'units', 'pixels'); p = get(gca, 'position'); set(gca, 'units', 'normalized');        
        ur = [p(1)+p(3), p(2)+p(4)];
        s = 32;
        ll = ur - s;
        hStim = axes('units', 'pixels', 'position', [ll, s, s]);
        set(hStim, 'units', 'normalized');
        frmId = find(frameStimIds == idx(ind),1);
        imagesc(getfrm(frmId)); set(hStim , 'xtick', [], 'ytick', []); colormap('gray');
        
%         ylim([0 600])
         ind = ind + 1;
    end
end
matchAxes('Y', hs(:))