
nTrials = 1;

nSpks = [5, 25, 100, 200, 500];
% nSpks = [5, 25, 100, 200, 500];

Nnspks = length(nSpks);

H = zeros(nTrials, Nnspks);

for ni = 1:Nnspks
    for ti = 1:nTrials
        R = zeros(36,10);
        stim_ids = randi(360,1,nSpks(ni));
        [uStim, stimCount] = uniqueCount(stim_ids);
        R(uStim) = stimCount;
        r_entropy = sum ( entropy (R(:) / sum(R(:))));
        
        H(ti,ni) = r_entropy;        
        figure(ni); 
        imagesc(R); title(sprintf('Nspikes = %d. H = %.2f', nSpks(ni),r_entropy)); %         
        3;
    end
end
figure(100); clf;

plot( nSpks, mean(H,1), '.-');
drawHorizontalLine(log2(360), 'linestyle', ':');
ylabel('H')
xlabel('# of spikes')

% figure(1);
% imagesc(R);
