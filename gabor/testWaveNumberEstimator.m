function testWaveNumberEstimator

    k = 2.5;
    x = 0:.1:5;
    y = sin(k*x);
    figure(22); 
    plot(x,y, 'bo-'); hold on
    
    [maxfreq, amp, allFreqs, powers] = findStrongestFrequencies(x, y, 1);
    figure(23);
    plot(allFreqs,powers, '.-', maxfreq, amp, 'ro');
    k_est = 2*pi*maxfreq;
    
    figure(22);
    plot(x,sin(k_est*x), 'go-');
    hold off
    


end