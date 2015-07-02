function testShiftSpikeDistance

    sp1 = sort(10*randn(50,1)+50);
    sp2 = sort(10*randn(50,1)+50);

    shifts = -25:25;
    ds12 = zeros(size(shifts));
    ds21 = zeros(size(shifts));
    for si = 1:length(shifts);
        ds12(si) = spikeShiftDistance(sp1, sp2 + shifts(si), [], [], 1);
        ds21(si) = spikeShiftDistance(sp2 + shifts(si), sp1, [], [], 2);
    end
        
%     idx = indmin([ max(ds12, ds21) ]);
    figure(10); clf; hold on;
    plot(shifts, ds12, 'bo', shifts, ds21, 'go');
    plot(shifts, smooth(ds12,3), 'b.:', shifts, smooth(ds21,3), 'g.:');
    3;
    hold off
    


end
