function testShiftOriTuningCurve(dirPref)
    nOri = 36;
    oris = linspace(0,360, nOri+1);    
    oris = oris(1:nOri);
    dOri = diff(oris(1:2));

%     dirPref = 203;
    sigma = 38;
    y = doubleGaussOri(20, 10, sigma, 3, oris, dirPref);
    
    figure(1);
    subplot(2,1,1); cla;
    plot(oris, y, 'o-');
    set(gca, 'xtick', [0 : 45 : 360]);
    xlim([0 360])
    
    drawVerticalLine(dirPref);
    ori_start = -90;
    [oris_shifted, otc_shifted] = shiftOriTuningCurve(oris, y, dirPref, ori_start);
    ori_shft_itp = linspace(oris_shifted(1), oris_shifted(end), nOri*50);
    otc_shifted_itp = interp1(oris_shifted, otc_shifted, ori_shft_itp, 'spline');
    
    subplot(2,1,2); cla;
    plot(oris_shifted, otc_shifted, 'bo-'); hold on;
    plot(ori_shft_itp, otc_shifted_itp, 'g-', 'markersize', 2);    
    oris_lim = [ori_start(1), 360+ori_start];
    set(gca, 'xtick', oris_lim(1):45:oris_lim(2));
    xlim(oris_lim);
    drawVerticalLine(0);
    drawVerticalLine(ori_shft_itp(indmax(otc_shifted_itp)), 'color', 'r');
    
    
    3;

end
