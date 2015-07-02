function testGetCosineOverDC

    nPh = 50;
    ph = linspace(0, 2*pi, nPh+1); ph = ph(1:nPh);

    T = 2*pi;
    Fn = 1;
    
%     y1 = rectified(cos(ph));
    y1 = (-cos(2*ph))+.2;
    figure(7); clf;
    plot(ph, y1, 'o');
    xlim([0 T]);
    [C1, DC1] = getF1oDC(ph, y1, T, 2, 'cos');
    fprintf('C1 = %.2f, DC1 = %.2f C/D = %.2f\n', C1, DC1, C1/DC1);
    [F1, DC] = getF1oDC(ph, y1, T, 2);
    fprintf('F1 = %.2f, DC = %.2f F1/DC = %.2f\n', F1, DC, F1/DC);
   

end
   
