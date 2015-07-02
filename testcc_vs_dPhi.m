function testcc_vs_dPhi
    B = 3000;
    N = 40;
    
    ccs = zeros(1,B);
    dphis = zeros(1,B);
    
    r = 1;
    s = 1;
    phis = linspace(0,360, N+1); phis = phis(1:N);
    for b = 1:B
        r1 = r*randn(1,N)*r + s*sin(2*pi/N*[0:N-1] + rand*2*pi);
        r2 = r*randn(1,N)*r + s*sin(2*pi/N*[0:N-1] + rand*2*pi);
        ccs(b) = pearsonR(r1, r2);
        dphis(b) = abs( deltaPhi(phis, r1, r2) );
    end
    
    
    figure(88);
    plot(ccs, dphis, '.');
    
    axis ij;
    ylim([0 180]);
    set(gca, 'ytick', [0:45:180]);
    
end                