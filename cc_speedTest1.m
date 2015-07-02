function cc_speedTest1

    dims = [64, 64];
    mid1_odd = randn(dims);
    mid1_even = randn(dims);
    mid1_all = (mid1_odd + mid1_even)/2;
    mid1_oe = cat(3, mid1_odd, mid1_even);
    
    mid2_odd = randn(dims);
    mid2_even = randn(dims);
    mid2_all = (mid2_odd + mid2_even)/2;
    mid2_oe = cat(3, mid2_odd, mid2_even);
    
    S.mid1_odd = mid1_odd;
    S.mid2_odd = mid2_odd;

    S.mid1_oe = mid1_oe;
    S.mid2_oe = mid2_oe;
    
    B = 10000;
    tic;   for b = 1:B,   cc = pearsonR(mid1_odd, mid2_odd);              end;     t1 = toc;
    tic;   for b = 1:B,   cc = pearsonR(mid1_oe(:,:,1), mid2_oe(:,:,1));  end;     t2 = toc;
    tic;   for b = 1:B,   cc = pearsonR(S.mid1_odd, S.mid2_odd);              end;     t3 = toc;
    tic;   for b = 1:B,   cc = pearsonR(S.mid1_oe(:,:,1), S.mid2_oe(:,:,1));  end;     t4 = toc;
    
    fprintf('[%.2f, %.2f,  %.2f]  \n', t2/t1, t3/t1, t4/t1);
        

end