% test whether cc vs cc_p depends on skewness of vector values s(no)

mu1 = 5; sig1 = .2; 
mu2 = 10; sig2 = 1; 


n1 = 8;
n2 = 4;
B = 1000;
[ccs1, cc_ps1] = deal(zeros(1,B));
[ccs2, cc_ps2] = deal(zeros(1,B));
tic;
for b = 1:B
    
%     tc1a = exp(normrnd(mu1, sig1, n, 1));
%     tc1b = exp(normrnd(mu1, sig1, n, 1));
    tc1a = normrnd(mu1, sig1, n1, 1);
    tc1b = normrnd(mu1, sig1, n1, 1);
    [ccs1(b), cc_ps1(b)] = corr(tc1a(:), tc1b(:), 'type', 'pearson', 'tail', 'gt');    

    tc2a = normrnd(mu2, sig2, n2, 1);
    tc2b = normrnd(mu2, sig2, n2, 1);
    [ccs2(b), cc_ps2(b)] = corr(tc2a(:), tc2b(:), 'type', 'pearson', 'tail', 'gt');        
end
toc;

figure(444); plot(ccs1, cc_ps1, 'b.', ccs2, cc_ps2, 'go');
figure(445); plot(ccs1, -log10(cc_ps1), 'b.', ccs2, -log10(cc_ps2), 'go');
