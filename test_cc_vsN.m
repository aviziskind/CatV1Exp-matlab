% test whether cc depends on sparseness of vectors

N = 50;
B = 50000;

ns = [1 2, 3, 4, 5, 6,  7, 8];%, 10, 15, 20];
% ns = [0];
ns = 1:N;


x2_sparse = 1;

doHists = 0;

ns(ns > N) = [];

nn = length(ns);

skw = zeros(1, nn);

X1_0 = zeros(N,B);
X2_0 = zeros(N,B);

binE = linspace(-1, 1, 40);
binC = binEdge2cent(binE);
for ni = 1:nn
    n = ns(ni);
    X1= X1_0;
    
    i_nz1 = ord( rand(N, B) );
    i_nz1 = i_nz1(1:n, :);
    
    
    for j = 1:B
        X1(i_nz1(:,j), j) = rand(n, 1);
    end

    if x2_sparse
        X2= X2_0;
        i_nz2 = ord( rand(N, B) );
        i_nz2 = i_nz2(1:n, :);    
        for j = 1:B
            X2(i_nz2(:,j), j) = rand(n, 1);
        end
        
    else
        x2 = rand(N, B);        
    end
    
    ccs = pearsonR_v(X1, X2);    
    ccs(ccs > 1) = 1;
    ccs(ccs < -1) = -1;
    
    if doHists
        figure(300+ns(ni)+x2_sparse*1000);
        binV = histcnt(ccs, binE); binV = binV / (sum(binV) * diff(binE(1:2)));
        bar(binC, binV, 1);
        title(sprintf('CCs: sparseness %.0f %%', n/N*100));
        xlabel('cc');
        xlim([-1, 1]);
        ylim([0, max(binV)*1.1]);
    end
    skw(ni) = skewness(ccs);
    
    fprintf('n = %d. m = %.2f, std = %.2f, skew = %.3f\n', n, mean(ccs), std(ccs), skw(ni))
    3;
    
end

% figure(14); 
% plot(ns, skw, 'bo-');

3;

return;
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
