
%%
v_s = cat(1, Si.shuffVal{:}); 
j_s = [Si.shuffJackStd];

v = Si.val;
j = Si.jackStd;

%%
v_s = v_s(:); j_s = j_s(:);
v = v(:); j = j(:);
%%
% jackstd_th = [.05 : .05 : 1.5];

jackstd_th = [1:1:20];

n = length(jackstd_th);
all_medians_sh = zeros(1,n);
all_means_sh = zeros(1,n);
all_medians = zeros(1,n);
all_means = zeros(1,n);
n_pr_sh = zeros(1,n);
n_pr = zeros(1,n);
tic;
for i = 1:n
    fprintf('*');
    idx_use_sh = j_s < jackstd_th(i);
    all_medians_sh(i) = nanmedian(v_s(idx_use_sh));
    all_means_sh(i) = nanmean(v_s(idx_use_sh));
    n_pr_sh(i) = nnz(idx_use_sh);
    
    idx_use = j < jackstd_th(i);
    all_medians(i) = nanmedian(v(idx_use));
    all_means(i) = nanmean(v(idx_use));

    n_pr(i) = nnz(idx_use);
    
end
%%
% figure(10); clf;
subplot(1,2,2); cla; hold on; 
plot(jackstd_th, all_medians, 'bo-')
plot(jackstd_th, all_medians_sh, 'ro-')

plot(jackstd_th, all_means, 'bs:')
plot(jackstd_th, all_means_sh, 'rs:')
legend({'Median (WS)', 'Median (Ctrl)', 'Mean (WS)', 'Mean (ctrl'}, 'location', 'best')
xlabel('Jackknife std error threshold'); ylabel('Median / Mean');
title('Drifting grating dphi');
%%
figure(11);
subplot(1,2,2); hold on; 
plot(jackstd_th, n_pr, 'bo-')
plot(jackstd_th, n_pr_sh, 'ro-')
legend({'Within-site', 'Control'});
xlabel('Jackknife std error threshold'); ylabel('N pairs');
title('Drifting grating dphi');
set(gca, 'yscale', 'log')

% toc;


