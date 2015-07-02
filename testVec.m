

nrep = 1000;

D = 8;
n = 10;
a = rand(D,n);
a(2,:) = 3*a(2,:);
b = bsxfun(@minus, a, mean(a,1));


r = sum(b,2);
r_shuff = zeros(D,nrep);
for rep_i = 1:nrep
    b_shuff = shuffle(b,1);
    r_shuff(:,rep_i) = sum(b_shuff,2);
end

figure(1); clf;
% quiver3(zeros(1,n), zeros(1,n), zeros(1,n), b(1,:), b(2,:), b(3,:), 'linewidth', 1);
hold on;
% quiver3(zeros(1,nrep), zeros(1,nrep), zeros(1,nrep), r_shuff(1,:), r_shuff(2,:), r_shuff(3,:), 'linewidth', 2, 'color', 'g');
% quiver3(zeros(1,nrep), zeros(1,nrep), zeros(1,nrep), r(1,:), r(2,:), r(3,:), 'linewidth', 2, 'color', 'r');
hold off;
axis equal

r_hat = r / norm(r);
r_hat_m = repmat(r_hat, 1, n);
projs = dot(b, r_hat_m);

% figure(2); clf;
% hist(projs);
% m = mean(projs);
% [h,pval] = ttest(projs);
% title( sprintf('m = %3.2f.  h = %d. p = %2.4g', m, h, pval));

figure(3); clf;
mag_rshuffs = normV(r_shuff, 1);
mag_r = norm(r);
hist(mag_rshuffs);
drawVerticalLine(mag_r, 'color', 'r');
drawVerticalLine(mean(mag_rshuffs), 'color', 'g');
drawVerticalLine(median(mag_rshuffs), 'color', 'g', 'linestyle', ':');
m = mean(mag_rshuffs);
p = mag_r/mean(mag_rshuffs);

% [h,pval] = ttest(mag_rshuffs, mag_r);
% title( sprintf('h = %d. p = %2.4g', h, pval));
title( sprintf('p = %2.4g', p));




