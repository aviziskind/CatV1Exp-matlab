
N = 64;
x = 1:N;
tc_orig = gaussSmooth(rand(N,1), 5);
tc_noisy = tc_orig + randn(N,1)*0;
figure(1); clf; hold on;
plot(x, tc_orig, ':');
plot(x, tc_noisy, 'g.');

% figure(2);
% [ks1, p1] = powerSpectrum(x,tc_orig);
% [ks2, p2] = powerSpectrum(x,tc_noisy);
% plot(ks1,p1, 'b', ks2,p2, 'g')

k_hp = 5;b_hp = .05;
% tc_fm_orig = fermiSmooth(tc_orig, k_hp, b_hp);
tc_fm_noisy1 = fermiSmooth(tc_noisy, k_hp, b_hp);
tc_fm_noisy2 = fermiSmooth(tc_noisy, k_hp, b_hp, 1);

% plot(x, tc_fm_orig, 'g:');
plot(x, tc_fm_noisy1, 'r');
plot(x, tc_fm_noisy2, 'g:');


% 1/(1+exp(-(k_hp-|k|)/Bhp))
% 
% B_hp = .05 k_hp
