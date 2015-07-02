function compareSmoothingMethods

figure(1); clf; hold on;
h_orig = plot(0,0);
h_noisy = plot(0,0, 'b.:');

h_g = plot(0,0, 'g:');
h_f = plot(0,0, 'r:');
h_a = plot(0,0, 'k:');
h_ax = gca;

ext = @(x) [x(:); x(end)+diff(x(1:2))];
wrp = @(x) [x(:); x(1)];

% 
% figure(2); clf; hold on
% h_cc_w_hp = plot(0,0, '.-');


function updatePlots(show, nx, w, k_hp, b_hp_frac, n_av, r_state, eta)
    vis = arrayfun(@(tf) iff(tf, 'on', 'off'), show, 'un', 0);
%     N = round(10.^nx);
    N = round(nx);
    rand('state', r_state);
    randn('state', r_state);
    x = 0:N-1;
    tc_orig = gaussSmooth(rand(N,1), N/8, [], 1);
    tc_noisy = tc_orig + randn(N,1)/eta;

    b_hp = k_hp / b_hp_frac;
    
    tc_gauss = gaussSmooth(tc_noisy, w, [], 1);
    tc_fermi = fermiSmooth(tc_noisy, k_hp, b_hp, 1);
    tc_alias = alias(tc_noisy, n_av);
    av_x = linspace(0, N-1, n_av);
    
    set(h_orig, 'xdata', ext(x), 'ydata', wrp(tc_orig), 'visible', vis{1});
    set(h_noisy, 'xdata', ext(x), 'ydata', wrp(tc_noisy), 'visible', vis{2});
    set(h_g, 'xdata', ext(x), 'ydata', wrp(tc_gauss), 'linewidth', 2, 'visible', vis{3});
    set(h_f, 'xdata', ext(x), 'ydata', wrp(tc_fermi), 'linewidth', 2, 'visible', vis{4});
    set(h_a, 'xdata', av_x, 'ydata', tc_alias, 'linewidth', 2, 'visible', vis{5});
    set(h_ax, 'xlim', [0, N]);
    
    ws = [0:.02:1];
    ks = zeros(size(ws));
%     for wi = 1:length(ws)
%         tc_gauss = gaussSmooth(tc_noisy, ws(wi), [], 1);
%         ks(wi) = fminsearch(   @(k) norm( [tc_gauss - fermiSmooth(tc_noisy, k, k/b_hp_frac, 1)])   , 5/(ws(wi)+.1) );
%     
%     	set(h_cc_w_hp, 'xdata', ws, 'ydata', ks);
%     end
    
end


show = { 'show', repmat({[false, true]},1,5), true(5,1), [], {'orig', 'noisy', 'gauss', 'fermi', 'alias'} };
%     showDepVars = { {{}, {}}


args = { show, {'nx', [10:10:150]}, {'w', [0:.5:100], 1}, ...
        {'k_hp', [0:.1:100], 3}, {'b_hp_frac', [1:1000], 100}, {'n_av', [1:30], 1}, {'r_state', [1:20]}, {'eta', [5:100], 30} };

manipulate(@updatePlots, args, 'FigId', 10)



end