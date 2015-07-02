%%
M = [0.023, -0.011, 0.026];
S = [0.235, 0.196, .212];
sem = S/sqrt(615);

figure(55); clf;
MM = diag(M);

h1 = bar(1:3, MM, .6, 'stacked');
set(h1(1), 'faceColor', 'm')
set(h1(2), 'faceColor', 'r')
set(h1(3), 'faceColor', 'b');
ylabel('<MID-cc>');

hold on;
for i = 1:3
h2(i) = errorbar(i, M(i), sem(i), 'k.');
end
h_ax = gca;
set(h_ax, 'xtick', 1:3, 'xticklabel', {'All - All', 'Odd - Even', 'Odd - Odd'})
axis(axis);

i = 1;
set([h1, h2], 'visible', 'off')
set([h1(1:i), h2(1:i)], 'visible', 'on')
    
