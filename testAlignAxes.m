h = zeros(3,3);
for i = 1:9
    h(i) = subplot(3,3,i);
    bar(1:5,rand(1,5))
end
h = h';


for i = 1:9
    p = get(h(i), 'position');
    p(1:4) = p(1:4) + randn(1,4)*.01;
    set(h(i), 'position', p)
end