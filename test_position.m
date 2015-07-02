op = get(gca, 'outerposition');
p = get(gca, 'position');
t = get(gca, 'tightInset');
M = [-1,0,0,0;
     0,-1,0,0;
     1,0,1,0;
     0,1,0,1];
t2 = (M*t')';
annotation('rectangle', get(gca, 'outerPosition'), 'color', 'y')
annotation('rectangle', p, 'color', 'g')
pt = [p(1)-t(1), p(2)-t(2), p(3)+t(1)+t(3), p(4)+t(2)+t(4)];
pt2 = p + t2;
annotation('rectangle', pt2, 'color', 'r')

pt3 = op-t2;
annotation('rectangle', pt3, 'color', 'k')