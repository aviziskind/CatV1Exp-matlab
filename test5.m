Gids1 = [461   462   570   588   589];
cellIds = {[0 3 6 8 9], [0 1 3 5 6 8 9], [0 2 3 4 6], [0 1 5 7 8 9], [0 1 3 4 5 8 9]};
figure(2); clf;

c2 = cell(1,30);
c_ind = 1;
for Gid_i = 1:5
    for cellId = cellIds{Gid_i}
        c2{c_ind} = calculatePSTH_STAs_OSP_ForOneCell(Gids1(Gid_i),cellId);
        fprintf('\n');
        subplot(10,3,c_ind);
        imageOSP(c2{c_ind}.OSP.R, 'pref:ori', 'SPO', 'noLabels', 'noTicks');
        ylabel(sprintf('%d:%d', Gids1(Gid_i), cellId));
        c_ind = c_ind +1;
    end    
end