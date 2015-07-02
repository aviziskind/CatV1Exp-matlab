N = length(allOSPs);
F1s = zeros(1,N);
j = 1;
for i = 1:length(allOSPs)
    OSP = allOSPs(i);
    ph = deg2rad(OSP.ph);
%     [ori_i, sp_i] = dealV(OSP.oriSp_maxR_av);
    [ori_i, sp_i] = dealV([randi(30),randi(10)]);
    R = OSP.R(ori_i, sp_i,:);
    if nnz(R)>=0
        R = squeeze(R);    
        R = R(randperm(length(R)));
        F1s(j) = getF1phase(ph, R, 2*pi);
        j = j+1;
    else
        3;
    end
end
F1s(j:end) = [];
figure(55);
hist(F1s, 100);