function testspmd

%     matlabpool;
    spmd
        fprintf('Here is #%d\n', labindex)
        eig(rand(1000));
        fprintf('Here is #%d\n', labindex)
    end
%     matlabpool close;

end