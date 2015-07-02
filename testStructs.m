function testStructs
    a = rand(1, 1000000);
    b = randn(1, 1000000);
    S_rec = struct('a', num2cell(a), 'b', num2cell(b))
    S_grp = struct('a', a, 'b', b)
    
    tic;
    a1 = [S_rec(:).a];
    toc;
    tic;
    a2 = [S_rec.a];
    toc;
    
    
    
end