% function test3_cc_dot
N = 8;

rnk = @(x) ord(x);

for b = 1:1000
    x1 = randn(1,N);
    y1 = randn(1,N);

    shftR = @(x) x([end, 1:end-1]);

    ccs = zeros(1,N);
    dots = zeros(1,N);

    ccs(1) = pearsonR(x1, y1);
    dots(1) = dot(x1, y1);
    for i = 1:N-1
        y1 = shftR(y1); 
        ccs(i+1) = pearsonR(x1, y1);
        dots(i+1) = dot(x1, y1);        
    end
    [tmp, indmax_cc] = max(ccs);
    [tmp, indmax_dots] = max(dots);

%     if indmax_cc ~= indmax_dots
%     d1 = dot(x1,y1) > dot(x2,y2);
%     c1 = pearsonR(x1, y1) > pearsonR(x2, y2);
%     if xor(d1, c1)
%         plot(1:N, x1, 'bo-', 1:N, y1, 'gs-');
    if ~isequal(rnk(ccs), rnk(dots))
        3;
    end
end

