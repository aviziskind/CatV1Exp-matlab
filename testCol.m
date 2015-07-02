Ncol = 200;
Nbr = 200;
cols0 = jet(Ncol);

A = zeros(Ncol, Nbr, 3);
B = zeros(Ncol, Nbr, 3);
for i = 1:Ncol
    curCol = cols0(i,:);
    for j = 1:Nbr
        A(i,j,:) = min(curCol * [1-(j-1)/Nbr], 1);
        B(i,j,:) = min([curCol + .1] * [1.3 + (j-1)/Nbr], 1);
    end
end

figure(3); imagesc(A);
figure(4); imagesc(B);
