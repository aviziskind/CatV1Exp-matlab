function testIndSubConversion
   n = 25;

   for ind = 1:n^2
       [y,x] = ind2sub([n,n], ind);
    
       subs = mtxInd2sub(ind, n);
       if (subs(1) ~= y) || (subs(2) ~= x)
           error('1');
       end
       ind2 = mtxSub2ind(subs, n);
       if (ind2 ~= ind)
           error('2');
       end
   end
    
end