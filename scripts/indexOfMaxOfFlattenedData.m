function i = indexOfMaxOfFlattenedData(x)

    dbug = false;

     a = max(x);
     N = length(x);
     indices = find(x==a);
     
     if dbug
         figure(103);
         plot(x, 'o'); hold on;
         plot(x(indices), 'go'); hold off
     end
 
     if length(indices) == 1
         i = indices;
     else
         numConsecutiveGroups = nnz(  diff(indices) > 1 ) + 1;
         if (numConsecutiveGroups == 1)
             i = round(mean(indices));
         elseif (numConsecutiveGroups == 2) && (indices(1) == 1) && (indices(end) == N)
             n = find(diff(indices) > 1) + 1;
             indices(n:end) = indices(n:end)-N;
             i = mod(round(mean(indices)),N);
         else
             % warning('more than 2 groups'); %#ok<WNTAG>
             i = round(mean(indices));
             if (x(i)) ~= a
                 i = indices(1);
             end
         end
     end

end