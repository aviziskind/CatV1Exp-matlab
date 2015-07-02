function xi = quadCenter(y)

    halfGrandSum = sum(y)/2;
    sequentialSum = 0;
    for xi = 1:length(y)
        sequentialSum = sequentialSum + y(xi);
        if (sequentialSum >= halfGrandSum)
            break;
        end
    end
    
end