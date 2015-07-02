function testUniqueInts

    i_max = 1e4;
    S = {randi(i_max, 1, 100), randi(i_max, 1, 200), {randi(i_max, 1, 100), randi(i_max, 1, 100)}, {{randi(i_max, 1, 100), randi(i_max, 1, 100)}, {randi(i_max, 1, 100), randi(i_max, 1, 100)}}};

    uS = uniqueInts(S);
    uS2 = uniqueInts2(S);

    assert(isequal(uS(:), uS2(:)));
    3;

end