G = [1 2 9;
    7 5 6;
    4 8 3];

    [myMax, pos1] = max(G);
    [Max,pos2] = max(myMax);
    col = pos2;
    row = pos1(col);