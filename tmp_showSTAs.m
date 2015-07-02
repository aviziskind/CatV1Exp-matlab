    idx = [112   113   128   131   135   142   144   159   160   162   164   176   178   179   181   182   183   184   185   186   187   189   190   192   193   194   195   196   197   198   199   200   201   202   203   204   205   206   207   208   209   210   211   212];
    load('cellsGroups_GLFcuw8_movie_fg');
    grps = movieGroups_fg(idx);
   for gi = 1:length(grps)
       Gid = grps(gi).Gid;
       cellIds = grps(gi).cellIds;
       cellIds = cellIds(cellIds > 0);
       for ci = 1:length(cellIds)
           calculatePSTH_STAs_OSP_ForOneCell(Gid, cellIds(ci));
           3;
       end
   end
