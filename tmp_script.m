% doCalculationsForClustersInGroup(getAllGids('f'));
% clear all;
% doCalculationsForClustersInGroup(getAllGids('s'));
% clear all;
% doCalculationsForClustersInGroup(getAllGids('o'));
% clear all;

% load tmplist.mat
% curMatchDB(1);
% curGroupingType(1);
% doCalculationsForAllCells(Gid_list_o, cellId_list_o)
% doCalculationsForAllCells(Gid_list_s, cellId_list_s)
% doCalculationsForAllCells(Gid_list_f, cellId_list_f)

curMatchDB(0);
doCalculationsForClustersInGroup(getAllGids('o'));
clear all;
doCalculationsForClustersInGroup(getAllGids('s'));
clear all;
doCalculationsForClustersInGroup(getAllGids('f'));
clear all;

curMatchDB(1);
doCalculationsForClustersInGroup(getAllGids('o'));
clear all;
doCalculationsForClustersInGroup(getAllGids('s'));
clear all;
doCalculationsForClustersInGroup(getAllGids('f'));
clear all;



