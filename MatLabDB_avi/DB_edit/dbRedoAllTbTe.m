function dbRedoAllTbTe


         

    MY_METHOD = 1; ORIG_METHOD = 2;
    method = MY_METHOD;
    
%     DidsNeedTbTeFixing = ...
%         {[1622 1675 1676 1707 1708 1711 1716 1717 1719 1744 1749 1755 1760 1766 1768 1770 1781 1790 1796 1806 1985 2070 2362 2639 2849 2996 3006 3007 3029 3310 3326], ...
%          [246, 622 634 703 742 763 1221 1471 1476 1484 1742 1860 1943 2052], ...
%          [76 191 337 419 467 474 661 729 1072 1075 1309 1344 1469 1517 1594 1673 1677 1680 1697 1703 1710 1714 1726 1754 1761 1773 1786 1787 1789 1794 1808 1833 1892 1918 1998 2339 ...
%                 2366 2395 2445 2449 2460 2621 2702 2726 2813 2855 2887 2913 2922 3008 3026 3028 3054 3259 3274 3659 3747 3768], ...
%          [1307 1308 1312 1313 1332]};
    
errorDids = [500 621 627 629 633 636 644 647 652 656 667 672 676 678 680 683 686 688 692 693 696 699 702 705 708 711 712  ...
714 715 718 721 728 732 734 735 738 740 745 746 751 761 764 767 770 771 773 777 780 781 784 786 788 793 796 ...
798 799 802 806 809 820 823 826 831 834 1010 1016 1017 1019 1023 1024 1026 1029 1031 1035 1036 1041 1043 ... 
1053 1057 1060 1063 1064 1067 1082 1090 1096 1103 1112 1113 1204 1205 1208 1403 1405 1415 1416 1419 1420 ...
1421 1428 1429 1432 1435 1436 1437 1441 1444 1445 1446 1453 1456 1463  1464 1469 1474 1477 1493 1497 1499 ...
1501 1510 1511 1513 1530 1532 1533 1547 1552 1553 1554 1561 1566 1569 1594 1599 1602 1618 1620 1728 1710  ...
1739 1743 1759 1750 1747 1748 1802 1803 1805 1751 1758 1769 1773 1775 1786 1788 1791 1793 1795];     
    DidsNeedTbTeFixing = errorDids;

        
    hnd = dbOpenExpDb;
    stimTypes = {'movie', 'noise', 'grating', 'mseq'};
    stimTypesToDo = [4];
    %gather all correct GroupIds, DatafileIds;

%     diary('allTbTe.log');
    
    for si = 1:length(stimTypesToDo)
        stim = stimTypesToDo(si);
        stimType = stimTypes{stim};
        disp([' *** STARTING ' upper(stimType) ' STIMULI ***']);
%         Dids = dbGetStimulusDids(stimType);        
        Dids = DidsNeedTbTeFixing; %{stim};
        
        progressBar('init=', length(Dids));
        for i = 1:length(Dids)
            progressBar(i);
            fprintf('(%s) Did = %d : ', outOf(i, length(Dids)), Dids(i));
            if (method == ORIG_METHOD)
                edbMakeTbTe(hnd, Dids(i));
            elseif (method == MY_METHOD)
                dbMakeTbTe(hnd, Dids(i));
            end
            fprintf('\n');
%             disp(['completed ' outOf(i, length(Dids)) ]);
        end
        
    end
    
    diary('off');
%     disp('Dids not done : ');
%     disp(DidsToComeBackTo);


end


