
% After pruning:

% -------------------------------------------
% 6 mseq groups (out of total of 24 in orig table)
mseqGids = [1025 1046 1047 1050 1051 1052, [1547]];

% N Values
% ----- ----------------
% 17 : db:noSpikes 
% 1 : db:noSyncs 
% 6 : ok 

% -------------------------------------------
% 266 noise groups (out of total of 358 in orig table)
noiseGids = [1 2 3 4 17 25 26 32 33 35 36 38 44 53 62 69 75 87 88 91 112 118 121 131 133 136 148 ...
 159 172 175 195 196 197 198 199 204 216 217 235 236 237 238 259 260 280 281 333 334 440 446 467 471 472 475 477 478 489 498 499 501 512 524 533 ...
 536 571 577 585 590 606 609 611 618 621 627 629 630 632 637 639 641 642 645 646 651 653 654 657 658 660 661 666 669 671 673 677 680 681 683 684 ...
 685 686 689 693 694 696 697 699 702 704 705 706 707 708 711 712 715 717 718 720 721 722 725 727 729 731 733 737 738 739 740 744 746 747 748 749 ...
 750 751 753 754 756 758 759 760 762 764 766 767 769 770 773 775 776 777 778 779 780 781 783 785 786 787 789 805 821 832 836 840 841 845 846 847 ...
 850 857 863 864 865 866 867 868 870 872 874 876 877 878 879 880 881 882 883 885 887 888 889 892 893 894 898 899 900 901 902 903 904 905 907 911 ...
 912 913 914 915 916 1003 1006 1011 1015 1016 1019 1030 1031 1038 1041 1045 1120 1133 1137 1139 1142 1165 1205 1211 1212 1222 1225 1228 1231 1235 1463 1467 1481 1594 1595 1711 ...
 1729 1739 1830 1879 1881 1883 1897 1943 1945 1947 1989 1991 1993 2030 2076 2078 2080 2111 2112 2122 2123 2124 4189];


%    N       Values
%  -----    ----------------
%     3  :  db:MultipleMissedFrames 
%     1  :  db:flawedSyncs 
%    92  :  db:noSpikes 
%   200  :  ok 
%    52  :  stim:interrupted 
%    10  :  warning 



% Grating 
gratingGids = [ ...
 5 6 7 8 9 10 11 12 13 14 15 16 18 19 20 21 22 23 24 27 28 29 30 31 34 37 39 ...
 40 41 42 43 45 47 48 49 50 51 52 54 55 56 57 58 59 60 ...
 61 63 64 65 66 67 68 70 71 72 73 74 76 77 78 79 80 81 ...
 82 83 84 85 86 89 90 92 93 94 95 96 97 98 99 100 101 102 ...
 103 104 105 106 107 108 109 110 111 113 114 115 116 117 119 120 122 123 ...
 124 125 126 127 128 129 130 132 135 137 138 139 140 141 142 143 144 145 ...
 146 147 149 150 151 152 153 154 155 156 157 158 160 161 162 163 164 165 ...
 166 167 168 169 170 171 173 174 176 177 178 179 180 181 182 183 184 185 ...
 186 187 188 189 190 191 192 193 194 200 201 202 203 205 206 207 208 209 ...
 210 211 212 213 214 215 218 219 220 221 222 223 224 225 226 227 228 229 ...
 230 231 232 233 234 239 240 241 242 243 244 245 246 247 248 249 250 251 ...
 252 253 254 255 256 257 258 261 262 263 264 265 266 267 268 269 270 271 ...
 272 273 274 275 276 277 278 279 282 283 284 285 286 287 288 289 290 291 ...
 292 293 294 295 296 297 298 299 300 301 302 303 304 305 306 307 308 309 ...
 310 311 312 313 314 315 316 317 318 319 320 321 322 323 324 325 326 327 ...
 328 329 330 331 332 335 336 337 338 339 340 341 342 421 422 423 424 425 ...
 426 427 428 429 430 431 432 433 434 435 436 437 438 439 441 442 443 444 ...
 445 447 448 449 450 451 452 453 454 455 457 458 459 460 461 462 463 464 ...
 465 466 468 469 473 474 476 479 480 481 482 483 484 485 486 487 488 490 ...
 491 492 493 494 495 496 497 502 503 504 505 506 507 508 509 510 511 513 ...
 514 515 516 517 518 519 520 521 522 523 525 526 527 528 529 530 531 532 ...
 534 535 537 538 539 540 541 542 543 544 545 546 547 548 549 550 551 552 ...
 553 554 555 556 557 558 559 560 561 562 563 564 565 566 567 568 569 570 ...
 572 573 574 575 576 578 579 580 581 582 583 584 586 587 588 589 595 596 ...
 597 598 599 600 601 602 605 607 608 610 612 613 614 615 616 617 619 620 ...
 623 624 625 626 628 631 633 634 640 644 647 648 649 650 652 659 662 663 ...
 664 665 667 668 670 672 674 675 676 678 679 682 687 695 698 700 701 703 ...
 709 710 713 714 716 719 723 724 726 728 730 732 734 735 736 741 742 743 ...
 745 752 755 757 761 763 765 768 771 772 774 782 784 788 790 803 807 808 ...
 810 813 814 817 818 820 822 826 829 830 831 833 834 835 837 838 839 842 ...
 843 844 848 849 851 852 853 854 855 856 858 859 860 861 862 869 871 873 ...
 875 884 886 895 906 908 917 918 919 920 921 922 923 924 925 926 927 928 ...
 929 931 932 933 1002 1004 1005 1007 1008 1009 1010 1012 1013 1014 1017 1018 1022 1023 ...
 1026 1027 1028 1029 1032 1033 1034 1035 1036 1037 1039 1040 1042 1044 1048 1101 1102 1103 ...
 1104 1105 1107 1108 1109 1110 1111 1112 1113 1114 1115 1116 1117 1118 1119 1121 1122 1123 ...
 1124 1125 1126 1127 1128 1129 1130 1131 1132 1134 1135 1136 1138 1140 1141 1143 1144 1145 ...
 1146 1147 1148 1149 1150 1151 1152 1153 1154 1155 1156 1157 1158 1159 1160 1162 1163 1164 ...
 1166 1167 1168 1169 1170 1171 1172 1173 1174 1175 1176 1192 1197 1200 1208 1209 1213 1216 ...
 1218 1219 1234 1238 1376 1377 1394 1398 1400 1401 1417 1422 1427 1428 1435 1438 1439 1440 ...
 1451 1453 1459 1460 1464 1468 1472 1473 1475 1476 1479 1483 1484 1486 1490 1492 1494 1498 ...
 1500 1503 1508 1511 1512 1513 1514 1516 1518 1519 1520 1523 1527 1528 1530 1533 1538 1539 ...
 1540 1556 1557 1559 1560 1583 1584 1585 1586 1589 1590 1597 1598 1599 1600 1601 1602 1607 ...
 1608 1610 1612 1623 1624 1626 1642 1656 1660 1678 1685 1695 1703 1705 1713 1715 1717 1721 ...
 1723 1733 1743 1747 1756 1768 1770 1772 1774 1780 1788 1802 1808 1818 1836 1838 1847 1859 ...
 1887 1889 1899 1901 1905 1917 1919 1923 1941 1949 1951 1953 1959 1961 1965 1967 1997 1998 ...
 2010 2014 2016 2032 2034 2036 2042 2046 2058 2082 2084 2087 2089 2090 2091 2092 2093 2094 ...
 2097 2098 2099 2102 2106 2107 2108 2109 2116 2119 2120 2127 2128 2129 2130 2133 2134 2135 ...
 2136 2140 2142 2143 2144 2146 2147 2151 2152 2154 2155 2156 2157 2160 2161 2163 2164 2167 ...
 2168 2171 2172 2173 2174 2175 2176 2177 2180 2181 2182 2183 2184 2185 2186 2187 2189 2190 ...
 2196 2198 2200 2202 2204 2206 2218 2230 2242 2252 2254 2264 2282 2294 2296 2306 2316 2318 ...
 2320 2322 2324 2326 2328 2330 2332 2342 2360 2378 2396 2398 2414 2424 2426 2436 2438 2450 ...
 2452 2466 2468 2472 2474 2488 2490 2492 2494 2508 2510 2512 2526 2530 2532 2536 2550 2552 ...
 2553 2558 2559 2564 2565 2568 2585 2587 2589 2591 2599 2603 2607 2611 2615 2627 2629 2637 ...
 2645 2647 2653 2663 2667 2669 2691 2693 2699 2713 2719 2727 2729 2739 2743 2757 2769 2779 ...
 2789 2795 2797 2807 2813 2819 2825 2831 2837 2843 2849 2863 2867 2877 2881 2893 2897 2905 ...
 2909 2915 2925 2939 2951 2955 2959 2963 2965 2973 2983 2985 2991 2997 3011 3017 3021 3023 ...
 3033 3037 3039 3045 3053 3055 3063 3069 3075 3081 3083 3089 3095 3103 3105 3121 3123 3127 ...
 3129 3135 3141 3143 3149 3161 3167 3169 3171 3183 3199 3205 3211 3227 3233 3243 3249 3255 ...
 3277 3287 3295 3301 3307 3309 3319 3337 3343 3349 3355 3369 3375 3381 3387 3393 3399 3407 ...
 3413 3421 3427 3451 3459 3465 3473 3485 3491 3497 3503 3507 3509 3521 3543 3561 3573 3585 ...
 3593 3599 3605 3607 3613 3619 3627 3629 3635 3641 3651 3657 3663 3677 3683 3689 3707 3717 ...
 3723 3733 3769 3773 3795 3799 3805 3813 3825 3829 3833 3839 3851 3853 3859 3865 3875 3881 ...
 3887 3893 3921 3927 3935 3947 3953 3959 3965 3971 3977 3985 3993 3995 4001 4007 4015 4021 ...
 4027 4035 4047 4049 4057 4061 4069 4081 4087 4093 4103 4111 4115 4125 4131 4137 4147 4153 ...
 4159 4187 4200 4201 4204 4205 4206 4207 4210 4211 4212 4213 4214 4215 4216 4217 4218 4219 ...
 4221 4222 4223 4224 4226 4228 4230 4232 4234 4236 4238 4240 4242 4244 4246 4252 4256 4270 ...
 4276 4282 4290 4310 4316 4350 4358 4368 4378 4386 4390 4392 4396 4402 4412 4422 4428 4448 ...
 4454 4456 4464 4472 4478 4484 4490 4496 4500 4504 4510 4512 4516 4524 4536 4540 4542 4550 ...
 4568 4572 4578 4584 4588 4592 4594 4610 4616 4618 4624 4648 4650 4652 4666 4674 4690 4694 ...
 4698 4704 4710 4720 4728 4734 4742 4754 4756 4760 4762 4764 4766 4770 4772 4776 4780 4782 ...
 4786 4790 4792 4800 4804 4812 4818 4822 4830 4834 4844 4858 4870 4876 4880 4886 4896 4902 ...
 4906 4908 4916 4920 4924 4926 4928 4932 4934 4940 4944 4948 4958 4960 4974 4978 4980 4986 ...
 4988 4992 4994 4996 5000 5002 5008 5010 5014 5016 5038 5046 5050 5052 5078 5082 5092 5096 ...
 5104 5106 5110 5116 5118 5124 5128 5130 5134 5136 5142 5144 5148 5154 5156 5160 5162 5166 ...
 5168 5172 5174 5178 5180 5182 5186 5188 5210 5212 5216 5218 5222 5224 5228 5232 5234 5242 ...
 5244 5252 5254 5258 5260 5268 5272 5274 5278 5280 5284 5286 5290 5292 5296 5298 5302 5306 ...
 5334 5336 ];

% N Values
% ----- ----------------
% 6 : db:badSyncs 
% 3 : db:flawedSyncs 
% 1299 : db:noSpikes 
% 5 : db:noSyncs 
% 16 : db:syncNframesMismatch 
% 976 : ok 
% 293 : stim:interrupted 
% 80 : warning 


% Movie cellGroups
movieGids = [...
 1193 1194 1195 1196 1198 1199 1201 1202 1203 1204 1206 1207 1210 1214 1215 1217 1220 1221 ...
 1223 1224 1226 1227 1229 1230 1232 1233 1236 1237 1396 1397 1432 1433 1436 1441 1442 1443 ...
 1444 1462 1465 1469 1474 1477 1480 1485 1491 1493 1495 1499 1506 1510 1515 1517 1521 1524 ...
 1531 1534 1535 1536 1628 1644 1646 1648 1652 1662 1664 1666 1668 1697 1699 1701 1707 1709 ...
 1725 1727 1735 1737 1741 1749 1751 1753 1758 1790 1792 1794 1810 1812 1814 1816 1820 1822 ...
 1824 1826 1828 1832 1834 1853 1855 1857 1861 1863 1867 1875 1877 1885 1891 1893 1895 1925 ...
 1927 1929 1931 1933 1935 1937 1939 1969 1973 1975 1977 1979 1981 1983 1985 1987 2000 2008 ...
 2018 2020 2022 2024 2026 2028 2060 2062 2064 2066 2068 2070 2072 2074 2110 2113 2114 2118 ...
 2121 2125 2126 2131 2132 2138 2139 2145 2153 2158 2159 2162 2165 2166 2170 2178 2179 2191 ...
 2214 2216 2220 2224 2226 2228 2234 2236 2238 2240 2244 2246 2248 2250 2256 2258 2260 2262 ...
 2268 2270 2272 2274 2276 2278 2284 2286 2288 2290 2292 2298 2300 2302 2304 2308 2310 2312 ...
 2314 2334 2336 2338 2340 2344 2346 2348 2350 2352 2354 2356 2358 2362 2368 2370 2372 2374 ...
 2376 2416 2418 2420 2422 2428 2430 2432 2434 2440 2442 2444 2446 2448 2454 2456 2458 2460 ...
 2462 2464 2476 2478 2480 2482 2484 2486 2496 2498 2500 2502 2504 2506 2514 2516 2518 2520 ...
 2522 2524 2538 2540 2542 2544 2546 2548 2595 2597 2601 2617 2619 2621 2625 2631 2635 2649 ...
 2655 2657 2661 2665 2671 2673 2675 2677 2695 2697 2701 2715 2717 2721 2723 2731 2733 2735 ...
 2737 2741 2745 2747 2749 2759 2761 2771 2773 2775 2777 2785 2787 2791 2793 2799 2801 2803 ...
 2805 2815 2817 2821 2827 2829 2833 2839 2841 2845 2851 2853 2855 2861 2871 2875 2883 2885 ...
 2887 2891 2903 2911 2913 2917 2919 2921 2923 2953 2957 2975 2977 2979 2981 2987 2989 2993 ...
 2995 3013 3019 3025 3027 3035 3041 3043 3047 3049 3051 3057 3059 3061 3065 3067 3071 3073 ...
 3077 3079 3085 3087 3091 3093 3097 3099 3107 3115 3117 3125 3137 3139 3145 3147 3151 3157 ...
 3159 3163 3165 3173 3185 3201 3203 3207 3209 3213 3215 3217 3219 3225 3229 3231 3235 3237 ...
 3239 3245 3247 3251 3253 3257 3259 3261 3263 3279 3281 3289 3293 3297 3299 3311 3315 3317 ...
 3339 3341 3345 3347 3351 3353 3357 3359 3361 3363 3365 3367 3371 3373 3377 3379 3383 3385 ...
 3389 3391 3395 3397 3401 3403 3409 3411 3415 3429 3431 3433 3435 3437 3439 3441 3453 3455 ...
 3461 3463 3467 3469 3471 3487 3489 3493 3495 3499 3501 3505 3511 3513 3595 3597 3601 3603 ...
 3609 3611 3615 3617 3621 3623 3625 3631 3633 3637 3639 3643 3645 3647 3649 3653 3655 3659 ...
 3661 3667 3669 3671 3673 3679 3681 3685 3687 3691 3693 3695 3709 3711 3719 3721 3725 3727 ...
 3729 3731 3775 3801 3803 3807 3809 3815 3817 3819 3821 3827 3831 3835 3837 3841 3843 3845 ...
 3847 3855 3857 3861 3863 3869 3873 3877 3879 3883 3885 3889 3891 3895 3923 3925 3931 3933 ...
 3937 3939 3941 3945 3949 3951 3955 3957 3961 3963 3967 3969 3973 3975 3979 3981 3983 3987 ...
 3991 3997 3999 4003 4005 4009 4011 4013 4017 4019 4023 4025 4029 4031 4033 4037 4041 4043 ...
 4045 4051 4053 4055 4065 4067 4071 4073 4075 4083 4085 4089 4091 4095 4097 4099 4101 4105 ...
 4107 4109 4113 4117 4119 4121 4127 4129 4133 4135 4139 4141 4143 4145 4149 4151 4155 4157 ...
 4161 4163 4165 4191 4254 4274 4278 4284 4288 4292 4294 4296 4300 4302 4304 4306 4308 4312 ...
 4314 4318 4320 4322 4324 4326 4328 4330 4344 4356 4362 4388 4398 4400 4404 4406 4408 4410 ...
 4424 4426 4430 4434 4436 4438 4440 4442 4444 4446 4450 4452 4460 4462 4466 4468 4470 4474 ...
 4476 4480 4482 4486 4488 4492 4494 4498 4502 4506 4508 4514 4518 4520 4522 4538 4544 4546 ...
 4548 4566 4570 4580 4590 4596 4614 4620 4622 4626 4628 4634 4636 4638 4640 4642 4646 4670 ...
 4672 4676 4680 4682 4684 4686 4688 4692 4696 4702 4706 4708 4712 4714 4716 4718 4722 4726 ...
 4730 4732 4738 4744 4746 4758 4768 4774 4778 4784 4788 4794 4796 4798 4806 4808 4810 4824 ...
 4826 4828 4832 4836 4838 4840 4842 4846 4848 4864 4866 4868 4872 4874 4878 4882 4884 4888 ...
 4890 4892 4898 4904 4910 4918 4930 4936 4946 4950 4976 4990 4998 5012 5048 5084 5086 5100 ...
 5112 5114 5126 5132 5146 5158 5164 5170 5176 5184 5214 5226 5236 5248 5256 5262 5270 5276 ...
 5282 5288 5294 5300 5308 5312 ];
% N Values
% ----- ----------------
% 17 : db:SinglePhaseMovie 
% 38 : db:flawedSyncs 
% 1177 : db:noSpikes 
% 2 : db:noSyncs 
% 14 : db:syncNframesMismatch 
% 706 : ok 
% 2 : stim:interrupted 
% 19 : warning 

% mseq:  Unique Dids: 7;  unique Gids: 7;
% noise:  Unique Dids: 259;  unique Gids: 266;
% grating: Unique Dids: 1320  unique Gids: 1379
% movie:  Unique Dids: 798  unique Gids: 798

total D : 2384
total G : 2450




