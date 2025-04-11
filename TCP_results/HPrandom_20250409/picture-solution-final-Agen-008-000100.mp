prologues := 3;
outputformat   := "svg";
outputtemplate := "%j.%o";
input mpcolornames;
beginfig(8);
u := 20cm;
path pols[];
pair p[];
R := .2063373359u;
D := 2*R;

pols[1] = 
(       -0.3695000000u,       -0.1938000000u)--
(       -0.4071000000u,       -0.1956000000u)--
(       -0.4383000000u,       -0.2155000000u)--
(       -0.4460000000u,       -0.2505000000u)--
(       -0.4110000000u,       -0.2998000000u)--
(       -0.3668000000u,       -0.3348000000u)--
(       -0.3695000000u,       -0.1938000000u)--
cycle;
pols[2] = 
(       -0.1408000000u,       -0.2635000000u)--
(       -0.1603000000u,       -0.2842000000u)--
(       -0.1702000000u,       -0.3218000000u)--
(       -0.1153000000u,       -0.4217000000u)--
(       -0.0670000000u,       -0.4295000000u)--
(       -0.1408000000u,       -0.2635000000u)--
cycle;
pols[3] = 
(        0.1748000000u,       -0.4230000000u)--
(        0.1990000000u,       -0.3919000000u)--
(        0.2067000000u,       -0.3576000000u)--
(        0.1838000000u,       -0.3687000000u)--
(        0.1748000000u,       -0.4230000000u)--
cycle;
pols[4] = 
(        0.0637000000u,       -0.4269000000u)--
(        0.0896000000u,       -0.4606000000u)--
(        0.1311000000u,       -0.4814000000u)--
(        0.1571000000u,       -0.4736000000u)--
(        0.1705000000u,       -0.4489000000u)--
(        0.1748000000u,       -0.4230000000u)--
(        0.1838000000u,       -0.3687000000u)--
(        0.0637000000u,       -0.4269000000u)--
cycle;
pols[5] = 
(        0.2414000000u,       -0.3080000000u)--
(        0.2790000000u,       -0.3365000000u)--
(        0.3010000000u,       -0.3365000000u)--
(        0.3259000000u,       -0.3094000000u)--
(        0.2414000000u,       -0.3080000000u)--
cycle;
pols[6] = 
(        0.3361000000u,       -0.2985000000u)--
(        0.3879000000u,       -0.2777000000u)--
(        0.3966000000u,       -0.2628000000u)--
(        0.3700000000u,       -0.2622000000u)--
(        0.3361000000u,       -0.2985000000u)--
cycle;
pols[7] = 
(        0.2180000000u,       -0.3076000000u)--
(        0.2414000000u,       -0.3080000000u)--
(        0.3259000000u,       -0.3094000000u)--
(        0.3361000000u,       -0.2985000000u)--
(        0.3700000000u,       -0.2622000000u)--
(        0.2289000000u,       -0.2591000000u)--
(        0.2180000000u,       -0.3076000000u)--
cycle;
pols[8] = 
(       -0.0297000000u,       -0.4113000000u)--
(        0.0352000000u,       -0.4407000000u)--
(        0.0637000000u,       -0.4269000000u)--
(        0.2067000000u,       -0.3576000000u)--
(        0.2180000000u,       -0.3076000000u)--
(        0.2289000000u,       -0.2591000000u)--
(        0.1797000000u,       -0.2580000000u)--
(       -0.0297000000u,       -0.4113000000u)--
cycle;
pols[9] = 
(       -0.1464000000u,       -0.2509000000u)--
(       -0.1408000000u,       -0.2635000000u)--
(       -0.0670000000u,       -0.4295000000u)--
(       -0.0569000000u,       -0.4312000000u)--
(       -0.0297000000u,       -0.4113000000u)--
(        0.1797000000u,       -0.2580000000u)--
(       -0.1464000000u,       -0.2509000000u)--
cycle;
pols[10] = 
(        0.4243000000u,       -0.2159000000u)--
(        0.4697000000u,       -0.2375000000u)--
(        0.4969000000u,       -0.2367000000u)--
(        0.5228000000u,       -0.2220000000u)--
(        0.5358000000u,       -0.1830000000u)--
(        0.5250000000u,       -0.1454000000u)--
(        0.4831000000u,       -0.1163000000u)--
(        0.4243000000u,       -0.2159000000u)--
cycle;
pols[11] = 
(        0.4035000000u,       -0.0961000000u)--
(        0.3620000000u,       -0.0339000000u)--
(        0.3316000000u,       -0.0066000000u)--
(        0.2824000000u,       -0.0299000000u)--
(        0.4035000000u,       -0.0961000000u)--
cycle;
pols[12] = 
(        0.4320000000u,       -0.1117000000u)--
(        0.4035000000u,       -0.0961000000u)--
(        0.2824000000u,       -0.0299000000u)--
(        0.1033000000u,       -0.1147000000u)--
(        0.4320000000u,       -0.1117000000u)--
cycle;
pols[13] = 
(       -0.1827000000u,       -0.2501000000u)--
(       -0.1464000000u,       -0.2509000000u)--
(        0.3966000000u,       -0.2628000000u)--
(        0.4243000000u,       -0.2159000000u)--
(        0.4831000000u,       -0.1163000000u)--
(        0.4761000000u,       -0.1113000000u)--
(        0.4320000000u,       -0.1117000000u)--
(        0.1033000000u,       -0.1147000000u)--
(       -0.1827000000u,       -0.2501000000u)--
cycle;
pols[14] = 
(        0.2314000000u,        0.1736000000u)--
(        0.2466000000u,        0.2592000000u)--
(        0.2453000000u,        0.2645000000u)--
(        0.2110000000u,        0.2280000000u)--
(        0.2314000000u,        0.1736000000u)--
cycle;
pols[15] = 
(        0.2440000000u,        0.1399000000u)--
(        0.2314000000u,        0.1736000000u)--
(        0.2110000000u,        0.2280000000u)--
(        0.1462000000u,        0.1592000000u)--
(        0.1686000000u,        0.1391000000u)--
(        0.2440000000u,        0.1399000000u)--
cycle;
pols[16] = 
(        0.3041000000u,        0.0180000000u)--
(        0.3326000000u,        0.0686000000u)--
(        0.3235000000u,        0.1088000000u)--
(        0.2842000000u,        0.1403000000u)--
(        0.2440000000u,        0.1399000000u)--
(        0.1686000000u,        0.1391000000u)--
(        0.3041000000u,        0.0180000000u)--
cycle;
pols[17] = 
(       -0.2839000000u,       -0.2980000000u)--
(       -0.1827000000u,       -0.2501000000u)--
(        0.3316000000u,       -0.0066000000u)--
(        0.3041000000u,        0.0180000000u)--
(        0.1462000000u,        0.1592000000u)--
(       -0.2839000000u,       -0.2980000000u)--
cycle;
pols[18] = 
(        0.2276000000u,        0.3435000000u)--
(        0.2959000000u,        0.3678000000u)--
(        0.3348000000u,        0.3941000000u)--
(        0.3560000000u,        0.4395000000u)--
(        0.3508000000u,        0.4849000000u)--
(        0.3387000000u,        0.5096000000u)--
(        0.3023000000u,        0.5385000000u)--
(        0.2427000000u,        0.5459000000u)--
(        0.2064000000u,        0.5368000000u)--
(        0.1877000000u,        0.5212000000u)--
(        0.2276000000u,        0.3435000000u)--
cycle;
pols[19] = 
(       -0.0751000000u,        0.5083000000u)--
(       -0.0984000000u,        0.5148000000u)--
(       -0.1220000000u,        0.4966000000u)--
(       -0.1126000000u,        0.4952000000u)--
(       -0.0751000000u,        0.5083000000u)--
cycle;
pols[20] = 
(        0.0507000000u,        0.4724000000u)--
(       -0.0466000000u,        0.5182000000u)--
(       -0.0751000000u,        0.5083000000u)--
(       -0.1126000000u,        0.4952000000u)--
(        0.0507000000u,        0.4724000000u)--
cycle;
pols[21] = 
(       -0.2061000000u,        0.4642000000u)--
(       -0.2826000000u,        0.4590000000u)--
(       -0.4015000000u,        0.4266000000u)--
(       -0.4538000000u,        0.3786000000u)--
(       -0.4544000000u,        0.3761000000u)--
(       -0.2061000000u,        0.4642000000u)--
cycle;
pols[22] = 
(       -0.1283000000u,        0.4918000000u)--
(       -0.2061000000u,        0.4642000000u)--
(       -0.4544000000u,        0.3761000000u)--
(       -0.4629000000u,        0.3384000000u)--
(       -0.4581000000u,        0.3085000000u)--
(       -0.4216000000u,        0.2683000000u)--
(       -0.1283000000u,        0.4918000000u)--
cycle;
pols[23] = 
(       -0.3371000000u,        0.1140000000u)--
(       -0.3760000000u,        0.1205000000u)--
(       -0.4305000000u,        0.1084000000u)--
(       -0.4586000000u,        0.0880000000u)--
(       -0.4642000000u,        0.0543000000u)--
(       -0.4279000000u,        0.0150000000u)--
(       -0.4178000000u,        0.0132000000u)--
(       -0.3410000000u,        0.0777000000u)--
(       -0.3371000000u,        0.1140000000u)--
cycle;
pols[24] = 
(       -0.3111000000u,        0.2346000000u)--
(       -0.3293000000u,        0.1866000000u)--
(       -0.3371000000u,        0.1140000000u)--
(       -0.3410000000u,        0.0777000000u)--
(       -0.2665000000u,        0.1403000000u)--
(       -0.3111000000u,        0.2346000000u)--
cycle;
pols[25] = 
(       -0.3215000000u,        0.2566000000u)--
(       -0.3111000000u,        0.2346000000u)--
(       -0.2665000000u,        0.1403000000u)--
(       -0.0899000000u,        0.2884000000u)--
(       -0.3215000000u,        0.2566000000u)--
cycle;
pols[26] = 
(        0.1182000000u,        0.4629000000u)--
(        0.0507000000u,        0.4724000000u)--
(       -0.1220000000u,        0.4966000000u)--
(       -0.1283000000u,        0.4918000000u)--
(       -0.4216000000u,        0.2683000000u)--
(       -0.4188000000u,        0.2653000000u)--
(       -0.3682000000u,        0.2502000000u)--
(       -0.3215000000u,        0.2566000000u)--
(       -0.0899000000u,        0.2884000000u)--
(        0.1182000000u,        0.4629000000u)--
cycle;
pols[27] = 
(       -0.0504000000u,       -0.0497000000u)--
(        0.2453000000u,        0.2645000000u)--
(        0.2276000000u,        0.3435000000u)--
(        0.1877000000u,        0.5212000000u)--
(        0.1182000000u,        0.4629000000u)--
(       -0.4178000000u,        0.0132000000u)--
(       -0.3098000000u,       -0.0053000000u)--
(       -0.0504000000u,       -0.0497000000u)--
cycle;
pols[28] = 
(       -0.1116000000u,       -0.1147000000u)--
(       -0.0504000000u,       -0.0497000000u)--
(       -0.3098000000u,       -0.0053000000u)--
(       -0.2839000000u,       -0.0196000000u)--
(       -0.1116000000u,       -0.1147000000u)--
cycle;
pols[29] = 
(       -0.3699000000u,       -0.1722000000u)--
(       -0.3695000000u,       -0.1938000000u)--
(       -0.3668000000u,       -0.3348000000u)--
(       -0.3630000000u,       -0.3378000000u)--
(       -0.3189000000u,       -0.3352000000u)--
(       -0.2839000000u,       -0.2980000000u)--
(       -0.1116000000u,       -0.1147000000u)--
(       -0.2839000000u,       -0.0196000000u)--
(       -0.3626000000u,       -0.1428000000u)--
(       -0.3699000000u,       -0.1722000000u)--
cycle;
fill pols[1] withcolor 0.9Dandelion;
fill pols[2] withcolor 0.9Dandelion;
fill pols[3] withcolor 0.9Dandelion;
fill pols[4] withcolor 0.9Dandelion;
fill pols[5] withcolor 0.9Dandelion;
fill pols[6] withcolor 0.9Dandelion;
fill pols[7] withcolor 0.9Dandelion;
fill pols[8] withcolor 0.9Dandelion;
fill pols[9] withcolor 0.9Dandelion;
fill pols[10] withcolor 0.9Dandelion;
fill pols[11] withcolor 0.9Dandelion;
fill pols[12] withcolor 0.9Dandelion;
fill pols[13] withcolor 0.9Dandelion;
fill pols[14] withcolor 0.9Dandelion;
fill pols[15] withcolor 0.9Dandelion;
fill pols[16] withcolor 0.9Dandelion;
fill pols[17] withcolor 0.9Dandelion;
fill pols[18] withcolor 0.9Dandelion;
fill pols[19] withcolor 0.9Dandelion;
fill pols[20] withcolor 0.9Dandelion;
fill pols[21] withcolor 0.9Dandelion;
fill pols[22] withcolor 0.9Dandelion;
fill pols[23] withcolor 0.9Dandelion;
fill pols[24] withcolor 0.9Dandelion;
fill pols[25] withcolor 0.9Dandelion;
fill pols[26] withcolor 0.9Dandelion;
fill pols[27] withcolor 0.9Dandelion;
fill pols[28] withcolor 0.9Dandelion;
fill pols[29] withcolor 0.9Dandelion;
draw (       -0.3668000000u,       -0.3348000000u)--(       -0.3695000000u,       -0.1938000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.0670000000u,       -0.4295000000u)--(       -0.1408000000u,       -0.2635000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.2067000000u,       -0.3576000000u)--(        0.1838000000u,       -0.3687000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1838000000u,       -0.3687000000u)--(        0.1748000000u,       -0.4230000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1748000000u,       -0.4230000000u)--(        0.1838000000u,       -0.3687000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1838000000u,       -0.3687000000u)--(        0.0637000000u,       -0.4269000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.3259000000u,       -0.3094000000u)--(        0.2414000000u,       -0.3080000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.3966000000u,       -0.2628000000u)--(        0.3700000000u,       -0.2622000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.3700000000u,       -0.2622000000u)--(        0.3361000000u,       -0.2985000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.2414000000u,       -0.3080000000u)--(        0.3259000000u,       -0.3094000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.3361000000u,       -0.2985000000u)--(        0.3700000000u,       -0.2622000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.3700000000u,       -0.2622000000u)--(        0.2289000000u,       -0.2591000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.2289000000u,       -0.2591000000u)--(        0.2180000000u,       -0.3076000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.0637000000u,       -0.4269000000u)--(        0.2067000000u,       -0.3576000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.2180000000u,       -0.3076000000u)--(        0.2289000000u,       -0.2591000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.2289000000u,       -0.2591000000u)--(        0.1797000000u,       -0.2580000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1797000000u,       -0.2580000000u)--(       -0.0297000000u,       -0.4113000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.1408000000u,       -0.2635000000u)--(       -0.0670000000u,       -0.4295000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.0297000000u,       -0.4113000000u)--(        0.1797000000u,       -0.2580000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1797000000u,       -0.2580000000u)--(       -0.1464000000u,       -0.2509000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.4831000000u,       -0.1163000000u)--(        0.4243000000u,       -0.2159000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.3316000000u,       -0.0066000000u)--(        0.2824000000u,       -0.0299000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.2824000000u,       -0.0299000000u)--(        0.4035000000u,       -0.0961000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.4035000000u,       -0.0961000000u)--(        0.2824000000u,       -0.0299000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.2824000000u,       -0.0299000000u)--(        0.1033000000u,       -0.1147000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1033000000u,       -0.1147000000u)--(        0.4320000000u,       -0.1117000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.1464000000u,       -0.2509000000u)--(        0.3966000000u,       -0.2628000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.4243000000u,       -0.2159000000u)--(        0.4831000000u,       -0.1163000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.4320000000u,       -0.1117000000u)--(        0.1033000000u,       -0.1147000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1033000000u,       -0.1147000000u)--(       -0.1827000000u,       -0.2501000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.2453000000u,        0.2645000000u)--(        0.2110000000u,        0.2280000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.2110000000u,        0.2280000000u)--(        0.2314000000u,        0.1736000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.2314000000u,        0.1736000000u)--(        0.2110000000u,        0.2280000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.2110000000u,        0.2280000000u)--(        0.1462000000u,        0.1592000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1462000000u,        0.1592000000u)--(        0.1686000000u,        0.1391000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1686000000u,        0.1391000000u)--(        0.2440000000u,        0.1399000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.2440000000u,        0.1399000000u)--(        0.1686000000u,        0.1391000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1686000000u,        0.1391000000u)--(        0.3041000000u,        0.0180000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.1827000000u,       -0.2501000000u)--(        0.3316000000u,       -0.0066000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.3041000000u,        0.0180000000u)--(        0.1462000000u,        0.1592000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1462000000u,        0.1592000000u)--(       -0.2839000000u,       -0.2980000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1877000000u,        0.5212000000u)--(        0.2276000000u,        0.3435000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.1220000000u,        0.4966000000u)--(       -0.1126000000u,        0.4952000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.1126000000u,        0.4952000000u)--(       -0.0751000000u,        0.5083000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.0751000000u,        0.5083000000u)--(       -0.1126000000u,        0.4952000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.1126000000u,        0.4952000000u)--(        0.0507000000u,        0.4724000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.4544000000u,        0.3761000000u)--(       -0.2061000000u,        0.4642000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.2061000000u,        0.4642000000u)--(       -0.4544000000u,        0.3761000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.4216000000u,        0.2683000000u)--(       -0.1283000000u,        0.4918000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.4178000000u,        0.0132000000u)--(       -0.3410000000u,        0.0777000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.3410000000u,        0.0777000000u)--(       -0.3371000000u,        0.1140000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.3371000000u,        0.1140000000u)--(       -0.3410000000u,        0.0777000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.3410000000u,        0.0777000000u)--(       -0.2665000000u,        0.1403000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.2665000000u,        0.1403000000u)--(       -0.3111000000u,        0.2346000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.3111000000u,        0.2346000000u)--(       -0.2665000000u,        0.1403000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.2665000000u,        0.1403000000u)--(       -0.0899000000u,        0.2884000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.0899000000u,        0.2884000000u)--(       -0.3215000000u,        0.2566000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.0507000000u,        0.4724000000u)--(       -0.1220000000u,        0.4966000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.1283000000u,        0.4918000000u)--(       -0.4216000000u,        0.2683000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.3215000000u,        0.2566000000u)--(       -0.0899000000u,        0.2884000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.0899000000u,        0.2884000000u)--(        0.1182000000u,        0.4629000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.0504000000u,       -0.0497000000u)--(        0.2453000000u,        0.2645000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.2276000000u,        0.3435000000u)--(        0.1877000000u,        0.5212000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1182000000u,        0.4629000000u)--(       -0.4178000000u,        0.0132000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.3098000000u,       -0.0053000000u)--(       -0.0504000000u,       -0.0497000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.1116000000u,       -0.1147000000u)--(       -0.0504000000u,       -0.0497000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.0504000000u,       -0.0497000000u)--(       -0.3098000000u,       -0.0053000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.2839000000u,       -0.0196000000u)--(       -0.1116000000u,       -0.1147000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.3695000000u,       -0.1938000000u)--(       -0.3668000000u,       -0.3348000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.2839000000u,       -0.2980000000u)--(       -0.1116000000u,       -0.1147000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.1116000000u,       -0.1147000000u)--(       -0.2839000000u,       -0.0196000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.3695000000u,       -0.1938000000u)--(       -0.4071000000u,       -0.1956000000u) withpen pencircle scaled 1.0;
draw (       -0.4071000000u,       -0.1956000000u)--(       -0.4383000000u,       -0.2155000000u) withpen pencircle scaled 1.0;
draw (       -0.4383000000u,       -0.2155000000u)--(       -0.4460000000u,       -0.2505000000u) withpen pencircle scaled 1.0;
draw (       -0.4460000000u,       -0.2505000000u)--(       -0.4110000000u,       -0.2998000000u) withpen pencircle scaled 1.0;
draw (       -0.4110000000u,       -0.2998000000u)--(       -0.3668000000u,       -0.3348000000u) withpen pencircle scaled 1.0;
draw (       -0.1408000000u,       -0.2635000000u)--(       -0.1603000000u,       -0.2842000000u) withpen pencircle scaled 1.0;
draw (       -0.1603000000u,       -0.2842000000u)--(       -0.1702000000u,       -0.3218000000u) withpen pencircle scaled 1.0;
draw (       -0.1702000000u,       -0.3218000000u)--(       -0.1153000000u,       -0.4217000000u) withpen pencircle scaled 1.0;
draw (       -0.1153000000u,       -0.4217000000u)--(       -0.0670000000u,       -0.4295000000u) withpen pencircle scaled 1.0;
draw (        0.1748000000u,       -0.4230000000u)--(        0.1990000000u,       -0.3919000000u) withpen pencircle scaled 1.0;
draw (        0.1990000000u,       -0.3919000000u)--(        0.2067000000u,       -0.3576000000u) withpen pencircle scaled 1.0;
draw (        0.0637000000u,       -0.4269000000u)--(        0.0896000000u,       -0.4606000000u) withpen pencircle scaled 1.0;
draw (        0.0896000000u,       -0.4606000000u)--(        0.1311000000u,       -0.4814000000u) withpen pencircle scaled 1.0;
draw (        0.1311000000u,       -0.4814000000u)--(        0.1571000000u,       -0.4736000000u) withpen pencircle scaled 1.0;
draw (        0.1571000000u,       -0.4736000000u)--(        0.1705000000u,       -0.4489000000u) withpen pencircle scaled 1.0;
draw (        0.1705000000u,       -0.4489000000u)--(        0.1748000000u,       -0.4230000000u) withpen pencircle scaled 1.0;
draw (        0.2414000000u,       -0.3080000000u)--(        0.2790000000u,       -0.3365000000u) withpen pencircle scaled 1.0;
draw (        0.2790000000u,       -0.3365000000u)--(        0.3010000000u,       -0.3365000000u) withpen pencircle scaled 1.0;
draw (        0.3010000000u,       -0.3365000000u)--(        0.3259000000u,       -0.3094000000u) withpen pencircle scaled 1.0;
draw (        0.3361000000u,       -0.2985000000u)--(        0.3879000000u,       -0.2777000000u) withpen pencircle scaled 1.0;
draw (        0.3879000000u,       -0.2777000000u)--(        0.3966000000u,       -0.2628000000u) withpen pencircle scaled 1.0;
draw (        0.2180000000u,       -0.3076000000u)--(        0.2414000000u,       -0.3080000000u) withpen pencircle scaled 1.0;
draw (        0.3259000000u,       -0.3094000000u)--(        0.3361000000u,       -0.2985000000u) withpen pencircle scaled 1.0;
draw (       -0.0297000000u,       -0.4113000000u)--(        0.0352000000u,       -0.4407000000u) withpen pencircle scaled 1.0;
draw (        0.0352000000u,       -0.4407000000u)--(        0.0637000000u,       -0.4269000000u) withpen pencircle scaled 1.0;
draw (        0.2067000000u,       -0.3576000000u)--(        0.2180000000u,       -0.3076000000u) withpen pencircle scaled 1.0;
draw (       -0.1464000000u,       -0.2509000000u)--(       -0.1408000000u,       -0.2635000000u) withpen pencircle scaled 1.0;
draw (       -0.0670000000u,       -0.4295000000u)--(       -0.0569000000u,       -0.4312000000u) withpen pencircle scaled 1.0;
draw (       -0.0569000000u,       -0.4312000000u)--(       -0.0297000000u,       -0.4113000000u) withpen pencircle scaled 1.0;
draw (        0.4243000000u,       -0.2159000000u)--(        0.4697000000u,       -0.2375000000u) withpen pencircle scaled 1.0;
draw (        0.4697000000u,       -0.2375000000u)--(        0.4969000000u,       -0.2367000000u) withpen pencircle scaled 1.0;
draw (        0.4969000000u,       -0.2367000000u)--(        0.5228000000u,       -0.2220000000u) withpen pencircle scaled 1.0;
draw (        0.5228000000u,       -0.2220000000u)--(        0.5358000000u,       -0.1830000000u) withpen pencircle scaled 1.0;
draw (        0.5358000000u,       -0.1830000000u)--(        0.5250000000u,       -0.1454000000u) withpen pencircle scaled 1.0;
draw (        0.5250000000u,       -0.1454000000u)--(        0.4831000000u,       -0.1163000000u) withpen pencircle scaled 1.0;
draw (        0.4035000000u,       -0.0961000000u)--(        0.3620000000u,       -0.0339000000u) withpen pencircle scaled 1.0;
draw (        0.3620000000u,       -0.0339000000u)--(        0.3316000000u,       -0.0066000000u) withpen pencircle scaled 1.0;
draw (        0.4320000000u,       -0.1117000000u)--(        0.4035000000u,       -0.0961000000u) withpen pencircle scaled 1.0;
draw (       -0.1827000000u,       -0.2501000000u)--(       -0.1464000000u,       -0.2509000000u) withpen pencircle scaled 1.0;
draw (        0.3966000000u,       -0.2628000000u)--(        0.4243000000u,       -0.2159000000u) withpen pencircle scaled 1.0;
draw (        0.4831000000u,       -0.1163000000u)--(        0.4761000000u,       -0.1113000000u) withpen pencircle scaled 1.0;
draw (        0.4761000000u,       -0.1113000000u)--(        0.4320000000u,       -0.1117000000u) withpen pencircle scaled 1.0;
draw (        0.2314000000u,        0.1736000000u)--(        0.2466000000u,        0.2592000000u) withpen pencircle scaled 1.0;
draw (        0.2466000000u,        0.2592000000u)--(        0.2453000000u,        0.2645000000u) withpen pencircle scaled 1.0;
draw (        0.2440000000u,        0.1399000000u)--(        0.2314000000u,        0.1736000000u) withpen pencircle scaled 1.0;
draw (        0.3041000000u,        0.0180000000u)--(        0.3326000000u,        0.0686000000u) withpen pencircle scaled 1.0;
draw (        0.3326000000u,        0.0686000000u)--(        0.3235000000u,        0.1088000000u) withpen pencircle scaled 1.0;
draw (        0.3235000000u,        0.1088000000u)--(        0.2842000000u,        0.1403000000u) withpen pencircle scaled 1.0;
draw (        0.2842000000u,        0.1403000000u)--(        0.2440000000u,        0.1399000000u) withpen pencircle scaled 1.0;
draw (       -0.2839000000u,       -0.2980000000u)--(       -0.1827000000u,       -0.2501000000u) withpen pencircle scaled 1.0;
draw (        0.3316000000u,       -0.0066000000u)--(        0.3041000000u,        0.0180000000u) withpen pencircle scaled 1.0;
draw (        0.2276000000u,        0.3435000000u)--(        0.2959000000u,        0.3678000000u) withpen pencircle scaled 1.0;
draw (        0.2959000000u,        0.3678000000u)--(        0.3348000000u,        0.3941000000u) withpen pencircle scaled 1.0;
draw (        0.3348000000u,        0.3941000000u)--(        0.3560000000u,        0.4395000000u) withpen pencircle scaled 1.0;
draw (        0.3560000000u,        0.4395000000u)--(        0.3508000000u,        0.4849000000u) withpen pencircle scaled 1.0;
draw (        0.3508000000u,        0.4849000000u)--(        0.3387000000u,        0.5096000000u) withpen pencircle scaled 1.0;
draw (        0.3387000000u,        0.5096000000u)--(        0.3023000000u,        0.5385000000u) withpen pencircle scaled 1.0;
draw (        0.3023000000u,        0.5385000000u)--(        0.2427000000u,        0.5459000000u) withpen pencircle scaled 1.0;
draw (        0.2427000000u,        0.5459000000u)--(        0.2064000000u,        0.5368000000u) withpen pencircle scaled 1.0;
draw (        0.2064000000u,        0.5368000000u)--(        0.1877000000u,        0.5212000000u) withpen pencircle scaled 1.0;
draw (       -0.0751000000u,        0.5083000000u)--(       -0.0984000000u,        0.5148000000u) withpen pencircle scaled 1.0;
draw (       -0.0984000000u,        0.5148000000u)--(       -0.1220000000u,        0.4966000000u) withpen pencircle scaled 1.0;
draw (        0.0507000000u,        0.4724000000u)--(       -0.0466000000u,        0.5182000000u) withpen pencircle scaled 1.0;
draw (       -0.0466000000u,        0.5182000000u)--(       -0.0751000000u,        0.5083000000u) withpen pencircle scaled 1.0;
draw (       -0.2061000000u,        0.4642000000u)--(       -0.2826000000u,        0.4590000000u) withpen pencircle scaled 1.0;
draw (       -0.2826000000u,        0.4590000000u)--(       -0.4015000000u,        0.4266000000u) withpen pencircle scaled 1.0;
draw (       -0.4015000000u,        0.4266000000u)--(       -0.4538000000u,        0.3786000000u) withpen pencircle scaled 1.0;
draw (       -0.4538000000u,        0.3786000000u)--(       -0.4544000000u,        0.3761000000u) withpen pencircle scaled 1.0;
draw (       -0.1283000000u,        0.4918000000u)--(       -0.2061000000u,        0.4642000000u) withpen pencircle scaled 1.0;
draw (       -0.4544000000u,        0.3761000000u)--(       -0.4629000000u,        0.3384000000u) withpen pencircle scaled 1.0;
draw (       -0.4629000000u,        0.3384000000u)--(       -0.4581000000u,        0.3085000000u) withpen pencircle scaled 1.0;
draw (       -0.4581000000u,        0.3085000000u)--(       -0.4216000000u,        0.2683000000u) withpen pencircle scaled 1.0;
draw (       -0.3371000000u,        0.1140000000u)--(       -0.3760000000u,        0.1205000000u) withpen pencircle scaled 1.0;
draw (       -0.3760000000u,        0.1205000000u)--(       -0.4305000000u,        0.1084000000u) withpen pencircle scaled 1.0;
draw (       -0.4305000000u,        0.1084000000u)--(       -0.4586000000u,        0.0880000000u) withpen pencircle scaled 1.0;
draw (       -0.4586000000u,        0.0880000000u)--(       -0.4642000000u,        0.0543000000u) withpen pencircle scaled 1.0;
draw (       -0.4642000000u,        0.0543000000u)--(       -0.4279000000u,        0.0150000000u) withpen pencircle scaled 1.0;
draw (       -0.4279000000u,        0.0150000000u)--(       -0.4178000000u,        0.0132000000u) withpen pencircle scaled 1.0;
draw (       -0.3111000000u,        0.2346000000u)--(       -0.3293000000u,        0.1866000000u) withpen pencircle scaled 1.0;
draw (       -0.3293000000u,        0.1866000000u)--(       -0.3371000000u,        0.1140000000u) withpen pencircle scaled 1.0;
draw (       -0.3215000000u,        0.2566000000u)--(       -0.3111000000u,        0.2346000000u) withpen pencircle scaled 1.0;
draw (        0.1182000000u,        0.4629000000u)--(        0.0507000000u,        0.4724000000u) withpen pencircle scaled 1.0;
draw (       -0.1220000000u,        0.4966000000u)--(       -0.1283000000u,        0.4918000000u) withpen pencircle scaled 1.0;
draw (       -0.4216000000u,        0.2683000000u)--(       -0.4188000000u,        0.2653000000u) withpen pencircle scaled 1.0;
draw (       -0.4188000000u,        0.2653000000u)--(       -0.3682000000u,        0.2502000000u) withpen pencircle scaled 1.0;
draw (       -0.3682000000u,        0.2502000000u)--(       -0.3215000000u,        0.2566000000u) withpen pencircle scaled 1.0;
draw (        0.2453000000u,        0.2645000000u)--(        0.2276000000u,        0.3435000000u) withpen pencircle scaled 1.0;
draw (        0.1877000000u,        0.5212000000u)--(        0.1182000000u,        0.4629000000u) withpen pencircle scaled 1.0;
draw (       -0.4178000000u,        0.0132000000u)--(       -0.3098000000u,       -0.0053000000u) withpen pencircle scaled 1.0;
draw (       -0.3098000000u,       -0.0053000000u)--(       -0.2839000000u,       -0.0196000000u) withpen pencircle scaled 1.0;
draw (       -0.3699000000u,       -0.1722000000u)--(       -0.3695000000u,       -0.1938000000u) withpen pencircle scaled 1.0;
draw (       -0.3668000000u,       -0.3348000000u)--(       -0.3630000000u,       -0.3378000000u) withpen pencircle scaled 1.0;
draw (       -0.3630000000u,       -0.3378000000u)--(       -0.3189000000u,       -0.3352000000u) withpen pencircle scaled 1.0;
draw (       -0.3189000000u,       -0.3352000000u)--(       -0.2839000000u,       -0.2980000000u) withpen pencircle scaled 1.0;
draw (       -0.2839000000u,       -0.0196000000u)--(       -0.3626000000u,       -0.1428000000u) withpen pencircle scaled 1.0;
draw (       -0.3626000000u,       -0.1428000000u)--(       -0.3699000000u,       -0.1722000000u) withpen pencircle scaled 1.0;

x1 = .3317165748u;
y1 = -.1525697484u;
draw fullcircle scaled D shifted z1 withcolor RoyalBlue withpen pencircle scaled 1.0;

x2 = -.0898266834u;
y2 = .0040553762u;
draw fullcircle scaled D shifted z2 withcolor RoyalBlue withpen pencircle scaled 1.0;

x3 = -.1632029328u;
y3 = .3628813280u;
draw fullcircle scaled D shifted z3 withcolor RoyalBlue withpen pencircle scaled 1.0;

x4 = -.3329716402u;
y4 = -.1336561778u;
draw fullcircle scaled D shifted z4 withcolor RoyalBlue withpen pencircle scaled 1.0;

x5 = .1672151002u;
y5 = .3948440071u;
draw fullcircle scaled D shifted z5 withcolor RoyalBlue withpen pencircle scaled 1.0;

x6 = .1268489942u;
y6 = .0842089497u;
draw fullcircle scaled D shifted z6 withcolor RoyalBlue withpen pencircle scaled 1.0;

x7 = -.4017156699u;
y7 = .2383062811u;
draw fullcircle scaled D shifted z7 withcolor RoyalBlue withpen pencircle scaled 1.0;

x8 = .0355961500u;
y8 = -.3068175203u;
draw fullcircle scaled D shifted z8 withcolor RoyalBlue withpen pencircle scaled 1.0;
draw (        0.3317165748u,       -0.1525697484u)--(        0.5358000000u,       -0.1830000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3317165748u,       -0.1525697484u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3317165748u,       -0.1525697484u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5358000000u,       -0.1830000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5358000000u,       -0.1830000000u) withcolor Yellow;
draw (        0.3317165748u,       -0.1525697484u)--(        0.1276373551u,       -0.1221268801u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1268489942u,        0.0842089497u)--(        0.1276373551u,       -0.1221268801u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3317165748u,       -0.1525697484u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3317165748u,       -0.1525697484u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1268489942u,        0.0842089497u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1268489942u,        0.0842089497u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1276373551u,       -0.1221268801u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1276373551u,       -0.1221268801u) withcolor MidnightBlue;
draw (        0.3317165748u,       -0.1525697484u)--(        0.1276346793u,       -0.1221448226u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0355961500u,       -0.3068175203u)--(        0.1276346793u,       -0.1221448226u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3317165748u,       -0.1525697484u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3317165748u,       -0.1525697484u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0355961500u,       -0.3068175203u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0355961500u,       -0.3068175203u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1276346793u,       -0.1221448226u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1276346793u,       -0.1221448226u) withcolor MidnightBlue;
draw (       -0.0898266834u,        0.0040553762u)--(       -0.2862270865u,        0.0673165638u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.3329716402u,       -0.1336561778u)--(       -0.2862270865u,        0.0673165638u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.0898266834u,        0.0040553762u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0898266834u,        0.0040553762u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.3329716402u,       -0.1336561778u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.3329716402u,       -0.1336561778u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.2862270865u,        0.0673165638u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.2862270865u,        0.0673165638u) withcolor MidnightBlue;
draw (       -0.0898266834u,        0.0040553762u)--(       -0.2862270863u,        0.0673165644u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.4017156699u,        0.2383062811u)--(       -0.2862270863u,        0.0673165644u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.0898266834u,        0.0040553762u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0898266834u,        0.0040553762u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.4017156699u,        0.2383062811u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.4017156699u,        0.2383062811u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.2862270863u,        0.0673165644u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.2862270863u,        0.0673165644u) withcolor MidnightBlue;
draw (       -0.1632029328u,        0.3628813280u)--(       -0.0097964136u,        0.5008724003u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1672151002u,        0.3948440071u)--(       -0.0097964136u,        0.5008724003u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.1632029328u,        0.3628813280u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1632029328u,        0.3628813280u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1672151002u,        0.3948440071u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1672151002u,        0.3948440071u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0097964136u,        0.5008724003u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0097964136u,        0.5008724003u) withcolor Green;
draw (       -0.1632029328u,        0.3628813280u)--(        0.0138085809u,        0.2568529348u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1672151002u,        0.3948440071u)--(        0.0138085809u,        0.2568529348u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.1632029328u,        0.3628813280u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1632029328u,        0.3628813280u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1672151002u,        0.3948440071u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1672151002u,        0.3948440071u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0138085809u,        0.2568529348u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0138085809u,        0.2568529348u) withcolor MidnightBlue;
draw (       -0.1632029328u,        0.3628813280u)--(        0.0137829193u,        0.2568101052u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1268489942u,        0.0842089497u)--(        0.0137829193u,        0.2568101052u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.1632029328u,        0.3628813280u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1632029328u,        0.3628813280u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1268489942u,        0.0842089497u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1268489942u,        0.0842089497u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0137829193u,        0.2568101052u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0137829193u,        0.2568101052u) withcolor MidnightBlue;
draw (       -0.3329716402u,       -0.1336561778u)--(       -0.3630000000u,       -0.3378000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.3329716402u,       -0.1336561778u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.3329716402u,       -0.1336561778u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.3630000000u,       -0.3378000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.3630000000u,       -0.3378000000u) withcolor Yellow;
draw (       -0.3329716402u,       -0.1336561778u)--(       -0.2862270870u,        0.0673165639u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.4017156699u,        0.2383062811u)--(       -0.2862270870u,        0.0673165639u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.3329716402u,       -0.1336561778u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.3329716402u,       -0.1336561778u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.4017156699u,        0.2383062811u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.4017156699u,        0.2383062811u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.2862270870u,        0.0673165639u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.2862270870u,        0.0673165639u) withcolor MidnightBlue;
draw (        0.1672151002u,        0.3948440071u)--(        0.3387000000u,        0.5096000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1672151002u,        0.3948440071u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1672151002u,        0.3948440071u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3387000000u,        0.5096000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3387000000u,        0.5096000000u) withcolor Yellow;
draw (        0.1672151002u,        0.3948440071u)--(        0.0138233127u,        0.2568365590u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1268489942u,        0.0842089497u)--(        0.0138233127u,        0.2568365590u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1672151002u,        0.3948440071u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1672151002u,        0.3948440071u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1268489942u,        0.0842089497u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1268489942u,        0.0842089497u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0138233127u,        0.2568365590u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0138233127u,        0.2568365590u) withcolor MidnightBlue;
draw (        0.1268489942u,        0.0842089497u)--(        0.3326000000u,        0.0686000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1268489942u,        0.0842089497u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1268489942u,        0.0842089497u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3326000000u,        0.0686000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3326000000u,        0.0686000000u) withcolor Yellow;
draw (        0.1268489942u,        0.0842089497u)--(        0.1275989568u,       -0.1221270232u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0355961500u,       -0.3068175203u)--(        0.1275989568u,       -0.1221270232u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1268489942u,        0.0842089497u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1268489942u,        0.0842089497u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0355961500u,       -0.3068175203u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0355961500u,       -0.3068175203u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1275989568u,       -0.1221270232u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1275989568u,       -0.1221270232u) withcolor MidnightBlue;
draw (        0.0355961500u,       -0.3068175203u)--(       -0.1702000000u,       -0.3218000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0355961500u,       -0.3068175203u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0355961500u,       -0.3068175203u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.1702000000u,       -0.3218000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1702000000u,       -0.3218000000u) withcolor Yellow;
draw (        0.0355961500u,       -0.3068175203u)--(        0.1571000000u,       -0.4736000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0355961500u,       -0.3068175203u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0355961500u,       -0.3068175203u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1571000000u,       -0.4736000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1571000000u,       -0.4736000000u) withcolor Yellow;
% r =         0.2063373359
endfig;
end;
