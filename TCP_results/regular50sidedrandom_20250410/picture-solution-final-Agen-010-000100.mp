prologues := 3;
outputformat   := "svg";
outputtemplate := "%j.%o";
input mpcolornames;
beginfig(10);
u := 20cm;
path pols[];
pair p[];
R := .1971924491u;
D := 2*R;

pols[1] = 
(       -0.0314000000u,       -0.4990000000u)--
(        0.0314000000u,       -0.4990000000u)--
(        0.0937000000u,       -0.4911000000u)--
(        0.1545000000u,       -0.4755000000u)--
(        0.2129000000u,       -0.4524000000u)--
(        0.2679000000u,       -0.4222000000u)--
(        0.3187000000u,       -0.3853000000u)--
(        0.3645000000u,       -0.3423000000u)--
(        0.4045000000u,       -0.2939000000u)--
(        0.4382000000u,       -0.2409000000u)--
(        0.4649000000u,       -0.1841000000u)--
(        0.4843000000u,       -0.1243000000u)--
(        0.4961000000u,       -0.0627000000u)--
(        0.5000000000u,        0.0000000000u)--
(        0.4961000000u,        0.0627000000u)--
(        0.4843000000u,        0.1243000000u)--
(        0.4649000000u,        0.1841000000u)--
(        0.4382000000u,        0.2409000000u)--
(        0.4045000000u,        0.2939000000u)--
(        0.3645000000u,        0.3423000000u)--
(        0.3187000000u,        0.3853000000u)--
(        0.2679000000u,        0.4222000000u)--
(        0.2129000000u,        0.4524000000u)--
(        0.1545000000u,        0.4755000000u)--
(        0.0937000000u,        0.4911000000u)--
(        0.0314000000u,        0.4990000000u)--
(       -0.0314000000u,        0.4990000000u)--
(       -0.0937000000u,        0.4911000000u)--
(       -0.1545000000u,        0.4755000000u)--
(       -0.2129000000u,        0.4524000000u)--
(       -0.2679000000u,        0.4222000000u)--
(       -0.3187000000u,        0.3853000000u)--
(       -0.3645000000u,        0.3423000000u)--
(       -0.4045000000u,        0.2939000000u)--
(       -0.4382000000u,        0.2409000000u)--
(       -0.4649000000u,        0.1841000000u)--
(       -0.4843000000u,        0.1243000000u)--
(       -0.4961000000u,        0.0627000000u)--
(       -0.5000000000u,       -0.0000000000u)--
(       -0.4961000000u,       -0.0627000000u)--
(       -0.4843000000u,       -0.1243000000u)--
(       -0.4649000000u,       -0.1841000000u)--
(       -0.4382000000u,       -0.2409000000u)--
(       -0.4045000000u,       -0.2939000000u)--
(       -0.3645000000u,       -0.3423000000u)--
(       -0.3187000000u,       -0.3853000000u)--
(       -0.2679000000u,       -0.4222000000u)--
(       -0.2129000000u,       -0.4524000000u)--
(       -0.1545000000u,       -0.4755000000u)--
(       -0.0937000000u,       -0.4911000000u)--
(       -0.0314000000u,       -0.4990000000u)--
cycle;
fill pols[1] withcolor 0.9Dandelion;
draw (       -0.0314000000u,       -0.4990000000u)--(        0.0314000000u,       -0.4990000000u) withpen pencircle scaled 1.0;
draw (        0.0314000000u,       -0.4990000000u)--(        0.0937000000u,       -0.4911000000u) withpen pencircle scaled 1.0;
draw (        0.0937000000u,       -0.4911000000u)--(        0.1545000000u,       -0.4755000000u) withpen pencircle scaled 1.0;
draw (        0.1545000000u,       -0.4755000000u)--(        0.2129000000u,       -0.4524000000u) withpen pencircle scaled 1.0;
draw (        0.2129000000u,       -0.4524000000u)--(        0.2679000000u,       -0.4222000000u) withpen pencircle scaled 1.0;
draw (        0.2679000000u,       -0.4222000000u)--(        0.3187000000u,       -0.3853000000u) withpen pencircle scaled 1.0;
draw (        0.3187000000u,       -0.3853000000u)--(        0.3645000000u,       -0.3423000000u) withpen pencircle scaled 1.0;
draw (        0.3645000000u,       -0.3423000000u)--(        0.4045000000u,       -0.2939000000u) withpen pencircle scaled 1.0;
draw (        0.4045000000u,       -0.2939000000u)--(        0.4382000000u,       -0.2409000000u) withpen pencircle scaled 1.0;
draw (        0.4382000000u,       -0.2409000000u)--(        0.4649000000u,       -0.1841000000u) withpen pencircle scaled 1.0;
draw (        0.4649000000u,       -0.1841000000u)--(        0.4843000000u,       -0.1243000000u) withpen pencircle scaled 1.0;
draw (        0.4843000000u,       -0.1243000000u)--(        0.4961000000u,       -0.0627000000u) withpen pencircle scaled 1.0;
draw (        0.4961000000u,       -0.0627000000u)--(        0.5000000000u,        0.0000000000u) withpen pencircle scaled 1.0;
draw (        0.5000000000u,        0.0000000000u)--(        0.4961000000u,        0.0627000000u) withpen pencircle scaled 1.0;
draw (        0.4961000000u,        0.0627000000u)--(        0.4843000000u,        0.1243000000u) withpen pencircle scaled 1.0;
draw (        0.4843000000u,        0.1243000000u)--(        0.4649000000u,        0.1841000000u) withpen pencircle scaled 1.0;
draw (        0.4649000000u,        0.1841000000u)--(        0.4382000000u,        0.2409000000u) withpen pencircle scaled 1.0;
draw (        0.4382000000u,        0.2409000000u)--(        0.4045000000u,        0.2939000000u) withpen pencircle scaled 1.0;
draw (        0.4045000000u,        0.2939000000u)--(        0.3645000000u,        0.3423000000u) withpen pencircle scaled 1.0;
draw (        0.3645000000u,        0.3423000000u)--(        0.3187000000u,        0.3853000000u) withpen pencircle scaled 1.0;
draw (        0.3187000000u,        0.3853000000u)--(        0.2679000000u,        0.4222000000u) withpen pencircle scaled 1.0;
draw (        0.2679000000u,        0.4222000000u)--(        0.2129000000u,        0.4524000000u) withpen pencircle scaled 1.0;
draw (        0.2129000000u,        0.4524000000u)--(        0.1545000000u,        0.4755000000u) withpen pencircle scaled 1.0;
draw (        0.1545000000u,        0.4755000000u)--(        0.0937000000u,        0.4911000000u) withpen pencircle scaled 1.0;
draw (        0.0937000000u,        0.4911000000u)--(        0.0314000000u,        0.4990000000u) withpen pencircle scaled 1.0;
draw (        0.0314000000u,        0.4990000000u)--(       -0.0314000000u,        0.4990000000u) withpen pencircle scaled 1.0;
draw (       -0.0314000000u,        0.4990000000u)--(       -0.0937000000u,        0.4911000000u) withpen pencircle scaled 1.0;
draw (       -0.0937000000u,        0.4911000000u)--(       -0.1545000000u,        0.4755000000u) withpen pencircle scaled 1.0;
draw (       -0.1545000000u,        0.4755000000u)--(       -0.2129000000u,        0.4524000000u) withpen pencircle scaled 1.0;
draw (       -0.2129000000u,        0.4524000000u)--(       -0.2679000000u,        0.4222000000u) withpen pencircle scaled 1.0;
draw (       -0.2679000000u,        0.4222000000u)--(       -0.3187000000u,        0.3853000000u) withpen pencircle scaled 1.0;
draw (       -0.3187000000u,        0.3853000000u)--(       -0.3645000000u,        0.3423000000u) withpen pencircle scaled 1.0;
draw (       -0.3645000000u,        0.3423000000u)--(       -0.4045000000u,        0.2939000000u) withpen pencircle scaled 1.0;
draw (       -0.4045000000u,        0.2939000000u)--(       -0.4382000000u,        0.2409000000u) withpen pencircle scaled 1.0;
draw (       -0.4382000000u,        0.2409000000u)--(       -0.4649000000u,        0.1841000000u) withpen pencircle scaled 1.0;
draw (       -0.4649000000u,        0.1841000000u)--(       -0.4843000000u,        0.1243000000u) withpen pencircle scaled 1.0;
draw (       -0.4843000000u,        0.1243000000u)--(       -0.4961000000u,        0.0627000000u) withpen pencircle scaled 1.0;
draw (       -0.4961000000u,        0.0627000000u)--(       -0.5000000000u,       -0.0000000000u) withpen pencircle scaled 1.0;
draw (       -0.5000000000u,       -0.0000000000u)--(       -0.4961000000u,       -0.0627000000u) withpen pencircle scaled 1.0;
draw (       -0.4961000000u,       -0.0627000000u)--(       -0.4843000000u,       -0.1243000000u) withpen pencircle scaled 1.0;
draw (       -0.4843000000u,       -0.1243000000u)--(       -0.4649000000u,       -0.1841000000u) withpen pencircle scaled 1.0;
draw (       -0.4649000000u,       -0.1841000000u)--(       -0.4382000000u,       -0.2409000000u) withpen pencircle scaled 1.0;
draw (       -0.4382000000u,       -0.2409000000u)--(       -0.4045000000u,       -0.2939000000u) withpen pencircle scaled 1.0;
draw (       -0.4045000000u,       -0.2939000000u)--(       -0.3645000000u,       -0.3423000000u) withpen pencircle scaled 1.0;
draw (       -0.3645000000u,       -0.3423000000u)--(       -0.3187000000u,       -0.3853000000u) withpen pencircle scaled 1.0;
draw (       -0.3187000000u,       -0.3853000000u)--(       -0.2679000000u,       -0.4222000000u) withpen pencircle scaled 1.0;
draw (       -0.2679000000u,       -0.4222000000u)--(       -0.2129000000u,       -0.4524000000u) withpen pencircle scaled 1.0;
draw (       -0.2129000000u,       -0.4524000000u)--(       -0.1545000000u,       -0.4755000000u) withpen pencircle scaled 1.0;
draw (       -0.1545000000u,       -0.4755000000u)--(       -0.0937000000u,       -0.4911000000u) withpen pencircle scaled 1.0;
draw (       -0.0937000000u,       -0.4911000000u)--(       -0.0314000000u,       -0.4990000000u) withpen pencircle scaled 1.0;

x1 = .3698278570u;
y1 = .0011127168u;
draw fullcircle scaled D shifted z1 withcolor RoyalBlue withpen pencircle scaled 1.0;

x2 = -.1808083746u;
y2 = -.3227676418u;
draw fullcircle scaled D shifted z2 withcolor RoyalBlue withpen pencircle scaled 1.0;

x3 = -.3476924275u;
y3 = -.1311159168u;
draw fullcircle scaled D shifted z3 withcolor RoyalBlue withpen pencircle scaled 1.0;

x4 = -.1873314678u;
y4 = .3212991463u;
draw fullcircle scaled D shifted z4 withcolor RoyalBlue withpen pencircle scaled 1.0;

x5 = -.0015788330u;
y5 = .0013298825u;
draw fullcircle scaled D shifted z5 withcolor RoyalBlue withpen pencircle scaled 1.0;

x6 = .2832805056u;
y6 = .2380333054u;
draw fullcircle scaled D shifted z6 withcolor RoyalBlue withpen pencircle scaled 1.0;

x7 = .0637217653u;
y7 = .3660854723u;
draw fullcircle scaled D shifted z7 withcolor RoyalBlue withpen pencircle scaled 1.0;

x8 = .0681383199u;
y8 = -.3634805494u;
draw fullcircle scaled D shifted z8 withcolor RoyalBlue withpen pencircle scaled 1.0;

x9 = .2840958080u;
y9 = -.2350714925u;
draw fullcircle scaled D shifted z9 withcolor RoyalBlue withpen pencircle scaled 1.0;

x10 = -.3507022702u;
y10 = .1238003361u;
draw fullcircle scaled D shifted z10 withcolor RoyalBlue withpen pencircle scaled 1.0;
draw (        0.3698278570u,        0.0011127168u)--(        0.1840857305u,       -0.0651045350u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.0015788330u,        0.0013298825u)--(        0.1840857305u,       -0.0651045350u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3698278570u,        0.0011127168u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3698278570u,        0.0011127168u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0015788330u,        0.0013298825u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0015788330u,        0.0013298825u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1840857305u,       -0.0651045350u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1840857305u,       -0.0651045350u) withcolor MidnightBlue;
draw (        0.3698278570u,        0.0011127168u)--(        0.1841632935u,        0.0675471343u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.0015788330u,        0.0013298825u)--(        0.1841632935u,        0.0675471343u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3698278570u,        0.0011127168u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3698278570u,        0.0011127168u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0015788330u,        0.0013298825u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0015788330u,        0.0013298825u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1841632935u,        0.0675471343u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1841632935u,        0.0675471343u) withcolor MidnightBlue;
draw (        0.3698278570u,        0.0011127168u)--(        0.1841675456u,        0.0675590166u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2832805056u,        0.2380333054u)--(        0.1841675456u,        0.0675590166u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3698278570u,        0.0011127168u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3698278570u,        0.0011127168u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2832805056u,        0.2380333054u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2832805056u,        0.2380333054u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1841675456u,        0.0675590166u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1841675456u,        0.0675590166u) withcolor MidnightBlue;
draw (        0.3698278570u,        0.0011127168u)--(        0.1840908635u,       -0.0651189313u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2840958080u,       -0.2350714925u)--(        0.1840908635u,       -0.0651189313u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3698278570u,        0.0011127168u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3698278570u,        0.0011127168u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2840958080u,       -0.2350714925u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2840958080u,       -0.2350714925u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1840908635u,       -0.0651189313u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1840908635u,       -0.0651189313u) withcolor MidnightBlue;
draw (       -0.1808083746u,       -0.3227676418u)--(       -0.1505259714u,       -0.1279142661u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.3476924275u,       -0.1311159168u)--(       -0.1505259714u,       -0.1279142661u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.1808083746u,       -0.3227676418u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1808083746u,       -0.3227676418u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.3476924275u,       -0.1311159168u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.3476924275u,       -0.1311159168u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.1505259714u,       -0.1279142661u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1505259714u,       -0.1279142661u) withcolor MidnightBlue;
draw (       -0.1808083746u,       -0.3227676418u)--(       -0.1505087493u,       -0.1279169434u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.0015788330u,        0.0013298825u)--(       -0.1505087493u,       -0.1279169434u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.1808083746u,       -0.3227676418u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1808083746u,       -0.3227676418u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0015788330u,        0.0013298825u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0015788330u,        0.0013298825u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.1505087493u,       -0.1279169434u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1505087493u,       -0.1279169434u) withcolor MidnightBlue;
draw (       -0.1808083746u,       -0.3227676418u)--(       -0.0318784583u,       -0.1935208159u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.0015788330u,        0.0013298825u)--(       -0.0318784583u,       -0.1935208159u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.1808083746u,       -0.3227676418u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1808083746u,       -0.3227676418u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0015788330u,        0.0013298825u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0015788330u,        0.0013298825u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0318784583u,       -0.1935208159u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0318784583u,       -0.1935208159u) withcolor MidnightBlue;
draw (       -0.1808083746u,       -0.3227676418u)--(       -0.0318703381u,       -0.1935301734u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0681383199u,       -0.3634805494u)--(       -0.0318703381u,       -0.1935301734u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.1808083746u,       -0.3227676418u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1808083746u,       -0.3227676418u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0681383199u,       -0.3634805494u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0681383199u,       -0.3634805494u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0318703381u,       -0.1935301734u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0318703381u,       -0.1935301734u) withcolor MidnightBlue;
draw (       -0.3476924275u,       -0.1311159168u)--(       -0.1987450041u,       -0.0018892669u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.0015788330u,        0.0013298825u)--(       -0.1987450041u,       -0.0018892669u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.3476924275u,       -0.1311159168u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.3476924275u,       -0.1311159168u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0015788330u,        0.0013298825u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0015788330u,        0.0013298825u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.1987450041u,       -0.0018892669u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1987450041u,       -0.0018892669u) withcolor MidnightBlue;
draw (       -0.3476924275u,       -0.1311159168u)--(       -0.1505262564u,       -0.1278967674u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.0015788330u,        0.0013298825u)--(       -0.1505262564u,       -0.1278967674u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.3476924275u,       -0.1311159168u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.3476924275u,       -0.1311159168u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0015788330u,        0.0013298825u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0015788330u,        0.0013298825u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.1505262564u,       -0.1278967674u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1505262564u,       -0.1278967674u) withcolor MidnightBlue;
draw (       -0.3476924275u,       -0.1311159168u)--(       -0.1987517846u,       -0.0018814521u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.3507022702u,        0.1238003361u)--(       -0.1987517846u,       -0.0018814521u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.3476924275u,       -0.1311159168u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.3476924275u,       -0.1311159168u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.3507022702u,        0.1238003361u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.3507022702u,        0.1238003361u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.1987517846u,       -0.0018814521u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1987517846u,       -0.0018814521u) withcolor MidnightBlue;
draw (       -0.1873314678u,        0.3212991463u)--(       -0.0353942785u,        0.1956012844u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.0015788330u,        0.0013298825u)--(       -0.0353942785u,        0.1956012844u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.1873314678u,        0.3212991463u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1873314678u,        0.3212991463u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0015788330u,        0.0013298825u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0015788330u,        0.0013298825u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0353942785u,        0.1956012844u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0353942785u,        0.1956012844u) withcolor MidnightBlue;
draw (       -0.1873314678u,        0.3212991463u)--(       -0.1535160223u,        0.1270277444u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.0015788330u,        0.0013298825u)--(       -0.1535160223u,        0.1270277444u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.1873314678u,        0.3212991463u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1873314678u,        0.3212991463u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0015788330u,        0.0013298825u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0015788330u,        0.0013298825u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.1535160223u,        0.1270277444u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1535160223u,        0.1270277444u) withcolor MidnightBlue;
draw (       -0.1873314678u,        0.3212991463u)--(       -0.0353877474u,        0.1956091793u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0637217653u,        0.3660854723u)--(       -0.0353877474u,        0.1956091793u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.1873314678u,        0.3212991463u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1873314678u,        0.3212991463u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0637217653u,        0.3660854723u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0637217653u,        0.3660854723u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0353877474u,        0.1956091793u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0353877474u,        0.1956091793u) withcolor MidnightBlue;
draw (       -0.1873314678u,        0.3212991463u)--(       -0.1535361767u,        0.1270242373u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.3507022702u,        0.1238003361u)--(       -0.1535361767u,        0.1270242373u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.1873314678u,        0.3212991463u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1873314678u,        0.3212991463u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.3507022702u,        0.1238003361u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.3507022702u,        0.1238003361u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.1535361767u,        0.1270242373u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1535361767u,        0.1270242373u) withcolor MidnightBlue;
draw (       -0.0015788330u,        0.0013298825u)--(        0.0975448360u,        0.1717979447u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2832805056u,        0.2380333054u)--(        0.0975448360u,        0.1717979447u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.0015788330u,        0.0013298825u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0015788330u,        0.0013298825u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2832805056u,        0.2380333054u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2832805056u,        0.2380333054u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0975448360u,        0.1717979447u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0975448360u,        0.1717979447u) withcolor MidnightBlue;
draw (       -0.0015788330u,        0.0013298825u)--(        0.1841568366u,        0.0675652433u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2832805056u,        0.2380333054u)--(        0.1841568366u,        0.0675652433u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.0015788330u,        0.0013298825u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0015788330u,        0.0013298825u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2832805056u,        0.2380333054u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2832805056u,        0.2380333054u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1841568366u,        0.0675652433u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1841568366u,        0.0675652433u) withcolor MidnightBlue;
draw (       -0.0015788330u,        0.0013298825u)--(       -0.0353788002u,        0.1956039780u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0637217653u,        0.3660854723u)--(       -0.0353788002u,        0.1956039780u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.0015788330u,        0.0013298825u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0015788330u,        0.0013298825u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0637217653u,        0.3660854723u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0637217653u,        0.3660854723u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0353788002u,        0.1956039780u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0353788002u,        0.1956039780u) withcolor MidnightBlue;
draw (       -0.0015788330u,        0.0013298825u)--(        0.0975217326u,        0.1718113767u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0637217653u,        0.3660854723u)--(        0.0975217326u,        0.1718113767u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.0015788330u,        0.0013298825u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0015788330u,        0.0013298825u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0637217653u,        0.3660854723u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0637217653u,        0.3660854723u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0975217326u,        0.1718113767u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0975217326u,        0.1718113767u) withcolor MidnightBlue;
draw (       -0.0015788330u,        0.0013298825u)--(        0.0984189425u,       -0.1686268970u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0681383199u,       -0.3634805494u)--(        0.0984189425u,       -0.1686268970u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.0015788330u,        0.0013298825u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0015788330u,        0.0013298825u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0681383199u,       -0.3634805494u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0681383199u,       -0.3634805494u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0984189425u,       -0.1686268970u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0984189425u,       -0.1686268970u) withcolor MidnightBlue;
draw (       -0.0015788330u,        0.0013298825u)--(       -0.0318594555u,       -0.1935237700u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0681383199u,       -0.3634805494u)--(       -0.0318594555u,       -0.1935237700u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.0015788330u,        0.0013298825u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0015788330u,        0.0013298825u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0681383199u,       -0.3634805494u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0681383199u,       -0.3634805494u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0318594555u,       -0.1935237700u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0318594555u,       -0.1935237700u) withcolor MidnightBlue;
draw (       -0.0015788330u,        0.0013298825u)--(        0.1840778342u,       -0.0651265987u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2840958080u,       -0.2350714925u)--(        0.1840778342u,       -0.0651265987u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.0015788330u,        0.0013298825u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0015788330u,        0.0013298825u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2840958080u,       -0.2350714925u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2840958080u,       -0.2350714925u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1840778342u,       -0.0651265987u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1840778342u,       -0.0651265987u) withcolor MidnightBlue;
draw (       -0.0015788330u,        0.0013298825u)--(        0.0984391408u,       -0.1686150112u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2840958080u,       -0.2350714925u)--(        0.0984391408u,       -0.1686150112u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.0015788330u,        0.0013298825u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0015788330u,        0.0013298825u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2840958080u,       -0.2350714925u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2840958080u,       -0.2350714925u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0984391408u,       -0.1686150112u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0984391408u,       -0.1686150112u) withcolor MidnightBlue;
draw (       -0.0015788330u,        0.0013298825u)--(       -0.1987452599u,       -0.0018735633u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.3507022702u,        0.1238003361u)--(       -0.1987452599u,       -0.0018735633u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.0015788330u,        0.0013298825u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0015788330u,        0.0013298825u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.3507022702u,        0.1238003361u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.3507022702u,        0.1238003361u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.1987452599u,       -0.0018735633u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1987452599u,       -0.0018735633u) withcolor MidnightBlue;
draw (       -0.0015788330u,        0.0013298825u)--(       -0.1535358433u,        0.1270037819u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.3507022702u,        0.1238003361u)--(       -0.1535358433u,        0.1270037819u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.0015788330u,        0.0013298825u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0015788330u,        0.0013298825u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.3507022702u,        0.1238003361u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.3507022702u,        0.1238003361u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.1535358433u,        0.1270037819u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1535358433u,        0.1270037819u) withcolor MidnightBlue;
draw (        0.2832805056u,        0.2380333054u)--(        0.0975389767u,        0.1718143777u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0637217653u,        0.3660854723u)--(        0.0975389767u,        0.1718143777u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2832805056u,        0.2380333054u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2832805056u,        0.2380333054u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0637217653u,        0.3660854723u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0637217653u,        0.3660854723u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0975389767u,        0.1718143777u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0975389767u,        0.1718143777u) withcolor MidnightBlue;
draw (        0.0681383199u,       -0.3634805494u)--(        0.0984340464u,       -0.1686292448u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2840958080u,       -0.2350714925u)--(        0.0984340464u,       -0.1686292448u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0681383199u,       -0.3634805494u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0681383199u,       -0.3634805494u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2840958080u,       -0.2350714925u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2840958080u,       -0.2350714925u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0984340464u,       -0.1686292448u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0984340464u,       -0.1686292448u) withcolor MidnightBlue;
% r =         0.1971924491
endfig;
end;
