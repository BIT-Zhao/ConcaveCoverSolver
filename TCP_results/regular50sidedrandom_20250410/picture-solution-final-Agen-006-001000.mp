prologues := 3;
outputformat   := "svg";
outputtemplate := "%j.%o";
input mpcolornames;
beginfig(6);
u := 20cm;
path pols[];
pair p[];
R := .2775544345u;
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

x1 = -.4140386174u;
y1 = -.0314637699u;
draw fullcircle scaled D shifted z1 withcolor RoyalBlue withpen pencircle scaled 1.0;

x2 = -.1609663359u;
y2 = .2037092007u;
draw fullcircle scaled D shifted z2 withcolor RoyalBlue withpen pencircle scaled 1.0;

x3 = -.1152466009u;
y3 = -.3172547395u;
draw fullcircle scaled D shifted z3 withcolor RoyalBlue withpen pencircle scaled 1.0;

x4 = .1729297645u;
y4 = -.2200628625u;
draw fullcircle scaled D shifted z4 withcolor RoyalBlue withpen pencircle scaled 1.0;

x5 = .1278134009u;
y5 = .3946693862u;
draw fullcircle scaled D shifted z5 withcolor RoyalBlue withpen pencircle scaled 1.0;

x6 = .3348373034u;
y6 = .0375206994u;
draw fullcircle scaled D shifted z6 withcolor RoyalBlue withpen pencircle scaled 1.0;
draw (       -0.1609663359u,        0.2037092007u)--(       -0.0453874153u,       -0.0486357560u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.1152466009u,       -0.3172547395u)--(       -0.0453874153u,       -0.0486357560u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.1609663359u,        0.2037092007u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1609663359u,        0.2037092007u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.1152466009u,       -0.3172547395u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1152466009u,       -0.3172547395u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0453874153u,       -0.0486357560u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0453874153u,       -0.0486357560u) withcolor MidnightBlue;
draw (       -0.1609663359u,        0.2037092007u)--(        0.0573036900u,        0.0322604915u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1729297645u,       -0.2200628625u)--(        0.0573036900u,        0.0322604915u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.1609663359u,        0.2037092007u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1609663359u,        0.2037092007u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1729297645u,       -0.2200628625u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1729297645u,       -0.2200628625u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0573036900u,        0.0322604915u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0573036900u,        0.0322604915u) withcolor MidnightBlue;
draw (       -0.1609663359u,        0.2037092007u)--(       -0.0453402614u,       -0.0486141533u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1729297645u,       -0.2200628625u)--(       -0.0453402614u,       -0.0486141533u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.1609663359u,        0.2037092007u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1609663359u,        0.2037092007u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1729297645u,       -0.2200628625u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1729297645u,       -0.2200628625u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0453402614u,       -0.0486141533u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0453402614u,       -0.0486141533u) withcolor MidnightBlue;
draw (       -0.1609663359u,        0.2037092007u)--(        0.0573320372u,        0.0322965863u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3348373034u,        0.0375206994u)--(        0.0573320372u,        0.0322965863u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.1609663359u,        0.2037092007u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1609663359u,        0.2037092007u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3348373034u,        0.0375206994u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3348373034u,        0.0375206994u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0573320372u,        0.0322965863u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0573320372u,        0.0322965863u) withcolor MidnightBlue;
draw (       -0.1152466009u,       -0.3172547395u)--(       -0.0453623487u,       -0.0486422763u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1729297645u,       -0.2200628625u)--(       -0.0453623487u,       -0.0486422763u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.1152466009u,       -0.3172547395u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1152466009u,       -0.3172547395u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1729297645u,       -0.2200628625u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1729297645u,       -0.2200628625u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0453623487u,       -0.0486422763u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0453623487u,       -0.0486422763u) withcolor MidnightBlue;
draw (        0.1729297645u,       -0.2200628625u)--(        0.0573324694u,        0.0322736776u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3348373034u,        0.0375206994u)--(        0.0573324694u,        0.0322736776u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1729297645u,       -0.2200628625u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1729297645u,       -0.2200628625u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3348373034u,        0.0375206994u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3348373034u,        0.0375206994u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0573324694u,        0.0322736776u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0573324694u,        0.0322736776u) withcolor MidnightBlue;
% r =         0.2775544345
endfig;
end;
