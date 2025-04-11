prologues := 3;
outputformat   := "svg";
outputtemplate := "%j.%o";
input mpcolornames;
beginfig(1);
u := 20cm;
path pols[];
pair p[];
R := .5000408314u;
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

x1 = .0000000001u;
y1 = -.0000000001u;
draw fullcircle scaled D shifted z1 withcolor RoyalBlue withpen pencircle scaled 1.0;
draw (        0.0000000001u,       -0.0000000001u)--(       -0.0314000000u,       -0.4990000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0314000000u,       -0.4990000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0314000000u,       -0.4990000000u) withcolor Yellow;
draw (        0.0000000001u,       -0.0000000001u)--(        0.0314000000u,       -0.4990000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0314000000u,       -0.4990000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0314000000u,       -0.4990000000u) withcolor Yellow;
draw (        0.0000000001u,       -0.0000000001u)--(        0.0937000000u,       -0.4911000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0937000000u,       -0.4911000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0937000000u,       -0.4911000000u) withcolor Yellow;
draw (        0.0000000001u,       -0.0000000001u)--(        0.1545000000u,       -0.4755000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1545000000u,       -0.4755000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1545000000u,       -0.4755000000u) withcolor Yellow;
draw (        0.0000000001u,       -0.0000000001u)--(        0.2129000000u,       -0.4524000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2129000000u,       -0.4524000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2129000000u,       -0.4524000000u) withcolor Yellow;
draw (        0.0000000001u,       -0.0000000001u)--(        0.2679000000u,       -0.4222000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2679000000u,       -0.4222000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2679000000u,       -0.4222000000u) withcolor Yellow;
draw (        0.0000000001u,       -0.0000000001u)--(        0.3187000000u,       -0.3853000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3187000000u,       -0.3853000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3187000000u,       -0.3853000000u) withcolor Yellow;
draw (        0.0000000001u,       -0.0000000001u)--(        0.3645000000u,       -0.3423000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3645000000u,       -0.3423000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3645000000u,       -0.3423000000u) withcolor Yellow;
draw (        0.0000000001u,       -0.0000000001u)--(        0.4045000000u,       -0.2939000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4045000000u,       -0.2939000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4045000000u,       -0.2939000000u) withcolor Yellow;
draw (        0.0000000001u,       -0.0000000001u)--(        0.4382000000u,       -0.2409000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4382000000u,       -0.2409000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4382000000u,       -0.2409000000u) withcolor Yellow;
draw (        0.0000000001u,       -0.0000000001u)--(        0.4649000000u,       -0.1841000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4649000000u,       -0.1841000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4649000000u,       -0.1841000000u) withcolor Yellow;
draw (        0.0000000001u,       -0.0000000001u)--(        0.4843000000u,       -0.1243000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4843000000u,       -0.1243000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4843000000u,       -0.1243000000u) withcolor Yellow;
draw (        0.0000000001u,       -0.0000000001u)--(        0.4961000000u,       -0.0627000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4961000000u,       -0.0627000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4961000000u,       -0.0627000000u) withcolor Yellow;
draw (        0.0000000001u,       -0.0000000001u)--(        0.5000000000u,        0.0000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5000000000u,        0.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000000u,        0.0000000000u) withcolor Yellow;
draw (        0.0000000001u,       -0.0000000001u)--(        0.4961000000u,        0.0627000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4961000000u,        0.0627000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4961000000u,        0.0627000000u) withcolor Yellow;
draw (        0.0000000001u,       -0.0000000001u)--(        0.4843000000u,        0.1243000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4843000000u,        0.1243000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4843000000u,        0.1243000000u) withcolor Yellow;
draw (        0.0000000001u,       -0.0000000001u)--(        0.4649000000u,        0.1841000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4649000000u,        0.1841000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4649000000u,        0.1841000000u) withcolor Yellow;
draw (        0.0000000001u,       -0.0000000001u)--(        0.4382000000u,        0.2409000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4382000000u,        0.2409000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4382000000u,        0.2409000000u) withcolor Yellow;
draw (        0.0000000001u,       -0.0000000001u)--(        0.4045000000u,        0.2939000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4045000000u,        0.2939000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4045000000u,        0.2939000000u) withcolor Yellow;
draw (        0.0000000001u,       -0.0000000001u)--(        0.3645000000u,        0.3423000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3645000000u,        0.3423000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3645000000u,        0.3423000000u) withcolor Yellow;
draw (        0.0000000001u,       -0.0000000001u)--(        0.3187000000u,        0.3853000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3187000000u,        0.3853000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3187000000u,        0.3853000000u) withcolor Yellow;
draw (        0.0000000001u,       -0.0000000001u)--(        0.2679000000u,        0.4222000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2679000000u,        0.4222000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2679000000u,        0.4222000000u) withcolor Yellow;
draw (        0.0000000001u,       -0.0000000001u)--(        0.2129000000u,        0.4524000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2129000000u,        0.4524000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2129000000u,        0.4524000000u) withcolor Yellow;
draw (        0.0000000001u,       -0.0000000001u)--(        0.1545000000u,        0.4755000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1545000000u,        0.4755000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1545000000u,        0.4755000000u) withcolor Yellow;
draw (        0.0000000001u,       -0.0000000001u)--(        0.0937000000u,        0.4911000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0937000000u,        0.4911000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0937000000u,        0.4911000000u) withcolor Yellow;
draw (        0.0000000001u,       -0.0000000001u)--(        0.0314000000u,        0.4990000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0314000000u,        0.4990000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0314000000u,        0.4990000000u) withcolor Yellow;
draw (        0.0000000001u,       -0.0000000001u)--(       -0.0314000000u,        0.4990000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0314000000u,        0.4990000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0314000000u,        0.4990000000u) withcolor Yellow;
draw (        0.0000000001u,       -0.0000000001u)--(       -0.0937000000u,        0.4911000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0937000000u,        0.4911000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0937000000u,        0.4911000000u) withcolor Yellow;
draw (        0.0000000001u,       -0.0000000001u)--(       -0.1545000000u,        0.4755000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.1545000000u,        0.4755000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1545000000u,        0.4755000000u) withcolor Yellow;
draw (        0.0000000001u,       -0.0000000001u)--(       -0.2129000000u,        0.4524000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.2129000000u,        0.4524000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.2129000000u,        0.4524000000u) withcolor Yellow;
draw (        0.0000000001u,       -0.0000000001u)--(       -0.2679000000u,        0.4222000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.2679000000u,        0.4222000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.2679000000u,        0.4222000000u) withcolor Yellow;
draw (        0.0000000001u,       -0.0000000001u)--(       -0.3187000000u,        0.3853000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.3187000000u,        0.3853000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.3187000000u,        0.3853000000u) withcolor Yellow;
draw (        0.0000000001u,       -0.0000000001u)--(       -0.3645000000u,        0.3423000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.3645000000u,        0.3423000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.3645000000u,        0.3423000000u) withcolor Yellow;
draw (        0.0000000001u,       -0.0000000001u)--(       -0.4045000000u,        0.2939000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.4045000000u,        0.2939000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.4045000000u,        0.2939000000u) withcolor Yellow;
draw (        0.0000000001u,       -0.0000000001u)--(       -0.4382000000u,        0.2409000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.4382000000u,        0.2409000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.4382000000u,        0.2409000000u) withcolor Yellow;
draw (        0.0000000001u,       -0.0000000001u)--(       -0.4649000000u,        0.1841000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.4649000000u,        0.1841000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.4649000000u,        0.1841000000u) withcolor Yellow;
draw (        0.0000000001u,       -0.0000000001u)--(       -0.4843000000u,        0.1243000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.4843000000u,        0.1243000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.4843000000u,        0.1243000000u) withcolor Yellow;
draw (        0.0000000001u,       -0.0000000001u)--(       -0.4961000000u,        0.0627000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.4961000000u,        0.0627000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.4961000000u,        0.0627000000u) withcolor Yellow;
draw (        0.0000000001u,       -0.0000000001u)--(       -0.5000000000u,       -0.0000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.5000000000u,       -0.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.5000000000u,       -0.0000000000u) withcolor Yellow;
draw (        0.0000000001u,       -0.0000000001u)--(       -0.4961000000u,       -0.0627000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.4961000000u,       -0.0627000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.4961000000u,       -0.0627000000u) withcolor Yellow;
draw (        0.0000000001u,       -0.0000000001u)--(       -0.4843000000u,       -0.1243000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.4843000000u,       -0.1243000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.4843000000u,       -0.1243000000u) withcolor Yellow;
draw (        0.0000000001u,       -0.0000000001u)--(       -0.4649000000u,       -0.1841000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.4649000000u,       -0.1841000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.4649000000u,       -0.1841000000u) withcolor Yellow;
draw (        0.0000000001u,       -0.0000000001u)--(       -0.4382000000u,       -0.2409000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.4382000000u,       -0.2409000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.4382000000u,       -0.2409000000u) withcolor Yellow;
draw (        0.0000000001u,       -0.0000000001u)--(       -0.4045000000u,       -0.2939000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.4045000000u,       -0.2939000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.4045000000u,       -0.2939000000u) withcolor Yellow;
draw (        0.0000000001u,       -0.0000000001u)--(       -0.3645000000u,       -0.3423000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.3645000000u,       -0.3423000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.3645000000u,       -0.3423000000u) withcolor Yellow;
draw (        0.0000000001u,       -0.0000000001u)--(       -0.3187000000u,       -0.3853000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.3187000000u,       -0.3853000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.3187000000u,       -0.3853000000u) withcolor Yellow;
draw (        0.0000000001u,       -0.0000000001u)--(       -0.2679000000u,       -0.4222000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.2679000000u,       -0.4222000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.2679000000u,       -0.4222000000u) withcolor Yellow;
draw (        0.0000000001u,       -0.0000000001u)--(       -0.2129000000u,       -0.4524000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.2129000000u,       -0.4524000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.2129000000u,       -0.4524000000u) withcolor Yellow;
draw (        0.0000000001u,       -0.0000000001u)--(       -0.1545000000u,       -0.4755000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.1545000000u,       -0.4755000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1545000000u,       -0.4755000000u) withcolor Yellow;
draw (        0.0000000001u,       -0.0000000001u)--(       -0.0937000000u,       -0.4911000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0937000000u,       -0.4911000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0937000000u,       -0.4911000000u) withcolor Yellow;
% r =         0.5000408314
endfig;
end;
