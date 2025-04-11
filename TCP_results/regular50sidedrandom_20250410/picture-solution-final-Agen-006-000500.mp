prologues := 3;
outputformat   := "svg";
outputtemplate := "%j.%o";
input mpcolornames;
beginfig(6);
u := 20cm;
path pols[];
pair p[];
R := .2777406691u;
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

x1 = -.3864181453u;
y1 = .1507275470u;
draw fullcircle scaled D shifted z1 withcolor RoyalBlue withpen pencircle scaled 1.0;

x2 = .3864181454u;
y2 = -.1507275470u;
draw fullcircle scaled D shifted z2 withcolor RoyalBlue withpen pencircle scaled 1.0;

x3 = -.2379881319u;
y3 = -.2302122804u;
draw fullcircle scaled D shifted z3 withcolor RoyalBlue withpen pencircle scaled 1.0;

x4 = .0480942575u;
y4 = -.2651609585u;
draw fullcircle scaled D shifted z4 withcolor RoyalBlue withpen pencircle scaled 1.0;

x5 = -.0480942575u;
y5 = .2651609586u;
draw fullcircle scaled D shifted z5 withcolor RoyalBlue withpen pencircle scaled 1.0;

x6 = .2379881319u;
y6 = .2302122805u;
draw fullcircle scaled D shifted z6 withcolor RoyalBlue withpen pencircle scaled 1.0;
draw (       -0.2379881319u,       -0.2302122804u)--(       -0.0661556323u,       -0.0120071674u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0480942575u,       -0.2651609585u)--(       -0.0661556323u,       -0.0120071674u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.2379881319u,       -0.2302122804u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.2379881319u,       -0.2302122804u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0480942575u,       -0.2651609585u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0480942575u,       -0.2651609585u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0661556323u,       -0.0120071674u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0661556323u,       -0.0120071674u) withcolor MidnightBlue;
draw (       -0.2379881319u,       -0.2302122804u)--(       -0.0661768749u,       -0.0119904409u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.0480942575u,        0.2651609586u)--(       -0.0661768749u,       -0.0119904409u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.2379881319u,       -0.2302122804u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.2379881319u,       -0.2302122804u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0480942575u,        0.2651609586u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0480942575u,        0.2651609586u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0661768749u,       -0.0119904409u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0661768749u,       -0.0119904409u) withcolor MidnightBlue;
draw (        0.0480942575u,       -0.2651609585u)--(       -0.0661259211u,       -0.0119937606u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.0480942575u,        0.2651609586u)--(       -0.0661259211u,       -0.0119937606u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0480942575u,       -0.2651609585u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0480942575u,       -0.2651609585u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0480942575u,        0.2651609586u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0480942575u,        0.2651609586u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0661259211u,       -0.0119937606u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0661259211u,       -0.0119937606u) withcolor MidnightBlue;
draw (        0.0480942575u,       -0.2651609585u)--(        0.0661259211u,        0.0119937607u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.0480942575u,        0.2651609586u)--(        0.0661259211u,        0.0119937607u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0480942575u,       -0.2651609585u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0480942575u,       -0.2651609585u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0480942575u,        0.2651609586u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0480942575u,        0.2651609586u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0661259211u,        0.0119937607u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0661259211u,        0.0119937607u) withcolor MidnightBlue;
draw (        0.0480942575u,       -0.2651609585u)--(        0.0661768749u,        0.0119904410u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2379881319u,        0.2302122805u)--(        0.0661768749u,        0.0119904410u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0480942575u,       -0.2651609585u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0480942575u,       -0.2651609585u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2379881319u,        0.2302122805u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2379881319u,        0.2302122805u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0661768749u,        0.0119904410u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0661768749u,        0.0119904410u) withcolor MidnightBlue;
draw (       -0.0480942575u,        0.2651609586u)--(        0.0661556323u,        0.0120071674u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2379881319u,        0.2302122805u)--(        0.0661556323u,        0.0120071674u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.0480942575u,        0.2651609586u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0480942575u,        0.2651609586u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2379881319u,        0.2302122805u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2379881319u,        0.2302122805u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0661556323u,        0.0120071674u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0661556323u,        0.0120071674u) withcolor MidnightBlue;
% r =         0.2777406691
endfig;
end;
