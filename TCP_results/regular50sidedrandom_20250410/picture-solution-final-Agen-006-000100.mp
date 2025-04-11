prologues := 3;
outputformat   := "svg";
outputtemplate := "%j.%o";
input mpcolornames;
beginfig(6);
u := 20cm;
path pols[];
pair p[];
R := .2777429041u;
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

x1 = -.0194019836u;
y1 = -.2687803285u;
draw fullcircle scaled D shifted z1 withcolor RoyalBlue withpen pencircle scaled 1.0;

x2 = .0194019837u;
y2 = .2687803285u;
draw fullcircle scaled D shifted z2 withcolor RoyalBlue withpen pencircle scaled 1.0;

x3 = -.3367953289u;
y3 = .2420386949u;
draw fullcircle scaled D shifted z3 withcolor RoyalBlue withpen pencircle scaled 1.0;

x4 = .3367953288u;
y4 = -.2420386949u;
draw fullcircle scaled D shifted z4 withcolor RoyalBlue withpen pencircle scaled 1.0;

x5 = -.2877573159u;
y5 = -.1638424354u;
draw fullcircle scaled D shifted z5 withcolor RoyalBlue withpen pencircle scaled 1.0;

x6 = .2877573159u;
y6 = .1638424356u;
draw fullcircle scaled D shifted z6 withcolor RoyalBlue withpen pencircle scaled 1.0;
draw (       -0.0194019836u,       -0.2687803285u)--(       -0.0670699564u,        0.0048414636u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0194019837u,        0.2687803285u)--(       -0.0670699564u,        0.0048414636u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.0194019836u,       -0.2687803285u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0194019836u,       -0.2687803285u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0194019837u,        0.2687803285u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0194019837u,        0.2687803285u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0670699564u,        0.0048414636u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0670699564u,        0.0048414636u) withcolor MidnightBlue;
draw (       -0.0194019836u,       -0.2687803285u)--(        0.0670699565u,       -0.0048414637u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0194019837u,        0.2687803285u)--(        0.0670699565u,       -0.0048414637u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.0194019836u,       -0.2687803285u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0194019836u,       -0.2687803285u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0194019837u,        0.2687803285u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0194019837u,        0.2687803285u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0670699565u,       -0.0048414637u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0670699565u,       -0.0048414637u) withcolor MidnightBlue;
draw (       -0.0194019836u,       -0.2687803285u)--(       -0.0671020431u,        0.0048358718u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.2877573159u,       -0.1638424354u)--(       -0.0671020431u,        0.0048358718u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.0194019836u,       -0.2687803285u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0194019836u,       -0.2687803285u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.2877573159u,       -0.1638424354u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.2877573159u,       -0.1638424354u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0671020431u,        0.0048358718u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0671020431u,        0.0048358718u) withcolor MidnightBlue;
draw (       -0.0194019836u,       -0.2687803285u)--(        0.0671184744u,       -0.0048573641u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2877573159u,        0.1638424356u)--(        0.0671184744u,       -0.0048573641u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.0194019836u,       -0.2687803285u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0194019836u,       -0.2687803285u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2877573159u,        0.1638424356u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2877573159u,        0.1638424356u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0671184744u,       -0.0048573641u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0671184744u,       -0.0048573641u) withcolor MidnightBlue;
draw (        0.0194019837u,        0.2687803285u)--(       -0.0671184743u,        0.0048573641u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.2877573159u,       -0.1638424354u)--(       -0.0671184743u,        0.0048573641u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0194019837u,        0.2687803285u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0194019837u,        0.2687803285u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.2877573159u,       -0.1638424354u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.2877573159u,       -0.1638424354u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0671184743u,        0.0048573641u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0671184743u,        0.0048573641u) withcolor MidnightBlue;
draw (        0.0194019837u,        0.2687803285u)--(        0.0671020432u,       -0.0048358719u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2877573159u,        0.1638424356u)--(        0.0671020432u,       -0.0048358719u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0194019837u,        0.2687803285u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0194019837u,        0.2687803285u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2877573159u,        0.1638424356u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2877573159u,        0.1638424356u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0671020432u,       -0.0048358719u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0671020432u,       -0.0048358719u) withcolor MidnightBlue;
% r =         0.2777429041
endfig;
end;
