prologues := 3;
outputformat   := "svg";
outputtemplate := "%j.%o";
input mpcolornames;
beginfig(9);
u := 20cm;
path pols[];
pair p[];
R := .2068181445u;
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

x1 = .3228185987u;
y1 = .2051158071u;
draw fullcircle scaled D shifted z1 withcolor RoyalBlue withpen pencircle scaled 1.0;

x2 = .0827873482u;
y2 = .3726340955u;
draw fullcircle scaled D shifted z2 withcolor RoyalBlue withpen pencircle scaled 1.0;

x3 = .3725132248u;
y3 = -.0833465119u;
draw fullcircle scaled D shifted z3 withcolor RoyalBlue withpen pencircle scaled 1.0;

x4 = -.3725132248u;
y4 = .0833465119u;
draw fullcircle scaled D shifted z4 withcolor RoyalBlue withpen pencircle scaled 1.0;

x5 = .2052505981u;
y5 = -.3230139772u;
draw fullcircle scaled D shifted z5 withcolor RoyalBlue withpen pencircle scaled 1.0;

x6 = .0000000000u;
y6 = .0000000001u;
draw fullcircle scaled D shifted z6 withcolor RoyalBlue withpen pencircle scaled 1.0;

x7 = -.3228185986u;
y7 = -.2051158072u;
draw fullcircle scaled D shifted z7 withcolor RoyalBlue withpen pencircle scaled 1.0;

x8 = -.0827873481u;
y8 = -.3726340956u;
draw fullcircle scaled D shifted z8 withcolor RoyalBlue withpen pencircle scaled 1.0;

x9 = -.2052505982u;
y9 = .3230139773u;
draw fullcircle scaled D shifted z9 withcolor RoyalBlue withpen pencircle scaled 1.0;
draw (        0.3228185987u,        0.2051158071u)--(        0.1191709060u,        0.1690414100u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0827873482u,        0.3726340955u)--(        0.1191709060u,        0.1690414100u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3228185987u,        0.2051158071u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3228185987u,        0.2051158071u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0827873482u,        0.3726340955u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0827873482u,        0.3726340955u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1191709060u,        0.1690414100u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1191709060u,        0.1690414100u) withcolor MidnightBlue;
draw (        0.3228185987u,        0.2051158071u)--(        0.2036580989u,        0.0360758114u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3725132248u,       -0.0833465119u)--(        0.2036580989u,        0.0360758114u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3228185987u,        0.2051158071u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3228185987u,        0.2051158071u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3725132248u,       -0.0833465119u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3725132248u,       -0.0833465119u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2036580989u,        0.0360758114u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2036580989u,        0.0360758114u) withcolor MidnightBlue;
draw (        0.3228185987u,        0.2051158071u)--(        0.2036459212u,        0.0360843964u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0000000000u,        0.0000000001u)--(        0.2036459212u,        0.0360843964u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3228185987u,        0.2051158071u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3228185987u,        0.2051158071u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2036459212u,        0.0360843964u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2036459212u,        0.0360843964u) withcolor MidnightBlue;
draw (        0.3228185987u,        0.2051158071u)--(        0.1191726775u,        0.1690314108u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0000000000u,        0.0000000001u)--(        0.1191726775u,        0.1690314108u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3228185987u,        0.2051158071u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3228185987u,        0.2051158071u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1191726775u,        0.1690314108u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1191726775u,        0.1690314108u) withcolor MidnightBlue;
draw (        0.0827873482u,        0.3726340955u)--(        0.1191610037u,        0.1690396406u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0000000000u,        0.0000000001u)--(        0.1191610037u,        0.1690396406u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0827873482u,        0.3726340955u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0827873482u,        0.3726340955u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1191610037u,        0.1690396406u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1191610037u,        0.1690396406u) withcolor MidnightBlue;
draw (        0.0827873482u,        0.3726340955u)--(       -0.0363736555u,        0.2035944550u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0000000000u,        0.0000000001u)--(       -0.0363736555u,        0.2035944550u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0827873482u,        0.3726340955u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0827873482u,        0.3726340955u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0363736555u,        0.2035944550u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0363736555u,        0.2035944550u) withcolor MidnightBlue;
draw (        0.0827873482u,        0.3726340955u)--(       -0.0363868921u,        0.2036037866u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.2052505982u,        0.3230139773u)--(       -0.0363868921u,        0.2036037866u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0827873482u,        0.3726340955u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0827873482u,        0.3726340955u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.2052505982u,        0.3230139773u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.2052505982u,        0.3230139773u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0363868921u,        0.2036037866u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0363868921u,        0.2036037866u) withcolor MidnightBlue;
draw (        0.3725132248u,       -0.0833465119u)--(        0.1688656441u,       -0.1194215412u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2052505981u,       -0.3230139772u)--(        0.1688656441u,       -0.1194215412u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3725132248u,       -0.0833465119u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3725132248u,       -0.0833465119u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2052505981u,       -0.3230139772u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2052505981u,       -0.3230139772u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1688656441u,       -0.1194215412u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1688656441u,       -0.1194215412u) withcolor MidnightBlue;
draw (        0.3725132248u,       -0.0833465119u)--(        0.1688636482u,       -0.1194102724u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0000000000u,        0.0000000001u)--(        0.1688636482u,       -0.1194102724u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3725132248u,       -0.0833465119u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3725132248u,       -0.0833465119u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1688636482u,       -0.1194102724u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1688636482u,       -0.1194102724u) withcolor MidnightBlue;
draw (        0.3725132248u,       -0.0833465119u)--(        0.2036495766u,        0.0360637606u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0000000000u,        0.0000000001u)--(        0.2036495766u,        0.0360637606u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3725132248u,       -0.0833465119u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3725132248u,       -0.0833465119u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2036495766u,        0.0360637606u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2036495766u,        0.0360637606u) withcolor MidnightBlue;
draw (       -0.3725132248u,        0.0833465119u)--(       -0.1688636483u,        0.1194102726u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0000000000u,        0.0000000001u)--(       -0.1688636483u,        0.1194102726u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.3725132248u,        0.0833465119u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.3725132248u,        0.0833465119u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.1688636483u,        0.1194102726u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1688636483u,        0.1194102726u) withcolor MidnightBlue;
draw (       -0.3725132248u,        0.0833465119u)--(       -0.2036495766u,       -0.0360637606u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0000000000u,        0.0000000001u)--(       -0.2036495766u,       -0.0360637606u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.3725132248u,        0.0833465119u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.3725132248u,        0.0833465119u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.2036495766u,       -0.0360637606u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.2036495766u,       -0.0360637606u) withcolor MidnightBlue;
draw (       -0.3725132248u,        0.0833465119u)--(       -0.2036580988u,       -0.0360758114u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.3228185986u,       -0.2051158072u)--(       -0.2036580988u,       -0.0360758114u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.3725132248u,        0.0833465119u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.3725132248u,        0.0833465119u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.3228185986u,       -0.2051158072u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.3228185986u,       -0.2051158072u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.2036580988u,       -0.0360758114u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.2036580988u,       -0.0360758114u) withcolor MidnightBlue;
draw (       -0.3725132248u,        0.0833465119u)--(       -0.1688656441u,        0.1194215413u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.2052505982u,        0.3230139773u)--(       -0.1688656441u,        0.1194215413u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.3725132248u,        0.0833465119u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.3725132248u,        0.0833465119u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.2052505982u,        0.3230139773u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.2052505982u,        0.3230139773u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.1688656441u,        0.1194215413u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1688656441u,        0.1194215413u) withcolor MidnightBlue;
draw (        0.2052505981u,       -0.3230139772u)--(        0.0363963611u,       -0.2035903970u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0000000000u,        0.0000000001u)--(        0.0363963611u,       -0.2035903970u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2052505981u,       -0.3230139772u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2052505981u,       -0.3230139772u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0363963611u,       -0.2035903970u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0363963611u,       -0.2035903970u) withcolor MidnightBlue;
draw (        0.2052505981u,       -0.3230139772u)--(        0.1688542370u,       -0.1194235801u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0000000000u,        0.0000000001u)--(        0.1688542370u,       -0.1194235801u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2052505981u,       -0.3230139772u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2052505981u,       -0.3230139772u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1688542370u,       -0.1194235801u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1688542370u,       -0.1194235801u) withcolor MidnightBlue;
draw (        0.2052505981u,       -0.3230139772u)--(        0.0363868920u,       -0.2036037865u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.0827873481u,       -0.3726340956u)--(        0.0363868920u,       -0.2036037865u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2052505981u,       -0.3230139772u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2052505981u,       -0.3230139772u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0827873481u,       -0.3726340956u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0827873481u,       -0.3726340956u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0363868920u,       -0.2036037865u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0363868920u,       -0.2036037865u) withcolor MidnightBlue;
draw (        0.0000000000u,        0.0000000001u)--(       -0.1191726774u,       -0.1690314107u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.3228185986u,       -0.2051158072u)--(       -0.1191726774u,       -0.1690314107u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.3228185986u,       -0.2051158072u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.3228185986u,       -0.2051158072u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.1191726774u,       -0.1690314107u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1191726774u,       -0.1690314107u) withcolor MidnightBlue;
draw (        0.0000000000u,        0.0000000001u)--(       -0.2036459212u,       -0.0360843964u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.3228185986u,       -0.2051158072u)--(       -0.2036459212u,       -0.0360843964u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.3228185986u,       -0.2051158072u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.3228185986u,       -0.2051158072u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.2036459212u,       -0.0360843964u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.2036459212u,       -0.0360843964u) withcolor MidnightBlue;
draw (        0.0000000000u,        0.0000000001u)--(        0.0363736553u,       -0.2035944548u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.0827873481u,       -0.3726340956u)--(        0.0363736553u,       -0.2035944548u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0827873481u,       -0.3726340956u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0827873481u,       -0.3726340956u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0363736553u,       -0.2035944548u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0363736553u,       -0.2035944548u) withcolor MidnightBlue;
draw (        0.0000000000u,        0.0000000001u)--(       -0.1191610034u,       -0.1690396406u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.0827873481u,       -0.3726340956u)--(       -0.1191610034u,       -0.1690396406u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0827873481u,       -0.3726340956u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0827873481u,       -0.3726340956u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.1191610034u,       -0.1690396406u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1191610034u,       -0.1690396406u) withcolor MidnightBlue;
draw (        0.0000000000u,        0.0000000001u)--(       -0.1688542371u,        0.1194235803u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.2052505982u,        0.3230139773u)--(       -0.1688542371u,        0.1194235803u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.2052505982u,        0.3230139773u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.2052505982u,        0.3230139773u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.1688542371u,        0.1194235803u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1688542371u,        0.1194235803u) withcolor MidnightBlue;
draw (        0.0000000000u,        0.0000000001u)--(       -0.0363963611u,        0.2035903971u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.2052505982u,        0.3230139773u)--(       -0.0363963611u,        0.2035903971u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.2052505982u,        0.3230139773u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.2052505982u,        0.3230139773u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0363963611u,        0.2035903971u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0363963611u,        0.2035903971u) withcolor MidnightBlue;
draw (       -0.3228185986u,       -0.2051158072u)--(       -0.1191709059u,       -0.1690414100u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.0827873481u,       -0.3726340956u)--(       -0.1191709059u,       -0.1690414100u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.3228185986u,       -0.2051158072u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.3228185986u,       -0.2051158072u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0827873481u,       -0.3726340956u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0827873481u,       -0.3726340956u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.1191709059u,       -0.1690414100u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1191709059u,       -0.1690414100u) withcolor MidnightBlue;
% r =         0.2068181445
endfig;
end;
