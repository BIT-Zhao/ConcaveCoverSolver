prologues := 3;
outputformat   := "svg";
outputtemplate := "%j.%o";
input mpcolornames;
beginfig(7);
u := 20cm;
path pols[];
pair p[];
R := .2496308138u;
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

x1 = .2891270550u;
y1 = .3213515648u;
draw fullcircle scaled D shifted z1 withcolor RoyalBlue withpen pencircle scaled 1.0;

x2 = -.2891270549u;
y2 = -.3213515648u;
draw fullcircle scaled D shifted z2 withcolor RoyalBlue withpen pencircle scaled 1.0;

x3 = -.4228363987u;
y3 = .0899996148u;
draw fullcircle scaled D shifted z3 withcolor RoyalBlue withpen pencircle scaled 1.0;

x4 = .4228363986u;
y4 = -.0899996148u;
draw fullcircle scaled D shifted z4 withcolor RoyalBlue withpen pencircle scaled 1.0;

x5 = .1337093440u;
y5 = -.4113511807u;
draw fullcircle scaled D shifted z5 withcolor RoyalBlue withpen pencircle scaled 1.0;

x6 = -.1337093440u;
y6 = .4113511808u;
draw fullcircle scaled D shifted z6 withcolor RoyalBlue withpen pencircle scaled 1.0;

x7 = .0000000001u;
y7 = .0000000001u;
draw fullcircle scaled D shifted z7 withcolor RoyalBlue withpen pencircle scaled 1.0;
draw (        0.2891270550u,        0.3213515648u)--(        0.2374144407u,        0.0771357816u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4228363986u,       -0.0899996148u)--(        0.2374144407u,        0.0771357816u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2891270550u,        0.3213515648u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2891270550u,        0.3213515648u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4228363986u,       -0.0899996148u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4228363986u,       -0.0899996148u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2374144407u,        0.0771357816u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2374144407u,        0.0771357816u) withcolor MidnightBlue;
draw (        0.2891270550u,        0.3213515648u)--(        0.0517126174u,        0.2442157883u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.1337093440u,        0.4113511808u)--(        0.0517126174u,        0.2442157883u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2891270550u,        0.3213515648u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2891270550u,        0.3213515648u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.1337093440u,        0.4113511808u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1337093440u,        0.4113511808u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0517126174u,        0.2442157883u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0517126174u,        0.2442157883u) withcolor MidnightBlue;
draw (        0.2891270550u,        0.3213515648u)--(        0.2374144357u,        0.0771357827u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0000000001u,        0.0000000001u)--(        0.2374144357u,        0.0771357827u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2891270550u,        0.3213515648u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2891270550u,        0.3213515648u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,        0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,        0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2374144357u,        0.0771357827u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2374144357u,        0.0771357827u) withcolor MidnightBlue;
draw (        0.2891270550u,        0.3213515648u)--(        0.0517126194u,        0.2442157821u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0000000001u,        0.0000000001u)--(        0.0517126194u,        0.2442157821u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2891270550u,        0.3213515648u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2891270550u,        0.3213515648u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,        0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,        0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0517126194u,        0.2442157821u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0517126194u,        0.2442157821u) withcolor MidnightBlue;
draw (       -0.2891270549u,       -0.3213515648u)--(       -0.2374144408u,       -0.0771357817u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.4228363987u,        0.0899996148u)--(       -0.2374144408u,       -0.0771357817u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.2891270549u,       -0.3213515648u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.2891270549u,       -0.3213515648u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.4228363987u,        0.0899996148u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.4228363987u,        0.0899996148u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.2374144408u,       -0.0771357817u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.2374144408u,       -0.0771357817u) withcolor MidnightBlue;
draw (       -0.2891270549u,       -0.3213515648u)--(       -0.0517126174u,       -0.2442157881u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1337093440u,       -0.4113511807u)--(       -0.0517126174u,       -0.2442157881u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.2891270549u,       -0.3213515648u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.2891270549u,       -0.3213515648u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1337093440u,       -0.4113511807u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1337093440u,       -0.4113511807u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0517126174u,       -0.2442157881u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0517126174u,       -0.2442157881u) withcolor MidnightBlue;
draw (       -0.2891270549u,       -0.3213515648u)--(       -0.2374144355u,       -0.0771357828u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0000000001u,        0.0000000001u)--(       -0.2374144355u,       -0.0771357828u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.2891270549u,       -0.3213515648u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.2891270549u,       -0.3213515648u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,        0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,        0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.2374144355u,       -0.0771357828u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.2374144355u,       -0.0771357828u) withcolor MidnightBlue;
draw (       -0.2891270549u,       -0.3213515648u)--(       -0.0517126193u,       -0.2442157820u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0000000001u,        0.0000000001u)--(       -0.0517126193u,       -0.2442157820u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.2891270549u,       -0.3213515648u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.2891270549u,       -0.3213515648u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,        0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,        0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0517126193u,       -0.2442157820u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0517126193u,       -0.2442157820u) withcolor MidnightBlue;
draw (       -0.4228363987u,        0.0899996148u)--(       -0.1854219635u,        0.1671353988u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.1337093440u,        0.4113511808u)--(       -0.1854219635u,        0.1671353988u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.4228363987u,        0.0899996148u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.4228363987u,        0.0899996148u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.1337093440u,        0.4113511808u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1337093440u,        0.4113511808u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.1854219635u,        0.1671353988u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1854219635u,        0.1671353988u) withcolor MidnightBlue;
draw (       -0.4228363987u,        0.0899996148u)--(       -0.1854219614u,        0.1671353925u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0000000001u,        0.0000000001u)--(       -0.1854219614u,        0.1671353925u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.4228363987u,        0.0899996148u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.4228363987u,        0.0899996148u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,        0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,        0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.1854219614u,        0.1671353925u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1854219614u,        0.1671353925u) withcolor MidnightBlue;
draw (       -0.4228363987u,        0.0899996148u)--(       -0.2374144372u,       -0.0771357777u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0000000001u,        0.0000000001u)--(       -0.2374144372u,       -0.0771357777u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.4228363987u,        0.0899996148u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.4228363987u,        0.0899996148u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,        0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,        0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.2374144372u,       -0.0771357777u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.2374144372u,       -0.0771357777u) withcolor MidnightBlue;
draw (        0.4228363986u,       -0.0899996148u)--(        0.1854219634u,       -0.1671353986u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1337093440u,       -0.4113511807u)--(        0.1854219634u,       -0.1671353986u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4228363986u,       -0.0899996148u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4228363986u,       -0.0899996148u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1337093440u,       -0.4113511807u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1337093440u,       -0.4113511807u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1854219634u,       -0.1671353986u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1854219634u,       -0.1671353986u) withcolor MidnightBlue;
draw (        0.4228363986u,       -0.0899996148u)--(        0.1854219614u,       -0.1671353926u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0000000001u,        0.0000000001u)--(        0.1854219614u,       -0.1671353926u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4228363986u,       -0.0899996148u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4228363986u,       -0.0899996148u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,        0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,        0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1854219614u,       -0.1671353926u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1854219614u,       -0.1671353926u) withcolor MidnightBlue;
draw (        0.4228363986u,       -0.0899996148u)--(        0.2374144373u,        0.0771357778u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0000000001u,        0.0000000001u)--(        0.2374144373u,        0.0771357778u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4228363986u,       -0.0899996148u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4228363986u,       -0.0899996148u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,        0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,        0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2374144373u,        0.0771357778u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2374144373u,        0.0771357778u) withcolor MidnightBlue;
draw (        0.1337093440u,       -0.4113511807u)--(       -0.0517126131u,       -0.2442157833u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0000000001u,        0.0000000001u)--(       -0.0517126131u,       -0.2442157833u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1337093440u,       -0.4113511807u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1337093440u,       -0.4113511807u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,        0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,        0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0517126131u,       -0.2442157833u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0517126131u,       -0.2442157833u) withcolor MidnightBlue;
draw (        0.1337093440u,       -0.4113511807u)--(        0.1854219572u,       -0.1671353973u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0000000001u,        0.0000000001u)--(        0.1854219572u,       -0.1671353973u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1337093440u,       -0.4113511807u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1337093440u,       -0.4113511807u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,        0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,        0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1854219572u,       -0.1671353973u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1854219572u,       -0.1671353973u) withcolor MidnightBlue;
draw (       -0.1337093440u,        0.4113511808u)--(        0.0517126130u,        0.2442157835u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0000000001u,        0.0000000001u)--(        0.0517126130u,        0.2442157835u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.1337093440u,        0.4113511808u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1337093440u,        0.4113511808u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,        0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,        0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0517126130u,        0.2442157835u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0517126130u,        0.2442157835u) withcolor MidnightBlue;
draw (       -0.1337093440u,        0.4113511808u)--(       -0.1854219570u,        0.1671353974u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0000000001u,        0.0000000001u)--(       -0.1854219570u,        0.1671353974u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.1337093440u,        0.4113511808u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1337093440u,        0.4113511808u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,        0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,        0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.1854219570u,        0.1671353974u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1854219570u,        0.1671353974u) withcolor MidnightBlue;
% r =         0.2496308138
endfig;
end;
