prologues := 3;
outputformat   := "svg";
outputtemplate := "%j.%o";
input mpcolornames;
beginfig(8);
u := 20cm;
path pols[];
pair p[];
R := .2222096427u;
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

x1 = -.0148493382u;
y1 = -.3995725322u;
draw fullcircle scaled D shifted z1 withcolor RoyalBlue withpen pencircle scaled 1.0;

x2 = .3047471603u;
y2 = -.2618258475u;
draw fullcircle scaled D shifted z2 withcolor RoyalBlue withpen pencircle scaled 1.0;

x3 = -.1610712004u;
y3 = .3660115443u;
draw fullcircle scaled D shifted z3 withcolor RoyalBlue withpen pencircle scaled 1.0;

x4 = .3956253256u;
y4 = .0755392888u;
draw fullcircle scaled D shifted z4 withcolor RoyalBlue withpen pencircle scaled 1.0;

x5 = -.3202808067u;
y5 = -.2368702462u;
draw fullcircle scaled D shifted z5 withcolor RoyalBlue withpen pencircle scaled 1.0;

x6 = .0030826349u;
y6 = .0006174802u;
draw fullcircle scaled D shifted z6 withcolor RoyalBlue withpen pencircle scaled 1.0;

x7 = -.3850538302u;
y7 = .1021713466u;
draw fullcircle scaled D shifted z7 withcolor RoyalBlue withpen pencircle scaled 1.0;

x8 = .1868431572u;
y8 = .3557186063u;
draw fullcircle scaled D shifted z8 withcolor RoyalBlue withpen pencircle scaled 1.0;
draw (       -0.0148493382u,       -0.3995725322u)--(        0.0902504087u,       -0.2037890587u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3047471603u,       -0.2618258475u)--(        0.0902504087u,       -0.2037890587u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.0148493382u,       -0.3995725322u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0148493382u,       -0.3995725322u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3047471603u,       -0.2618258475u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3047471603u,       -0.2618258475u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0902504087u,       -0.2037890587u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0902504087u,       -0.2037890587u) withcolor MidnightBlue;
draw (       -0.0148493382u,       -0.3995725322u)--(       -0.1020182167u,       -0.1951741176u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.3202808067u,       -0.2368702462u)--(       -0.1020182167u,       -0.1951741176u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.0148493382u,       -0.3995725322u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0148493382u,       -0.3995725322u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.3202808067u,       -0.2368702462u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.3202808067u,       -0.2368702462u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.1020182167u,       -0.1951741176u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1020182167u,       -0.1951741176u) withcolor MidnightBlue;
draw (       -0.0148493382u,       -0.3995725322u)--(       -0.1020091667u,       -0.1951702583u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0030826349u,        0.0006174802u)--(       -0.1020091667u,       -0.1951702583u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.0148493382u,       -0.3995725322u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0148493382u,       -0.3995725322u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0030826349u,        0.0006174802u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0030826349u,        0.0006174802u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.1020091667u,       -0.1951702583u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1020091667u,       -0.1951702583u) withcolor MidnightBlue;
draw (       -0.0148493382u,       -0.3995725322u)--(        0.0902424634u,       -0.2037847938u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0030826349u,        0.0006174802u)--(        0.0902424634u,       -0.2037847938u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.0148493382u,       -0.3995725322u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0148493382u,       -0.3995725322u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0030826349u,        0.0006174802u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0030826349u,        0.0006174802u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0902424634u,       -0.2037847938u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0902424634u,       -0.2037847938u) withcolor MidnightBlue;
draw (        0.3047471603u,       -0.2618258475u)--(        0.2175863194u,       -0.0574240053u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3956253256u,        0.0755392888u)--(        0.2175863194u,       -0.0574240053u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3047471603u,       -0.2618258475u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3047471603u,       -0.2618258475u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3956253256u,        0.0755392888u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3956253256u,        0.0755392888u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2175863194u,       -0.0574240053u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2175863194u,       -0.0574240053u) withcolor MidnightBlue;
draw (        0.3047471603u,       -0.2618258475u)--(        0.0902527497u,       -0.2037804072u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0030826349u,        0.0006174802u)--(        0.0902527497u,       -0.2037804072u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3047471603u,       -0.2618258475u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3047471603u,       -0.2618258475u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0030826349u,        0.0006174802u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0030826349u,        0.0006174802u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0902527497u,       -0.2037804072u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0902527497u,       -0.2037804072u) withcolor MidnightBlue;
draw (        0.3047471603u,       -0.2618258475u)--(        0.2175770454u,       -0.0574279601u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0030826349u,        0.0006174802u)--(        0.2175770454u,       -0.0574279601u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3047471603u,       -0.2618258475u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3047471603u,       -0.2618258475u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0030826349u,        0.0006174802u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0030826349u,        0.0006174802u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2175770454u,       -0.0574279601u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2175770454u,       -0.0574279601u) withcolor MidnightBlue;
draw (       -0.1610712004u,        0.3660115443u)--(        0.0087943820u,        0.2227537025u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0030826349u,        0.0006174802u)--(        0.0087943820u,        0.2227537025u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.1610712004u,        0.3660115443u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1610712004u,        0.3660115443u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0030826349u,        0.0006174802u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0030826349u,        0.0006174802u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0087943820u,        0.2227537025u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0087943820u,        0.2227537025u) withcolor MidnightBlue;
draw (       -0.1610712004u,        0.3660115443u)--(       -0.1667829475u,        0.1438753220u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0030826349u,        0.0006174802u)--(       -0.1667829475u,        0.1438753220u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.1610712004u,        0.3660115443u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1610712004u,        0.3660115443u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0030826349u,        0.0006174802u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0030826349u,        0.0006174802u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.1667829475u,        0.1438753220u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1667829475u,        0.1438753220u) withcolor MidnightBlue;
draw (       -0.1610712004u,        0.3660115443u)--(       -0.1667927878u,        0.1438755752u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.3850538302u,        0.1021713466u)--(       -0.1667927878u,        0.1438755752u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.1610712004u,        0.3660115443u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1610712004u,        0.3660115443u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.3850538302u,        0.1021713466u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.3850538302u,        0.1021713466u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.1667927878u,        0.1438755752u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1667927878u,        0.1438755752u) withcolor MidnightBlue;
draw (       -0.1610712004u,        0.3660115443u)--(        0.0088002006u,        0.2227606021u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1868431572u,        0.3557186063u)--(        0.0088002006u,        0.2227606021u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.1610712004u,        0.3660115443u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1610712004u,        0.3660115443u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1868431572u,        0.3557186063u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1868431572u,        0.3557186063u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0088002006u,        0.2227606021u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0088002006u,        0.2227606021u) withcolor MidnightBlue;
draw (        0.3956253256u,        0.0755392888u)--(        0.2175802975u,       -0.0574159415u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0030826349u,        0.0006174802u)--(        0.2175802975u,       -0.0574159415u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3956253256u,        0.0755392888u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3956253256u,        0.0755392888u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0030826349u,        0.0006174802u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0030826349u,        0.0006174802u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2175802975u,       -0.0574159415u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2175802975u,       -0.0574159415u) withcolor MidnightBlue;
draw (        0.3956253256u,        0.0755392888u)--(        0.1811276630u,        0.1335727105u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0030826349u,        0.0006174802u)--(        0.1811276630u,        0.1335727105u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3956253256u,        0.0755392888u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3956253256u,        0.0755392888u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0030826349u,        0.0006174802u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0030826349u,        0.0006174802u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1811276630u,        0.1335727105u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1811276630u,        0.1335727105u) withcolor MidnightBlue;
draw (        0.3956253256u,        0.0755392888u)--(        0.1811302883u,        0.1335824129u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1868431572u,        0.3557186063u)--(        0.1811302883u,        0.1335824129u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3956253256u,        0.0755392888u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3956253256u,        0.0755392888u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1868431572u,        0.3557186063u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1868431572u,        0.3557186063u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1811302883u,        0.1335824129u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1811302883u,        0.1335824129u) withcolor MidnightBlue;
draw (       -0.3202808067u,       -0.2368702462u)--(       -0.2151780992u,       -0.0410883620u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0030826349u,        0.0006174802u)--(       -0.2151780992u,       -0.0410883620u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.3202808067u,       -0.2368702462u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.3202808067u,       -0.2368702462u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0030826349u,        0.0006174802u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0030826349u,        0.0006174802u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.2151780992u,       -0.0410883620u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.2151780992u,       -0.0410883620u) withcolor MidnightBlue;
draw (       -0.3202808067u,       -0.2368702462u)--(       -0.1020200726u,       -0.1951644040u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0030826349u,        0.0006174802u)--(       -0.1020200726u,       -0.1951644040u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.3202808067u,       -0.2368702462u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.3202808067u,       -0.2368702462u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0030826349u,        0.0006174802u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0030826349u,        0.0006174802u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.1020200726u,       -0.1951644040u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1020200726u,       -0.1951644040u) withcolor MidnightBlue;
draw (       -0.3202808067u,       -0.2368702462u)--(       -0.2151861689u,       -0.0410840301u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.3850538302u,        0.1021713466u)--(       -0.2151861689u,       -0.0410840301u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.3202808067u,       -0.2368702462u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.3202808067u,       -0.2368702462u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.3850538302u,        0.1021713466u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.3850538302u,        0.1021713466u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.2151861689u,       -0.0410840301u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.2151861689u,       -0.0410840301u) withcolor MidnightBlue;
draw (        0.0030826349u,        0.0006174802u)--(       -0.2151802645u,       -0.0410770285u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.3850538302u,        0.1021713466u)--(       -0.2151802645u,       -0.0410770285u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0030826349u,        0.0006174802u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0030826349u,        0.0006174802u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.3850538302u,        0.1021713466u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.3850538302u,        0.1021713466u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.2151802645u,       -0.0410770285u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.2151802645u,       -0.0410770285u) withcolor MidnightBlue;
draw (        0.0030826349u,        0.0006174802u)--(       -0.1667909308u,        0.1438658553u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.3850538302u,        0.1021713466u)--(       -0.1667909308u,        0.1438658553u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0030826349u,        0.0006174802u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0030826349u,        0.0006174802u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.3850538302u,        0.1021713466u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.3850538302u,        0.1021713466u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.1667909308u,        0.1438658553u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1667909308u,        0.1438658553u) withcolor MidnightBlue;
draw (        0.0030826349u,        0.0006174802u)--(        0.0088055682u,        0.2227534146u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1868431572u,        0.3557186063u)--(        0.0088055682u,        0.2227534146u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0030826349u,        0.0006174802u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0030826349u,        0.0006174802u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1868431572u,        0.3557186063u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1868431572u,        0.3557186063u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0088055682u,        0.2227534146u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0088055682u,        0.2227534146u) withcolor MidnightBlue;
draw (        0.0030826349u,        0.0006174802u)--(        0.1811202238u,        0.1335826720u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1868431572u,        0.3557186063u)--(        0.1811202238u,        0.1335826720u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0030826349u,        0.0006174802u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0030826349u,        0.0006174802u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1868431572u,        0.3557186063u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1868431572u,        0.3557186063u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1811202238u,        0.1335826720u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1811202238u,        0.1335826720u) withcolor MidnightBlue;
% r =         0.2222096427
endfig;
end;
