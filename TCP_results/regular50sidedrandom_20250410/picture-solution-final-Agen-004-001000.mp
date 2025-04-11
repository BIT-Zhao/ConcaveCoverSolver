prologues := 3;
outputformat   := "svg";
outputtemplate := "%j.%o";
input mpcolornames;
beginfig(4);
u := 20cm;
path pols[];
pair p[];
R := .3529842570u;
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

x1 = .1701603353u;
y1 = -.3092625172u;
draw fullcircle scaled D shifted z1 withcolor RoyalBlue withpen pencircle scaled 1.0;

x2 = -.3092771974u;
y2 = -.1701336383u;
draw fullcircle scaled D shifted z2 withcolor RoyalBlue withpen pencircle scaled 1.0;

x3 = .3092771974u;
y3 = .1701336383u;
draw fullcircle scaled D shifted z3 withcolor RoyalBlue withpen pencircle scaled 1.0;

x4 = -.1701603348u;
y4 = .3092625165u;
draw fullcircle scaled D shifted z4 withcolor RoyalBlue withpen pencircle scaled 1.0;
draw (        0.1701603353u,       -0.3092625172u)--(        0.0000000288u,        0.0000000828u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.3092771974u,       -0.1701336383u)--(        0.0000000288u,        0.0000000828u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1701603353u,       -0.3092625172u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1701603353u,       -0.3092625172u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.3092771974u,       -0.1701336383u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.3092771974u,       -0.1701336383u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000288u,        0.0000000828u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000288u,        0.0000000828u) withcolor MidnightBlue;
draw (        0.1701603353u,       -0.3092625172u)--(       -0.0000000853u,        0.0000000200u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3092771974u,        0.1701336383u)--(       -0.0000000853u,        0.0000000200u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1701603353u,       -0.3092625172u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1701603353u,       -0.3092625172u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3092771974u,        0.1701336383u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3092771974u,        0.1701336383u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0000000853u,        0.0000000200u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0000000853u,        0.0000000200u) withcolor MidnightBlue;
draw (       -0.3092771974u,       -0.1701336383u)--(        0.0000000857u,       -0.0000000208u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.1701603348u,        0.3092625165u)--(        0.0000000857u,       -0.0000000208u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.3092771974u,       -0.1701336383u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.3092771974u,       -0.1701336383u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.1701603348u,        0.3092625165u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1701603348u,        0.3092625165u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000857u,       -0.0000000208u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000857u,       -0.0000000208u) withcolor MidnightBlue;
draw (        0.3092771974u,        0.1701336383u)--(       -0.0000000284u,       -0.0000000836u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.1701603348u,        0.3092625165u)--(       -0.0000000284u,       -0.0000000836u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3092771974u,        0.1701336383u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3092771974u,        0.1701336383u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.1701603348u,        0.3092625165u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1701603348u,        0.3092625165u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0000000284u,       -0.0000000836u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0000000284u,       -0.0000000836u) withcolor MidnightBlue;
% r =         0.3529842570
endfig;
end;
