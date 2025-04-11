prologues := 3;
outputformat   := "svg";
outputtemplate := "%j.%o";
input mpcolornames;
beginfig(3);
u := 20cm;
path pols[];
pair p[];
R := .5038736962u;
D := 2*R;

pols[1] = 
(       -0.5000000000u,       -0.5000000000u)--
(        0.5000000000u,       -0.5000000000u)--
(        0.5000000000u,        0.5000000000u)--
(       -0.5000000000u,        0.5000000000u)--
(       -0.5000000000u,       -0.5000000000u)--
cycle;
fill pols[1] withcolor 0.9Dandelion;
draw (       -0.5000000000u,       -0.5000000000u)--(        0.5000000000u,       -0.5000000000u) withpen pencircle scaled 1.0;
draw (        0.5000000000u,       -0.5000000000u)--(        0.5000000000u,        0.5000000000u) withpen pencircle scaled 1.0;
draw (        0.5000000000u,        0.5000000000u)--(       -0.5000000000u,        0.5000000000u) withpen pencircle scaled 1.0;
draw (       -0.5000000000u,        0.5000000000u)--(       -0.5000000000u,       -0.5000000000u) withpen pencircle scaled 1.0;

x1 = -.2500093583u;
y1 = -.0625147074u;
draw fullcircle scaled D shifted z1 withcolor RoyalBlue withpen pencircle scaled 1.0;

x2 = .2500093584u;
y2 = -.0625147073u;
draw fullcircle scaled D shifted z2 withcolor RoyalBlue withpen pencircle scaled 1.0;

x3 = -.0000000001u;
y3 = .4374906403u;
draw fullcircle scaled D shifted z3 withcolor RoyalBlue withpen pencircle scaled 1.0;
draw (       -0.2500093583u,       -0.0625147074u)--(       -0.5000000000u,       -0.5000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.2500093583u,       -0.0625147074u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.2500093583u,       -0.0625147074u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.5000000000u,       -0.5000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.5000000000u,       -0.5000000000u) withcolor Yellow;
draw (       -0.2500093583u,       -0.0625147074u)--(        0.0000000001u,       -0.4999893037u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2500093584u,       -0.0625147073u)--(        0.0000000001u,       -0.4999893037u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.2500093583u,       -0.0625147074u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.2500093583u,       -0.0625147074u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2500093584u,       -0.0625147073u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2500093584u,       -0.0625147073u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.4999893037u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.4999893037u) withcolor Green;
draw (       -0.2500093583u,       -0.0625147074u)--(       -0.4999812814u,        0.3749812805u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.0000000001u,        0.4374906403u)--(       -0.4999812814u,        0.3749812805u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.2500093583u,       -0.0625147074u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.2500093583u,       -0.0625147074u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0000000001u,        0.4374906403u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0000000001u,        0.4374906403u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.4999812814u,        0.3749812805u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.4999812814u,        0.3749812805u) withcolor Green;
draw (        0.2500093584u,       -0.0625147073u)--(        0.5000000000u,       -0.5000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2500093584u,       -0.0625147073u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2500093584u,       -0.0625147073u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5000000000u,       -0.5000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000000u,       -0.5000000000u) withcolor Yellow;
draw (        0.2500093584u,       -0.0625147073u)--(        0.4999812813u,        0.3749812807u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.0000000001u,        0.4374906403u)--(        0.4999812813u,        0.3749812807u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2500093584u,       -0.0625147073u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2500093584u,       -0.0625147073u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0000000001u,        0.4374906403u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0000000001u,        0.4374906403u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4999812813u,        0.3749812807u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4999812813u,        0.3749812807u) withcolor Green;
draw (       -0.0000000001u,        0.4374906403u)--(        0.5000000000u,        0.5000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.0000000001u,        0.4374906403u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0000000001u,        0.4374906403u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5000000000u,        0.5000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000000u,        0.5000000000u) withcolor Yellow;
draw (       -0.0000000001u,        0.4374906403u)--(       -0.5000000000u,        0.5000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.0000000001u,        0.4374906403u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0000000001u,        0.4374906403u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.5000000000u,        0.5000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.5000000000u,        0.5000000000u) withcolor Yellow;
% r =         0.5038736962
endfig;
end;
