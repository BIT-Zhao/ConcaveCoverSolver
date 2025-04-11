prologues := 3;
outputformat   := "svg";
outputtemplate := "%j.%o";
input mpcolornames;
beginfig(1);
u := 20cm;
path pols[];
pair p[];
R := .7070568509u;
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

x1 = .0000000001u;
y1 = -.0000000001u;
draw fullcircle scaled D shifted z1 withcolor RoyalBlue withpen pencircle scaled 1.0;
draw (        0.0000000001u,       -0.0000000001u)--(       -0.5000000000u,       -0.5000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.5000000000u,       -0.5000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.5000000000u,       -0.5000000000u) withcolor Yellow;
draw (        0.0000000001u,       -0.0000000001u)--(        0.5000000000u,       -0.5000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5000000000u,       -0.5000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000000u,       -0.5000000000u) withcolor Yellow;
draw (        0.0000000001u,       -0.0000000001u)--(        0.5000000000u,        0.5000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5000000000u,        0.5000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000000u,        0.5000000000u) withcolor Yellow;
draw (        0.0000000001u,       -0.0000000001u)--(       -0.5000000000u,        0.5000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,       -0.0000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.5000000000u,        0.5000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.5000000000u,        0.5000000000u) withcolor Yellow;
% r =         0.7070568509
endfig;
end;
