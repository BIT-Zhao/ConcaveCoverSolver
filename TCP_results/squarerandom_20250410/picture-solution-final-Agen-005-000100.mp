prologues := 3;
outputformat   := "svg";
outputtemplate := "%j.%o";
input mpcolornames;
beginfig(5);
u := 20cm;
path pols[];
pair p[];
R := .3261408863u;
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

x1 = .2500066394u;
y1 = .2905173613u;
draw fullcircle scaled D shifted z1 withcolor RoyalBlue withpen pencircle scaled 1.0;

x2 = -.3517515381u;
y2 = -.2094999604u;
draw fullcircle scaled D shifted z2 withcolor RoyalBlue withpen pencircle scaled 1.0;

x3 = .3517515381u;
y3 = -.2094999604u;
draw fullcircle scaled D shifted z3 withcolor RoyalBlue withpen pencircle scaled 1.0;

x4 = -.2500066392u;
y4 = .2905173614u;
draw fullcircle scaled D shifted z4 withcolor RoyalBlue withpen pencircle scaled 1.0;

x5 = -.0000000001u;
y5 = -.2450903120u;
draw fullcircle scaled D shifted z5 withcolor RoyalBlue withpen pencircle scaled 1.0;
draw (        0.2500066394u,        0.2905173613u)--(        0.5000000000u,        0.5000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2500066394u,        0.2905173613u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2500066394u,        0.2905173613u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5000000000u,        0.5000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000000u,        0.5000000000u) withcolor Yellow;
draw (        0.2500066394u,        0.2905173613u)--(        0.4999631588u,        0.0810188741u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3517515381u,       -0.2094999604u)--(        0.4999631588u,        0.0810188741u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2500066394u,        0.2905173613u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2500066394u,        0.2905173613u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3517515381u,       -0.2094999604u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3517515381u,       -0.2094999604u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4999631588u,        0.0810188741u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4999631588u,        0.0810188741u) withcolor Green;
draw (        0.2500066394u,        0.2905173613u)--(        0.0000000000u,        0.0810786876u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.2500066392u,        0.2905173614u)--(        0.0000000000u,        0.0810786876u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2500066394u,        0.2905173613u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2500066394u,        0.2905173613u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.2500066392u,        0.2905173614u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.2500066392u,        0.2905173614u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0810786876u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0810786876u) withcolor MidnightBlue;
draw (        0.2500066394u,        0.2905173613u)--(        0.0000000001u,        0.4999560350u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.2500066392u,        0.2905173614u)--(        0.0000000001u,        0.4999560350u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2500066394u,        0.2905173613u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2500066394u,        0.2905173613u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.2500066392u,        0.2905173614u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.2500066392u,        0.2905173614u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000001u,        0.4999560350u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000001u,        0.4999560350u) withcolor Green;
draw (        0.2500066394u,        0.2905173613u)--(        0.0000235549u,        0.0810505734u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.0000000001u,       -0.2450903120u)--(        0.0000235549u,        0.0810505734u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2500066394u,        0.2905173613u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2500066394u,        0.2905173613u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0000000001u,       -0.2450903120u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0000000001u,       -0.2450903120u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000235549u,        0.0810505734u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000235549u,        0.0810505734u) withcolor MidnightBlue;
draw (       -0.3517515381u,       -0.2094999604u)--(       -0.5000000000u,       -0.5000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.3517515381u,       -0.2094999604u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.3517515381u,       -0.2094999604u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.5000000000u,       -0.5000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.5000000000u,       -0.5000000000u) withcolor Yellow;
draw (       -0.3517515381u,       -0.2094999604u)--(       -0.4999631587u,        0.0810188742u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.2500066392u,        0.2905173614u)--(       -0.4999631587u,        0.0810188742u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.3517515381u,       -0.2094999604u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.3517515381u,       -0.2094999604u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.2500066392u,        0.2905173614u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.2500066392u,        0.2905173614u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.4999631587u,        0.0810188742u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.4999631587u,        0.0810188742u) withcolor Green;
draw (       -0.3517515381u,       -0.2094999604u)--(       -0.2034662465u,       -0.4999811991u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.0000000001u,       -0.2450903120u)--(       -0.2034662465u,       -0.4999811991u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.3517515381u,       -0.2094999604u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.3517515381u,       -0.2094999604u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0000000001u,       -0.2450903120u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0000000001u,       -0.2450903120u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.2034662465u,       -0.4999811991u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.2034662465u,       -0.4999811991u) withcolor Green;
draw (        0.3517515381u,       -0.2094999604u)--(        0.5000000000u,       -0.5000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3517515381u,       -0.2094999604u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3517515381u,       -0.2094999604u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5000000000u,       -0.5000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000000u,       -0.5000000000u) withcolor Yellow;
draw (        0.3517515381u,       -0.2094999604u)--(        0.2034662464u,       -0.4999811991u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.0000000001u,       -0.2450903120u)--(        0.2034662464u,       -0.4999811991u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3517515381u,       -0.2094999604u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3517515381u,       -0.2094999604u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0000000001u,       -0.2450903120u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0000000001u,       -0.2450903120u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2034662464u,       -0.4999811991u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2034662464u,       -0.4999811991u) withcolor Green;
draw (       -0.2500066392u,        0.2905173614u)--(       -0.5000000000u,        0.5000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.2500066392u,        0.2905173614u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.2500066392u,        0.2905173614u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.5000000000u,        0.5000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.5000000000u,        0.5000000000u) withcolor Yellow;
draw (       -0.2500066392u,        0.2905173614u)--(       -0.0000235548u,        0.0810505734u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.0000000001u,       -0.2450903120u)--(       -0.0000235548u,        0.0810505734u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.2500066392u,        0.2905173614u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.2500066392u,        0.2905173614u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0000000001u,       -0.2450903120u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0000000001u,       -0.2450903120u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0000235548u,        0.0810505734u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0000235548u,        0.0810505734u) withcolor MidnightBlue;
% r =         0.3261408863
endfig;
end;
