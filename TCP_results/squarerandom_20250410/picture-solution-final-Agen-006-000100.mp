prologues := 3;
outputformat   := "svg";
outputtemplate := "%j.%o";
input mpcolornames;
beginfig(6);
u := 20cm;
path pols[];
pair p[];
R := .2987124703u;
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

x1 = -.0402288121u;
y1 = .2589747241u;
draw fullcircle scaled D shifted z1 withcolor RoyalBlue withpen pencircle scaled 1.0;

x2 = .0402288122u;
y2 = -.2589747242u;
draw fullcircle scaled D shifted z2 withcolor RoyalBlue withpen pencircle scaled 1.0;

x3 = .3583691394u;
y3 = -.2369984019u;
draw fullcircle scaled D shifted z3 withcolor RoyalBlue withpen pencircle scaled 1.0;

x4 = -.3583691394u;
y4 = .2369984019u;
draw fullcircle scaled D shifted z4 withcolor RoyalBlue withpen pencircle scaled 1.0;

x5 = .3181290669u;
y5 = .2630094094u;
draw fullcircle scaled D shifted z5 withcolor RoyalBlue withpen pencircle scaled 1.0;

x6 = -.3181290670u;
y6 = -.2630094093u;
draw fullcircle scaled D shifted z6 withcolor RoyalBlue withpen pencircle scaled 1.0;
draw (       -0.0402288121u,        0.2589747241u)--(        0.1416292772u,        0.0220005160u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0402288122u,       -0.2589747242u)--(        0.1416292772u,        0.0220005160u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.0402288121u,        0.2589747241u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0402288121u,        0.2589747241u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0402288122u,       -0.2589747242u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0402288122u,       -0.2589747242u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1416292772u,        0.0220005160u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1416292772u,        0.0220005160u) withcolor MidnightBlue;
draw (       -0.0402288121u,        0.2589747241u)--(       -0.1416292771u,       -0.0220005161u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0402288122u,       -0.2589747242u)--(       -0.1416292771u,       -0.0220005161u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.0402288121u,        0.2589747241u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0402288121u,        0.2589747241u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0402288122u,       -0.2589747242u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0402288122u,       -0.2589747242u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.1416292771u,       -0.0220005161u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1416292771u,       -0.0220005161u) withcolor MidnightBlue;
draw (       -0.0402288121u,        0.2589747241u)--(       -0.2167062285u,        0.4999827353u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.3583691394u,        0.2369984019u)--(       -0.2167062285u,        0.4999827353u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.0402288121u,        0.2589747241u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0402288121u,        0.2589747241u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.3583691394u,        0.2369984019u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.3583691394u,        0.2369984019u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.2167062285u,        0.4999827353u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.2167062285u,        0.4999827353u) withcolor Green;
draw (       -0.0402288121u,        0.2589747241u)--(        0.1362594651u,        0.4999747821u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3181290669u,        0.2630094094u)--(        0.1362594651u,        0.4999747821u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.0402288121u,        0.2589747241u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0402288121u,        0.2589747241u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3181290669u,        0.2630094094u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3181290669u,        0.2630094094u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1362594651u,        0.4999747821u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1362594651u,        0.4999747821u) withcolor Green;
draw (       -0.0402288121u,        0.2589747241u)--(        0.1416407898u,        0.0220093514u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3181290669u,        0.2630094094u)--(        0.1416407898u,        0.0220093514u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.0402288121u,        0.2589747241u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0402288121u,        0.2589747241u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3181290669u,        0.2630094094u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3181290669u,        0.2630094094u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1416407898u,        0.0220093514u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1416407898u,        0.0220093514u) withcolor MidnightBlue;
draw (       -0.0402288121u,        0.2589747241u)--(       -0.1416757851u,       -0.0219837276u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.3181290670u,       -0.2630094093u)--(       -0.1416757851u,       -0.0219837276u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.0402288121u,        0.2589747241u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0402288121u,        0.2589747241u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.3181290670u,       -0.2630094093u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.3181290670u,       -0.2630094093u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.1416757851u,       -0.0219837276u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1416757851u,       -0.0219837276u) withcolor MidnightBlue;
draw (        0.0402288122u,       -0.2589747242u)--(        0.2167062286u,       -0.4999827353u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3583691394u,       -0.2369984019u)--(        0.2167062286u,       -0.4999827353u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0402288122u,       -0.2589747242u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0402288122u,       -0.2589747242u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3583691394u,       -0.2369984019u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3583691394u,       -0.2369984019u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2167062286u,       -0.4999827353u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2167062286u,       -0.4999827353u) withcolor Green;
draw (        0.0402288122u,       -0.2589747242u)--(        0.1416757855u,        0.0219837274u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3181290669u,        0.2630094094u)--(        0.1416757855u,        0.0219837274u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0402288122u,       -0.2589747242u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0402288122u,       -0.2589747242u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3181290669u,        0.2630094094u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3181290669u,        0.2630094094u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1416757855u,        0.0219837274u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1416757855u,        0.0219837274u) withcolor MidnightBlue;
draw (        0.0402288122u,       -0.2589747242u)--(       -0.1362594652u,       -0.4999747820u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.3181290670u,       -0.2630094093u)--(       -0.1362594652u,       -0.4999747820u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0402288122u,       -0.2589747242u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0402288122u,       -0.2589747242u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.3181290670u,       -0.2630094093u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.3181290670u,       -0.2630094093u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.1362594652u,       -0.4999747820u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1362594652u,       -0.4999747820u) withcolor Green;
draw (        0.0402288122u,       -0.2589747242u)--(       -0.1416407896u,       -0.0220093515u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.3181290670u,       -0.2630094093u)--(       -0.1416407896u,       -0.0220093515u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0402288122u,       -0.2589747242u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0402288122u,       -0.2589747242u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.3181290670u,       -0.2630094093u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.3181290670u,       -0.2630094093u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.1416407896u,       -0.0220093515u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1416407896u,       -0.0220093515u) withcolor MidnightBlue;
draw (        0.3583691394u,       -0.2369984019u)--(        0.5000000000u,       -0.5000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3583691394u,       -0.2369984019u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3583691394u,       -0.2369984019u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5000000000u,       -0.5000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000000u,       -0.5000000000u) withcolor Yellow;
draw (        0.3583691394u,       -0.2369984019u)--(        0.4999679402u,        0.0260204559u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3181290669u,        0.2630094094u)--(        0.4999679402u,        0.0260204559u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3583691394u,       -0.2369984019u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3583691394u,       -0.2369984019u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3181290669u,        0.2630094094u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3181290669u,        0.2630094094u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4999679402u,        0.0260204559u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4999679402u,        0.0260204559u) withcolor Green;
draw (       -0.3583691394u,        0.2369984019u)--(       -0.5000000000u,        0.5000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.3583691394u,        0.2369984019u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.3583691394u,        0.2369984019u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.5000000000u,        0.5000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.5000000000u,        0.5000000000u) withcolor Yellow;
draw (       -0.3583691394u,        0.2369984019u)--(       -0.4999679403u,       -0.0260204558u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.3181290670u,       -0.2630094093u)--(       -0.4999679403u,       -0.0260204558u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.3583691394u,        0.2369984019u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.3583691394u,        0.2369984019u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.3181290670u,       -0.2630094093u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.3181290670u,       -0.2630094093u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.4999679403u,       -0.0260204558u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.4999679403u,       -0.0260204558u) withcolor Green;
draw (        0.3181290669u,        0.2630094094u)--(        0.5000000000u,        0.5000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3181290669u,        0.2630094094u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3181290669u,        0.2630094094u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5000000000u,        0.5000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000000u,        0.5000000000u) withcolor Yellow;
draw (       -0.3181290670u,       -0.2630094093u)--(       -0.5000000000u,       -0.5000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.3181290670u,       -0.2630094093u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.3181290670u,       -0.2630094093u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.5000000000u,       -0.5000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.5000000000u,       -0.5000000000u) withcolor Yellow;
% r =         0.2987124703
endfig;
end;
