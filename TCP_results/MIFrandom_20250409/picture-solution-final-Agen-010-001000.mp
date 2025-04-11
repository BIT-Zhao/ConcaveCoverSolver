prologues := 3;
outputformat   := "svg";
outputtemplate := "%j.%o";
input mpcolornames;
beginfig(10);
u := 20cm;
path pols[];
pair p[];
R := .1668792477u;
D := 2*R;

pols[1] = 
(       -0.3816000000u,       -0.1780000000u)--
(       -0.3816000000u,       -0.3433000000u)--
(       -0.2109000000u,       -0.3433000000u)--
(       -0.2109000000u,       -0.1780000000u)--
(       -0.3816000000u,       -0.1780000000u)--
cycle;
pols[2] = 
(       -0.0483000000u,       -0.5059000000u)--
(       -0.0483000000u,       -0.5086000000u)--
(        0.1198000000u,       -0.5086000000u)--
(        0.1198000000u,       -0.5059000000u)--
(       -0.0483000000u,       -0.5059000000u)--
cycle;
pols[3] = 
(        0.2851000000u,       -0.3406000000u)--
(        0.2878000000u,       -0.3406000000u)--
(        0.2878000000u,       -0.1780000000u)--
(        0.2851000000u,       -0.1780000000u)--
(        0.2851000000u,       -0.3406000000u)--
cycle;
pols[4] = 
(        0.1198000000u,       -0.3433000000u)--
(        0.2851000000u,       -0.3433000000u)--
(        0.2851000000u,       -0.3406000000u)--
(        0.1198000000u,       -0.3433000000u)--
cycle;
pols[5] = 
(        0.2851000000u,       -0.3406000000u)--
(        0.2851000000u,       -0.1780000000u)--
(        0.1198000000u,       -0.1780000000u)--
(        0.1198000000u,       -0.3433000000u)--
(        0.2851000000u,       -0.3406000000u)--
cycle;
pols[6] = 
(       -0.0510000000u,       -0.1780000000u)--
(       -0.0510000000u,       -0.5059000000u)--
(       -0.0483000000u,       -0.5059000000u)--
(       -0.0510000000u,       -0.1780000000u)--
cycle;
pols[7] = 
(       -0.0483000000u,       -0.5059000000u)--
(        0.1198000000u,       -0.5059000000u)--
(        0.1198000000u,       -0.3433000000u)--
(       -0.0483000000u,       -0.5059000000u)--
cycle;
pols[8] = 
(        0.1198000000u,       -0.3433000000u)--
(        0.1198000000u,       -0.1780000000u)--
(       -0.0510000000u,       -0.1780000000u)--
(       -0.0483000000u,       -0.5059000000u)--
(        0.1198000000u,       -0.3433000000u)--
cycle;
pols[9] = 
(       -0.5470000000u,       -0.1752000000u)--
(       -0.5469000000u,       -0.1780000000u)--
(       -0.3816000000u,       -0.1780000000u)--
(       -0.5470000000u,       -0.1752000000u)--
cycle;
pols[10] = 
(        0.2878000000u,       -0.1780000000u)--
(        0.2878000000u,       -0.1752000000u)--
(        0.2851000000u,       -0.1752000000u)--
(        0.2878000000u,       -0.1780000000u)--
cycle;
pols[11] = 
(       -0.5470000000u,       -0.1752000000u)--
(       -0.3816000000u,       -0.1780000000u)--
(       -0.2109000000u,       -0.1780000000u)--
(       -0.5470000000u,       -0.1752000000u)--
cycle;
pols[12] = 
(       -0.5470000000u,       -0.1752000000u)--
(       -0.2109000000u,       -0.1780000000u)--
(       -0.0510000000u,       -0.1780000000u)--
(       -0.5470000000u,       -0.1752000000u)--
cycle;
pols[13] = 
(       -0.0510000000u,       -0.1780000000u)--
(        0.2878000000u,       -0.1780000000u)--
(        0.2851000000u,       -0.1752000000u)--
(       -0.5470000000u,       -0.1752000000u)--
(       -0.0510000000u,       -0.1780000000u)--
cycle;
pols[14] = 
(        0.2851000000u,       -0.1752000000u)--
(        0.2851000000u,       -0.1725000000u)--
(        0.1198000000u,       -0.1725000000u)--
(        0.2851000000u,       -0.1752000000u)--
cycle;
pols[15] = 
(        0.1198000000u,       -0.1725000000u)--
(       -0.5470000000u,       -0.1725000000u)--
(       -0.5470000000u,       -0.1752000000u)--
(        0.2851000000u,       -0.1752000000u)--
(        0.1198000000u,       -0.1725000000u)--
cycle;
pols[16] = 
(        0.4504000000u,       -0.0099000000u)--
(        0.4531000000u,       -0.0099000000u)--
(        0.4531000000u,        0.1608000000u)--
(        0.4504000000u,        0.1608000000u)--
(        0.4504000000u,       -0.0099000000u)--
cycle;
pols[17] = 
(        0.2851000000u,        0.3261000000u)--
(        0.2851000000u,        0.3288000000u)--
(        0.1198000000u,        0.3288000000u)--
(        0.1198000000u,        0.3261000000u)--
(        0.2851000000u,        0.3261000000u)--
cycle;
pols[18] = 
(        0.2878000000u,        0.1608000000u)--
(        0.2878000000u,        0.3261000000u)--
(        0.2851000000u,        0.3261000000u)--
(        0.2878000000u,        0.1608000000u)--
cycle;
pols[19] = 
(        0.2851000000u,        0.3261000000u)--
(        0.1198000000u,        0.3261000000u)--
(        0.1198000000u,        0.1608000000u)--
(        0.2878000000u,        0.1608000000u)--
(        0.2851000000u,        0.3261000000u)--
cycle;
pols[20] = 
(        0.1198000000u,       -0.0126000000u)--
(        0.4504000000u,       -0.0126000000u)--
(        0.4504000000u,       -0.0099000000u)--
(        0.1198000000u,       -0.0126000000u)--
cycle;
pols[21] = 
(        0.4504000000u,       -0.0099000000u)--
(        0.4504000000u,        0.1608000000u)--
(        0.2878000000u,        0.1608000000u)--
(        0.4504000000u,       -0.0099000000u)--
cycle;
pols[22] = 
(        0.2878000000u,        0.1608000000u)--
(        0.1198000000u,        0.1608000000u)--
(        0.1198000000u,       -0.0126000000u)--
(        0.4504000000u,       -0.0099000000u)--
(        0.2878000000u,        0.1608000000u)--
cycle;
pols[23] = 
(        0.1170000000u,       -0.1725000000u)--
(        0.1198000000u,       -0.1725000000u)--
(        0.1198000000u,       -0.0126000000u)--
(        0.1170000000u,       -0.1725000000u)--
cycle;
pols[24] = 
(        0.1198000000u,        0.3288000000u)--
(        0.1170000000u,        0.3288000000u)--
(        0.1170000000u,        0.1608000000u)--
(        0.1198000000u,        0.3288000000u)--
cycle;
pols[25] = 
(        0.1170000000u,        0.1608000000u)--
(        0.1170000000u,       -0.1725000000u)--
(        0.1198000000u,       -0.0126000000u)--
(        0.1198000000u,        0.3288000000u)--
(        0.1170000000u,        0.1608000000u)--
cycle;
pols[26] = 
(       -0.3789000000u,        0.3261000000u)--
(       -0.3816000000u,        0.3261000000u)--
(       -0.3816000000u,        0.1608000000u)--
(       -0.3789000000u,        0.1608000000u)--
(       -0.3789000000u,        0.3261000000u)--
cycle;
pols[27] = 
(       -0.2163000000u,        0.3288000000u)--
(       -0.3789000000u,        0.3288000000u)--
(       -0.3789000000u,        0.3261000000u)--
(       -0.2163000000u,        0.3288000000u)--
cycle;
pols[28] = 
(       -0.3789000000u,        0.3261000000u)--
(       -0.3789000000u,        0.1608000000u)--
(       -0.2163000000u,        0.1608000000u)--
(       -0.2163000000u,        0.3288000000u)--
(       -0.3789000000u,        0.3261000000u)--
cycle;
pols[29] = 
(       -0.2136000000u,        0.4914000000u)--
(       -0.2163000000u,        0.4914000000u)--
(       -0.2163000000u,        0.3288000000u)--
(       -0.2136000000u,        0.4914000000u)--
cycle;
pols[30] = 
(       -0.2163000000u,        0.3288000000u)--
(       -0.2163000000u,        0.1608000000u)--
(       -0.2136000000u,        0.1608000000u)--
(       -0.2136000000u,        0.4914000000u)--
(       -0.2163000000u,        0.3288000000u)--
cycle;
pols[31] = 
(       -0.0456000000u,        0.4941000000u)--
(       -0.2136000000u,        0.4941000000u)--
(       -0.2136000000u,        0.4914000000u)--
(       -0.0456000000u,        0.4941000000u)--
cycle;
pols[32] = 
(       -0.2136000000u,        0.4914000000u)--
(       -0.2136000000u,        0.1608000000u)--
(       -0.0456000000u,        0.1608000000u)--
(       -0.0456000000u,        0.4941000000u)--
(       -0.2136000000u,        0.4914000000u)--
cycle;
pols[33] = 
(        0.1170000000u,        0.1581000000u)--
(        0.1170000000u,        0.1608000000u)--
(       -0.0456000000u,        0.1608000000u)--
(        0.1170000000u,        0.1581000000u)--
cycle;
pols[34] = 
(       -0.3816000000u,        0.1608000000u)--
(       -0.3816000000u,        0.1581000000u)--
(       -0.2163000000u,        0.1581000000u)--
(       -0.3816000000u,        0.1608000000u)--
cycle;
pols[35] = 
(       -0.2163000000u,        0.1581000000u)--
(        0.1170000000u,        0.1581000000u)--
(       -0.0456000000u,        0.1608000000u)--
(       -0.3816000000u,        0.1608000000u)--
(       -0.2163000000u,        0.1581000000u)--
cycle;
pols[36] = 
(        0.1170000000u,        0.1581000000u)--
(       -0.2163000000u,        0.1581000000u)--
(       -0.2163000000u,       -0.0072000000u)--
(        0.1170000000u,        0.1581000000u)--
cycle;
pols[37] = 
(       -0.2163000000u,       -0.0072000000u)--
(       -0.2163000000u,       -0.1725000000u)--
(        0.1170000000u,       -0.1725000000u)--
(        0.1170000000u,        0.1581000000u)--
(       -0.2163000000u,       -0.0072000000u)--
cycle;
pols[38] = 
(       -0.5469000000u,       -0.0072000000u)--
(       -0.5470000000u,       -0.1725000000u)--
(       -0.2163000000u,       -0.1725000000u)--
(       -0.2163000000u,       -0.0072000000u)--
(       -0.5469000000u,       -0.0072000000u)--
cycle;
fill pols[1] withcolor 0.9Dandelion;
fill pols[2] withcolor 0.9Dandelion;
fill pols[3] withcolor 0.9Dandelion;
fill pols[4] withcolor 0.9Dandelion;
fill pols[5] withcolor 0.9Dandelion;
fill pols[6] withcolor 0.9Dandelion;
fill pols[7] withcolor 0.9Dandelion;
fill pols[8] withcolor 0.9Dandelion;
fill pols[9] withcolor 0.9Dandelion;
fill pols[10] withcolor 0.9Dandelion;
fill pols[11] withcolor 0.9Dandelion;
fill pols[12] withcolor 0.9Dandelion;
fill pols[13] withcolor 0.9Dandelion;
fill pols[14] withcolor 0.9Dandelion;
fill pols[15] withcolor 0.9Dandelion;
fill pols[16] withcolor 0.9Dandelion;
fill pols[17] withcolor 0.9Dandelion;
fill pols[18] withcolor 0.9Dandelion;
fill pols[19] withcolor 0.9Dandelion;
fill pols[20] withcolor 0.9Dandelion;
fill pols[21] withcolor 0.9Dandelion;
fill pols[22] withcolor 0.9Dandelion;
fill pols[23] withcolor 0.9Dandelion;
fill pols[24] withcolor 0.9Dandelion;
fill pols[25] withcolor 0.9Dandelion;
fill pols[26] withcolor 0.9Dandelion;
fill pols[27] withcolor 0.9Dandelion;
fill pols[28] withcolor 0.9Dandelion;
fill pols[29] withcolor 0.9Dandelion;
fill pols[30] withcolor 0.9Dandelion;
fill pols[31] withcolor 0.9Dandelion;
fill pols[32] withcolor 0.9Dandelion;
fill pols[33] withcolor 0.9Dandelion;
fill pols[34] withcolor 0.9Dandelion;
fill pols[35] withcolor 0.9Dandelion;
fill pols[36] withcolor 0.9Dandelion;
fill pols[37] withcolor 0.9Dandelion;
fill pols[38] withcolor 0.9Dandelion;
draw (       -0.2109000000u,       -0.1780000000u)--(       -0.3816000000u,       -0.1780000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1198000000u,       -0.5059000000u)--(       -0.0483000000u,       -0.5059000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.2878000000u,       -0.1780000000u)--(        0.2851000000u,       -0.1780000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.2851000000u,       -0.1780000000u)--(        0.2851000000u,       -0.3406000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.2851000000u,       -0.3406000000u)--(        0.1198000000u,       -0.3433000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.2851000000u,       -0.3406000000u)--(        0.2851000000u,       -0.1780000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.2851000000u,       -0.1780000000u)--(        0.1198000000u,       -0.1780000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1198000000u,       -0.1780000000u)--(        0.1198000000u,       -0.3433000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1198000000u,       -0.3433000000u)--(        0.2851000000u,       -0.3406000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.0483000000u,       -0.5059000000u)--(       -0.0510000000u,       -0.1780000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.0483000000u,       -0.5059000000u)--(        0.1198000000u,       -0.5059000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1198000000u,       -0.3433000000u)--(       -0.0483000000u,       -0.5059000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1198000000u,       -0.3433000000u)--(        0.1198000000u,       -0.1780000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1198000000u,       -0.1780000000u)--(       -0.0510000000u,       -0.1780000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.0510000000u,       -0.1780000000u)--(       -0.0483000000u,       -0.5059000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.0483000000u,       -0.5059000000u)--(        0.1198000000u,       -0.3433000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.3816000000u,       -0.1780000000u)--(       -0.5470000000u,       -0.1752000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.2851000000u,       -0.1752000000u)--(        0.2878000000u,       -0.1780000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.5470000000u,       -0.1752000000u)--(       -0.3816000000u,       -0.1780000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.3816000000u,       -0.1780000000u)--(       -0.2109000000u,       -0.1780000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.2109000000u,       -0.1780000000u)--(       -0.5470000000u,       -0.1752000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.5470000000u,       -0.1752000000u)--(       -0.2109000000u,       -0.1780000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.0510000000u,       -0.1780000000u)--(       -0.5470000000u,       -0.1752000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.0510000000u,       -0.1780000000u)--(        0.2878000000u,       -0.1780000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.2878000000u,       -0.1780000000u)--(        0.2851000000u,       -0.1752000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.2851000000u,       -0.1752000000u)--(       -0.5470000000u,       -0.1752000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.5470000000u,       -0.1752000000u)--(       -0.0510000000u,       -0.1780000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1198000000u,       -0.1725000000u)--(        0.2851000000u,       -0.1752000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1198000000u,       -0.1725000000u)--(       -0.5470000000u,       -0.1725000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.5470000000u,       -0.1752000000u)--(        0.2851000000u,       -0.1752000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.2851000000u,       -0.1752000000u)--(        0.1198000000u,       -0.1725000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.4504000000u,        0.1608000000u)--(        0.4504000000u,       -0.0099000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1198000000u,        0.3288000000u)--(        0.1198000000u,        0.3261000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1198000000u,        0.3261000000u)--(        0.2851000000u,        0.3261000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.2851000000u,        0.3261000000u)--(        0.2878000000u,        0.1608000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.2851000000u,        0.3261000000u)--(        0.1198000000u,        0.3261000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1198000000u,        0.3261000000u)--(        0.1198000000u,        0.1608000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1198000000u,        0.1608000000u)--(        0.2878000000u,        0.1608000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.2878000000u,        0.1608000000u)--(        0.2851000000u,        0.3261000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.4504000000u,       -0.0099000000u)--(        0.1198000000u,       -0.0126000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.4504000000u,       -0.0099000000u)--(        0.4504000000u,        0.1608000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.2878000000u,        0.1608000000u)--(        0.4504000000u,       -0.0099000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.2878000000u,        0.1608000000u)--(        0.1198000000u,        0.1608000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1198000000u,        0.1608000000u)--(        0.1198000000u,       -0.0126000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1198000000u,       -0.0126000000u)--(        0.4504000000u,       -0.0099000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.4504000000u,       -0.0099000000u)--(        0.2878000000u,        0.1608000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1170000000u,       -0.1725000000u)--(        0.1198000000u,       -0.1725000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1198000000u,       -0.0126000000u)--(        0.1170000000u,       -0.1725000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1170000000u,        0.1608000000u)--(        0.1198000000u,        0.3288000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1170000000u,        0.1608000000u)--(        0.1170000000u,       -0.1725000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1170000000u,       -0.1725000000u)--(        0.1198000000u,       -0.0126000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1198000000u,       -0.0126000000u)--(        0.1198000000u,        0.3288000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1198000000u,        0.3288000000u)--(        0.1170000000u,        0.1608000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.3816000000u,        0.1608000000u)--(       -0.3789000000u,        0.1608000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.3789000000u,        0.1608000000u)--(       -0.3789000000u,        0.3261000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.3789000000u,        0.3261000000u)--(       -0.2163000000u,        0.3288000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.3789000000u,        0.3261000000u)--(       -0.3789000000u,        0.1608000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.3789000000u,        0.1608000000u)--(       -0.2163000000u,        0.1608000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.2163000000u,        0.1608000000u)--(       -0.2163000000u,        0.3288000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.2163000000u,        0.3288000000u)--(       -0.3789000000u,        0.3261000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.2163000000u,        0.3288000000u)--(       -0.2136000000u,        0.4914000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.2163000000u,        0.3288000000u)--(       -0.2163000000u,        0.1608000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.2163000000u,        0.1608000000u)--(       -0.2136000000u,        0.1608000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.2136000000u,        0.1608000000u)--(       -0.2136000000u,        0.4914000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.2136000000u,        0.4914000000u)--(       -0.2163000000u,        0.3288000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.2136000000u,        0.4914000000u)--(       -0.0456000000u,        0.4941000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.2136000000u,        0.4914000000u)--(       -0.2136000000u,        0.1608000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.2136000000u,        0.1608000000u)--(       -0.0456000000u,        0.1608000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.0456000000u,        0.4941000000u)--(       -0.2136000000u,        0.4914000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1170000000u,        0.1581000000u)--(        0.1170000000u,        0.1608000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.0456000000u,        0.1608000000u)--(        0.1170000000u,        0.1581000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.2163000000u,        0.1581000000u)--(       -0.3816000000u,        0.1608000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.2163000000u,        0.1581000000u)--(        0.1170000000u,        0.1581000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1170000000u,        0.1581000000u)--(       -0.0456000000u,        0.1608000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.0456000000u,        0.1608000000u)--(       -0.3816000000u,        0.1608000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.3816000000u,        0.1608000000u)--(       -0.2163000000u,        0.1581000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1170000000u,        0.1581000000u)--(       -0.2163000000u,        0.1581000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.2163000000u,       -0.0072000000u)--(        0.1170000000u,        0.1581000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.2163000000u,       -0.0072000000u)--(       -0.2163000000u,       -0.1725000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.2163000000u,       -0.1725000000u)--(        0.1170000000u,       -0.1725000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1170000000u,       -0.1725000000u)--(        0.1170000000u,        0.1581000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1170000000u,        0.1581000000u)--(       -0.2163000000u,       -0.0072000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.5470000000u,       -0.1725000000u)--(       -0.2163000000u,       -0.1725000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.2163000000u,       -0.1725000000u)--(       -0.2163000000u,       -0.0072000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (       -0.3816000000u,       -0.1780000000u)--(       -0.3816000000u,       -0.3433000000u) withpen pencircle scaled 1.0;
draw (       -0.3816000000u,       -0.3433000000u)--(       -0.2109000000u,       -0.3433000000u) withpen pencircle scaled 1.0;
draw (       -0.2109000000u,       -0.3433000000u)--(       -0.2109000000u,       -0.1780000000u) withpen pencircle scaled 1.0;
draw (       -0.0483000000u,       -0.5059000000u)--(       -0.0483000000u,       -0.5086000000u) withpen pencircle scaled 1.0;
draw (       -0.0483000000u,       -0.5086000000u)--(        0.1198000000u,       -0.5086000000u) withpen pencircle scaled 1.0;
draw (        0.1198000000u,       -0.5086000000u)--(        0.1198000000u,       -0.5059000000u) withpen pencircle scaled 1.0;
draw (        0.2851000000u,       -0.3406000000u)--(        0.2878000000u,       -0.3406000000u) withpen pencircle scaled 1.0;
draw (        0.2878000000u,       -0.3406000000u)--(        0.2878000000u,       -0.1780000000u) withpen pencircle scaled 1.0;
draw (        0.1198000000u,       -0.3433000000u)--(        0.2851000000u,       -0.3433000000u) withpen pencircle scaled 1.0;
draw (        0.2851000000u,       -0.3433000000u)--(        0.2851000000u,       -0.3406000000u) withpen pencircle scaled 1.0;
draw (       -0.0510000000u,       -0.1780000000u)--(       -0.0510000000u,       -0.5059000000u) withpen pencircle scaled 1.0;
draw (       -0.0510000000u,       -0.5059000000u)--(       -0.0483000000u,       -0.5059000000u) withpen pencircle scaled 1.0;
draw (        0.1198000000u,       -0.5059000000u)--(        0.1198000000u,       -0.3433000000u) withpen pencircle scaled 1.0;
draw (       -0.5470000000u,       -0.1752000000u)--(       -0.5469000000u,       -0.1780000000u) withpen pencircle scaled 1.0;
draw (       -0.5469000000u,       -0.1780000000u)--(       -0.3816000000u,       -0.1780000000u) withpen pencircle scaled 1.0;
draw (        0.2878000000u,       -0.1780000000u)--(        0.2878000000u,       -0.1752000000u) withpen pencircle scaled 1.0;
draw (        0.2878000000u,       -0.1752000000u)--(        0.2851000000u,       -0.1752000000u) withpen pencircle scaled 1.0;
draw (       -0.2109000000u,       -0.1780000000u)--(       -0.0510000000u,       -0.1780000000u) withpen pencircle scaled 1.0;
draw (        0.2851000000u,       -0.1752000000u)--(        0.2851000000u,       -0.1725000000u) withpen pencircle scaled 1.0;
draw (        0.2851000000u,       -0.1725000000u)--(        0.1198000000u,       -0.1725000000u) withpen pencircle scaled 1.0;
draw (       -0.5470000000u,       -0.1725000000u)--(       -0.5470000000u,       -0.1752000000u) withpen pencircle scaled 1.0;
draw (        0.4504000000u,       -0.0099000000u)--(        0.4531000000u,       -0.0099000000u) withpen pencircle scaled 1.0;
draw (        0.4531000000u,       -0.0099000000u)--(        0.4531000000u,        0.1608000000u) withpen pencircle scaled 1.0;
draw (        0.4531000000u,        0.1608000000u)--(        0.4504000000u,        0.1608000000u) withpen pencircle scaled 1.0;
draw (        0.2851000000u,        0.3261000000u)--(        0.2851000000u,        0.3288000000u) withpen pencircle scaled 1.0;
draw (        0.2851000000u,        0.3288000000u)--(        0.1198000000u,        0.3288000000u) withpen pencircle scaled 1.0;
draw (        0.2878000000u,        0.1608000000u)--(        0.2878000000u,        0.3261000000u) withpen pencircle scaled 1.0;
draw (        0.2878000000u,        0.3261000000u)--(        0.2851000000u,        0.3261000000u) withpen pencircle scaled 1.0;
draw (        0.1198000000u,       -0.0126000000u)--(        0.4504000000u,       -0.0126000000u) withpen pencircle scaled 1.0;
draw (        0.4504000000u,       -0.0126000000u)--(        0.4504000000u,       -0.0099000000u) withpen pencircle scaled 1.0;
draw (        0.4504000000u,        0.1608000000u)--(        0.2878000000u,        0.1608000000u) withpen pencircle scaled 1.0;
draw (        0.1198000000u,       -0.1725000000u)--(        0.1198000000u,       -0.0126000000u) withpen pencircle scaled 1.0;
draw (        0.1198000000u,        0.3288000000u)--(        0.1170000000u,        0.3288000000u) withpen pencircle scaled 1.0;
draw (        0.1170000000u,        0.3288000000u)--(        0.1170000000u,        0.1608000000u) withpen pencircle scaled 1.0;
draw (       -0.3789000000u,        0.3261000000u)--(       -0.3816000000u,        0.3261000000u) withpen pencircle scaled 1.0;
draw (       -0.3816000000u,        0.3261000000u)--(       -0.3816000000u,        0.1608000000u) withpen pencircle scaled 1.0;
draw (       -0.2163000000u,        0.3288000000u)--(       -0.3789000000u,        0.3288000000u) withpen pencircle scaled 1.0;
draw (       -0.3789000000u,        0.3288000000u)--(       -0.3789000000u,        0.3261000000u) withpen pencircle scaled 1.0;
draw (       -0.2136000000u,        0.4914000000u)--(       -0.2163000000u,        0.4914000000u) withpen pencircle scaled 1.0;
draw (       -0.2163000000u,        0.4914000000u)--(       -0.2163000000u,        0.3288000000u) withpen pencircle scaled 1.0;
draw (       -0.0456000000u,        0.4941000000u)--(       -0.2136000000u,        0.4941000000u) withpen pencircle scaled 1.0;
draw (       -0.2136000000u,        0.4941000000u)--(       -0.2136000000u,        0.4914000000u) withpen pencircle scaled 1.0;
draw (       -0.0456000000u,        0.1608000000u)--(       -0.0456000000u,        0.4941000000u) withpen pencircle scaled 1.0;
draw (        0.1170000000u,        0.1608000000u)--(       -0.0456000000u,        0.1608000000u) withpen pencircle scaled 1.0;
draw (       -0.3816000000u,        0.1608000000u)--(       -0.3816000000u,        0.1581000000u) withpen pencircle scaled 1.0;
draw (       -0.3816000000u,        0.1581000000u)--(       -0.2163000000u,        0.1581000000u) withpen pencircle scaled 1.0;
draw (       -0.2163000000u,        0.1581000000u)--(       -0.2163000000u,       -0.0072000000u) withpen pencircle scaled 1.0;
draw (       -0.5469000000u,       -0.0072000000u)--(       -0.5470000000u,       -0.1725000000u) withpen pencircle scaled 1.0;
draw (       -0.2163000000u,       -0.0072000000u)--(       -0.5469000000u,       -0.0072000000u) withpen pencircle scaled 1.0;

x1 = .0183326335u;
y1 = .0614979716u;
draw fullcircle scaled D shifted z1 withcolor RoyalBlue withpen pencircle scaled 1.0;

x2 = .1606840943u;
y2 = -.2324312440u;
draw fullcircle scaled D shifted z2 withcolor RoyalBlue withpen pencircle scaled 1.0;

x3 = -.2948673863u;
y3 = -.2218777103u;
draw fullcircle scaled D shifted z3 withcolor RoyalBlue withpen pencircle scaled 1.0;

x4 = -.2553787574u;
y4 = .2168949154u;
draw fullcircle scaled D shifted z4 withcolor RoyalBlue withpen pencircle scaled 1.0;

x5 = -.4035222289u;
y5 = -.0925902722u;
draw fullcircle scaled D shifted z5 withcolor RoyalBlue withpen pencircle scaled 1.0;

x6 = -.1056146552u;
y6 = -.0702527012u;
draw fullcircle scaled D shifted z6 withcolor RoyalBlue withpen pencircle scaled 1.0;

x7 = .3104884939u;
y7 = .0741000000u;
draw fullcircle scaled D shifted z7 withcolor RoyalBlue withpen pencircle scaled 1.0;

x8 = -.1295908509u;
y8 = .3498902273u;
draw fullcircle scaled D shifted z8 withcolor RoyalBlue withpen pencircle scaled 1.0;

x9 = .0357437096u;
y9 = -.3644236080u;
draw fullcircle scaled D shifted z9 withcolor RoyalBlue withpen pencircle scaled 1.0;

x10 = .2195771750u;
y10 = .2709894837u;
draw fullcircle scaled D shifted z10 withcolor RoyalBlue withpen pencircle scaled 1.0;
draw (        0.0183326335u,        0.0614979716u)--(        0.0585384441u,       -0.1004655308u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1606840943u,       -0.2324312440u)--(        0.0585384441u,       -0.1004655308u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0183326335u,        0.0614979716u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0183326335u,        0.0614979716u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1606840943u,       -0.2324312440u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1606840943u,       -0.2324312440u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0585384441u,       -0.1004655308u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0585384441u,       -0.1004655308u) withcolor MidnightBlue;
draw (        0.0183326335u,        0.0614979716u)--(       -0.0911132091u,        0.1874753142u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.2553787574u,        0.2168949154u)--(       -0.0911132091u,        0.1874753142u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0183326335u,        0.0614979716u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0183326335u,        0.0614979716u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.2553787574u,        0.2168949154u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.2553787574u,        0.2168949154u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0911132091u,        0.1874753142u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0911132091u,        0.1874753142u) withcolor MidnightBlue;
draw (        0.0183326335u,        0.0614979716u)--(        0.0585053368u,       -0.1004737457u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.1056146552u,       -0.0702527012u)--(        0.0585053368u,       -0.1004737457u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0183326335u,        0.0614979716u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0183326335u,        0.0614979716u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.1056146552u,       -0.0702527012u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1056146552u,       -0.0702527012u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0585053368u,       -0.1004737457u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0585053368u,       -0.1004737457u) withcolor MidnightBlue;
draw (        0.0183326335u,        0.0614979716u)--(       -0.0910622463u,        0.1875195713u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.1295908509u,        0.3498902273u)--(       -0.0910622463u,        0.1875195713u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0183326335u,        0.0614979716u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0183326335u,        0.0614979716u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.1295908509u,        0.3498902273u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1295908509u,        0.3498902273u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0910622463u,        0.1875195713u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0910622463u,        0.1875195713u) withcolor MidnightBlue;
draw (        0.1606840943u,       -0.2324312440u)--(        0.2878000000u,       -0.3406000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1606840943u,       -0.2324312440u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1606840943u,       -0.2324312440u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2878000000u,       -0.3406000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2878000000u,       -0.3406000000u) withcolor Yellow;
draw (        0.1606840943u,       -0.2324312440u)--(       -0.0034321554u,       -0.2021898834u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.1056146552u,       -0.0702527012u)--(       -0.0034321554u,       -0.2021898834u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1606840943u,       -0.2324312440u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1606840943u,       -0.2324312440u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.1056146552u,       -0.0702527012u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1056146552u,       -0.0702527012u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0034321554u,       -0.2021898834u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0034321554u,       -0.2021898834u) withcolor MidnightBlue;
draw (        0.1606840943u,       -0.2324312440u)--(        0.0585015945u,       -0.1004940618u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.1056146552u,       -0.0702527012u)--(        0.0585015945u,       -0.1004940618u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1606840943u,       -0.2324312440u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1606840943u,       -0.2324312440u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.1056146552u,       -0.0702527012u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1056146552u,       -0.0702527012u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0585015945u,       -0.1004940618u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0585015945u,       -0.1004940618u) withcolor MidnightBlue;
draw (        0.1606840943u,       -0.2324312440u)--(       -0.0034356294u,       -0.2022087421u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0357437096u,       -0.3644236080u)--(       -0.0034356294u,       -0.2022087421u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1606840943u,       -0.2324312440u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1606840943u,       -0.2324312440u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0357437096u,       -0.3644236080u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0357437096u,       -0.3644236080u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0034356294u,       -0.2022087421u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0034356294u,       -0.2022087421u) withcolor MidnightBlue;
draw (       -0.2553787574u,        0.2168949154u)--(       -0.3816000000u,        0.3261000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.2553787574u,        0.2168949154u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.2553787574u,        0.2168949154u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.3816000000u,        0.3261000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.3816000000u,        0.3261000000u) withcolor Yellow;
draw (       -0.2553787574u,        0.2168949154u)--(       -0.0911071953u,        0.1875089120u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.1295908509u,        0.3498902273u)--(       -0.0911071953u,        0.1875089120u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.2553787574u,        0.2168949154u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.2553787574u,        0.2168949154u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.1295908509u,        0.3498902273u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1295908509u,        0.3498902273u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0911071953u,        0.1875089120u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0911071953u,        0.1875089120u) withcolor MidnightBlue;
draw (       -0.4035222289u,       -0.0925902722u)--(       -0.5469000000u,       -0.1780000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.4035222289u,       -0.0925902722u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.4035222289u,       -0.0925902722u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.5469000000u,       -0.1780000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.5469000000u,       -0.1780000000u) withcolor Yellow;
draw (       -0.4035222289u,       -0.0925902722u)--(       -0.5469000000u,       -0.0072000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.4035222289u,       -0.0925902722u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.4035222289u,       -0.0925902722u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.5469000000u,       -0.0072000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.5469000000u,       -0.0072000000u) withcolor Yellow;
draw (       -0.4035222289u,       -0.0925902722u)--(       -0.2601321317u,       -0.0072207032u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.1056146552u,       -0.0702527012u)--(       -0.2601321317u,       -0.0072207032u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.4035222289u,       -0.0925902722u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.4035222289u,       -0.0925902722u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.1056146552u,       -0.0702527012u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1056146552u,       -0.0702527012u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.2601321317u,       -0.0072207032u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.2601321317u,       -0.0072207032u) withcolor Green;
draw (       -0.1056146552u,       -0.0702527012u)--(       -0.0034659852u,       -0.2022160769u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0357437096u,       -0.3644236080u)--(       -0.0034659852u,       -0.2022160769u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.1056146552u,       -0.0702527012u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1056146552u,       -0.0702527012u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0357437096u,       -0.3644236080u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0357437096u,       -0.3644236080u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0034659852u,       -0.2022160769u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0034659852u,       -0.2022160769u) withcolor MidnightBlue;
draw (        0.3104884939u,        0.0741000000u)--(        0.4531000000u,        0.1608000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3104884939u,        0.0741000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3104884939u,        0.0741000000u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4531000000u,        0.1608000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4531000000u,        0.1608000000u) withcolor Yellow;
draw (       -0.1295908509u,        0.3498902273u)--(       -0.0456000000u,        0.4941000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.1295908509u,        0.3498902273u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1295908509u,        0.3498902273u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0456000000u,        0.4941000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0456000000u,        0.4941000000u) withcolor Yellow;
draw (       -0.1295908509u,        0.3498902273u)--(       -0.2136000000u,        0.4941000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.1295908509u,        0.3498902273u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1295908509u,        0.3498902273u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.2136000000u,        0.4941000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.2136000000u,        0.4941000000u) withcolor Yellow;
draw (        0.0357437096u,       -0.3644236080u)--(       -0.0483000000u,       -0.5086000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0357437096u,       -0.3644236080u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0357437096u,       -0.3644236080u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0483000000u,       -0.5086000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0483000000u,       -0.5086000000u) withcolor Yellow;
draw (        0.0357437096u,       -0.3644236080u)--(        0.1198000000u,       -0.5086000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0357437096u,       -0.3644236080u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0357437096u,       -0.3644236080u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1198000000u,       -0.5086000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1198000000u,       -0.5086000000u) withcolor Yellow;
% r =         0.1668792477
endfig;
end;
