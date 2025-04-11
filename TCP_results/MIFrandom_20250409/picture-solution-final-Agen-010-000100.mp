prologues := 3;
outputformat   := "svg";
outputtemplate := "%j.%o";
input mpcolornames;
beginfig(10);
u := 20cm;
path pols[];
pair p[];
R := .1668793128u;
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

x1 = .3104885285u;
y1 = .0740999999u;
draw fullcircle scaled D shifted z1 withcolor RoyalBlue withpen pencircle scaled 1.0;

x2 = -.4035221824u;
y2 = -.0925903216u;
draw fullcircle scaled D shifted z2 withcolor RoyalBlue withpen pencircle scaled 1.0;

x3 = .2385323142u;
y3 = .2259305859u;
draw fullcircle scaled D shifted z3 withcolor RoyalBlue withpen pencircle scaled 1.0;

x4 = .0357437414u;
y4 = -.3644235775u;
draw fullcircle scaled D shifted z4 withcolor RoyalBlue withpen pencircle scaled 1.0;

x5 = -.2328283159u;
y5 = -.2896953748u;
draw fullcircle scaled D shifted z5 withcolor RoyalBlue withpen pencircle scaled 1.0;

x6 = .0183327558u;
y6 = .0614980435u;
draw fullcircle scaled D shifted z6 withcolor RoyalBlue withpen pencircle scaled 1.0;

x7 = .1606841788u;
y7 = -.2324312843u;
draw fullcircle scaled D shifted z7 withcolor RoyalBlue withpen pencircle scaled 1.0;

x8 = -.1295908973u;
y8 = .3498902182u;
draw fullcircle scaled D shifted z8 withcolor RoyalBlue withpen pencircle scaled 1.0;

x9 = -.1056145309u;
y9 = -.0702525382u;
draw fullcircle scaled D shifted z9 withcolor RoyalBlue withpen pencircle scaled 1.0;

x10 = -.2553787511u;
y10 = .2168950289u;
draw fullcircle scaled D shifted z10 withcolor RoyalBlue withpen pencircle scaled 1.0;
draw (        0.3104885285u,        0.0740999999u)--(        0.4531000000u,        0.1608000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3104885285u,        0.0740999999u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3104885285u,        0.0740999999u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4531000000u,        0.1608000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4531000000u,        0.1608000000u) withcolor Yellow;
draw (       -0.4035221824u,       -0.0925903216u)--(       -0.5469000000u,       -0.1780000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.4035221824u,       -0.0925903216u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.4035221824u,       -0.0925903216u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.5469000000u,       -0.1780000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.5469000000u,       -0.1780000000u) withcolor Yellow;
draw (       -0.4035221824u,       -0.0925903216u)--(       -0.5469000000u,       -0.0072000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.4035221824u,       -0.0925903216u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.4035221824u,       -0.0925903216u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.5469000000u,       -0.0072000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.5469000000u,       -0.0072000000u) withcolor Yellow;
draw (       -0.4035221824u,       -0.0925903216u)--(       -0.2601321013u,       -0.0072205981u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.1056145309u,       -0.0702525382u)--(       -0.2601321013u,       -0.0072205981u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.4035221824u,       -0.0925903216u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.4035221824u,       -0.0925903216u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.1056145309u,       -0.0702525382u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1056145309u,       -0.0702525382u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.2601321013u,       -0.0072205981u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.2601321013u,       -0.0072205981u) withcolor Green;
draw (        0.0357437414u,       -0.3644235775u)--(       -0.0483000000u,       -0.5086000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0357437414u,       -0.3644235775u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0357437414u,       -0.3644235775u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0483000000u,       -0.5086000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0483000000u,       -0.5086000000u) withcolor Yellow;
draw (        0.0357437414u,       -0.3644235775u)--(        0.1198000000u,       -0.5086000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0357437414u,       -0.3644235775u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0357437414u,       -0.3644235775u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1198000000u,       -0.5086000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1198000000u,       -0.5086000000u) withcolor Yellow;
draw (        0.0357437414u,       -0.3644235775u)--(       -0.0034355852u,       -0.2022086417u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1606841788u,       -0.2324312843u)--(       -0.0034355852u,       -0.2022086417u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0357437414u,       -0.3644235775u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0357437414u,       -0.3644235775u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1606841788u,       -0.2324312843u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1606841788u,       -0.2324312843u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0034355852u,       -0.2022086417u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0034355852u,       -0.2022086417u) withcolor MidnightBlue;
draw (        0.0357437414u,       -0.3644235775u)--(       -0.0034657872u,       -0.2022159393u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.1056145309u,       -0.0702525382u)--(       -0.0034657872u,       -0.2022159393u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0357437414u,       -0.3644235775u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0357437414u,       -0.3644235775u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.1056145309u,       -0.0702525382u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1056145309u,       -0.0702525382u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0034657872u,       -0.2022159393u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0034657872u,       -0.2022159393u) withcolor MidnightBlue;
draw (        0.0183327558u,        0.0614980435u)--(        0.0585384400u,       -0.1004655574u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1606841788u,       -0.2324312843u)--(        0.0585384400u,       -0.1004655574u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0183327558u,        0.0614980435u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0183327558u,        0.0614980435u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1606841788u,       -0.2324312843u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1606841788u,       -0.2324312843u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0585384400u,       -0.1004655574u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0585384400u,       -0.1004655574u) withcolor MidnightBlue;
draw (        0.0183327558u,        0.0614980435u)--(       -0.0910624318u,        0.1875194622u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.1295908973u,        0.3498902182u)--(       -0.0910624318u,        0.1875194622u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0183327558u,        0.0614980435u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0183327558u,        0.0614980435u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.1295908973u,        0.3498902182u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1295908973u,        0.3498902182u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0910624318u,        0.1875194622u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0910624318u,        0.1875194622u) withcolor MidnightBlue;
draw (        0.0183327558u,        0.0614980435u)--(        0.0585055001u,       -0.1004737308u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.1056145309u,       -0.0702525382u)--(        0.0585055001u,       -0.1004737308u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0183327558u,        0.0614980435u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0183327558u,        0.0614980435u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.1056145309u,       -0.0702525382u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1056145309u,       -0.0702525382u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0585055001u,       -0.1004737308u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0585055001u,       -0.1004737308u) withcolor MidnightBlue;
draw (        0.0183327558u,        0.0614980435u)--(       -0.0911131363u,        0.1874754294u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.2553787511u,        0.2168950289u)--(       -0.0911131363u,        0.1874754294u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0183327558u,        0.0614980435u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0183327558u,        0.0614980435u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.2553787511u,        0.2168950289u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.2553787511u,        0.2168950289u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0911131363u,        0.1874754294u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0911131363u,        0.1874754294u) withcolor MidnightBlue;
draw (        0.1606841788u,       -0.2324312843u)--(        0.2878000000u,       -0.3406000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1606841788u,       -0.2324312843u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1606841788u,       -0.2324312843u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2878000000u,       -0.3406000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2878000000u,       -0.3406000000u) withcolor Yellow;
draw (        0.1606841788u,       -0.2324312843u)--(       -0.0034321288u,       -0.2021898785u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.1056145309u,       -0.0702525382u)--(       -0.0034321288u,       -0.2021898785u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1606841788u,       -0.2324312843u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1606841788u,       -0.2324312843u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.1056145309u,       -0.0702525382u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1056145309u,       -0.0702525382u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0034321288u,       -0.2021898785u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0034321288u,       -0.2021898785u) withcolor MidnightBlue;
draw (        0.1606841788u,       -0.2324312843u)--(        0.0585017767u,       -0.1004939440u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.1056145309u,       -0.0702525382u)--(        0.0585017767u,       -0.1004939440u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1606841788u,       -0.2324312843u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1606841788u,       -0.2324312843u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.1056145309u,       -0.0702525382u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1056145309u,       -0.0702525382u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0585017767u,       -0.1004939440u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0585017767u,       -0.1004939440u) withcolor MidnightBlue;
draw (       -0.1295908973u,        0.3498902182u)--(       -0.0456000000u,        0.4941000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.1295908973u,        0.3498902182u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1295908973u,        0.3498902182u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0456000000u,        0.4941000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0456000000u,        0.4941000000u) withcolor Yellow;
draw (       -0.1295908973u,        0.3498902182u)--(       -0.2136000000u,        0.4941000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.1295908973u,        0.3498902182u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1295908973u,        0.3498902182u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.2136000000u,        0.4941000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.2136000000u,        0.4941000000u) withcolor Yellow;
draw (       -0.1295908973u,        0.3498902182u)--(       -0.0911071530u,        0.1875088570u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.2553787511u,        0.2168950289u)--(       -0.0911071530u,        0.1875088570u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.1295908973u,        0.3498902182u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1295908973u,        0.3498902182u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.2553787511u,        0.2168950289u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.2553787511u,        0.2168950289u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0911071530u,        0.1875088570u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0911071530u,        0.1875088570u) withcolor MidnightBlue;
draw (       -0.2553787511u,        0.2168950289u)--(       -0.3816000000u,        0.3261000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.2553787511u,        0.2168950289u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.2553787511u,        0.2168950289u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.3816000000u,        0.3261000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.3816000000u,        0.3261000000u) withcolor Yellow;
% r =         0.1668793128
endfig;
end;
