prologues := 3;
outputformat   := "svg";
outputtemplate := "%j.%o";
input mpcolornames;
beginfig(8);
u := 20cm;
path pols[];
pair p[];
R := .2041845240u;
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

x1 = .1184000001u;
y1 = -.3919036668u;
draw fullcircle scaled D shifted z1 withcolor RoyalBlue withpen pencircle scaled 1.0;

x2 = .0962252592u;
y2 = -.2072626634u;
draw fullcircle scaled D shifted z2 withcolor RoyalBlue withpen pencircle scaled 1.0;

x3 = -.2106469655u;
y3 = .3738114510u;
draw fullcircle scaled D shifted z3 withcolor RoyalBlue withpen pencircle scaled 1.0;

x4 = .3189838238u;
y4 = .1436735242u;
draw fullcircle scaled D shifted z4 withcolor RoyalBlue withpen pencircle scaled 1.0;

x5 = -.3789314635u;
y5 = -.2272213821u;
draw fullcircle scaled D shifted z5 withcolor RoyalBlue withpen pencircle scaled 1.0;

x6 = -.4591503219u;
y6 = .0732194052u;
draw fullcircle scaled D shifted z6 withcolor RoyalBlue withpen pencircle scaled 1.0;

x7 = -.1782399033u;
y7 = .0129234810u;
draw fullcircle scaled D shifted z7 withcolor RoyalBlue withpen pencircle scaled 1.0;

x8 = .1014304136u;
y8 = .1725260599u;
draw fullcircle scaled D shifted z8 withcolor RoyalBlue withpen pencircle scaled 1.0;
draw (        0.1184000001u,       -0.3919036668u)--(       -0.0510000000u,       -0.5059000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1184000001u,       -0.3919036668u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1184000001u,       -0.3919036668u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0510000000u,       -0.5059000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0510000000u,       -0.5059000000u) withcolor Yellow;
draw (        0.1184000001u,       -0.3919036668u)--(        0.2877993974u,       -0.2779073334u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0962252592u,       -0.2072626634u)--(        0.2877993974u,       -0.2779073334u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1184000001u,       -0.3919036668u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1184000001u,       -0.3919036668u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0962252592u,       -0.2072626634u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0962252592u,       -0.2072626634u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2877993974u,       -0.2779073334u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2877993974u,       -0.2779073334u) withcolor Green;
draw (        0.0962252592u,       -0.2072626634u)--(       -0.1058512840u,       -0.1779985125u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.1782399033u,        0.0129234810u)--(       -0.1058512840u,       -0.1779985125u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0962252592u,       -0.2072626634u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0962252592u,       -0.2072626634u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.1782399033u,        0.0129234810u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1782399033u,        0.0129234810u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.1058512840u,       -0.1779985125u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1058512840u,       -0.1779985125u) withcolor Green;
draw (        0.0962252592u,       -0.2072626634u)--(        0.0238366399u,       -0.0163406698u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.1782399033u,        0.0129234810u)--(        0.0238366399u,       -0.0163406698u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0962252592u,       -0.2072626634u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0962252592u,       -0.2072626634u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.1782399033u,        0.0129234810u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1782399033u,        0.0129234810u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0238366399u,       -0.0163406698u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0238366399u,       -0.0163406698u) withcolor MidnightBlue;
draw (        0.0962252592u,       -0.2072626634u)--(        0.0238370270u,       -0.0163405231u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1014304136u,        0.1725260599u)--(        0.0238370270u,       -0.0163405231u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0962252592u,       -0.2072626634u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0962252592u,       -0.2072626634u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1014304136u,        0.1725260599u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1014304136u,        0.1725260599u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0238370270u,       -0.0163405231u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0238370270u,       -0.0163405231u) withcolor MidnightBlue;
draw (       -0.2106469655u,        0.3738114510u)--(       -0.0456000000u,        0.4941000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.2106469655u,        0.3738114510u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.2106469655u,        0.3738114510u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.0456000000u,        0.4941000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.0456000000u,        0.4941000000u) withcolor Yellow;
draw (       -0.2106469655u,        0.3738114510u)--(       -0.2882462599u,        0.1849472953u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.4591503219u,        0.0732194052u)--(       -0.2882462599u,        0.1849472953u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.2106469655u,        0.3738114510u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.2106469655u,        0.3738114510u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.4591503219u,        0.0732194052u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.4591503219u,        0.0732194052u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.2882462599u,        0.1849472953u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.2882462599u,        0.1849472953u) withcolor MidnightBlue;
draw (       -0.2106469655u,        0.3738114510u)--(       -0.1006468944u,        0.2017902192u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.1782399033u,        0.0129234810u)--(       -0.1006468944u,        0.2017902192u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.2106469655u,        0.3738114510u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.2106469655u,        0.3738114510u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.1782399033u,        0.0129234810u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1782399033u,        0.0129234810u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.1006468944u,        0.2017902192u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1006468944u,        0.2017902192u) withcolor MidnightBlue;
draw (       -0.2106469655u,        0.3738114510u)--(       -0.2882399744u,        0.1849447128u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.1782399033u,        0.0129234810u)--(       -0.2882399744u,        0.1849447128u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.2106469655u,        0.3738114510u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.2106469655u,        0.3738114510u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.1782399033u,        0.0129234810u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1782399033u,        0.0129234810u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.2882399744u,        0.1849447128u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.2882399744u,        0.1849447128u) withcolor MidnightBlue;
draw (       -0.2106469655u,        0.3738114510u)--(       -0.1006460502u,        0.2017907591u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1014304136u,        0.1725260599u)--(       -0.1006460502u,        0.2017907591u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.2106469655u,        0.3738114510u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.2106469655u,        0.3738114510u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1014304136u,        0.1725260599u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1014304136u,        0.1725260599u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.1006460502u,        0.2017907591u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1006460502u,        0.2017907591u) withcolor MidnightBlue;
draw (        0.3189838238u,        0.1436735242u)--(        0.4504000000u,       -0.0126000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3189838238u,        0.1436735242u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3189838238u,        0.1436735242u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4504000000u,       -0.0126000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4504000000u,       -0.0126000000u) withcolor Yellow;
draw (        0.3189838238u,        0.1436735242u)--(        0.1875684768u,       -0.0125997186u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1014304136u,        0.1725260599u)--(        0.1875684768u,       -0.0125997186u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3189838238u,        0.1436735242u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3189838238u,        0.1436735242u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1014304136u,        0.1725260599u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1014304136u,        0.1725260599u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1875684768u,       -0.0125997186u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1875684768u,       -0.0125997186u) withcolor Green;
draw (        0.3189838238u,        0.1436735242u)--(        0.2328457606u,        0.3287993027u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1014304136u,        0.1725260599u)--(        0.2328457606u,        0.3287993027u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3189838238u,        0.1436735242u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3189838238u,        0.1436735242u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1014304136u,        0.1725260599u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1014304136u,        0.1725260599u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2328457606u,        0.3287993027u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2328457606u,        0.3287993027u) withcolor Green;
draw (       -0.3789314635u,       -0.2272213821u)--(       -0.2109000000u,       -0.3433000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.3789314635u,       -0.2272213821u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.3789314635u,       -0.2272213821u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.2109000000u,       -0.3433000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.2109000000u,       -0.3433000000u) withcolor Yellow;
draw (       -0.4591503219u,        0.0732194052u)--(       -0.2882432164u,        0.1849426396u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (       -0.1782399033u,        0.0129234810u)--(       -0.2882432164u,        0.1849426396u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.4591503219u,        0.0732194052u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.4591503219u,        0.0732194052u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.1782399033u,        0.0129234810u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1782399033u,        0.0129234810u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.2882432164u,        0.1849426396u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.2882432164u,        0.1849426396u) withcolor MidnightBlue;
draw (       -0.1782399033u,        0.0129234810u)--(       -0.1006461714u,        0.2017899222u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1014304136u,        0.1725260599u)--(       -0.1006461714u,        0.2017899222u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.1782399033u,        0.0129234810u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1782399033u,        0.0129234810u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1014304136u,        0.1725260599u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1014304136u,        0.1725260599u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (       -0.1006461714u,        0.2017899222u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1006461714u,        0.2017899222u) withcolor MidnightBlue;
draw (       -0.1782399033u,        0.0129234810u)--(        0.0238366817u,       -0.0163403812u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1014304136u,        0.1725260599u)--(        0.0238366817u,       -0.0163403812u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (       -0.1782399033u,        0.0129234810u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (       -0.1782399033u,        0.0129234810u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1014304136u,        0.1725260599u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1014304136u,        0.1725260599u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0238366817u,       -0.0163403812u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0238366817u,       -0.0163403812u) withcolor MidnightBlue;
% r =         0.2041845240
endfig;
end;
