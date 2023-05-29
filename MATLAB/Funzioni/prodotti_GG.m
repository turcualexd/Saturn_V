clear, clc, close all;

ch4     = 0.15220;
co      = 0.28364;
co2     = 0.06268;
c2h4    = 0.00003;
c2h6    = 0.00011;
h2      = 0.05018;
h2o     = 0.09705;
c       = 0.35410;

tot1 = ch4 + co + co2 + c2h4 + c2h6 + h2 + h2o + c;
tot2 = tot1 - c;
r = tot1/tot2;

ch4 = ch4 * r
co = co * r
co2 = co2 * r
c2h4 = c2h4 * r
c2h6 = c2h6 * r
h2 = h2 * r
h2o = h2o * r

test = tot1 - (ch4 + co + co2 + c2h4 + c2h6 + h2 + h2o)