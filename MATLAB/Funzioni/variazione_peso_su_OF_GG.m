clear; clc; close all;
x=CEA('problem','tp','equilibrium','o/f', [0.3:0.1:5],'case','CEAM-TP-GG-F1','p,bar', 67.57,'t(k)', 1062,'reactants','fuel','RP-1','t(k)',293.15,'oxid','O2(L)','t(k)',90.15,'output','massf','transport','end','screen');


plot(x.input.problem.mixture.value, x.output.mw)



