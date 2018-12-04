%madras validation of model
a = 40.5;
b = 2.39;
t = 98;
l = 14;
s = 1;
m = 18000;
f = 24; 
inf = 17; %from Sackett and Mundt 2005
size = 60;
n = 89;
[pt,ut,vt,remt,dist] = epidemic_SLIR_loc(n, t, l, s, 5, m, f, .062, 1,inf , 1,a,b,size); 
save('madcheck_block','pt','ut','vt','remt','dist');
