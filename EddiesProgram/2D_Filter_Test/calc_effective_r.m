function [start,stop] = calc_effective_r(x,y,eps)
window = x(y< 1- eps && y > eps);
R_eff = window(end) - window(1):