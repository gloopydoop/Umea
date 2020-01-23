function [nmma,mmma,iter,xmin,xmax,fval,dfdx, ...
    low,upp,a0mma,amma,cmma,dmma] = initmma(param)

if nargin<1
    error('Error in initmma.. Requires one input argument');
end

nmma = param.nel;
mmma = 2; %compliance f_0 note put penalty in f_1

cmma = ones(mmma,1);
dmma = zeros(mmma,1);
    
iter = 0;
xmin = zeros(nmma,1);
xmax = ones(nmma,1);
low = xmin;
upp = xmax;
a0mma = 1;
amma = zeros(mmma,1);

fval = zeros(mmma,1);
dfdx = zeros(mmma,nmma);
