function [xf,B] = filterRect(x,r,param)

filterBoxSize = 1+2*r;

B = ones(filterBoxSize,filterBoxSize);

xf = filter2(B,reshape(x,param.nely,param.nelx));
xf = xf(:);
