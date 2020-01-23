function [xf,B] = filterOct(x,r,param)

p = round( (2*r+1)/(2*sqrt(2)+2) );
nDiag = 1+r-p;
filterBoxSize = 1+2*r; %2*nDiag+nVert;

B = ones(filterBoxSize,filterBoxSize);
for nn = 1:nDiag-1
    B(nn,1:nDiag-nn) = 0;
    B(filterBoxSize+1-nn,1:nDiag-nn) = 0;
    B(nn,filterBoxSize-(0:nDiag-nn-1)) = 0;
    B(filterBoxSize+1-nn,filterBoxSize-(0:nDiag-nn-1)) = 0;
end
xf = filter2(B,reshape(x,param.nely,param.nelx));
xf = xf(:);