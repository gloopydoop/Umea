function [xf,B] = filterCirc(x,r,param)
line_BCs = reshape(1-param.nullel,param.nely,param.nelx);
% HARRY BEING SNEAKY AGAIN!
if param.boundary_treatment == 1
    pad = param.rFactor + 1;
    sizex = param.nelx +2*pad;
    sizey = param.nely +2*pad;
elseif param.boundary_treatment == 2
    x(find(line_BCs == 0)) = 0;
    sizex = param.nelx;
    sizey = param.nely;
end


filterBoxSize = 1+2*floor(r); %2*nDiag+nVert;

B = ones(filterBoxSize,filterBoxSize);
mid = floor(r)+1;
for nn=1:filterBoxSize
    for mm = 1:filterBoxSize
        if norm([(nn-mid), (mm-mid)],2)>=r
            B(nn,mm) = 0;
        end
    end
end

xf = filter2(B,reshape(x,sizey,sizex));

% again being sneaky
xf = xf(:);

if param.boundary_treatment == 2
    xf(find(line_BCs == 0)) = 0;
end