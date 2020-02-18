function [xf,B] = filterCirc(x,r,param)

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

if param.fill_expansion == true
    xf = filter2(B,reshape(x,param.nely_exp,param.nelx_exp));
end

xf = xf(:);