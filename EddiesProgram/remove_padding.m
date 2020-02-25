function x = remove_padding(xf,param)
pad = param.rFactor + 1;
if param.boundary_treatment == 1
    xf = reshape(xf,param.sizey,param.sizex);
    x = xf(pad+1:end-pad,pad+1:end-pad);
    x = x(:);
else
    x = xf;
end