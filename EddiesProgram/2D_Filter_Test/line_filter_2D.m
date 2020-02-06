function [xf,B] = line_filter_2D(x,r)

filterBoxSize = 1+2*r;
B = ones(filterBoxSize);
xf = imfilter(x,B);
xf = xf(:);
