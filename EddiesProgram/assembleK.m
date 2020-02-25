function [K] = assembleK(param,alpha)
% ASSEMBLE Assembles stiffness matrix for heat conduction problem
%
%    [K] = ASSEMBLE(PARAM) assembles the stiffness matrix K for the problem
%    size and element order defined by PARAM, with the domain filled with high
%    conductivity material...
%
%    [K] = ASSEMBLE(PARAM,ALPHA) assembles the stiffness matrix K for the
%    problem size and element order defined by PARAM, with the domain filled
%    with material having conductivity defined by ALPHA.
%
%    NOTE essentially no error control is performed....

if nargin<2
    if nargin<1
        error('Error in assembleK; this function requires one or two input arguments.');
    end
    alpha = ones(param.nel,1);
end

ipos = zeros(1,param.nel*(param.order+1)^4);
jpos = zeros(size(ipos));
Kval = zeros(size(ipos));

nList = zeros(1,(param.order+1)^2);
for ii=0:param.order
    for jj=0:param.order
        nList(1+jj+(param.order+1)*ii) = jj + ii*param.nny;
    end
end
niList = reshape(repmat(nList,1,(param.order+1)^2), ...
    1,(param.order+1)^4);
njList = reshape(repmat(nList,(param.order+1)^2,1),  ...
    1,(param.order+1)^4);

[~,Kel] = elMatD(param.order);
Kelv = reshape(Kel,1,(param.order+1)^4);
for nn=1:param.nelx
    for mm=1:param.nely
        % HARRY: you are only considering things in the BC
        if param.nullel(param.nelx,param.nely) == 0
            el = (nn-1)*param.nely+mm;
            pos = ((el-1)*(param.order+1)^4 + 1) : (el*(param.order+1)^4);
            stpos = (nn-1)*param.order*param.nny + (mm-1)*param.order+1;
            ipos(pos) = stpos + niList;
            jpos(pos) = stpos + njList;
            Kval(pos) = alpha(el)*Kelv;
        end
    end
end
K = sparse(ipos,jpos,Kval);