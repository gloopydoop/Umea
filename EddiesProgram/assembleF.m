function [f] = assembleF(param)
% ASSEMBLEF Assembles matrix for heat conduction problem
%
%    [F] = ASSEMBLE(PARAM) assembles the right hand side vector F
%    corresponding to the applied forces (the applied heating) for the
%    problem size and element order defined by PARAM.

if nargin<1
    error('Error in assembleF; this function requires one input argument.');
end

ipos = zeros(1,param.nel*(param.order+1)^4);
jpos = zeros(size(ipos));
Mval = zeros(size(ipos));
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

[Mel,~] = elMat(param.order);
Melv = reshape(Mel,1,(param.order+1)^4);
for nn=1:param.nelx
    for mm=1:param.nely
        % HARRY: you are only considering things in the BC
        if param.nullel(param.nelx,param.nely) == 0
            el = (nn-1)*param.nely+mm;
            pos = ((el-1)*(param.order+1)^4 + 1) : (el*(param.order+1)^4);
            stpos = (nn-1)*param.order*param.nny + (mm-1)*param.order+1;
            ipos(pos) = stpos + niList;
            jpos(pos) = stpos + njList;
            Mval(pos) = Melv;
        end
    end
end
M = sparse(ipos,jpos,Mval);

%NEED TO FIX FORCE FOR THIS PROBLEM...
f = M*ones(param.nn,1)/param.nel;
% This is a cheeky Harry edit: