function [obj,dcda] = getobj(param,u,f)
%GETOBJ computes the value and gradient of the objective function 
%
%    [OBJ,GRAD] = GETOBJ(PARAM,U,F) computes OBJ, the value of the
%    objective function as well as GRAD, the gradient of the objective
%    function with respect to elementwise changes in the design. The input
%    parameters are: PARAM, a structure that holds the parameters defining
%    the problem; U, the solution of the state equation; and F, the forcing
%    of the sytem.
%

if nargin<3
    error('Error in getobj; this function requires three input arguments.');
end

nList = zeros(1,(param.order+1)^2);
for ii=0:param.order
    for jj=0:param.order
        nList(1+jj+(param.order+1)*ii) = jj + ii*param.nny;
    end
end

[~,Kel] = elMatD(param.order);
obj = f'*u;
if nargout>1
    dcda = zeros(param.nel,1);
    for nn=1:param.nelx
        for mm=1:param.nely
            el = (nn-1)*param.nely+mm;
            stpos = (nn-1)*param.order*param.nny + (mm-1)*param.order+1;
            dcda(el) = -transpose(u(stpos+nList))*Kel*u(stpos+nList);
        end
    end
end
