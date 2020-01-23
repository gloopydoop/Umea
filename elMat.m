function [Mel,Kel]=elMat(o,h)
% ELMAT computes element matrices for square elements
%
%    [MEL,KEL] = ELMAT computes the mass MEL and stiffless KEL element
%    matrices for a first order element on a unit size element.
%
%    [MEL,KEL] = ELMAT(ORDER) computes the mass MEL and stiffless KEL 
%    element matrices for an element of order ORDER on a unit size
%    element.
%
%    [MEL,KEL] = ELMAT(ORDER,H) computes the mass MEL and stiffless KEL 
%    element matrices for an element of order ORDER on and size 
%    length H.
%
%    NOTE no error control is performed....
  
%Initial step is to compute the element matrices for the element
%with side 2 centered at the origin.
%
%Basis function are of the form
%[sum(n=0..o)a_nx^n] * [sum(m=0..o) a_m y^n]
 
if nargin<2
    h = 1;
    if (nargin<1)
        o=1;
    end
end

%  Gausspoints =[-sqrt(3/5) 0 sqrt(3/5)];
%  Gaussweights=[5/9 8/9 5/9];

Gausspoints =[-1/21*sqrt(245+14*sqrt(70)) -1/21*sqrt(245-14*sqrt(70)) ...
    0 ...
    1/21*sqrt(245-14*sqrt(70)) 1/21*sqrt(245+14*sqrt(70))];
Gaussweights=[1/900*(322-13*sqrt(70)) 1/900*(322+13*sqrt(70)) 128/225 ...
    1/900*(322+13*sqrt(70)) 1/900*(322-13*sqrt(70))];


xx=linspace(-1,1,o+1)';


for nn=0:o
    A(:,nn+1)=xx.^nn;
end
coeffs = zeros(o+1,o+1);
for nn=1:o+1
    b=zeros(o+1,1);
    b(nn)=1;
    coeffs(:,nn)=A\b;
end

%values contains the function values of the basis function at the
%gauss points.
%values(i,j) = the value of function i at point j...
tmpvalues =zeros(o+1,length(Gausspoints));
tmpdiff   =zeros(o+1,length(Gausspoints));
for nn=1:o+1
    for point=1:length(Gausspoints)
        x=Gausspoints(point);
        for mm=1:o+1
            tmpvalues(nn,point)=tmpvalues(nn,point)+coeffs(mm,nn)*x^(mm-1);
        end
        for mm=1:o
            tmpdiff(nn,point)=tmpdiff(nn,point)+coeffs(mm+1,nn)*mm*x^(mm-1);
        end
    end
end
values = zeros((o+1)^2,length(Gausspoints)^2);
diffx = zeros((o+1)^2,length(Gausspoints)^2);
for nnx=1:o+1
    for nny=1:o+1
        nn=nny+(nnx-1)*(o+1);
        for mmx=1:length(Gausspoints)
            for mmy=1:length(Gausspoints)
                mm=mmy+(mmx-1)*length(Gausspoints);
                values(nn,mm)=tmpvalues(nnx,mmx)*tmpvalues(nny,mmy);
                diffx(nn,mm) =tmpdiff(nnx,mmx)*tmpvalues(nny,mmy);
            end
        end
    end
end

Mel=zeros((o+1)^2);
Kxel=zeros((o+1)^2);
Kyel=zeros((o+1)^2);
for nn=1:(o+1)^2
    for mm=1:(o+1)^2
        for ii=1:length(Gausspoints)
            for jj=1:length(Gausspoints)
                Mel(nn,mm)=Mel(nn,mm)+Gaussweights(ii)*Gaussweights(jj)*...
                    values(nn,ii+length(Gausspoints)*(jj-1))*...
                    values(mm,ii+length(Gausspoints)*(jj-1));
                Kxel(nn,mm)=Kxel(nn,mm)+Gaussweights(ii)*Gaussweights(jj)*...
                    diffx(nn,ii+length(Gausspoints)*(jj-1))*...
                    diffx(mm,ii+length(Gausspoints)*(jj-1));
            end
        end
    end
end
perm=zeros((o+1)^2,1);
for nn=1:o+1
    for mm=1:o+1
        perm(nn+(mm-1)*(o+1))=mm+(nn-1)*(o+1);
    end
end
Kyel(perm,perm)=Kxel;

Kel = Kxel + Kyel;
Mel  = (h^2/4)*Mel;

