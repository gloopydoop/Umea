function [Mel,Kel] = elMatD(o,h)
% ELMATD computes element matrices for square elements
%
%    [MEL,KEL] = ELMATD computes the mass MEL and stiffless KEL element
%    matrices for a first order element on a unit size element.
%
%    [MEL,KEL] = ELMATD(ORDER) computes the mass MEL and stiffless KEL 
%    element matrices for an element of order ORDER on a unit size
%    element.
%
%    [MEL,KEL] = ELMATD(ORDER,H) computes the mass MEL and stiffless KEL 
%    element matrices for an element of order ORDER on and size 
%    length H.
%

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
        o=2;
    end
end
    
    
if o>=3
    [Mel,Kel] = elMat(o,h);
    return;
end
if o==1
    Mel = [4 2 2 1; 2 4 1 2; 2 1 4 2; 1 2 2 4]*h^2/36;
    Kel = [4 -1 -1 -2; -1 4 -2 -1; -1 -2 4 -1; -2 -1 -1 4]/6;
end
if o==2
    Mel = ...
       [16     8    -4     8     4    -2    -4    -2     1
         8    64     8     4    32     4    -2   -16    -2
        -4     8    16    -2     4     8     1    -2    -4
         8     4    -2    64    32   -16     8     4    -2
         4    32     4    32   256    32     4    32     4
        -2     4     8   -16    32    64    -2     4     8
        -4    -2     1     8     4    -2    16     8    -4
        -2   -16    -2     4    32     4     8    64     8
         1    -2    -4    -2     4     8    -4     8    16]*h^2/900;
    Kel = ...
       [56   -18    -3   -18   -32    10    -3    10    -2
       -18   176   -18   -32   -96   -32    10     0    10
        -3   -18    56    10   -32   -18    -2    10    -3
       -18   -32    10   176   -96     0   -18   -32    10
       -32   -96   -32   -96   512   -96   -32   -96   -32
        10   -32   -18     0   -96   176    10   -32   -18
        -3    10    -2   -18   -32    10    56   -18    -3
        10     0    10   -32   -96   -32   -18   176   -18
        -2    10    -3    10   -32   -18    -3   -18    56]/90;
end