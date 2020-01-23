%OCTAGONALFILTER two-dimensinoal filtering with octagonal neighbourhood
%
%   AF = OCTAGONALFILTER(A,R) computes AF, a moving sum over an octagonal 
%   shaped neighbourhood of radius R for a two dimensional input matrix A.
%   More precisely, the neighbours of position (i,j) are all positions
%   (k,l) such that norm([k-i,l-j],inf)<=R and norm([k-i,l-j],1)<=R+Q,
%   where the constant Q = round((sqrt(2)*(R+1)-1)/(sqrt(2)+2)). That is,
%   the octagon approximates a circle with radius R. The dimensions of the
%   output AF is the same as the dimensions of the input A.
%   
%   AF = OCTAGONALFILTER(A,R,S) computes AF, a moving sum over an octagonal 
%   shaped neighbourhood of radius R for a two dimensional input matrix A.
%   More precisely, the neighbours of position (i,j) are all positions
%   (k,l) such that norm([k-i,l-j],inf)<=R and norm([k-i,l-j],1)<=R+Q,
%   where the constant Q = round((sqrt(2)*(R+1)-1)/(sqrt(2)+2)). That is,
%   the octagon approximates a circle with radius R. The input matrix A is
%   stored in column (row) major format and the input parameter S is the
%   major-stride. That is, S is the distance between columns (rows). The
%   dimensions of the output AF is NM x 1.
%
%   To compile the underlying C-code issue the command
%   mex octagonalfilter.c -lmwblas -largeArrayDims
%   in the MATLAB prompt.
%