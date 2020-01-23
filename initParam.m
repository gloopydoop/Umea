function [param] = initParam(meshsize,rFactor,octFilter)
% INITPARAM Creates a structure containing geometric and disctetization
% information for topopt - heat conduction.
%
%    [PARAM] = INITPARAM returns information of the setup when computing on
%    a square discretized by 256 by 256 biquadratic elements.
%
%    [PARAM] = INITPARAM(MESHSIZE) returns information of the setup when
%    computing on a square discretized by N by N bilinear elements,
%    where N = 256*MESHSIZE.
%
%    [PARAM] = INITINFO(MESHSIZE,RFACTOR) returns information of the setup
%    when computing on a square discretized by N by N bilinear elements,
%    where N = 256*MESHSIZE. RFACTOR is used to determing the radius R of
%    each filter component R = RFACTOR*MESHSIZE.
%
%    [PARAM] = INITINFO(MESHSIZE,ORDER,OCTFILTER) returns information of
%    the setup when computing on a square discretized by N by N bilinear
%    elements, where N = 256*MESHSIZE. RFACTOR is used to determing the 
%    radius R of each filter component R = RFACTOR*MESHSIZE. OCFILTER is a
%    binary variable that specifies of the filtering should be performed
%    using octagonal neighborhoods (OCTFITLER=TRUE) of circular
%    neighborhoods (OCTFILTER=FALSE).
%
%    NOTE componentwise basis function...
%
%    NOTE no error control is performed....
if nargin<3
    octFilter = 0;
    if nargin<2
        rFactor = 2;
        if nargin<1
            meshsize = 1;
        end
    end
end

param.meshsize = meshsize;   %meshsize
param.rFactor = rFactor;     %factor used to determing filter radius...
param.octFilter = octFilter; %if true then use octagonal neighborhoods else
                             %use circular neighborhoods for the filtering.
param.order = 1;             %element order in the FE discretization

% THIS IS ME BEING CHEEKY!
param.nelx = 16*meshsize;
param.nely = 16*meshsize;
param.nel = param.nelx*param.nely;

%number of nodes in the left pipe
param.nnx = param.nelx*param.order+1;
param.nny = param.nely*param.order+1;
param.nn  = param.nnx*param.nny;

%% harry edit
%Fix 1/8 of the left hand boundary...%(zero-Dirichlet condition)
%Note offset by one due to 1-based indexing of arrays
%param.fixedNodes = 1+((7*param.nely)/16:(9*param.nely)/16);
%param.comNodes = setdiff(1:param.nn,param.fixedNodes);

%param.nUnknowns = length(param.comNodes);

%Work in matrix form first
BCs = zeros(param.nely,param.nelx);
BCs(param.nelx/2,:) = 1;

param.fixedNodes = BCsToIndex(BCs);

param.fixedNodes = param.nely*8 + 1+((7*param.nely)/16:(9*param.nely)/16);


param.comNodes = setdiff(1:param.nn,param.fixedNodes);

param.nUnknowns = length(param.comNodes);

param.weakMaterial = 0.001;


param.ocSmooth = 1/2;  %OC smoothing parameter (default = 1/2)
penal1 = 1:0.5:3;
alpha1 = 10*ones(size(penal1));
alpha2 = 10.^(0.5:-0.5:-8);
penal2 = 3*ones(size(alpha2));
param.penal = [penal1 penal2];  %SIMP penalty paramter (default = 3)
param.alpha = [alpha1 alpha2];

%param.alpha = [10 alpha2];
%param.penal = ones(size(param.alpha));

% %If needed add explicit penalty to cover difference between open and close
% %filter...
param.coPower = 2;
%param.coPenal = 0;


param.saveresults = false;
param.plotDesign = true;

param.plotmod = 10;
