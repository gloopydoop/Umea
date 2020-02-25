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
param.nelx = 128*meshsize;
param.nely = 128*meshsize;
param.nel = param.nelx*param.nely;

%number of nodes in the left pipe
param.nnx = param.nelx*param.order+1;
param.nny = param.nely*param.order+1;
param.nn  = param.nnx*param.nny;

% HARRY. FILT_TYPE = 1 FOR OPEN. 2 FOR CLOSE
param.filt_type = 1;

% HARRY. 
% 1) for padding with edge values (edge values = 0)
% 2) for ignoring boundaires 
% 3) Potentially "extend" according to Wikipedia
param.boundary_treatment = 1;


% this is my sneaky padding:
if param.boundary_treatment == 1
    param.pad = param.rFactor + 1;
    pad = param.pad;
    param.sizex = param.nelx +2*pad;
    param.sizey = param.nely +2*pad;
else
    param.sizex = param.nelx;
    param.sizey = param.nely;
end


%% harry edit
%Fix 1/8 of the left hand boundary...%(zero-Dirichlet condition)
%Note offset by one due to 1-based indexing of arrays
%param.fixedNodes = 1+((7*param.nely)/16:(9*param.nely)/16);
%param.comNodes = setdiff(1:param.nn,param.fixedNodes);

%param.nUnknowns = length(param.comNodes);

%Work in matrix form first
param.BCs = zeros(param.nny,param.nnx);
param.BCs(:,param.nnx) = 1;
param.BCs(param.nny,:) = 1;
% %define square as: square_ratio
 param.square_ratio = 0.3;
 edge = floor(param.nny*(1-param.square_ratio)/2);
% 
 param.BCs(edge:param.nny-edge+1,edge:param.nny-edge+1) = 1;

param.nullel = zeros(param.nelx,param.nely);
param.nullel(edge:param.nelx-edge+1,edge:param.nely-edge+1) = 1;
param.fixedEl = BCsToIndex(param.nullel);
param.comEl = setdiff(1:param.nel,param.fixedEl);


param.fixedNodes = BCsToIndex(param.BCs);

%param.fixedNodes = param.nely*8 + 1+((7*param.nely)/16:(9*param.nely)/16);



param.comNodes = setdiff(1:param.nn,param.fixedNodes);

param.nUnknowns = length(param.comNodes);

param.weakMaterial = 0.001;


param.ocSmooth = 1/2;  %OC smoothing parameter (default = 1/2)

%% My changes to this
penal1 = 1:0.5:3;
alpha1 = 10*ones(size(penal1));
alpha2 = 10.^(1:-1:-5);
penal2 = 3*ones(size(alpha2));
param.penal = [penal1,penal2];  %SIMP penalty paramter (default = 3)
param.alpha = [alpha1,alpha2];

% penal1 = 1:0.5:3;
% alpha1 = 10*ones(size(penal1));
% alpha2 = 10.^(0.5:-0.5:-8);
% penal2 = 3*ones(size(alpha2));
% param.penal = [penal1 penal2];  %SIMP penalty paramter (default = 3)
% param.alpha = [alpha1 alpha2];

%param.alpha = [10 alpha2];
%param.penal = ones(size(param.alpha));

% %If needed add explicit penalty to cover difference between open and close
% %filter...
param.coPower = 2;
%param.coPenal = 0;


param.saveresults = false;
param.plotDesign = false;

param.plotmod = 20;
