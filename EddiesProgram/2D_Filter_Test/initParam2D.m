function [param] = initParam2D(meshsize)
param.nelx = meshsize;
param.nely = 1;
param.nel = param.nelx*param.nely;
param.nnx = param.nelx+1;

param.ocSmooth = 1/2;  %OC smoothing parameter (default = 1/2)
penal1 = 1:0.5:3;
alpha1 = 10*ones(size(penal1));
alpha2 = 10.^(0.5:-0.5:-8);
penal2 = 3*ones(size(alpha2));
param.penal = [penal1 penal2];  %SIMP penalty paramter (default = 3)
param.alpha = [alpha1 alpha2];


param.coPower = 2;
%param.coPenal = 0;


param.saveresults = false;
param.plotDesign = true;

param.plotmod = 10;
