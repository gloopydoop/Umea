function relObj = plotResults(meshsize,rFactor,volfrac)
%PLOTRESULTS plots results from the heat conduction problem
%
%    PLOTRESULTS(MESHSIZE,RFACTOR) plots the resulting designs obtained by
%    solving the minimal heat compliance problem using the input arguments
%    MESHSIZE and RFACTOR. For more information about the parameter
%    settings, please see HELP INITPARAM.
%    
%    The default values of the input parameters are:
%    MESHSIZE = 1
%    RFACTOR = 2
if nargin<3
    volfrac = 0.5;
    if nargin<2
        rFactor = 2;
        if nargin<1
            meshsize = 1;
        end
    end
end

closedCirc = false;
volfrac = round(100*volfrac);

folderName = sprintf('./ResultsC');
if meshsize < 10 
    nametype = 1;
    filename = sprintf('%s/Mesh%1dV%02dR%03dfinal.mat',...
        folderName,meshsize,volfrac,round(10*rFactor));
elseif meshsize < 100
    nametype = 2;
    meshsize = meshsize - 10;
    filename = sprintf('%s/sMesh%1dV%02dR%03dfinal.mat',...
        folderName,meshsize,volfrac,round(10*rFactor));
else
    nametype = 3;
    meshsize = meshsize - 100;
    filename = sprintf('%s/dfMesh%1dR%03dfinal.mat',...
        folderName,meshsize,round(10*rFactor));
end
load(filename)

if nametype <3
    filterRadius = meshsize.*rFactor;
    %Ugly hack to make domains correct irrespective of open or close currently
    %implementeted in filterCirc...
    if closedCirc
        filterRadius = filterRadius + eps(filterRadius);
    else
        filterRadius = filterRadius - eps(filterRadius);
    end
    [~,B] = filterCirc(rho(:),filterRadius,param);
    
    xB = zeros(param.nely,32*meshsize);
    empty = size(xB,2)-size(B,2);
    ind = round(empty/2)+(1:size(B,2));
    xB(ind,ind) = B;
else
    filterRadius = meshsize.*rFactor;
    [~,Bo] = filterOct(rho(:),filterRadius,param);
    [~,Br] = filterRect(rho(:),filterRadius,param);
    xBo = zeros(param.nely,32*meshsize);
    empty = size(xBo,2)-size(Bo,2);
    ind = round(empty/2)+(1:size(Bo,2));
    xBo(ind+20*meshsize,ind) = Bo;
    xBr = zeros(param.nely,32*meshsize);
    xBr(ind-4*meshsize,ind) = Br;
end

rhofc = rhof{2};
%rhofc = rhof{1}.^3; %SIMP with penalty 3... -> physical design
%rhop = rhof{1}.^3;

hh = figure(1000*meshsize + round(10*rFactor));
set(hh,'Position',[2500 200 1000 800]);
if nametype < 3
    imagesc([reshape(1-rhofc,param.nely,param.nelx) xB]);
else
    imagesc([reshape(1-rhofc,param.nely,param.nelx) xBr+xBo]);    
end
axis equal;axis off;colormap('gray');

if nametype < 3
    A = zeros(param.nely,param.nelx+size(xB,2),3);
    A(:,:,1) = [reshape(1-rhofc,param.nely,param.nelx) 1-1*xB];
    A(:,:,2) = [reshape(1-rhofc,param.nely,param.nelx) 1-1*xB];
    A(:,:,3) = [reshape(1-rhofc,param.nely,param.nelx) 1-0*xB];
else
    A = zeros(param.nely,param.nelx+size(xBo,2),3);
    A(:,:,1) = [reshape(1-rhofc,param.nely,param.nelx) 1-0*xBo-1*xBr];
    A(:,:,2) = [reshape(1-rhofc,param.nely,param.nelx) 1-1*xBo-1*xBr];
    A(:,:,3) = [reshape(1-rhofc,param.nely,param.nelx) 1-1*xBo-0*xBr];
end
if nametype == 1
    imwrite(A,sprintf('heatComp_M%dV%2dR%03d.png',meshsize,volfrac,round(10*rFactor)));
elseif nametype == 2
    imwrite(A,sprintf('heatComp_sM%dV%2dR%03d.png',meshsize,volfrac,round(10*rFactor)));
else
    imwrite(A,sprintf('heatComp_dfM%dV%2dR%03d.png',meshsize,volfrac,round(10*rFactor)));    
end

fprintf('Mesh: %d rFactor: %4.1f   MND: %7.5f %%   MOCD: %8.6f %% fraqDiff5e-1: %7.5f %%   Vol: %5.3f\n',...
    meshsize,rFactor,mean(rhofc.*(1-rhofc)*400),...
    100*norm(rhof{2}-rhof{1},1)/param.nel,...
    100*sum(abs(rhof{2}-rhof{1})>0.5)/param.nel, ...
    mean(rhof{2}));

fprintf('                        MNDc: %8.2e   MNDo: %8.2e    MOCD: %8.2e   numDiff5e-1: %d\n',...
    mean(rhofc.*(1-rhofc)*4),...
    mean(rhof{1}.*(1-rhof{1})*4),...
    norm(rhof{2}-rhof{1},1)/param.nel,...
    sum(abs(rhof{2}-rhof{1})>0.5));


fnUPen = sprintf('ResultsC/vMesh%1dV%02dfinal.mat',meshsize,volfrac);
temp = load(fnUPen);
relObj = objC/temp.obj;
fprintf('\t       RelObj: %6.4f            NumOCD: %d\n',relObj,sum(abs(rhof{2}-rhof{1})>0.5));

if nametype == 3 
    ind = round(param.nely/2-size(Br,2)/2)+(1:size(Br,2));
    if meshsize == 4 && rFactor == 4
        temp = reshape(1-rhofc,param.nely,param.nelx);
        halfSize = 64;
        temp = temp(param.nely/2+(-halfSize:halfSize),564+(-halfSize:halfSize));
        imagesc(temp)
        B = ones(2*halfSize+1,2*halfSize+1,3);
        B(:,:,1) = temp;B(:,:,2) = temp;B(:,:,3) = temp;
        Btemp = B;
        filterRadius = rFactor*meshsize;
        ind = 1+halfSize+(-filterRadius:filterRadius);
        B(ind,filterRadius+ind,1) = 1-Br;
        B(ind,filterRadius+ind,2) = 1-Br;
        imwrite(B,sprintf('heatComp_dfM%dV%2dR%03dCloseupR.png',meshsize,volfrac,round(10*rFactor)));    
        B = Btemp;
        B(ind,filterRadius+ind,2) = 1-Bo;
        B(ind,filterRadius+ind,3) = 1-Bo;
        imwrite(B,sprintf('heatComp_dfM%dV%2dR%03dCloseupO.png',meshsize,volfrac,round(10*rFactor)));            
    elseif meshsize == 8 && rFactor == 4
        temp = reshape(1-rhofc,param.nely,param.nelx);
        halfSize = 128;
        temp = temp(param.nely/2+(-halfSize:halfSize),1124+(-halfSize:halfSize));
        B = ones(2*halfSize+1,2*halfSize+1,3);
        B(:,:,1) = temp;B(:,:,2) = temp;B(:,:,3) = temp;
        Btemp = B;
        filterRadius = rFactor*meshsize;
        ind = 1+halfSize+(-filterRadius:filterRadius);
        B(ind,filterRadius+ind,1) = 1-Br;
        B(ind,filterRadius+ind,2) = 1-Br;
        imwrite(B,sprintf('heatComp_dfM%dV%2dR%03dCloseupR.png',meshsize,volfrac,round(10*rFactor)));    
        B = Btemp;
        B(ind,filterRadius+ind,2) = 1-Bo;
        B(ind,filterRadius+ind,3) = 1-Bo;
        imwrite(B,sprintf('heatComp_dfM%dV%2dR%03dCloseupO.png',meshsize,volfrac,round(10*rFactor)));            
%         temp(ind,filterRadius+ind) = 1-Br;
%         imagesc(temp);colormap('gray');        
%         keyboard
    end
%    imagesc(A)
end