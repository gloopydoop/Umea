function main(meshsize,rFactor)
%MAIN the main function for the heat conduction problem
%
%    MAIN(MESHSIZE,RFACTOR) solves the minimal heat compliance
%    problem by using the specified MESHSIZE and RFACTOR. For
%    more information about the parameter setting, please see HELP
%    INITPARAM.
%    
%    The default values of the input parameters are:
%    MESHSIZE = 1
%    RFACTOR = 2
if nargin<2
    rFactor = 2;
    if nargin<1
        meshsize = 1;
    end
end

octFilter = false;
volfrac = 0.5;
[param] = initParam(meshsize,rFactor,octFilter);
param.rFactor = 3;

filterParam = initFilter(param,param.alpha(1));
rho = volfrac*ones(param.nel,1);
rhof = cell(filterParam.numCascades);
for n = 1:filterParam.numCascades;
    rhof{n} = rho;
end
iter = 0;

f = assembleF(param);

for m = 1:length(param.penal)
    
    penal = param.penal(m);
    filterParam = initFilter(param,param.alpha(m));
    
    s = cell(filterParam.Nmax,filterParam.numCascades);

    change = 1.;
    inneriter = 0;

    rho = rhof{2};
    for n = 1:filterParam.numCascades;
        s{1,n} = filterParam.cascade{n}.G{1}(...
            filterParam.cascade{n}.f{1}(rho(:)))./filterParam.cascade{n}.Ni{1};
        for k = 2:filterParam.cascade{n}.N
            s{k,n} = filterParam.cascade{n}.G{k}(filterParam.cascade{n}.f{k}(...
                filterParam.cascade{n}.g{k-1}(s{k-1,n})) )./ ...
                filterParam.cascade{n}.Ni{k};
        end
        rhof{n} = filterParam.cascade{n}.g{filterParam.cascade{n}.N}( ...
            s{filterParam.cascade{n}.N,n});
        rhof{n} = min(1,max(rhof{n},0));
    end
    
    while change > 0.01 && inneriter<50
        tic;
        %   THIS IS WHERE I ADD MY SNEAKY CHANGE
        rhop = param.weakMaterial + (1-param.weakMaterial)*rhof{2}.^penal;
        [K] = assembleK(param,rhop);
        % Harry edit I think this was a typo
        %u = zeros(param.nUnknowns,1);
        u = zeros(param.nn,1);
        u(param.comNodes) = K(param.comNodes,param.comNodes)\f(param.comNodes);
        
        if param.saveresults && mod(iter,10)==0
            fn = sprintf('./Trace/trace_M%1d_R%03d_I%05d.mat',...
                meshsize,round(10*rFactor),iter);
            save(fn,'iter','param','penal','rho','rhof','rhop')
        end
        
        [c,dcdap] = getobj(param,u,f);
        pen = sum((rhof{2}-rhof{1}).^param.coPower);
        
        iter = iter + 1;
        inneriter = inneriter + 1;
        if penal>1
            dc = (1-param.weakMaterial)*penal*(rhof{1}.^(penal-1)).*dcdap;
        else
            dc = (1-param.weakMaterial)*dcdap;
        end
        dv = ones(param.nel,1);
                
        
        %% MODIFICATION OF SENSITIVITIES
        %Filter 1: open - use for compliance...
        for k = filterParam.cascade{1}.N:-1:2
            dc(:) = filterParam.cascade{1}.df{k}( ...
                    filterParam.cascade{1}.g{k-1}(s{k-1,1})).* ...
                filterParam.cascade{1}.GT{k}(dc(:).* ...
                    filterParam.cascade{1}.dg{k}(s{k,1})./ ...
                filterParam.cascade{1}.Ni{k});
        end
        dc(:) = filterParam.cascade{1}.df{1}(rho(:)) .* ...
                filterParam.cascade{1}.GT{1}(dc(:).* ...
                    filterParam.cascade{1}.dg{1}(s{1,1})./ ...
                filterParam.cascade{1}.Ni{1});
        %Filter 2: close - use for volume...    
        for k = filterParam.cascade{2}.N:-1:2
            dv(:) = filterParam.cascade{2}.df{k}( ...
                    filterParam.cascade{2}.g{k-1}(s{k-1,2})).* ...
                filterParam.cascade{2}.GT{k}(dv(:).* ...
                    filterParam.cascade{2}.dg{k}(s{k,1})./ ...
                filterParam.cascade{2}.Ni{k});
        end
        dv(:) = filterParam.cascade{2}.df{1}(rho(:)) .* ...
                filterParam.cascade{2}.GT{1}(dv(:).* ...
                    filterParam.cascade{2}.dg{1}(s{1,2})./ ...
                filterParam.cascade{2}.Ni{1});

        %%OC UPDATE
        l1 = 0; l2 = 1e16;
        move = 0.2;
        while l2-l1 > l2*1e-12;
            lmid = 0.5*(l2+l1);
            
            Bvec = -dc./dv; %gain per volume change...
            scaleX = (Bvec/lmid).^param.ocSmooth;
            rhonew = max(0,max(rho-move,min(1,min(rho+move,...
                 rho.*scaleX))));            
            
            %% FILTERING OF DESIGN VARIABLES
            for n = 1:filterParam.numCascades
                s{1,n} = filterParam.cascade{n}.G{1}(...
                    filterParam.cascade{n}.f{1}(rhonew(:)))./filterParam.cascade{n}.Ni{1};
                for k = 2:filterParam.cascade{n}.N
                    s{k,n} = filterParam.cascade{n}.G{k}(filterParam.cascade{n}.f{k}(...
                        filterParam.cascade{n}.g{k-1}(s{k-1,n})) )./ ...
                        filterParam.cascade{n}.Ni{k};
                end
                rhof{n} = filterParam.cascade{n}.g{filterParam.cascade{n}.N}( ...
                    s{filterParam.cascade{n}.N,n});
                rhof{n} = min(1,max(rhof{n},0));
            end
            if mean(rhof{2}) > volfrac
                l1 = lmid;
            else
                l2 = lmid;
            end            
            
        end
        change = max(abs(rho-rhonew));
        t2 = tic;        
        if param.plotDesign
            if mod(inneriter,param.plotmod) == 0
                figure(10+meshsize);
                subplot(2,2,1);
                imagesc(reshape(1-rho,param.nely,param.nelx));
                axis equal;axis off;
% plot filter -------------------------------------------------------------
                r = param.rFactor;
                filterBoxSize = 1+2*floor(r); %2*nDiag+nVert;
                B = ones(filterBoxSize,filterBoxSize);
                mid = floor(r)+1;
                for nn=1:filterBoxSize
                    for mm = 1:filterBoxSize
                        if norm([(nn-mid), (mm-mid)],2)>=r
                            B(nn,mm) = 0;
                        end
                    end
                end
                background = zeros(param.nelx,param.nely);
                background(end - 2*filterBoxSize +1:end-filterBoxSize,end - 2*filterBoxSize +1:end-filterBoxSize) = B;
                green = cat(3, zeros(size(background)),ones(size(background)), zeros(size(background)));
                hold on 
                h = imshow(green); 
                hold off 
                set(h, 'AlphaData', background) 
%--------------------------------------------------------------------------
                caxis([0 1]);colormap(gray);drawnow;    
                subplot(2,2,2);
                imagesc(reshape(1-rhop,param.nely,param.nelx));
                axis equal;axis off;
                caxis([0 1]);colormap(gray);drawnow;
                subplot(2,2,3);
                imagesc(reshape(rhof{2}-rhof{1},param.nely,param.nelx));
                axis equal;axis off;colorbar;drawnow;
                subplot(2,2,4);
                imagesc(reshape(dc./dv,param.nely,param.nelx));
                axis equal;axis off;colorbar;drawnow;
                omega = ones(size(reshape(rhop,param.nely,param.nelx)));

            end
        end
        fprintf(1,'Iteration: %5d Penal: %3.1f Obj: %8.4e coDiff: %8.4e ',...
            iter,penal,c,pen/param.nel);
        fprintf(1,'Design: [%4.1e,%4.1e] Density: [%4.1e,%4.1e] ',...
            min(rho),max(rho),min(rhop),max(rhop));
        fprintf(1,'MND: %6.2f%% Vol: %5.3f change: %6.4f \n',...
            400*rhop'*(1-rhop)/param.nel,mean(rhof{2}),change);
        rho = rhonew;
        time2 = toc(t2);
        time1 = toc;
        disp(num2str(time2/time1*100))
    end
end



fn = sprintf('Mesh%1dR%03dfinal.mat',meshsize,round(10*rFactor));
save(fn,'filterParam','param','rho','rhof','rhop')
