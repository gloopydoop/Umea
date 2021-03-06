function [f0val,df0dx,fval,dfdx] =  HEATCONJ(meshsize,rFactor)
close all
clear all
clc
%MAIN the main function for the heat conduction problem
%
%    MAIN(MESHSIZE,RFACTOR) solves the minimal heat compliance
%    problem by using the specified MESHSIZE and RFACTOR. For
%    more information about the parameter setting, please see HELP
%    INITPARAM.
%    
%    The default values of the input parameters are:
%    MESHSIZE = 1;
%    RFACTOR = 2;
%if nargin<2
 %   rFactor = 5;
%    if nargin<1
        meshsize = 1;
%    end
%end

HarrysVolume = true;
octFilter = false;
volfrac = 0.5;

rFactor = 3;
%param.rFactor = rFactor ;
[param] = initParam(meshsize,rFactor,octFilter);

% This is if I want a movie------------------------------------------------
movie_plot = true;
% Axis requirements
movie_max_iterations = 50*length(param.alpha);
movie_frameskip = 10;
movie_max_obj = 0.7;
movie_fig_i = 1;
movie_figure_size = [100,100, 1100,300];
movie_filename = 'bigDog.gif';
%movie_obj = zeros(movie_max_iterations);
%movie_iter = zeros(movie_max_iterations);
movie_i = 1;
% -------------------------------------------------------------------------







if HarrysVolume == true
    volfrac = volfrac*(1-param.square_ratio^2);
end


filterParam = initFilter(param,param.alpha(1));

% This is where we have rho_omega
% one means it's in the domain
% zero means it's not in the domain
%--------------------------------------------------------------------------


rho = volfrac*ones(param.nel,1);
rhof = cell(filterParam.numCascades);

%--------------------------------------------------------------------------


for n = 1:filterParam.numCascades;
    rhof{n} = rho;
end
iter = 0;
f = assembleF(param);
c = 0;
for m = 1:length(param.penal)
    movie_iter_m(m) = iter;
    movie_J_m(m) = c;
    
    penal = param.penal(m);
    filterParam = initFilter(param,param.alpha(m));
    
    s = cell(filterParam.Nmax,filterParam.numCascades);

    change = 1.0;
    inneriter = 0;

    rho = rhof{param.filt_type};
    % Harry being sneaky!
    rho = start_padding(rho(:),param);
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
        rhof{n} = remove_padding(rhof{n},param);

    end
    rho = remove_padding(rho,param);
    
    while change > 0.01 && inneriter<50
        tic;
        %   THIS IS WHERE I ADD MY SNEAKY CHANGE
        %   THIS IS WHERE WE CHOSE OPEN OR CLOSE -> I chose open!!!
        rhop = param.weakMaterial + (1-param.weakMaterial)*rhof{param.filt_type}.^penal;
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
            dc = (1-param.weakMaterial)*penal*(rhof{param.filt_type}.^(penal-1)).*dcdap;
        else
            dc = (1-param.weakMaterial)*dcdap;
        end
        % HARRY being sneaky 
        dc = start_padding(dc,param);
        dv = start_padding(ones(param.nel,1),param);
        rho = start_padding(rho(:),param);
        %dv = reshape((1-param.nullel),param.nelx*param.nely,1);      
        
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
%        Filter 2: close - use for volume...    
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

        dv = remove_padding(dv,param);
        dc = remove_padding(dc,param);
        rho = remove_padding(rho,param);
        %%OC UPDATE
        l1 = 0; l2 = 1e16;
        move = 0.01;
        while l2-l1 > l2*1e-12;
            lmid = 0.5*(l2+l1);
            
            if HarrysVolume == true
                Bvec = -dc(param.comEl)./dv(param.comEl); %gain per volume change...
                scaleX = (Bvec/lmid).^param.ocSmooth;
                rhonew(param.comEl) = max(0,max(rho(param.comEl)-move,min(1,min(rho(param.comEl)+move,...
                     rho(param.comEl).*scaleX)))); 
                rhonew(param.fixedEl) = 0;
                 if size(rhonew,1) ~= size(rho,1);
                    rhonew = rhonew';
                 end
            else
                Bvec = -dc./dv; %gain per volume change...
                scaleX = (Bvec/lmid).^param.ocSmooth;
                rhonew = max(0,max(rho-move,min(1,min(rho+move,...
                     rho.*scaleX)))); 
            end
            
            %% FILTERING OF DESIGN VARIABLES
            rhonew = start_padding(rhonew(:),param);
            for n = 1:filterParam.numCascades
                s{1,n} = filterParam.cascade{n}.G{1}(...
                    filterParam.cascade{n}.f{1}(rhonew(:)))./filterParam.cascade{n}.Ni{1};
                for k = 2:filterParam.cascade{n}.N
                    s{k,n} = filterParam.cascade{n}.G{k}(filterParam.cascade{n}.f{k}(...
                        filterParam.cascade{n}.g{k-1}(s{k-1,n})) )./ ...
                        filterParam.cascade{n}.Ni{k};
                end
                % This is me being dodgy!
                rhof{n} = filterParam.cascade{n}.g{filterParam.cascade{n}.N}( ...
                    s{filterParam.cascade{n}.N,n});
                rhof{n} = min(1,max(rhof{n},0));
                rhof{n} = remove_padding(rhof{n},param);
            end
            rhonew = remove_padding(rhonew,param);
            % harry again!
            if mean(rhof{param.filt_type}(param.comEl)) > volfrac
                l1 = lmid;
            else
                l2 = lmid;
            end            
            
        end
        change = max(abs(rho-rhonew));
        t2 = tic; 
        % This is the movie stuff!!! --------------------------------------
        if movie_plot
            movie_obj(movie_i) = c;
            movie_iter(movie_i) = iter;
            movie_i = movie_i + 1;
            if mod(iter,movie_frameskip) == 0
                if movie_fig_i == 1
                    if param.boundary_treatment == 1
                        str_boundary = 'Padded boundaries';
                    else
                        str_boundary = 'Non-padded boundaries';
                    end
                    
                    if param.filt_type == 1
                        str_filt = 'Open';
                    else
                        str_filt = 'Close';
                    end
                    movie_fig = figure('position',movie_figure_size);
                                       
                end
                axis tight manual
                subplot(1,3,1)
                plot_filter(param, rho)
                
                subplot(1,3,2)
                plot(movie_iter,movie_obj,'k')
                axis([0, movie_max_iterations, 0, movie_max_obj])
                xlabel('Iteration')
                ylabel('Obj')
                grid on
                hold on
                if m > 1
                    plot(movie_iter_m(2:end),movie_J_m(2:end),'ro')
                    for mm = 2:m
                    text(movie_iter_m(mm),movie_J_m(mm),num2str(mm-1),...
                    'horizontalAlignment','center','verticalAlignment','top')
                end
                end
                hold off
                title([str_filt,' filter, ',str_boundary]);
                drawnow;
                
                subplot(1,3,3)
                semilogy(param.penal,param.alpha,'k-o','MarkerFaceColor','w')
                hold on
                semilogy(param.penal(1:m),param.alpha(1:m),'ko','MarkerFaceColor','r')
                for mm = 1:m
                text(param.penal(mm),param.alpha(mm),num2str(mm),...
                    'horizontalAlignment','center','verticalAlignment','top')
                end
                grid on
                axis([min(param.penal)-1, max(param.penal) + 1, min(param.alpha)/10, max(param.alpha)*10])
                xlabel('p')
                ylabel('\beta')
                drawnow;
                hold off
                frame = getframe(movie_fig); 
                im = frame2im(frame); 
                [imind,cm] = rgb2ind(im,256); 
                % Write to the GIF File
                if movie_fig_i == 1
                    imwrite(imind,cm,movie_filename,'gif', 'Loopcount',inf);
                else
                    imwrite(imind,cm,movie_filename,'gif','WriteMode','append');
                end
                movie_fig_i = movie_fig_i + 1;
                
            end
        end
        
        if param.plotDesign
            x_center = linspace(1/param.nelx/2,1-1/param.nelx/2,param.nelx);
            y_center = linspace(1/param.nely/2,1-1/param.nely/2,param.nely);
            x_node = linspace(0,1,param.nelx + 1);
            y_node = linspace(0,1,param.nely + 1);
            
            if mod(inneriter,param.plotmod) == 0
                figure(10+meshsize);
                subplot(2,2,1);
                imagesc(x_center,y_center,reshape(1-rho,param.nely,param.nelx));
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
                blue = cat(3, zeros(size(background)),zeros(size(background)), ones(size(background)));
                hold on 
                h = imagesc(x_center,y_center,blue); 
                hold off 
                set(h, 'AlphaData', background) 
                %blue = cat(3, zeros(size(param.BCs)),zeros(size(param.BCs)), ones(size(param.BCs)));
                hold on
                h = imagesc(x_center,y_center,blue); 
                 set(h, 'AlphaData', param.nullel/2) 
                
%--------------------------------------------------------------------------
                caxis([0 1]);colormap(gray);drawnow;    
                subplot(2,2,2);
                imagesc(reshape(1-rhop,param.nely,param.nelx));
                axis equal;axis off;
                caxis([0 1]);colormap(gray);drawnow;
                subplot(2,2,3);
                
                imagesc(x_center,y_center,reshape(rhof{2}-rhof{1},param.nely,param.nelx));
                %imagesc(x_center,y_center,reshape(filterParam.cascade{2}.Ni{k},param.nely,param.nelx));
                hold on
                h = imagesc(x_center,y_center,blue); 
                set(h, 'AlphaData', param.nullel/2) 
                hold off
                axis equal;axis off;colorbar;drawnow;
                
                subplot(2,2,4);
                imagesc(x_center,y_center,reshape(dc./dv,param.nely,param.nelx));
                omega = ones(size(reshape(rhop,param.nely,param.nelx)));
                hold on
                h = imagesc(x_center,y_center,blue); 
                set(h, 'AlphaData', param.nullel/2) 
                hold off
                axis equal;axis off;colorbar;drawnow;

            end
        end
        fprintf(1,'Iteration: %5d Penal: %3.1f Obj: %8.4e coDiff: %8.4e ',...
            iter,penal,c,pen/param.nel);
        fprintf(1,'Design: [%4.1e,%4.1e] Density: [%4.1e,%4.1e] ',...
            min(rho),max(rho),min(rhop),max(rhop));
        fprintf(1,'MND: %6.2f%% Vol: %5.3f change: %6.4f \n',...
            400*rhop'*(1-rhop)/param.nel,mean(rhof{param.filt_type}(param.comEl)),change);
        rho = rhonew;
        time2 = toc(t2);
        time1 = toc;
        disp(num2str(time2/time1*100))
    end
end



fn = sprintf('Mesh%1dR%03dfinal.mat',meshsize,round(10*rFactor));
save(fn,'filterParam','param','rho','rhof','rhop')
