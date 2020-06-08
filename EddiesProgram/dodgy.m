%function main(meshsize,rFactor)
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
        meshsize = 4;
%    end
%end

HarrysVolume = false;
octFilter = false;
volfrac = 0.5;

rFactor = 2;
%param.rFactor = rFactor ;
[param] = initParam(meshsize,rFactor,octFilter);

% This is if I want a movie------------------------------------------------
movie_plot = true;
% Axis requirements
movie_max_iterations = sum(param.maxiter);
movie_frameskip = 10;
movie_max_obj = 0.4;
movie_fig_i = 1;
movie_max_MND = 60;
final_plot = true;

movie_filename = 'bigDog.gif';
right_plot = true;
title_plot = false;
mid_plot = false;
numsubplots = 1;
if right_plot == true
numsubplots = numsubplots + 1;
end
if mid_plot == true
numsubplots = numsubplots + 1;
end

movie_figure_size = [100,100, numsubplots*400,300]*2;
% if right_plot == true
%     movie_figure_size = [100,100, 1400,300];
% else
%     movie_figure_size = [100,100, 700,300]*2;
% end
%movie_obj = zeros(movie_max_iterations);
%movie_iter = zeros(movie_max_iterations);
movie_i = 1;
% -------------------------------------------------------------------------



%param.plotDesign = true;


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
    
    obj_old = 1000;
    movie_iter_m(m) = iter;
    movie_J_m(m) = c;
    
    penal = param.penal(m);
    filterParam = initFilter(param,param.alpha(m));
    
    

    change = 1.0;
    inneriter = 0;


    s = cell(filterParam.Nmax,filterParam.numCascades);
    
    rho = rhof{2}; % THEY PICK THE CLOSE HERE
    rho = start_padding(rho(:),param,1);
    % Harry being sneaky!
    %if m == 1
    
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

    %end
    %rhof{param.filt_type} = remove_padding(rhof{param.filt_type},param);
    %rho = rhof{param.filt_type}; % THIS one is sketch
    
    end
    rho = remove_padding(rho,param);
    %end
    while change > 0.005 && inneriter<param.maxiter(m)
        tic;
        %   THIS IS WHERE I ADD MY SNEAKY CHANGE
        %   THIS IS WHERE WE CHOSE OPEN OR CLOSE -> I chose open!!!
        rhop = param.weakMaterial + (1-param.weakMaterial)*rhof{param.filt_type}.^penal;
        %rhop = param.weakMaterial + (1-param.weakMaterial)*rhof{1}.^penal;

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
            %dc = (1-param.weakMaterial)*penal*(rhof{1}.^(penal-1)).*dcdap;
        else
            dc = (1-param.weakMaterial)*dcdap;
        end
        % HARRY being sneaky 
        dc = start_padding(dc,param,1);
        dv = start_padding(ones(param.nel,1),param,0);
        dvo = start_padding(ones(param.nel,1),param,1);
        rho = start_padding(rho(:),param,0);
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
%      Filter 2: close - use for volume...    
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
            
%      Filter 3: open - use for volume...    
        for k = filterParam.cascade{1}.N:-1:2
            dvo(:) = filterParam.cascade{1}.df{k}( ...
                    filterParam.cascade{1}.g{k-1}(s{k-1,1})).* ...
                filterParam.cascade{1}.GT{k}(dvo(:).* ...
                    filterParam.cascade{1}.dg{k}(s{k,1})./ ...
                filterParam.cascade{1}.Ni{k});
        end
        dvo(:) = filterParam.cascade{1}.df{1}(rho(:)) .* ...
                filterParam.cascade{1}.GT{1}(dvo(:).* ...
                    filterParam.cascade{1}.dg{1}(s{1,1})./ ...
                filterParam.cascade{1}.Ni{1});
%      Filter 3: open - use for volume...    
        for k = filterParam.cascade{1}.N:-1:2
            dv(:) = filterParam.cascade{2}.df{k}( ...
                    filterParam.cascade{2}.g{k-1}(s{k-1,2})).* ...
                filterParam.cascade{2}.GT{2}(dv(:).* ...
                    filterParam.cascade{2}.dg{k}(s{k,2})./ ...
                filterParam.cascade{2}.Ni{k});
        end
        dv(:) = filterParam.cascade{2}.df{1}(rho(:)) .* ...
                filterParam.cascade{2}.GT{1}(dv(:).* ...
                    filterParam.cascade{2}.dg{1}(s{1,2})./ ...
                filterParam.cascade{2}.Ni{1});
        dv = remove_padding(dv,param);
        dvo = remove_padding(dvo,param);
        dc = remove_padding(dc,param);
        rho = remove_padding(rho,param);
        k_codiff = 0.001;
        %dc = dc + (dv - dvo).*k_codiff;
        %%OC UPDATE
        l1 = 0; l2 = 1e16;
        move = param.moves(m);
        obj_new = c;
        %dv = ones(param.nel,1);
        while l2-l1 > l2*1e-12;
            lmid = 0.5*(l2+l1);
                Bvec = -dc./dv; %gain per volume change...
                %Bvec = filterParam.minismooth(Bvec,filterParam.minirad(m));
                scaleX = (Bvec/lmid).^param.ocSmooth;
                rhonew = max(0,max(rho-move,min(1,min(rho+move,...
                     rho.*scaleX)))); 
            
            %% FILTERING OF DESIGN VARIABLES
            rhonew = start_padding(rhonew(:),param,1);
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
            %rhonew = rhof{param.filt_type};
            % harry again!
            if mean(rhof{2}(param.comEl)) > volfrac %we want the close here!!!
                l1 = lmid;
            else
                l2 = lmid;
            end            
            
        end
        
%         while obj_new >= obj_old && move > 1e-10
%         move = move/2;
%           while l2-l1 > l2*1e-12;
%             lmid = 0.5*(l2+l1);
%                 rhonew = max(0,max(rho-move,min(1,min(rho+move,...
%                      rho.*scaleX)))); 
%             
%             %% FILTERING OF DESIGN VARIABLES
%             rhonew = start_padding(rhonew(:),param);
%             for n = 1:filterParam.numCascades
%                 s{1,n} = filterParam.cascade{n}.G{1}(...
%                     filterParam.cascade{n}.f{1}(rhonew(:)))./filterParam.cascade{n}.Ni{1};
%                 for k = 2:filterParam.cascade{n}.N
%                     s{k,n} = filterParam.cascade{n}.G{k}(filterParam.cascade{n}.f{k}(...
%                         filterParam.cascade{n}.g{k-1}(s{k-1,n})) )./ ...
%                         filterParam.cascade{n}.Ni{k};
%                 end
%                 % This is me being dodgy!
%                 rhof{n} = filterParam.cascade{n}.g{filterParam.cascade{n}.N}( ...
%                     s{filterParam.cascade{n}.N,n});
%                 rhof{n} = min(1,max(rhof{n},0));
%                 rhof{n} = remove_padding(rhof{n},param);
%             end
%             %rhonew = remove_padding(rhonew,param);
%             %rhonew = rhof{param.filt_type};
%             % harry again!
%             if mean(rhof{2}(param.comEl)) > volfrac %we want the close here!!!
%                 l1 = lmid;
%             else
%                 l2 = lmid;
%             end  
%         end
%         rhop = param.weakMaterial + (1-param.weakMaterial)*rhof{param.filt_type}.^penal;
%         [K] = assembleK(param,rhop);
%         u = zeros(param.nn,1);
%         u(param.comNodes) = K(param.comNodes,param.comNodes)\f(param.comNodes);
%         [obj_new,~] = getobj(param,u,f);
%         disp(['Move: ',num2str(move),'  Obj: ',num2str(obj_new)]);
%         
%         end
        c = obj_new;
        obj_old = obj_new;
        %HARRY wtf?
        %rhonew = rhof{param.filt_type};
        change = max(abs(rho-rhonew)); 
        t2 = tic; 
        %% This is the movie stuff!!! --------------------------------------
        if movie_plot
            cur_plot = 1;
            movie_MND(movie_i) = 400*rhof{2}'*(1-rhof{2})/param.nel;
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
                    movie_fig = figure('position',movie_figure_size,'color','w');
                                       
                end
                axis tight manual
                
                
                subplot(1,numsubplots,cur_plot)
                plot_filter(param, rhof{2})
                cur_plot = cur_plot + 1;
                
                if mid_plot == true
                fig = subplot(1,numsubplots,cur_plot);
                
                left_color = [0,0,0];
                right_color = [0,87/256,159/256];
                yyaxis left
                plot(movie_iter,movie_obj,'color',left_color)
                xlim([0, movie_max_iterations])
                ylim([0, movie_max_obj])
                xlabel('Iteration')
                ylabel('Objective Function')
                grid on
                hold on
                if m > 1
                    plot(movie_iter_m(2:end)+5,movie_J_m(2:end),'ro')
                    for mm = 2:m
                    text(movie_iter_m(mm)+5,movie_J_m(mm),num2str(mm-1),...
                    'horizontalAlignment','center','verticalAlignment','top')
                end
                end
                hold off
                
                yyaxis right
                plot(movie_iter,movie_MND,'color',right_color,'linestyle','-');
                ylim([0,movie_max_MND])
                ylabel('MND (Measure of black and whiteness)')
                yticks = [get(gca,'ytick')]';
                percentsy = repmat('%', length(yticks),1);
                yticklabel = [num2str(yticks) percentsy];
                set(gca,'yticklabel',yticklabel) 
                ax = gca;
                ax.YAxis(1).Color = left_color;
                ax.YAxis(2).Color = right_color;
                
                if title_plot == true
                title([str_filt,' filter, ',str_boundary]);
                end
                drawnow;
                end
                
                
                if right_plot == true
                    subplot(1,numsubplots,cur_plot)
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
                end
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
            
            if mod(iter,param.plotmod) == 0
                figure(10+meshsize);
                subplot(2,3,1);
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
                subplot(2,3,4);
                imagesc(reshape(1-rhop,param.nely,param.nelx));
                axis equal;axis off;
                caxis([0 1]);colormap(gray);drawnow;
                subplot(2,3,5);
                
                imagesc(x_center,y_center,reshape(rhof{2}-rhof{1},param.nely,param.nelx));
                %imagesc(x_center,y_center,reshape(filterParam.cascade{2}.Ni{k},param.nely,param.nelx));
                hold on
                h = imagesc(x_center,y_center,blue); 
                set(h, 'AlphaData', param.nullel/2) 
                hold off
                axis equal;axis off;colorbar;drawnow;
                
                subplot(2,3,6);
                imagesc(x_center,y_center,reshape(dc./dv,param.nely,param.nelx));
                omega = ones(size(reshape(rhop,param.nely,param.nelx)));
                hold on
                h = imagesc(x_center,y_center,blue); 
                set(h, 'AlphaData', param.nullel/2) 
                hold off
                axis equal;axis off;colorbar;drawnow;

                
                subplot(2,3,2);
                imagesc(reshape(dc,param.nely,param.nelx));
                axis equal;axis off;
                colormap(gray);drawnow;
                subplot(2,3,5);
                
                
                subplot(2,3,3);
                imagesc(reshape(dv,param.nely,param.nelx));
                axis equal;axis off;
                colormap(gray);drawnow;
                subplot(2,3,5);
            end
        end
        fprintf(1,'Iteration: %5d Penal: %3.1f Obj: %8.4e coDiff: %8.4e ',...
            iter,penal,c,pen/param.nel);
        fprintf(1,'Design: [%4.1e,%4.1e] Density: [%4.1e,%4.1e] ',...
            min(rho),max(rho),min(rhop),max(rhop));
        fprintf(1,'MND: %6.2f%% Vol: %5.3f change: %6.4f \n',...
            400*rhof{2}'*(1-rhof{2})/param.nel,mean(rhof{2}(param.comEl)),change); %the volume is calculated on the close!!
        disp(num2str(param.alpha(m)))
        %%
        rho = rhonew;

        time2 = toc(t2);
        time1 = toc;
        %disp(num2str(time2/time1*100))
    end
end

if final_plot == true
    figure()
    plot_filter(param, rhof{2})
end


fn = sprintf('Mesh%1dR%03dfinal.mat',meshsize,round(10*rFactor));
save(fn,'filterParam','param','rho','rhof','rhop')
