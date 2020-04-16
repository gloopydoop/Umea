[f0val,df0dx,fval,dfdx] = HeatCompliance(xval,param);


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
move = 1;