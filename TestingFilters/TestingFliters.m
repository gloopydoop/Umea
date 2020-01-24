clear all
clc

%% 
filterParam.numCascades = 2;
filterParam.Nmax = 2; %maximum number of filters in any cascade
filterParam.cascade = cell(filterParam.numCascades,1);


N = 2;                                  %number of filters in each cascade
r = [1,1]*param.rFactor*param.meshsize; %radius in pixels...

%filter 1: harmonic open
filterParam.cascade{1}.N = N;
filterParam.cascade{1}.radius = r;

filterParam.cascade{1}.f = cell(N,1);
filterParam.cascade{1}.fMax = cell(N,1);
filterParam.cascade{1}.fMin = cell(N,1);
filterParam.cascade{1}.g = cell(N,1);
filterParam.cascade{1}.df = cell(N,1);
filterParam.cascade{1}.dg = cell(N,1);
filterParam.cascade{1}.G = cell(N,1);
filterParam.cascade{1}.Ni = cell(N,1); % Number of elements in each neighbourhood
    
filterParam.cascade{1}.f{1} = @(x)(1./(x+alpha)); % Erode
filterParam.cascade{1}.f{2} = @(x)(filterParam.cascade{1}.f{1}(1-x)); % Dilate

filterParam.cascade{1}.fMax{1} = filterParam.cascade{1}.f{1}(0);
filterParam.cascade{1}.fMax{2} = filterParam.cascade{1}.f{2}(1);

filterParam.cascade{1}.fMin{1} = filterParam.cascade{1}.f{1}(1);
filterParam.cascade{1}.fMin{2} = filterParam.cascade{1}.f{2}(0);

filterParam.cascade{1}.g{1} = @(x)(1./x-alpha);
filterParam.cascade{1}.g{2} = @(x)(1-filterParam.cascade{1}.g{1}(x));    

filterParam.cascade{1}.G{1} = @(x)(filterCirc(x,r(1),param));
        filterParam.cascade{1}.G{2} = @(x)(filterCirc(x,r(2),param));

%%


    for n = 1:filterParam.numCascades
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
    