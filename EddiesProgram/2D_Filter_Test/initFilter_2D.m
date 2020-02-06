function filterParam = initFilter(param,alpha,extraFix)
if nargin<3
    extraFix = false;
end

filterParam.numCascades = 2;
filterParam.Nmax = 2; %maximum number of filters in any cascade
filterParam.cascade = cell(filterParam.numCascades,1);


N = 2;                                  %number of filters in each cascade
r = [param.r,param.r]; %radius in pixels...

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


filterParam.cascade{1}.df{1} = @(x)(-1./(x+alpha).^2);
filterParam.cascade{1}.df{2} = @(x)(-1*filterParam.cascade{1}.df{1}(1-x));

filterParam.cascade{1}.g{1} = @(x)(1./x-alpha);
filterParam.cascade{1}.g{2} = @(x)(1-filterParam.cascade{1}.g{1}(x));    

filterParam.cascade{1}.dg{1} = @(x)(-1./x.^2);
filterParam.cascade{1}.dg{2} = @(x)(-filterParam.cascade{1}.dg{1}(x));
    
filterParam.cascade{1}.G{1} = @(x)(line_filter_2D(x,r(1)));
filterParam.cascade{1}.G{2} = @(x)(line_filter_2D(x,r(2)));

filterParam.cascade{1}.GT = filterParam.cascade{1}.G;    

oneVec = ones(param.nelx*param.nely,1);
for nn = 1:N
    filterParam.cascade{1}.Ni{nn} = filterParam.cascade{1}.G{nn}(oneVec);
end
