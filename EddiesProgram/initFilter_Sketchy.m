function filterParam = initFilter(param,alpha,extraFix)
if nargin<3
    extraFix = false;
end

linfilter = 1;                    % if I want just a linear filter
r = [1,1]*param.rFactor*param.meshsize; %radius in pixels...


if linfilter == 0

filterParam.numCascades = 2;
filterParam.Nmax = 2; %maximum number of filters in any cascade
filterParam.cascade = cell(filterParam.numCascades,1);

N = 2;                                  %number of filters in each cascade



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


% HARRY: PRIME x!!!


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



if extraFix
    filterParam.cascade{1}.G{1} = @(x)(filterRect(x,r(2),param));
    filterParam.cascade{1}.G{2} = @(x)(filterRect(x,r(2),param));
else
    if param.octFilter == 0
        filterParam.cascade{1}.G{1} = @(x)(filterCirc(x,r(1),param)); %IM MULTIPLYING BY BCS!
        filterParam.cascade{1}.G{2} = @(x)(filterCirc(x,r(2),param));
    else
        filterParam.cascade{1}.G{1} = @(x)(filterOct(x,r(1),param));
        filterParam.cascade{1}.G{2} = @(x)(filterOct(x,r(2),param));
    end
end
filterParam.cascade{1}.GT = filterParam.cascade{1}.G;    

% HARRY-----This is to account for interior section of rho
oneVec = ones(param.sizex*param.sizey,1);

% Now this has to be modified to INCLUDE the middle guys
for nn = 1:N
    filterParam.cascade{1}.Ni{nn} = filterParam.cascade{1}.G{nn}(oneVec);
    % you did this to remove NaN's but I'm not sure if that's correct
    %filterParam.cascade{1}.Ni{nn}(filterParam.cascade{1}.Ni{nn}==0) = 1;
end

%filter 2: harmonic close
filterParam.cascade{2}.N = N;
filterParam.cascade{2}.radius = r;

filterParam.cascade{2}.f = cell(N,1);
filterParam.cascade{2}.fMax = cell(N,1);
filterParam.cascade{2}.fMin = cell(N,1);
filterParam.cascade{2}.g = cell(N,1);
filterParam.cascade{2}.df = cell(N,1);
filterParam.cascade{2}.dg = cell(N,1);
filterParam.cascade{2}.G = cell(N,1);
filterParam.cascade{2}.Ni = cell(N,1);

filterParam.cascade{2}.f{1} = filterParam.cascade{1}.f{2};
filterParam.cascade{2}.f{2} = filterParam.cascade{1}.f{1};

filterParam.cascade{2}.fMax{1} = filterParam.cascade{1}.fMax{2}; 
filterParam.cascade{2}.fMax{2} = filterParam.cascade{1}.fMax{1};

filterParam.cascade{2}.fMin{1} = filterParam.cascade{1}.fMin{2};
filterParam.cascade{2}.fMin{2} = filterParam.cascade{1}.fMin{1};

filterParam.cascade{2}.df{1} = filterParam.cascade{1}.df{2};
filterParam.cascade{2}.df{2} = filterParam.cascade{1}.df{1};

filterParam.cascade{2}.g{1} = filterParam.cascade{1}.g{2};
filterParam.cascade{2}.g{2} = filterParam.cascade{1}.g{1};

filterParam.cascade{2}.dg{1} = filterParam.cascade{1}.dg{2};
filterParam.cascade{2}.dg{2} = filterParam.cascade{1}.dg{1};

if extraFix
        filterParam.cascade{2}.G{1} = @(x)(filterOct(x.*filterParam.line_BCs,r(1),param).*filterParam.line_BCs);
        filterParam.cascade{2}.G{2} = @(x)(filterOct(x.*filterParam.line_BCs,r(2),param).*filterParam.line_BCs);
else
    filterParam.cascade{2}.G{1} = filterParam.cascade{1}.G{2};
    filterParam.cascade{2}.G{2} = filterParam.cascade{1}.G{1};
end

filterParam.cascade{2}.GT = filterParam.cascade{2}.G;

if extraFix
    oneVec = ones(param.nelx*param.nely,1);
    for nn = 1:N
        filterParam.cascade{2}.Ni{nn} = filterParam.cascade{2}.G{nn}(oneVec);
    end
else
    filterParam.cascade{2}.Ni{1} = filterParam.cascade{1}.Ni{2};
    filterParam.cascade{2}.Ni{2} = filterParam.cascade{1}.Ni{1};
end

else
    
   
    
        filterParam.numCascades = 2;
filterParam.Nmax = 2; %maximum number of filters in any cascade
filterParam.cascade = cell(filterParam.numCascades,1);
        N = 2;
    for i = 1:2

filterParam.cascade{i}.N = N;
filterParam.cascade{i}.radius = r;

filterParam.cascade{i}.f = cell(N,1);
filterParam.cascade{i}.fMax = cell(N,1);
filterParam.cascade{i}.fMin = cell(N,1);
filterParam.cascade{i}.g = cell(N,1);
filterParam.cascade{i}.df = cell(N,1);
filterParam.cascade{i}.dg = cell(N,1);
filterParam.cascade{i}.G = cell(N,1);
filterParam.cascade{i}.Ni = cell(N,1); % Number of elements in each neighbourhood


% HARRY: PRIME x!!!


filterParam.cascade{i}.f{1} = @(x)(x); % Erode
filterParam.cascade{i}.f{2} = @(x)(filterParam.cascade{1}.f{1}(x)); % Dilate

filterParam.cascade{i}.fMax{1} = filterParam.cascade{1}.f{1}(0);
filterParam.cascade{i}.fMax{2} = filterParam.cascade{1}.f{2}(1);

filterParam.cascade{i}.fMin{1} = filterParam.cascade{1}.f{1}(1);
filterParam.cascade{i}.fMin{2} = filterParam.cascade{1}.f{2}(0);


filterParam.cascade{i}.df{1} = @(x)(1);
filterParam.cascade{i}.df{2} = @(x)(1);

filterParam.cascade{i}.g{1} = @(x)(x);
filterParam.cascade{i}.g{2} = @(x)(x);    

filterParam.cascade{i}.dg{1} = @(x)(1);
filterParam.cascade{i}.dg{2} = @(x)(1);



if extraFix
    filterParam.cascade{i}.G{1} = @(x)(filterRect(x,r(2),param));
    filterParam.cascade{i}.G{2} = @(x)(filterRect(x,r(2),param));
else
    if param.octFilter == 0
        filterParam.cascade{i}.G{1} = @(x)(filterCirc(x,r(1),param)); %IM MULTIPLYING BY BCS!
        filterParam.cascade{i}.G{2} = @(x)(filterCirc(x,r(2),param));
    else
        filterParam.cascade{i}.G{1} = @(x)(filterOct(x,r(1),param));
        filterParam.cascade{i}.G{2} = @(x)(filterOct(x,r(2),param));
    end
end
filterParam.cascade{i}.GT = filterParam.cascade{1}.G;    

% HARRY-----This is to account for interior section of rho
oneVec = ones(param.sizex*param.sizey,1);

% Now this has to be modified to INCLUDE the middle guys
for nn = 1:N
    filterParam.cascade{i}.Ni{nn} = filterParam.cascade{1}.G{nn}(oneVec);
    % you did this to remove NaN's but I'm not sure if that's correct
    %filterParam.cascade{1}.Ni{nn}(filterParam.cascade{1}.Ni{nn}==0) = 1;
end
    end

end