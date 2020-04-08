close all
clear all
R = 1;
L = 3;
dx = 0.01;
x = -L:dx:L;

rho_unf = ones(size(x));
rho_unf(x<0) = 0;
k = 200;
beta = 0.01;
logistic = @(x,k)(1./(1+exp(-k*x)));
inside = @(u,beta,k,R)(-1./(2*R*k*(1+beta))*(log(u) + 1/beta*log(1 + beta.*(1 + u))));
mid = @(x,beta,k,R)( inside(exp(-k*(x+R)),beta,k,R) - inside(exp(-k*(x-R)),beta,k,R));
fin = @(x,beta)(1./x - beta);

% filtery stuff
n = 1;
[param] = initParam2D(2*L/dx+1);
param.r = R/dx;
filterParam =  initFilter_2D(param,beta);

    s{1,n} = filterParam.cascade{1}.G{1}(...
        filterParam.cascade{1}.f{1}(rho_unf))./...
        filterParam.cascade{1}.Ni{1}; 
    s{2,n} = filterParam.cascade{1}.G{2}(filterParam.cascade{1}.f{2}(...
                    filterParam.cascade{1}.g{1}(s{1,n})) )./ ...
                    filterParam.cascade{1}.Ni{2};
    s{2,n} = filterParam.cascade{1}.g{2}(s{2,n});
rho_filt = s{2,n}

erode = fin(mid(x,beta,k,R),beta);

plot(x,logistic(x,k))
hold on
plot(x,erode)
plot(x,rho_filt)


%plot(x,rho_unf)