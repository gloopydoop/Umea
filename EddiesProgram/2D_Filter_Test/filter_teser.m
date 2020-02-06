clear all
clc
close all

dx = 0.01;
L = 1;
[param] = initParam2D(L/dx);

x = 0:dx:dx*(param.nelx-1);
rho = zeros(size(x));
rho(1:floor(param.nelx/2)) = 1;

param.r = 0.1/dx;
alphas = [0.5,0.1,0.05];
plot(x,rho,'k');
for n = 1:3
    filterParam =  initFilter_2D(param,alphas(n));
    s{1,n} = filterParam.cascade{1}.G{1}(...
        filterParam.cascade{1}.f{1}(rho))./...
        filterParam.cascade{1}.Ni{1}; 
    s{2,n} = filterParam.cascade{1}.G{2}(filterParam.cascade{1}.f{2}(...
                    filterParam.cascade{1}.g{1}(s{1,n})) )./ ...
                    filterParam.cascade{1}.Ni{2};
    s{2,n} = filterParam.cascade{1}.g{2}(s{2,n});
    hold on
    plot(x,s{2,n})
end
    
plot([x(floor(param.nelx/2)),x(floor(param.nelx/2)+param.r)],[-0.1,-0.1],'k')
text(mean([x(floor(param.nelx/2)),x(floor(param.nelx/2)+param.r)]),-0.05,'R','HorizontalAlignment','center')

plot([x(floor(param.nelx/2)-2*param.r),x(floor(param.nelx/2)+2*param.r)],[-0.2,-0.2],'k')
text(x(floor(param.nelx/2)),-0.15,'4 x R','HorizontalAlignment','center')

