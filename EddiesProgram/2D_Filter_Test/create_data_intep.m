clear all
clc
close all

dx = 0.0025;
L = 4;          % 4 * R
[param] = initParam2D(L/dx + 1);

x = -2:dx:2;
rho = zeros(size(x));
rho(1:floor(param.nelx/2)) = 1;

legstr = {};
r = 1;
for j = 1:length(r)
    
param.r = int32(r(j)/dx);
eps = [0.01, 0.005, 0.001];
alphas = 10.^(0:-0.1:-4.5);
%plot(x,rho,'k');
for m = 1:length(eps)
param.eps = eps(m);
for n = 1:length(alphas)
    disp(num2str(alphas(n)))
    filterParam =  initFilter_2D(param,alphas(n));
    s{1} = filterParam.cascade{1}.G{1}(...
        filterParam.cascade{1}.f{1}(rho))./...
        filterParam.cascade{1}.Ni{1}; 
    s{2} = filterParam.cascade{1}.G{2}(filterParam.cascade{1}.f{2}(...
                    filterParam.cascade{1}.g{1}(s{1})) )./ ...
                    filterParam.cascade{1}.Ni{2};
    s{2} = filterParam.cascade{1}.g{2}(s{2});
    %hold on
    %plot(x,s{2,n})
    window = x(s{2} < 1- param.eps & s{2} > param.eps);
    %plot([window(1), window(end)],[1+0.1*n,1+0.1*n])
    window_length(n,m,j) = (window(end)-window(1))/(4*double(param.r)*dx);
end
legstr = {legstr{:}, ['\epsilon = ',num2str(eps(m))]}
end
semilogx(alphas,window_length(:,:,j))
hold on

end




xlabel('\beta')
ylabel('Effective Smoothing (R''/4R)')
save('interp_data.mat','window_length','alphas','eps')
legend(legstr)