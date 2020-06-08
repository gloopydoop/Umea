clear all
clc
close all

 figure_width = 700;
figure_height = 600;
probe = figure();
probe.Position(3) = figure_width;
probe.Position(4) = figure_height;

recompute = true;
dx = 0.0005;
L = 1;
[param] = initParam2D(L/dx);

x = 0:dx:dx*(param.nelx-1);
rho = zeros(size(x));
rho(1:floor(param.nelx/2)) = 1;

legstr = {};
sym_map = {'o','+','*','.','>'};
c_map = lines;
r = 0.1:0.05:0.2;
eps = [0.1,0.05, 0.01, 0.005, 0.001];
alpha1 = linspace(0.001,0.01,20);
alpha2 = linspace(0.01,0.1,20);
alphas = [alpha1,alpha2(2:end)];

if recompute == true
for j = 1:length(r)
    
param.r = int32(r(j)/dx);
%plot(x,rho,'k');
for m = 1:length(eps)
param.eps = eps(m);
for n = 1:length(alphas)
    filterParam =  initFilter_2D(param,alphas(n));
    s{1,n} = filterParam.cascade{1}.G{1}(...
        filterParam.cascade{1}.f{1}(rho))./...
        filterParam.cascade{1}.Ni{1}; 
    s{2,n} = filterParam.cascade{1}.G{2}(filterParam.cascade{1}.f{2}(...
                    filterParam.cascade{1}.g{1}(s{1,n})) )./ ...
                    filterParam.cascade{1}.Ni{2};
    s{2,n} = filterParam.cascade{1}.g{2}(s{2,n});
    %hold on
    %plot(x,s{2,n})
    window = x(s{2,n} < 1- param.eps & s{2,n} > param.eps);
    %plot([window(1), window(end)],[1+0.1*n,1+0.1*n])
    window_length(n,m,j) = (window(end)-window(1))/(4*double(param.r)*dx);
end
end

end
save('interp_info.mat','window_length','alphas','eps','r')
else
    load('interp_info.mat')
end

for j = 1:length(r) 
for m = 1:length(eps)
plot(alphas,window_length(:,m,j),'color',c_map(m,:),'marker',sym_map{j})
hold on
if j == 1
legstr = {legstr{:}, ['\epsilon = ',num2str(100*eps(m)),'%']};
end
end
end

xlabel('\beta')
ylabel('Effective Smoothing ($\overline{R}$ /4R)','interpreter','latex')
legend(legstr,'Location','southeast')
ylim([0,1])
    
%plot([x(floor(param.nelx/2)),x(floor(param.nelx/2)+param.r)],[-0.1,-0.1],'k')
%text(mean([x(floor(param.nelx/2)),x(floor(param.nelx/2)+param.r)]),-0.05,'R','HorizontalAlignment','center')

%plot([x(floor(param.nelx/2)-2*param.r),x(floor(param.nelx/2)+2*param.r)],[-0.2,-0.2],'k')
%text(x(floor(param.nelx/2)),-0.15,'4 x R','HorizontalAlignment','center')


set(gcf, 'PaperUnits', 'normalized')
set(gcf,'renderer','Painters')
saveas(gcf,'interp_data.eps','epsc')