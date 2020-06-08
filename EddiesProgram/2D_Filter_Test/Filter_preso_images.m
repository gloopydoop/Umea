clear all
clc
close all

figure_width = 700;
figure_height = 400;
probe = figure();
probe.Position(3) = figure_width;
probe.Position(4) = figure_height;

dashes = false;

dx = 0.001;
L = 1;
[param] = initParam2D(L/dx);

x = 0:dx:dx*(param.nelx-1);
rho = zeros(size(x));
rho(1:floor(param.nelx/2)) = 1;

param.r = 0.1/dx;
param.eps = 0.1;
%alphas = [0.01,0.1,1];
alphas = [0.1];
plot(x,rho,'k','linewidth',2);
c_map = lines;

for n = 1:length(alphas)
    filterParam =  initFilter_2D(param,alphas(n));
    s{1,n} = filterParam.cascade{1}.G{1}(...
        filterParam.cascade{1}.f{1}(rho))./...
        filterParam.cascade{1}.Ni{1}; 
    s{2,n} = filterParam.cascade{1}.G{2}(filterParam.cascade{1}.f{2}(...
                    filterParam.cascade{1}.g{1}(s{1,n})) )./ ...
                    filterParam.cascade{1}.Ni{2};
    s{2,n} = filterParam.cascade{1}.g{2}(s{2,n});
    hold on
    filt(n) = plot(x,s{2,n},'color',c_map(n,:),'linewidth',2);
    window = x(s{2,n} < 1- param.eps & s{2,n} > param.eps);
    if dashes == true
        plot([window(1),window(1)],[s{2,n}(x==window(1)),1+0.1*n],'--k');
        plot([window(end),window(end)],[s{2,n}(x==window(end)),1+0.1*n],'--k');
        plot([window(1), window(end)],[1+0.1*n,1+0.1*n],'color',c_map(n,:),'linewidth',2)
    end
end
if dashes == true
    plot([0,x(floor(param.nelx/2))],[1-param.eps,1-param.eps],'--k')
    plot([x(floor(param.nelx/2)),x(end)],[param.eps,param.eps],'--k')
    plot([0.1,0.1],[1-param.eps,1],'k')
    text(0.11,1-param.eps/2,'\epsilon')
    plot([0.9,0.9],[0,param.eps],'k')
    text(0.91,param.eps/2,'\epsilon')

end
    
plot([x(floor(param.nelx/2)),x(floor(param.nelx/2)+param.r)],[-0.1,-0.1],'k')
text(mean([x(floor(param.nelx/2)),x(floor(param.nelx/2)+param.r)]),-0.05,'R','HorizontalAlignment','center')

plot([x(floor(param.nelx/2)-2*param.r),x(floor(param.nelx/2)+2*param.r)],[-0.2,-0.2],'k')
text(x(floor(param.nelx/2)),-0.15,'4 x R','HorizontalAlignment','center')

%legend([filt(1),filt(2),filt(3)],num2str(alphas))
xlabel('Physical length')
ylabel('\rho')
ylim([-0.3,1.2])

set(gcf, 'PaperUnits', 'normalized')
set(gcf,'renderer','Painters')
saveas(gcf,'Mcenter.eps','epsc')
