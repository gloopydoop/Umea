function plot_filter(param, rho)
r = param.rFactor;
x_center = linspace(1/param.nelx/2,1-1/param.nelx/2,param.nelx);
y_center = linspace(1/param.nely/2,1-1/param.nely/2,param.nely);
x_node = linspace(0,1,param.nelx + 1);
y_node = linspace(0,1,param.nely + 1);
imagesc(x_center,y_center,reshape(1-rho,param.nely,param.nelx));
axis([x_node(1), x_node(end), y_node(1), y_node(end)]);
axis off;


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
green = cat(3, zeros(size(background)),ones(size(background)), zeros(size(background)));
hold on 
h = imagesc(x_center,y_center,green); 
hold off 
set(h, 'AlphaData', background) 
%blue = cat(3, zeros(size(param.BCs)),zeros(size(param.BCs)), ones(size(param.BCs)));
hold on
h = imagesc(x_center,y_center,blue); 
 set(h, 'AlphaData', param.nullel) 
plot(x_node,y_node(end)*ones(size(x_node)),'b','LineWidth',3) 
plot(x_node(end)*ones(size(x_node)),y_node,'b','LineWidth',3) 
plot(x_node,y_node(1)*ones(size(x_node)),'r','LineWidth',3) 
plot(x_node(1)*ones(size(x_node)),y_node,'r','LineWidth',3) 
hold off
 caxis([0 1]);colormap(gray);drawnow;