clear all
clc
close all

R = 0.03;
lam = 1.1111111110;
A = .1111111111;
x = 0:0.01:lam;
y = x;
[X,Y] = meshgrid(x,y);

lam = 1.1111111110;
A = .1111111111;
y_eq = @(x,A,lam)(A - A*cos(2*pi*x/lam));

for i = 1:length(x)
    for j = 1:length(y)
        [dist,~,~,above] = perp_distance(y(j),x(i));
        if dist < R
            scaler = dist/R;
        else
            scaler = 0;
        end
        
        if above == true
            if scaler > 0
            Z(i,j) = 1/2 + scaler/2;
            else
            Z(i,j) = 1;
            end
        else
            if scaler > 0
            Z(i,j) = 1/2 - scaler/2;
            else
            Z(i,j) = 0;
            end
        end
        
        %if Z(i,j) > 6
        %    Z(i,j);
        %    Z(i,j) = 0; 
        %end
    end
end

a = 0.91;
b = 1;
[dist,x0,y0]= perp_distance(a,b);


D = @(x,A,lam,a,b)((b-A+A*cos(2*pi*x/lam)).^2 + (a-x).^2);
D(x,A,lam,a,b);
surf(X,Y,Z)
