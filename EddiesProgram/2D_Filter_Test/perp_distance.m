function [dist,x_0,y_0,abovebellow] = perp_distance(a,b)

lam = 1.1111111110;
A = .1111111111;



y = @(x,A,lam)(A - A*cos(2*pi*x/lam));
D = @(x,A,lam,a,b)((b-A+A*cos(2*pi*x/lam))^2 + (a-x)^2);
D_p = @(x,A,lam,a,b)(2*pi*A/lam*sin(2*pi*x/lam).*(A-A*cos(2*pi*x/lam)-b) - a +x);
D_pp = @(x,A,lam,a,b)(1 + A*(2*pi/lam)^2*cos(2*pi*x/lam).*(A-A*cos(2*pi*x/lam)-b) + (2*A*pi/lam*sin(2*pi*x/lam))^2  );

x_val = a;
tol = 0.000001;
residu = 10000;
max_iter = 50;
iter = 0;
while residu > tol && iter < max_iter && x_val > 0 && x_val < lam
    x_val = x_val - D_p(x_val,A,lam,a,b)/D_pp(x_val,A,lam,a,b);
    residu = abs(D_p(x_val,A,lam,a,b));
    iter = iter + 1;
end

x_val = max(x_val,0);
x_val = min(x_val,lam);

dist = D(x_val,A,lam,a,b); 
x_0 = x_val;
y_0 = y(x_0,A,lam);

abovebellow = false;
if b > y(a,A,lam)
    abovebellow = true;
end