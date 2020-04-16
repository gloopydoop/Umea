clear all
close all
clc

dx = 0.01;
L = 4;
[param] = initParam2D(L/dx+1);
x = -2:dx:2;
rho = zeros(size(x));
rho(x>0) = 1;

beta = 0.001;
R = 1/dx;
rho_f = filter_function(rho,beta,R,param);

k = 1;
rho_logistic = 1./(1+exp(-k*x/R));

modelfun = @(b,x)1./(1+exp(-b*x));

opts = statset('nlinfit');
opts.RobustWgtFun = 'bisquare';
beta = nlinfit(x,rho_f',modelfun,1);


plot(x,rho)
hold on
plot(x,rho_f)
plot(x,modelfun(beta,x))