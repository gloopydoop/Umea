close all
clear all
clc

x = -3:0.01:3;
% my parameters
k = 1000;
beta = 0.001;
R = 1;


logistic = @(x,k)(1./(1+exp(-k.*x)));
discrete = zeros(size(x));
discrete(x>0) = 1;

% cascade
erode_mid = @(u,k,beta,R)( -1./(2*R*k*(1+beta)).*(log(u) + 1/beta*log(1+beta*(1+u))));
x2u = @(x,k)(exp(-k*(x)));
u_fn = @(u,beta)(1 + beta*(1+u));
erode_2 = @(x,k,beta,R,u1,u2)( -1./(2*R*k*(1+beta)).*(-2*R*k+1/beta.*log(u1./u2)));

erode_2 = @(x,k,beta,R)( -1./(2*R*k*(1+beta)).*(-2*R*k+1/beta.*log((1 + beta*(1 + exp(-k*(x+R))))./(1 + beta*(1 + exp(-k*(x-R)))))));

erode_finv = @(x,beta)(1./x - beta);


erode = erode_2(x,k,beta,R);
erode = erode_finv(erode,beta);
%erode = erode_mid(exp(-k*(x+R)),k,beta,R) - erode_mid(exp(-k*(x-R)),k,beta,R);
%erode = erode_finv(erode,beta);

plot(x,discrete,'k')
hold on
plot(x,logistic(x,k),'r')
plot(x,erode)

figure()
logbit = @(x,k,beta,R)(log((1 + beta*(1 + exp(-k*(x+R))))./(1 + beta*(1 + exp(-k*(x-R))))));
plot(x,logbit(x,k,beta,R))