clc
clear 
close all

r = 128;
T = 1e-3;
h = 1e-3;
N = floor(T/h);
n = 128;
theta = 1e-2;

L = diag(ones(n-1,1),1) + diag(ones(n-1,1),-1) - 2*diag(ones(n,1)); 
L = n*n/(4*pi*pi)*L;
x = linspace(0,2*pi,n);
[X,Y] = ndgrid(x);

A0 = (exp(-tan(X).^2) + exp(-tan(Y).^2)).*sin(X).*sin(Y)./(1+ exp(abs(csc(-0.5*X))) + exp(abs(csc(-0.5*Y))));

f = @(A) theta*( L*A + A*L ) + A - A.*A.*A;

% 
%bug
tic
Y = heunBUG(A0,N,h,f,r);
toc

%dlr
tic
%Y = dlr(A0,f,r,h,T);
toc

%heun classic
tic
Ak = heun(A0,f,h,T);
toc

errinf = max(max(abs(Y-Ak)));
err = norm(Ak - Y,"fro");


%rk45 matlab
% tic
% fv = @(t,y) reshape(f(reshape(y,n,n)),[],1);
% [~,yk] = ode45(fv,[0 T],A0(:),odeset(RelTol=1e-10));
% Ak = reshape(yk(end,:),n,n);
% toc