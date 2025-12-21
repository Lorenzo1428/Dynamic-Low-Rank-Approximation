clc
clear
close all

r = 20;
T = 1;
h = 5e-4;
N = floor(T/h);
n = 128;
theta = 1e-5;

L = diag(ones(n-1,1),1) + diag(ones(n-1,1),-1) - 2*diag(ones(n,1)); 
L = n*n/(4*pi*pi)*L;
x = linspace(-pi,pi,n);
[X,Y] = ndgrid(x);

C = zeros(n);
for l = 1:11
    C = C + (10^(-l+1))*exp(-l*(X.^2+Y.^2));
end
nC = C/norm(C,"fro");

A0 = sin(X).*sin(Y);

f = @(A) L*A + A*L + theta*nC;

%dlr classic
tic
Y = dlr(A0,f,r,h,T);
toc

%bug
tic
%Yk = heunBUG(A0,N,h,f,r);
toc

%heun classic
tic
Ak = heun(A0,f,h,T);
toc

err = norm(Ak-Y,"fro");
errinf = norm(Ak(:)-Y(:),inf);



% %rk45 matlab
% tic
% fv = @(t,y) reshape(f(reshape(y,n,n)),[],1);
% [~,yk] = ode45(fv,[0 T],A0(:),odeset(RelTol=1e-10));
% Ak = reshape(yk(end,:),n,n);
% toc