clc
clear
close all

r = 5;
h = 5e-4;
T = 1;
N = floor(T/h);
n = 128;
theta = 0.01;

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
rho = max(abs(svd(L)));
cfl = 1/rho;

f = @(A) L*A + A*L + theta*nC;

tic
    Ak = heun(A0,f,h,T);
toc

r = [5 15 30 50 n];
for i = 1:length(r)
    tic
        Y_bug(:,:,i) = heunBUGLyapunov(A0,N,h,L,theta*nC,r(i));
    toc
    tic
        Y_dlr(:,:,i) = dlr(A0,f,r(i),h,T);
    toc
    err_bug(i) = norm(Ak - Y_bug(:,:,i),"fro");
    err_dlr(i) = norm(Ak - Y_dlr(:,:,i),"fro");
end
T = table(err_bug',err_dlr');
disp(T)
semilogy(r,err_bug)