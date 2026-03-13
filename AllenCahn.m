clc
clear 
close all

r = 5;
h = 5e-3;
T = 1;
N = floor(T/h);
n =  128;
theta = 0.01;

L = diag(ones(n-1,1),1) + diag(ones(n-1,1),-1) - 2*diag(ones(n,1)); 
L = n*n/(4*pi*pi)*L;
x = linspace(0,2*pi,n);
[X,Y] = ndgrid(x);

A0 = (exp(-tan(X).^2) + exp(-tan(Y).^2)).*sin(X).*sin(Y)./(1+ exp(abs(csc(-0.5*X))) + exp(abs(csc(-0.5*Y))));

f = @(A) theta*( L*A + A*L ) + A - A.*A.*A;

tic
    Ak = heun(A0,f,h,T);
toc

r = [5 15 30 50 n];
for i = 1:length(r)
    tic
        Y_bug(:,:,i) = heunBUGAllenCahn(A0,N,h,L,theta,r(i));
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



