clc
clear
close all

r = 5;
T = 1;
h = 5e-4;
N = floor(T/h);
n = 12;
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

nC = reshape(nC,[],1);
L = reshape(L,[],1);
A0 = reshape(A0,[],1);


f = @(t,A) L'*A + A'*L + theta*nC;

[t,Ak] =  ode45(f,[0 1],A0);
A = reshape(Ak(end,:),[n,n]);

