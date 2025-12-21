clc
clear
close all

n = 126;
theta = 0.1;
A = [1 2 3; 4 5 6];

K = ones(126,126,3);

D = diag([1,-2,1]);
L = D; 
for i = 1:n/3-1
    L = blkdiag(L,D);
end
L = n*n/(4*pi*pi)*L;
x = linspace(-pi,pi,n);
[X,Y] = ndgrid(x);

C = zeros(n);
for l = 1:11
    C = C + 10^(-l+1)*exp(-l*(X.^2+Y.^2));
end
nC = C/norm(C,"fro");
f = @(A) L*A + A*L + theta*nC;
f(K(:,:,1:2))







