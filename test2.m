clc
clear
close all

n = 256;
x = linspace(-pi,pi,n);
[X,Y] = ndgrid(x);
A = sin(X).*sin(Y);

S = svd(A)

