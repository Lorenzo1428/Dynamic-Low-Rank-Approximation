clc
clear
close all

n = 5;
x = linspace(-pi,pi,n);
[X,Y] = ndgrid(x);
A = sin(X).*sin(Y);

[U,S,V] = svd(A)

U*U'

