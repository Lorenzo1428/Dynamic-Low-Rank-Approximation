clc
clear
close all

A = rand(128);

[U,S,V] = svd(A);

A1 = U*S*V';

norm(A1-A,"fro")







