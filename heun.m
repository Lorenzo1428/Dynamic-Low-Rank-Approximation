function A = heun(A,F,h,T)
    Nt = floor(T/h);
    for n = 1:Nt-1
        A1 = A + h*F(A);
        A = A + +0.5*h*F(A) + 0.5*h*F(A1);
    end
end