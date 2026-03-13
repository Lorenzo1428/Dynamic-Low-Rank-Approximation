function A = heun(A,F,h,T)
    Nt = floor(T/h);
     %m = size(A,1);
     %Ak = zeros(m,m,Nt);
     %Ak(:,:,1) = A;
    for n = 1:Nt
        Fa = F(A);
        A1 = A + h*Fa;
        A = A + +0.5*h*Fa + 0.5*h*F(A1);

        %Ak(:,:,n+1) = A;
    end
end