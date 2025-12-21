function Ak1 = rk45(Ak,h,F,eps)
    b4 = [25/216 0 1408/2565 2197/4104 -1/5 0];
    b5 = [16/135 0 6656/12825 28561/56430 -9/50 2/55];
    A = [1/4 0 0 0 0 0;
        3/32 9/32 0 0 0 0;
        1932/2197 -7200/2197 7296/2197 0 0 0; 
        439/216 -8 3680/513 -845/4104 0 0;
        -8/27 2 -3544/2565 1859/4104 -11/40 0];
    K = zeros(size(Ak,1),size(Ak,2),6);

    K(:,:,1) = Ak;
    for i = 2:6
        K(:,:,i) = Ak + h*tensorprod(F(K(:,:,1:i-1)),A(i-1,1:i-1),3,2);
    end

    Fk = F(K);
    err = norm(tensorprod(Fk,b5 - b4,3,2),"fro");
    if err >= eps
        h = h*0.84*(eps/err)^0.25;
        for i = 2:6
            K(:,:,i) = Ak + h*tensorprod(Fk(:,:,1:i-1),A(i-1,1:i-1),3,2);
        end
        Ak1 = Ak + h*tensorprod(Fk,b5,3,2);
    else
        Ak1 = Ak + h*tensorprod(Fk,b5,3,2);
    end
end

