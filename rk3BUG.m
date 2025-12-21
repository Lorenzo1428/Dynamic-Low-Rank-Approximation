function Y = rk3BUG(Y,N,h,f,r)
    tol = 1e-15;
    [U,~,V] = svd(Y);
    U = U(:,1:r);
    V = V(:,1:r);

    for k = 1:N-1
        F1 = f(Y);
        U2 = orth([U F1*V],tol);
        V2 = orth([V F1'*U],tol);
        S2 = U2'*(Y + h*F1)*V2;
        [P2, T2, W2] = svd(S2);
        U2 = U2*P2(:,1:r);
        V2 = V2*W2(:,1:r);
        F2 = f(U2*T2(1:r,1:r)*V2');
        
        U3 = orth([U F1*V U2 F2*V2] ,tol);
        V3 = orth([V F1'*U V2 F2'*U2] ,tol);
        S3 = U3'*(Y + 0.25*h*(F1 + F2) )*V3;
        [P3, T3, W3] = svd(S3);
        U3 = U3*P3(:,1:r);
        V3 = V3*W3(:,1:r);

        Fk = f(U3*T3(1:r,1:r)*V3');
        Uk = orth([U F1*V U2 F2*V2 U3 Fk*V3],tol);
        Vk = orth([V F1'*U V2 F2'*U2 V3 Fk'*U3],tol);
        Sk = Uk'*(Y + h*( (1/6)*F1 + (1/6)*F2 + (2/3)*Fk ) )*Vk;
        [P, S, W] = svd(Sk);
        U = Uk*P(:,1:r);
        V = Vk*W(:,1:r);
        Y = U*S(1:r,1:r)*V';
    end
end