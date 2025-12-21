function Y = eulespBUG(Y,N,h,f,r)
    tol = 1e-12;
    [U,~,V] = svd(Y);
    U = U(:,1:r);
    V = V(:,1:r);
    for k = 0:N-1
        Fk = f(Y);
        Uk = orth([U Fk*V],tol);
        Vk = orth([V Fk'*U],tol);
        Sk = Uk'*(Y+h*Fk)*Vk;
        [P,S,W] = svd(Sk);
        U = Uk*P(:,1:r);
        V = Vk*W(:,1:r);
        Y = U*S(1:r,1:r)*V';
    end
end

