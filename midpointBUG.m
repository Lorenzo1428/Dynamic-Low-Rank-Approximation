function Y = midpointBUG(Y,N,h,f,r)
    tol = 1e-15;
    [U,~,V] = svd(Y);
    U = U(:,1:r);
    V = V(:,1:r);
    for k = 0:N-1 
        Fk1 = f(Y);
        Uk1 = orth([U Fk1*V],tol);
        Vk1 = orth([V Fk1'*U],tol);
        Sk1 = Uk1'*(Y + 0.5*h*Fk1)*Vk1;
        [P,S,W] = svd(Sk1);
        Uk = Uk1*P(:,1:r);
        Vk = Vk1*W(:,1:r);

        Yks = Uk*S(1:r,1:r)*Vk';
        Fks = f(Yks);
        Uks = orth([U Uk Fks*Vk],tol);
        Vks = orth([V Vk Fks'*Uk],tol);
        Sks = Uks'*(Y + h*Fks)*Vks;
        [P,S,W] = svd(Sks);
        U = Uks*P(:,1:r);
        V = Vks*W(:,1:r);
        Y = U*S(1:r,1:r)*V';
    end
end

