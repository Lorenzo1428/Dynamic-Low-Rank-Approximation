function  Y = heunBUG(Y,N,h,f,r)
    tol = 0;
    [U,S,V] = svd(Y);
    U = U(:,1:r);
    V = V(:,1:r);
    S = S(1:r,1:r);
    for k = 0:N-1 
        Y = U*S*V';
        Fk1 = f(Y);
        [Uk1,~] = qr([U Fk1*V]);
        [Vk1,~] = qr([V Fk1'*U]);
        Y1 = Y + h*Fk1;
        Sk1 = Uk1'*(Y1)*Vk1;
        [P,S,W] = svd(Sk1);
        Uk = Uk1*P(:,1:r);
        Vk = Vk1*W(:,1:r);

        Yks = Uk*S(1:r,1:r)*Vk';
        Fks = f(Yks);
        [Uks,~] = qr([U Fk1*V Uk Fks*Vk]);
        [Vks,~] = qr([V Fk1'*U Vk Fks'*Uk]);
        Sks = Uks'*(Y + 0.5*h*Fk1 + 0.5*h*Fks)*Vks;
        [P,S,W] = svd(Sks);
        U = Uks*P(:,1:r);
        V = Vks*W(:,1:r);
        S = S(1:r,1:r);
    end   
    Y = U*S(1:r,1:r)*V';
end 