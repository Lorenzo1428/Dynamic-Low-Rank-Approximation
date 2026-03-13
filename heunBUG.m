function  Y = heunBUG(Y,N,h,f,r)
    [U,S,V] = svd(Y);
    U = U(:,1:r);
    V = V(:,1:r);
    S = S(1:r,1:r);
    for k = 1:N 
        Y = U*S*V';
        Fk1 = f(Y);
        [Uk1,~] = qr([U Fk1*V],0);
        [Vk1,~] = qr([V Fk1'*U],0);
        Sk1 = Uk1'*(Y + h*Fk1)*Vk1;
        [P,S,W] = svd(Sk1,"econ");
        Uk = Uk1*P(:,1:r);
        Vk = Vk1*W(:,1:r);

        Yks = Uk*S(1:r,1:r)*Vk';
        Fks = f(Yks);
        [Uks,~] = qr([U Fk1*V Uk Fks*Vk],0);
        [Vks,~] = qr([V Fk1'*U Vk Fks'*Uk],0);
        Sks = Uks'*(Y + 0.5*h*Fk1 + 0.5*h*Fks)*Vks;
        [P,S,W] = svd(Sks,'econ');
        U = Uks*P(:,1:r);
        V = Vks*W(:,1:r);
        S = S(1:r,1:r);
    end   
    Y = U*S*V';
end 