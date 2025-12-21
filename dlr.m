function Y = dlr(Y,F,r,h,T)
    [U,S,V] = svd(Y);
    S = S(1:r, 1:r);
    U = U(:,1:r);
    V = V(:,1:r);
    Nt = floor(T/h);
    Im = eye(size(U,1));

    for n = 1:Nt-1
        Z = F(Y);
        S1 = U'*Z*V;
        U1s = (Im - U*U')*Z*V;
        V1s = (Im - V*V')*Z'*U;
        U1 = qr(U1s);
        V1 = qr(V1s);

        Y = Y + h*(U1*S*V' + U*S1*V' + U*S*V1');
        [U,S,V] = svd(Y);
        S = S(1:r,1:r);
        U = U(:,1:r);
        V = V(:,1:r);
    end
end