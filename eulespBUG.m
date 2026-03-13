function [Yk,S] = eulespBUG(Y,N,h,f,r)
    tol = 1e-12;
    [U,~,V] = svd(Y);
    U = U(:,1:r);
    V = V(:,1:r);
    n = size(Y,1);
    Yk = zeros(n,n,N);
    Yk(:,:,1) = Y;
    for k = 1:N
        Fk = f(Y);
        Uk = orth([U Fk*V],tol);
        Vk = orth([V Fk'*U],tol);
        Sk = Uk'*(Y+h*Fk)*Vk;
        [P,S,W] = svd(Sk);
        U = Uk*P(:,1:r);
        V = Vk*W(:,1:r);
        S = S(1:r,1:r);
        Y = U*S*V';

        Yk(:,:,k+1) = Y;
    end
end

