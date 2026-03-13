function [Y,S,U,V] = dlr(Y0,F,r,h,T)
    tol = 1e-14;
    [U,S,V] = svd(Y0);
    U = U(:,1:r);
    V = V(:,1:r);
    S = S(1:r,1:r);
    W = F(Y0);
    Nt = floor(T/h);
    I = eye(size(Y0,1));

    %Yk(:,:,1) = Y0;
    for n = 1:Nt
        Sinv = 1./diag(max(S,tol));
        U = U + h*(I - U*U')*W*V*Sinv;
        S = S + h*(U'*W*V); 
        V = V + h*(I - V*V')*W'*U*Sinv;

        %Yk(:,:,n+1) = U*S*V';
        W = F(U*S*V');
    end
    Y = U*S*V';
end