function A = eulesp(A,F,h,T)
    Nt = floor(T/h);
    % m = size(A,1);
    % Ak = zeros(m,m,Nt);
    % Ak(:,:,1) = A;
    for n = 1:Nt
        Fa = F(A);
        A = A + h*Fa;
        %Ak(:,:,n+1) = A;
    end
end