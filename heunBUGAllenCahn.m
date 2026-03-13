function Y = heunBUGAllenCahn(Y,N,h,L,theta,r) %versione heun BUG ottimizzata per Allen-Cahn
    [U,S,V] = svd(Y);
    U = U(:,1:r);
    V = V(:,1:r);
    S = S(1:r,1:r);
    % 
    % n = size(Y,1);
    % Yk = zeros(n,n,N);
    % Yk(:,:,1) = Y;
    for k = 1:N
       
        Us = U*S;
        Vs = V*S;
        Fu1 = theta*(L*Us + Us*(V'*L*V)) + Us - ((Us*V').*(Us*V').*(Us*V'))*V;
        Fv1 = theta*(Vs*(U'*L*U) + L*Vs) + Vs -  ((Us*V').*(Us*V').*(Us*V'))'*U;
        [Uk1,~] = qr([U Fu1],0);
        [Vk1,~] = qr([V Fv1],0);

        Uk1u = Uk1'*U;
        Vk1v = V'*Vk1;
        Fk1 = (Uk1'*L*U)*S*Vk1v;
        Fk2 = Uk1u*S*(V'*L*Vk1);
        Fk3 = Uk1u*S*Vk1v;
        Fk4 = Uk1'*((Us*V').*(Us*V').*(Us*V'))*Vk1;

        Sk1 = Uk1u*S*Vk1v + h*(theta*(Fk1 + Fk2) + Fk3 - Fk4);
        [P,S1,W] = svd(Sk1,"econ");
        Uk = Uk1*P(:,1:r);
        Vk = Vk1*W(:,1:r);
        S1 = S1(1:r,1:r);

        Uks = Uk*S1;
        Vks = Vk*S1;
        Fus = theta*(L*Uks + Uks*(Vk'*L*Vk)) + Uks - ((Uks*Vk').*(Uks*Vk').*(Uks*Vk'))*Vk;
        Fvs = theta*(Vks*(Uk'*L*Uk) + L*Vks) + Vks - ((Uks*Vk').*(Uks*Vk').*(Uks*Vk'))'*Uk;
 
        [Uks,~] = qr([U Fu1 Uk Fus],0);
        [Vks,~] = qr([V Fv1 Vk Fvs],0);

        Uksu = Uks'*U;
        Vksv = V'*Vks;
        Fk1 = (Uks'*L*U)*S*Vksv;
        Fk2 = Uksu*S*(V'*L*Vks);
        Fk3 = Uksu*S*Vksv;
        Fk4 = Uks'*((Us*V').*(Us*V').*(Us*V'))*Vks;
        
        Uksuk = Uks'*Uk;
        Vksvk = Vk'*Vks;
        Fks1 = (Uks'*L*Uk)*S1*Vksvk;
        Fks2 = Uksuk*S1*(Vk'*L*Vks);
        Fks3 = Uksuk*S1*Vksvk;
        Fks4 = Uks'*((Uk*S1*Vk').*(Uk*S1*Vk').*(Uk*S1*Vk'))*Vks;

        Sks = Uksu*S*Vksv + 0.5*h*(theta*(Fk1 + Fk2 + Fks1 + Fks2) + Fk3 + Fks3 - Fk4 - Fks4);
        [P,S,W] = svd(Sks,'econ');
        U = Uks*P(:,1:r);
        V = Vks*W(:,1:r);
        S = S(1:r,1:r);

        %Yk(:,:,k+1) = U*S*V';
    end   
    Y = U*S*V'; 
end