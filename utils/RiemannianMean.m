function [M, tC] = RiemannianMean(tC)

Np = size(tC, 3);
M  = mean(tC, 3);
kk=[];
for jj = 1 : Np
    if any(svd(tC(:,:,jj))<1e-8)
        disp(jj);
        kk(end+1)=jj;
        tC(:,:,jj)=tC(:,:,jj) + eye(size(tC,1))*0.01;
    end
end
    
for ii = 1 : 200
    A = M ^ (1/2);      %-- A = C^(1/2)
    B = A ^ (-1);       %-- B = C^(-1/2)
        
    S = zeros(size(M));
    for jj = 1 : Np
        C = tC(:,:,jj);
        S = S + A * logm(B * C * B) * A;
    end
    S = S / Np;
    
    M = A * expm(B * S * B) * A; 
    
    eps = norm(S, 'fro')
    if (eps < 1e-6)
        break;
    end
end

end