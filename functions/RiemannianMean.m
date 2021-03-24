function [M, tC] = RiemannianMean(tC)
Np=size(tC,3);
K = (size(tC,1 ));
M  = mean(tC, 3);
for jj = 1 : Np
    if any(svd(tC(:,:,jj))<1e-10)
        tC(:,:,jj)=tC(:,:,jj)+eye(K)*0.01;
    end
end
for ii = 1 : 200
    A = M ^ (1/2);      %-- A = C^(1/2)
    B = A ^ (-1);       %-- B = C^(-1/2)
    NN=0;    
    S = zeros(size(M));
    for jj = 1 : Np
        C = tC(:,:,jj);
        del = logm(B * C * B);
    
        S = S + A * del * A;
        NN=NN+1;
    end
    S = S / NN;
    
    M = A * expm(B * S * B) * A; 
    
    eps = norm(S, 'fro')
    if (eps < 1e-6)
        break;
    end
end

end